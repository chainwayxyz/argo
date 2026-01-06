//! Scalar mul [x]P where x is a scalar, P is a fixed curve point
use ark_ec::{
    AdditiveGroup, CurveConfig, CurveGroup, ScalarMul,
    scalar_mul::glv::GLVConfig,
    short_weierstrass::{Affine, Projective},
};
use ark_ff::Field;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{UniformRand, rand::rngs::StdRng};
use std::ops::Mul;

use crate::{
    garble::{Circuit, Nonce},
    mac::{
        IteratedMac, Random, RandomFromPhi,
        endo::{EndoMac, EndoMacKey, EndoScalar},
        group::{GroupMac, GroupMacKey},
    },
    utils::PRG,
};

pub(crate) type Input<C> = <C as CurveConfig>::ScalarField;
pub(crate) type InputMacKey<C> = GroupMacKey<<C as CurveConfig>::ScalarField>;
pub(crate) type InputMac<C> = GroupMac<<C as CurveConfig>::ScalarField>;

pub(crate) type Output<C> = Affine<C>;
pub(crate) type OutputMac<C, const KAPPA: usize> = IteratedMac<Affine<C>, EndoMac<C>, KAPPA>;
// { \phi_i, b • c_i }
pub(crate) type OutputMacKey<C, const KAPPA: usize> = [EndoMacKey<C>; KAPPA];

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub(crate) struct Table<C: GLVConfig, const KAPPA: usize> {
    point: Affine<C>,
    // S_i = [\phi_i]P + [a • c_i]P
    s: Vec<Affine<C>>,
    // T_i = [c_i]P
    t: Vec<Affine<C>>,
}

impl<C: GLVConfig + Clone, const KAPPA: usize> Table<C, KAPPA> {
    pub(crate) fn output_mac_key_from_phis(
        rng: &mut StdRng,
        point: Projective<C>,
        input_mac_key: &InputMacKey<C>,
        phis: &[EndoScalar<C>; KAPPA],
    ) -> OutputMacKey<C, KAPPA> {
        let mut randomised_ks = vec![C::ScalarField::ZERO; KAPPA];
        randomised_ks.iter_mut().for_each(|k| {
            let c = C::ScalarField::rand(rng);
            *k = input_mac_key.k * c
        });
        let k = point.batch_mul(&randomised_ks);

        let mut keys = Vec::with_capacity(KAPPA);
        phis.iter().zip(k.into_iter()).for_each(|(phi, k)| {
            keys.push(EndoMacKey {
                phi: phi.clone(),
                k: k.into(),
            })
        });

        keys.try_into().unwrap_or_else(|_| panic!())
    }
}

impl<C, const KAPPA: usize> Circuit<Input<C>, InputMac<C>, Output<C>, OutputMac<C, KAPPA>>
    for Table<C, KAPPA>
where
    C: GLVConfig + Clone,
    C::ScalarField: Random,
    EndoMacKey<C>: RandomFromPhi<Affine<C>, EndoMac<C>>,
{
    type Aux = (Affine<C>, [EndoMacKey<C>; KAPPA]);

    fn garble<P: PRG>(
        _parent_nonce: &mut Nonce,
        _prg_unseeded: &P::Unseeded,
        _rng: &mut StdRng,
        aux: Option<&Self::Aux>,
        input_mac_key: &InputMacKey<C>,
    ) -> (Self, OutputMacKey<C, KAPPA>) {
        let (point, output_mac_key) = aux.unwrap();
        let a = input_mac_key.phi;
        let b = input_mac_key.k;
        let b_inv = b.inverse().unwrap();
        let t: Vec<Projective<C>> = output_mac_key.iter().map(|key| key.k.mul(b_inv)).collect();
        let t = Projective::<C>::normalize_batch(&t);

        let s: Vec<Projective<C>> = output_mac_key
            .iter()
            .map(|key| &key.phi)
            .zip(t.iter())
            .map(|(phi, ci_p)| phi.endo(*point) - ci_p.mul(a))
            .collect();
        let s = Projective::<C>::normalize_batch(&s);

        (
            Self {
                point: *point,
                s,
                t,
            },
            output_mac_key.clone(),
        )
    }

    fn evaluate<P: PRG>(
        &self,
        _parent_nonce: &mut Nonce,
        _prg_unseeded: &P::Unseeded,
        input: &Input<C>,
        input_mac: &InputMac<C>,
    ) -> (Output<C>, OutputMac<C, KAPPA>) {
        let output_mac: Box<[EndoMac<C>; KAPPA]> = self
            .t
            .iter()
            .zip(self.s.iter())
            .map(|(t, s)| EndoMac::from((t.mul(**input_mac) + s.mul(input)).into_affine()))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap_or_else(|_| panic!());
        (self.point.mul(input).into_affine(), output_mac.into())
    }
}

#[test]
fn test_scalar_mul_fixed() {
    use crate::{garble::pairing::KAPPA_B, mac::MacKey};
    use aes::Aes128;
    use ark_bn254::{Fr, G1Projective, g1::Config};
    use ark_std::rand::SeedableRng;
    use std::ops::Mul;

    const SEED: u64 = 0;

    let mut rng = StdRng::seed_from_u64(SEED);
    let unseeded_prg = Aes128::init(&[]);

    let point = G1Projective::rand(&mut rng);
    let input = Fr::rand(&mut rng);

    let input_mac_key = <InputMacKey<Config> as MacKey<Fr, _>>::random_key(&mut rng);
    let phis: [EndoScalar<Config>; KAPPA_B] = (0..KAPPA_B)
        .map(|_| EndoScalar::random(&mut rng))
        .collect::<Vec<_>>()
        .try_into()
        .unwrap_or_else(|_| panic!());

    let output_mac_key = Table::output_mac_key_from_phis(&mut rng, point, &input_mac_key, &phis);
    let expected_output_mac = output_mac_key.encode(&point.mul(input).into_affine());

    let input_mac = input_mac_key.encode(&input);
    let circuit = Table::garble::<Aes128>(
        &mut Nonce::new(),
        &unseeded_prg,
        &mut rng,
        Some(&(point.into_affine(), output_mac_key)),
        &input_mac_key,
    )
    .0;
    let output_mac = circuit
        .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
        .1;

    assert!(output_mac == expected_output_mac);
}
