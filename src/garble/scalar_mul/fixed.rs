//! Scalar mul [x]P where x is a scalar, P is a fixed curve point
use ark_ec::{
    AdditiveGroup, AffineRepr, CurveConfig, CurveGroup, ScalarMul,
    scalar_mul::glv::GLVConfig,
    short_weierstrass::{Affine, Projective},
};
use ark_ff::Field;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{UniformRand, rand::rngs::StdRng};
use std::ops::Mul;

use crate::{
    garble::{Circuit, Nonce, bit_adaptor},
    mac::{
        IteratedMac, Random, RandomFromPhi,
        bit::{bytes::Bit16ByteMac, group::BitGroupMac},
        digit::{bytes, group},
        endo::{EndoMac, EndoMacKey, EndoScalar},
        group::{GroupMac, GroupMacKey},
    },
    utils::{
        PRG,
        bits::{FromBits, NumBits, ToBits},
        group::FromLeBytesModOrder,
    },
};

pub(crate) type Input<C> = <C as CurveConfig>::ScalarField;
pub(crate) type InputMacKey<C> = bytes::DigitMacKey<<C as CurveConfig>::ScalarField>;
pub(crate) type InputMac<C> = bytes::DigitMac<<C as CurveConfig>::ScalarField>;

pub(crate) type Output<C> = Affine<C>;
pub(crate) type OutputMac<C, const KAPPA: usize> = IteratedMac<Affine<C>, EndoMac<C>, KAPPA>;
// { \phi_i, b • c_i }
pub(crate) type OutputMacKey<C, const KAPPA: usize> = [EndoMacKey<C>; KAPPA];

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub(crate) struct Table<C: GLVConfig, const KAPPA: usize>
where
    C::ScalarField: Random,
    <C::ScalarField as Field>::BasePrimeField: FromLeBytesModOrder,
{
    adaptors: Vec<bit_adaptor::Table<Bit16ByteMac, BitGroupMac<C::ScalarField>>>,
    point: Affine<C>,
    // S_i = [\phi_i]P + [a • c_i]P
    s: Vec<Affine<C>>,
    // T_i = [c_i]P
    t: Vec<Affine<C>>,
}

impl<C: GLVConfig + Clone, const KAPPA: usize> Table<C, KAPPA>
where
    C::ScalarField: Random,
    <C::ScalarField as Field>::BasePrimeField: FromLeBytesModOrder,
{
    pub(crate) fn output_mac_key_from_phis(
        rng: &mut StdRng,
        point: Projective<C>,
        input_mac_key: &GroupMacKey<C::ScalarField>,
        phis: &[EndoScalar<C>; KAPPA],
    ) -> OutputMacKey<C, KAPPA> {
        let mut randomised_ks = vec![C::ScalarField::ZERO; KAPPA];
        randomised_ks.iter_mut().for_each(|k| {
            let c = C::ScalarField::rand(rng);
            *k = input_mac_key.k * c
        });
        let k = point.batch_mul(&randomised_ks);

        let mut keys: Vec<EndoMacKey<C>> = Vec::with_capacity(KAPPA);
        phis.iter().zip(k).for_each(|(phi, k)| {
            keys.push(EndoMacKey {
                phi: phi.clone(),
                k,
            })
        });

        keys.try_into().unwrap_or_else(|_| panic!())
    }
}

impl<C, const KAPPA: usize> Circuit<Input<C>, InputMac<C>, Output<C>, OutputMac<C, KAPPA>>
    for Table<C, KAPPA>
where
    C: GLVConfig + Clone,
    C::ScalarField: Random + FromBits<C::ScalarField, C::ScalarField> + ToBits + NumBits,
    <C::ScalarField as Field>::BasePrimeField: FromLeBytesModOrder,
    EndoMacKey<C>: RandomFromPhi<Affine<C>, EndoMac<C>>,
{
    type Aux = (Affine<C>, [EndoScalar<C>; KAPPA]);

    fn garble<P: PRG>(
        parent_nonce: &mut Nonce,
        prg_unseeded: &P::Unseeded,
        rng: &mut StdRng,
        aux: Option<&Self::Aux>,
        input_mac_key: &InputMacKey<C>,
    ) -> (Self, OutputMacKey<C, KAPPA>) {
        let (point, phis) = aux.unwrap();
        let mut nonce = parent_nonce.fork();

        // Adaptor tables for the scalar
        let phi = C::ScalarField::rand(rng);
        let (adaptors, input_mac_key): (Vec<_>, Vec<_>) = input_mac_key
            .iter()
            .map(|bit_mac_key| {
                bit_adaptor::Table::<Bit16ByteMac, BitGroupMac<C::ScalarField>>::garble::<P>(
                    &mut nonce,
                    prg_unseeded,
                    rng,
                    Some(&phi),
                    bit_mac_key,
                )
            })
            .collect();

        let input_mac_key: group::DigitMacKey<C::ScalarField> = input_mac_key.into();
        let input_mac_key: GroupMacKey<C::ScalarField> = (&input_mac_key).into();

        let a = input_mac_key.phi;
        let b = input_mac_key.k;
        let b_inv = b.inverse().unwrap();

        let output_mac_key =
            Self::output_mac_key_from_phis(rng, point.into_group(), &input_mac_key, phis);

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
                adaptors,
                point: *point,
                s,
                t,
            },
            output_mac_key.clone(),
        )
    }

    fn evaluate<P: PRG>(
        &self,
        parent_nonce: &mut Nonce,
        prg_unseeded: &P::Unseeded,
        input: &Input<C>,
        input_mac: &InputMac<C>,
    ) -> (Output<C>, OutputMac<C, KAPPA>) {
        let mut nonce = parent_nonce.fork();

        let input_mac: Vec<_> = self
            .adaptors
            .iter()
            .zip(input.to_bits().iter())
            .zip(input_mac.iter())
            .map(|((adaptor, input), input_mac)| {
                adaptor
                    .evaluate::<P>(&mut nonce, prg_unseeded, input, input_mac)
                    .1
            })
            .collect();

        let input_mac: group::DigitMac<C::ScalarField> = input_mac.into();
        let input_mac: GroupMac<C::ScalarField> = (&input_mac).into();

        let output_mac: Box<[EndoMac<C>; KAPPA]> = self
            .t
            .iter()
            .zip(self.s.iter())
            .map(|(t, s)| EndoMac::from((t.mul(*input_mac) + s.mul(input)).into_affine()))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap_or_else(|_| panic!());
        (self.point.mul(input).into_affine(), output_mac.into())
    }
}

#[test]
fn test_scalar_mul_fixed() {
    use crate::{garble::pairing::KAPPA_B, mac::MacKey, utils::prg::aes_ctr::clear_set};
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

    let input_mac = input_mac_key.encode(&input);
    let (circuit, output_mac_key) = Table::garble::<Aes128>(
        &mut Nonce::new(),
        &unseeded_prg,
        &mut rng,
        Some(&(point.into_affine(), phis)),
        &input_mac_key,
    );
    clear_set();

    let expected_output_mac = output_mac_key.encode(&point.mul(input).into_affine());

    let output_mac = circuit
        .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
        .1;
    clear_set();

    assert!(output_mac == expected_output_mac);
}
