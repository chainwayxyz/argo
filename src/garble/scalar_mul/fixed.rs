//! Scalar mul [x]P where x is a scalar, P is a fixed curve point
use ark_ec::{
    AdditiveGroup, CurveConfig, CurveGroup, PrimeGroup, ScalarMul,
    scalar_mul::glv::GLVConfig,
    short_weierstrass::{Affine, Projective},
};
use ark_ff::Field;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{UniformRand, rand::rngs::StdRng};
use std::{any::type_name, ops::Mul};

use crate::{
    garble::{Circuit, Nonce, bit_adaptor},
    mac::{
        IteratedMac, Random, RandomFromPhi,
        bit::{bytes::Bit16ByteMac, group::BitGroupMac},
        digit::{bytes, group},
        endo::{EndoMac, EndoMacKey, EndoScalar},
        group::GroupMacKey,
    },
    utils::{
        aes_hash::aes_ctr::AesCtr,
        bits::{FromBits, NumBits, ToBits},
        error::{EvalError, GarbleError},
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

#[derive(Debug, CanonicalSerialize, CanonicalDeserialize)]
pub(crate) struct Table<C: GLVConfig, const KAPPA: usize>
where
    C::ScalarField: Random,
    <C::ScalarField as Field>::BasePrimeField: FromLeBytesModOrder,
{
    pub(crate) adaptors: Vec<bit_adaptor::Table<AesCtr, Bit16ByteMac, BitGroupMac<C::ScalarField>>>,
    pub(crate) point: Affine<C>,
    // S_i = [\phi_i]P + [a • c_i]P
    pub(crate) s: [Affine<C>; KAPPA],
    // T_i = [c_i]P
    pub(crate) t: [Affine<C>; KAPPA],
}

impl<C: GLVConfig + Clone, const KAPPA: usize> Table<C, KAPPA>
where
    C::ScalarField: Random,
    <C::ScalarField as Field>::BasePrimeField: FromLeBytesModOrder,
{
    pub(crate) fn output_mac_key_from_phis(
        rng: &mut StdRng,
        input_mac_key: &GroupMacKey<C::ScalarField>,
        phis: &[EndoScalar<C>; KAPPA],
    ) -> OutputMacKey<C, KAPPA> {
        let point = Projective::<C>::generator();
        let mut k_scalars = vec![C::ScalarField::ZERO; KAPPA];
        k_scalars.iter_mut().for_each(|k| {
            let c = C::ScalarField::rand(rng);
            *k = input_mac_key.k * c
        });
        let k = point.batch_mul(&k_scalars);

        let mut keys: Vec<EndoMacKey<C>> = Vec::with_capacity(KAPPA);
        phis.iter()
            .zip(k_scalars)
            .zip(k)
            .for_each(|((phi, k_scalar), k)| {
                keys.push(EndoMacKey {
                    phi: phi.clone(),
                    k_scalar: Some(k_scalar),
                    k,
                })
            });

        // This unwrap cannot fail since `phi`, `k_scalars`, `k` are all
        // known to be length `KAPPA`.
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

    /// # Panics
    /// 
    /// Panics if reconstructed input_mac_key.k is C::ScalarField::ZERO.
    fn garble(
        parent_nonce: &mut Nonce,
        rng: &mut StdRng,
        aux: Option<&Self::Aux>,
        input_mac_key: &InputMacKey<C>,
    ) -> Result<(Self, OutputMacKey<C, KAPPA>), GarbleError> {
        let (point, phis) = aux.ok_or(GarbleError::AuxError(
            type_name::<Self::Aux>(),
            type_name::<Self>(),
        ))?;
        let mut nonce = parent_nonce.fork();

        // Adaptor tables for the scalar
        let phi = C::ScalarField::rand(rng);

        let mut adaptors = Vec::with_capacity(input_mac_key.len());
        let mut adapted_mac_key = Vec::with_capacity(input_mac_key.len());

        for bit_mac_key in input_mac_key.iter() {
            let (table, mac_key) = bit_adaptor::Table::<
                AesCtr,
                Bit16ByteMac,
                BitGroupMac<C::ScalarField>,
            >::garble(&mut nonce, rng, Some(&phi), bit_mac_key)?;
            adaptors.push(table);
            adapted_mac_key.push(mac_key);
        }

        let input_mac_key: group::DigitMacKey<C::ScalarField> = adapted_mac_key.into();
        let input_mac_key = input_mac_key
            .try_into_group_mac_key()
            .ok_or(GarbleError::DigitMacKeyReconstructionError)?;

        let a = input_mac_key.phi;
        let b = input_mac_key.k;
        // We assume that input_mac_key.k is nonzero.
        let b_inv = b.inverse().expect("Reconstructed input_mac_key.k must be nonzero.");

        let output_mac_key = Self::output_mac_key_from_phis(rng, &input_mac_key, phis);

        let t: Vec<Projective<C>> = output_mac_key.iter().map(|key| key.k.mul(b_inv)).collect();
        let t = Projective::<C>::normalize_batch(&t);

        let s: Vec<Projective<C>> = output_mac_key
            .iter()
            .map(|key| &key.phi)
            .zip(t.iter())
            .map(|(phi, ci_p)| phi.endo(*point) - ci_p.mul(a))
            .collect();
        let s = Projective::<C>::normalize_batch(&s);

        Ok((
            Self {
                adaptors,
                point: *point,
                // This unwrap cannot fail because `output_mac_key` is a fixed
                // array of `[EndoMacKey; KAPPA]`.
                s: s.try_into().unwrap(),
                // This unwrap cannot fail because `output_mac_key` is a fixed
                // array of `[EndoMacKey; KAPPA]`.
                t: t.try_into().unwrap(),
            },
            output_mac_key.clone(),
        ))
    }

    fn evaluate(
        &self,
        parent_nonce: &mut Nonce,
        input: &Input<C>,
        input_mac: &InputMac<C>,
    ) -> Result<(Output<C>, OutputMac<C, KAPPA>), EvalError> {
        let mut nonce = parent_nonce.fork();

        let mut adapted_input_mac = Vec::with_capacity(self.adaptors.len());

        for ((adaptor, input), input_mac) in self
            .adaptors
            .iter()
            .zip(input.to_bits().iter())
            .zip(input_mac.iter())
        {
            let mac = adaptor.evaluate(&mut nonce, input, input_mac)?.1;
            adapted_input_mac.push(mac);
        }

        let input_mac: group::DigitMac<C::ScalarField> = adapted_input_mac.into();
        let input_mac = input_mac.try_into_group_mac()?;

        let output_mac: Vec<_> = self
            .t
            .iter()
            .zip(self.s.iter())
            .map(|(t, s)| EndoMac::from((t.mul(*input_mac) + s.mul(input)).into_affine()))
            .collect();
        // This unwrap cannot fail because both `self.t`, `self.s` are fixed
        // arrays of length `KAPPA`.
        let output_mac: Box<[EndoMac<C>; KAPPA]> =
            output_mac.try_into().unwrap_or_else(|_| panic!());

        Ok((self.point.mul(input).into_affine(), output_mac.into()))
    }
}

#[test]
fn test_scalar_mul_fixed() {
    use crate::{
        garble::groth16_verifier::KAPPA_B, mac::MacKey, utils::aes_hash::aes_ctr::clear_set,
    };
    use ark_bn254::{Fr, G1Projective, g1::Config};
    use ark_std::rand::SeedableRng;
    use std::ops::Mul;

    const SEED: u64 = 0;

    let mut rng = StdRng::seed_from_u64(SEED);

    let point = G1Projective::rand(&mut rng);
    let input = Fr::rand(&mut rng);

    let input_mac_key = <InputMacKey<Config> as MacKey<Fr, _>>::random_key(&mut rng);
    let phis: [EndoScalar<Config>; KAPPA_B] = (0..KAPPA_B)
        .map(|_| EndoScalar::random(&mut rng))
        .collect::<Vec<_>>()
        .try_into()
        .unwrap_or_else(|_| panic!());

    let input_mac = input_mac_key.encode(&input);
    let (circuit, output_mac_key) = Table::garble(
        &mut Nonce::new(),
        &mut rng,
        Some(&(point.into_affine(), phis)),
        &input_mac_key,
    )
    .unwrap();
    clear_set();

    let expected_output_mac = output_mac_key.encode(&point.mul(input).into_affine());

    let output_mac = circuit
        .evaluate(&mut Nonce::new(), &input, &input_mac)
        .unwrap()
        .1;
    clear_set();

    assert!(output_mac == expected_output_mac);
}
