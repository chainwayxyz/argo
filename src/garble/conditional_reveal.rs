//! The conditional_reveal table stores the MAC of 1 `one_mac`, and a hash
//! of the secret `enc_s`.
//!
//! When evaluated on an input != one and its MAC, it returns the `secret`.

use std::marker::PhantomData;

use ark_bn254::{Fq, Fq12, g1::Config as G1Config};
use ark_ec::{AffineRepr, CurveGroup, PrimeGroup, scalar_mul::glv::GLVConfig};
use ark_ff::{
    AdditiveGroup, BigInteger, Field, Fp12, PrimeField, batch_inversion,
    field_hashers::{DefaultFieldHasher, HashToField},
};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

#[cfg(test)]
use crate::garble::groth16_verifier::reset_preprocessed_vk;
use crate::{
    garble::groth16_verifier::{self, PreprocessedVerifierKey},
    mac::{
        IdentityMac, IdentityMacKey, MacKey,
        endo::{ENDO_SIZE, EndoMacKey, EndoMultGroupMacKey, EndoScalar, IteratedEndoMultGroupMac},
    },
    utils::unwrap_single,
};

/// Input is a single P::TargetField
type Input = Fq12;
pub(crate) type InputMac<const KAPPA: usize> = IteratedEndoMultGroupMac<KAPPA>;
type InputMacKey<const KAPPA: usize> = [EndoMultGroupMacKey; KAPPA];

pub(crate) type Output = Fq;

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct Table<const KAPPA: usize> {
    pub(crate) one_mac: InputMac<KAPPA>,
    pub enc_s: Fq,
}

impl<const KAPPA: usize> Table<KAPPA> {
    pub(crate) fn groth16_input_mac_key(
        verifier_key: &groth16_verifier::VerifierKey,
        pairing_output_mac_key: &[EndoMultGroupMacKey; KAPPA],
        c_ec_mac_key: &[EndoMacKey<G1Config>; KAPPA],
        x_G_output_mac_key: &[EndoMacKey<G1Config>; KAPPA],
    ) -> InputMacKey<KAPPA> {
        let PreprocessedVerifierKey {
            h1_target_pow_2,
            h2_target_pow_2,
            ..
        } = verifier_key.preprocess();

        #[cfg(test)]
        reset_preprocessed_vk();

        let mut key = Vec::with_capacity(pairing_output_mac_key.len());

        // e(A,B) e(C, H_1) = c^{-1}
        pairing_output_mac_key
            .iter()
            .zip(c_ec_mac_key.iter())
            .zip(x_G_output_mac_key.iter())
            .for_each(|((AB, C), x_G)| {
                // All these MAC keys use the same phi
                let phi = C.phi.clone();

                // e(G, H1)^{C_k_scalar}
                let C_k_H1 =
                    Fq12::pow_with_table(&h1_target_pow_2, C.k_scalar.unwrap().into_bigint())
                        .unwrap();
                assert_eq!(
                    ark_bn254::G1Projective::generator() * C.k_scalar.unwrap(),
                    C.k.into_group()
                );

                // e(G, H2)^{xG_k_scalar}
                let xG_k_H2 =
                    Fq12::pow_with_table(&h2_target_pow_2, x_G.k_scalar.unwrap().into_bigint())
                        .unwrap();
                assert_eq!(
                    (ark_bn254::G1Projective::generator() * x_G.k_scalar.unwrap()).into_affine(),
                    x_G.k
                );

                // e(C.k, H1) e(xG.k, H2)
                let k = C_k_H1 * xG_k_H2;

                #[cfg(test)]
                {
                    use ark_bn254::Bn254;
                    use ark_ec::pairing::Pairing;

                    assert_eq!(C_k_H1, Bn254::pairing(C.k, verifier_key.H1).0);
                    assert_eq!(xG_k_H2, Bn254::pairing(x_G.k, verifier_key.H2).0);

                    assert_eq!(
                        k,
                        Bn254::multi_pairing(
                            vec![C.k, x_G.k],
                            vec![verifier_key.H1, verifier_key.H2]
                        )
                        .0
                    );
                }

                // c^{p^(2phi + 6)} = c^{-p^(2phi)}
                // This has the effect of multiplying the MAC'd value by c
                let c_term = phi.endo_frob(&verifier_key.c.inverse().unwrap());
                let k = k * AB.k * c_term;

                key.push(EndoMultGroupMacKey { phi, k });
            });

        key.try_into().unwrap_or_else(|_| panic!())
    }

    pub(crate) fn hash_packed_phis<I, C: GLVConfig>(phis: I) -> Fq
    where
        I: IntoIterator<Item = EndoScalar<C>>,
        I::IntoIter: DoubleEndedIterator,
    {
        let packed_phis = EndoScalar::<C>::pack(phis);
        let hasher = <DefaultFieldHasher<blake3::Hasher> as HashToField<Fq>>::new(&[]);
        hasher.hash_to_field::<1>(&packed_phis.into_bigint().to_bytes_le())[0]
    }
}

impl<const KAPPA: usize> Table<KAPPA> {
    pub(crate) fn garble(enc_s: Fq, input_mac_key: &InputMacKey<KAPPA>) -> (Self, IdentityMacKey) {
        let one_mac = input_mac_key.encode(&Fp12::ONE);

        (Self { one_mac, enc_s }, IdentityMacKey)
    }

    /// This method takes in four InputMac's. Each InputMac is a sequence of
    /// [EndoMac; KAPPA_B].
    ///
    /// For each position `i` in KAPPA_B, we check for a match with the one_mac
    /// in at least one of the InputMac's.
    ///
    /// If none of the InputMac's matches the one_mac at a given `i`, we can
    /// recover phi_i using one of the InputMac's.
    pub(crate) fn evaluate(
        &self,
        input: &Input,
        input_macs: &[InputMac<KAPPA>],
    ) -> (Output, IdentityMac<Output>) {
        let phis = EndoScalar::<G1Config>::iter();
        let input_exp = if *input == Fq12::ONE {
            vec![Fq12::ONE]
        } else {
            // Precompute the powers of input
            let mut input_exp = vec![Fq12::ZERO; ENDO_SIZE];

            input_exp
                .iter_mut()
                .zip(phis.iter())
                .for_each(|(input_exp, phi)| *input_exp = phi.endo_frob(input));

            input_exp
        };

        // If none of the input_macs is equal to one_mac, we can recover the `phi`s of
        // the input_mac_key in the following way:
        // - one_mac   = 1 ^ {phi_i} * k_i
        // - input_mac = input ^ {phi_i} * k_i
        // - input_mac / one_mac = input^{phi_i}
        let mut one_inv: Vec<_> = self.one_mac.iter().map(|v| **v).collect();
        batch_inversion(&mut one_inv);

        let input_mac_div_one = input_macs
            .iter()
            .map(|input_mac| {
                input_mac
                    .iter()
                    .zip(one_inv.iter())
                    .map(|(mac, one_inv)| **mac * one_inv)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        // For each input^{phi_i}, brute force the `phi`. This is efficient
        // because there are only seven possible values of `phi`.
        let mut recovered_phis = Vec::with_capacity(KAPPA);
        (0..KAPPA).for_each(|idx| {
            let find_position = |val: Fq12| input_exp.iter().position(|&x| x == val);
            let input_mac_div_one_slice =
                input_mac_div_one.iter().map(|mac| find_position(mac[idx]));
            let phi = unwrap_single(input_mac_div_one_slice);

            recovered_phis.push(phi);
        });

        let phis = if recovered_phis.iter().all(Option::is_some) {
            recovered_phis
                .iter()
                .map(|&idx| phis[idx.unwrap()].clone())
                .collect()
        } else {
            EndoScalar::<G1Config>::decode(&recovered_phis)
        };

        if *input == Fq12::ONE && phis.iter().all(|phi| matches!(phi, EndoScalar::Zero(_))) {
            return (Fq::ZERO, IdentityMac(PhantomData));
        }

        // Recover the secret using the phis.
        let s = self.enc_s - Self::hash_packed_phis(phis);
        (s, IdentityMac(PhantomData))
    }
}

#[test]
fn test_conditional_reveal_valid() {
    use crate::garble::groth16_verifier::KAPPA_B;
    use ark_std::{
        UniformRand,
        rand::{SeedableRng, rngs::StdRng},
    };

    const SEED: u64 = 0;

    let mut rng = StdRng::seed_from_u64(SEED);
    let secret = Fq::rand(&mut rng);

    let input = Fq12::ONE;
    let input_mac_key = InputMacKey::random_key(&mut rng);
    let input_mac = input_mac_key.encode(&input);

    let enc_s = secret
        + Table::<KAPPA_B>::hash_packed_phis(input_mac_key.iter().map(|key| key.phi.clone()));

    let table = Table::<KAPPA_B>::garble(enc_s, &input_mac_key).0;

    let output = table.evaluate(&input, &[input_mac]).0;

    assert_eq!(output, Fq::ZERO);
}

#[test]
fn test_conditional_reveal_invalid() {
    use crate::garble::groth16_verifier::KAPPA_B;
    use ark_std::{
        UniformRand,
        rand::{SeedableRng, rngs::StdRng},
    };

    const SEED: u64 = 0;

    let mut rng = StdRng::seed_from_u64(SEED);
    let secret = Fq::rand(&mut rng);

    let input = Input::rand(&mut rng);
    let input_mac_key = InputMacKey::random_key(&mut rng);
    let input_mac = input_mac_key.encode(&input);

    let enc_s = secret
        + Table::<KAPPA_B>::hash_packed_phis(input_mac_key.iter().map(|key| key.phi.clone()));

    let table = Table::<KAPPA_B>::garble(enc_s, &input_mac_key).0;

    let output = table.evaluate(&input, &[input_mac]).0;

    assert_eq!(output, secret);
}

#[test]
fn test_preprocessed_pow_two() {
    use crate::garble::groth16_verifier::{VerifierKey, reset_preprocessed_vk};
    use ark_bn254::{Bn254, Fr, G1Projective, G2Projective};
    use ark_ec::{PrimeGroup, pairing::Pairing};
    use ark_ff::UniformRand;
    use ark_std::rand::{SeedableRng, rngs::StdRng};

    const SEED: u64 = 0;

    let mut rng = StdRng::seed_from_u64(SEED);

    let verifier_key = VerifierKey {
        G1: G1Projective::rand(&mut rng),
        H1: G2Projective::rand(&mut rng),
        H2: G2Projective::rand(&mut rng),
        c: Bn254::pairing(G1Projective::rand(&mut rng), G2Projective::generator()).0,
    };

    let PreprocessedVerifierKey {
        h1_target_pow_2,
        h2_target_pow_2,
        ..
    } = verifier_key.preprocess();

    let scalar = Fr::rand(&mut rng);
    let point = G1Projective::generator() * scalar;

    let h1_pairing = Fq12::pow_with_table(&h1_target_pow_2, scalar.into_bigint()).unwrap();
    assert_eq!(h1_pairing, Bn254::pairing(point, verifier_key.H1).0);

    let h2_pairing = Fq12::pow_with_table(&h2_target_pow_2, scalar.into_bigint()).unwrap();
    assert_eq!(h2_pairing, Bn254::pairing(point, verifier_key.H2).0);

    reset_preprocessed_vk();
}
