//! Garbling table for a Groth16 verifier that commits to a secret `s`.
//! This table reveals `s` if it is evaluated on:
//! - (A, B, C) that does not satisfy e(A, B) e(C, H_1) = c^{-1}, where
//!   H_1, c are from the verifier key;
//! - (A, B, C) where any of the points is an invalid curve point.
//!
//! This implementation does not support Groth16 public inputs.

use core::panic;
use std::marker::PhantomData;

use aes::Aes128;
use ark_bn254::{Bn254, Fq, Fq2, Fq12, G1Affine, G1Projective, G2Affine, G2Projective, g1, g2};
use ark_ec::{AffineRepr, pairing::Pairing};
use ark_ff::Field;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    Zero,
    rand::{SeedableRng, rngs::StdRng},
};

use crate::{
    garble::{
        Circuit, Nonce,
        bivariate::{self, curve_membership},
        conditional_reveal, field_mac_to_ec_mac,
        g2_subgroup_check::check_g2_membership,
        pairing::{self, KAPPA_A, KAPPA_B},
    },
    mac::{
        IdentityMac, IdentityMacKey, Mac, MacKey, RandomFromPhi,
        endo::{EndoMacKey, EndoMultGroupMac},
    },
    utils::prg::PRG,
};

#[cfg(test)]
use crate::utils::prg::aes_ctr::clear_set;
#[cfg(test)]
use ark_ff::UniformRand;

/// Groth16 verifier key
#[derive(Debug, Clone, Copy, CanonicalSerialize, CanonicalDeserialize)]
pub struct VerifierKey {
    pub G1: G1Projective,
    pub H1: G2Projective,
    pub c: Fq12,
}

impl VerifierKey {
    /// Given the input (A, B, C), returns e(A, B) e(C, H_1) * c.
    /// Recall that if (A, B, C) is a valid Groth16 proof, it should satisfy
    ///                 e(A, B) e(C, H_1) = c^{-1}.
    /// In other words, this method should return Fq12::ONE on a valid proof.
    fn pairing(&self, input: &Input) -> Fq12 {
        Bn254::multi_pairing(
            vec![input.A.into_group(), input.C.into_group()],
            vec![input.B.into_group(), self.H1],
        )
        .0 * self.c
    }
}

/// Groth16 proof. This implementation does not support Groth16 public inputs.
#[derive(Debug, Clone, Copy)]
pub struct Input {
    pub A: G1Affine,
    pub B: G2Affine,
    pub C: G1Affine,
}

/// Input MAC key from the Groth16 proof.
#[derive(Clone)]
pub struct InputMacKey {
    A: bivariate::InputMacKey<Fq>,
    B: bivariate::InputMacKey<Fq2>,
    C: bivariate::InputMacKey<Fq>,
}

/// Input MAC of the Groth16 proof.
#[derive(Debug, Clone)]
pub struct InputMac {
    A: bivariate::InputMac<Fq>,
    B: bivariate::InputMac<Fq2>,
    C: bivariate::InputMac<Fq>,
}

type Output = conditional_reveal::Output;

impl MacKey<Input, InputMac> for InputMacKey {
    fn encode(&self, value: &Input) -> InputMac {
        InputMac {
            A: self.A.encode(&value.A.into()),
            B: self.B.encode(&value.B.into()),
            C: self.C.encode(&value.C.into()),
        }
    }
    #[cfg(test)]
    fn decode(&self, mac: &InputMac) -> Input {
        Input {
            A: self.A.decode(&mac.A).into(),
            B: self.B.decode(&mac.B).into(),
            C: self.C.decode(&mac.C).into(),
        }
    }
    fn random_key(rng: &mut StdRng) -> Self {
        Self {
            A: bivariate::InputMacKey::random_key(rng),
            B: bivariate::InputMacKey::random_key(rng),
            C: bivariate::InputMacKey::random_key(rng),
        }
    }
}

impl Mac<Input> for InputMac {
    type Key = InputMacKey;

    #[cfg(test)]
    fn random_input(rng: &mut StdRng) -> Input {
        use ark_ff::UniformRand;

        Input {
            A: G1Affine::rand(rng),
            B: G2Affine::rand(rng),
            C: G1Affine::rand(rng),
        }
    }
}

/// Given an RNG seed, a domain separation tag for the fixed-key AES hash, and
/// the verifier key, returns a garbling table and the input MAC key used
pub fn garble_groth16(rng_seed: u64, verifier_key: VerifierKey) -> (Table, InputMacKey) {
    let mut rng = StdRng::seed_from_u64(rng_seed);
    let mut nonce = Nonce::new();
    let unseeded_prg = Aes128::init(&[]);
    let input_mac_key = InputMacKey::random_key(&mut rng);

    let table = Table::garble::<Aes128>(
        &mut nonce,
        &unseeded_prg,
        &mut rng,
        Some(&verifier_key),
        &input_mac_key,
    )
    .0;

    #[cfg(test)]
    clear_set();

    (table, input_mac_key)
}

#[derive(CanonicalSerialize, CanonicalDeserialize)]
/// This table is structured as a collection of smaller modular tables.
pub struct Table {
    // Maps input MAC of (A.x, A.y) to iterated MAC of A projective coordinates
    a_ec_mac: field_mac_to_ec_mac::Table<g1::Config, KAPPA_A>,
    // Checks that A is a valid G1 curve point, leaking `secret` if not
    a_curve_check: curve_membership::Table<g1::Config>,
    // Maps input MAC of (B.x, B.y) to iterated MAC of B projective coordinates
    b_ec_mac: field_mac_to_ec_mac::Table<g2::Config, KAPPA_B>,
    // Checks that B is a valid G2 curve point, leaking `secret` if not
    b_curve_check: curve_membership::Table<g2::Config>,
    // Maps input MAC of (C.x, C.y) to iterated MAC of C projective coordinates
    c_ec_mac: field_mac_to_ec_mac::Table<g1::Config, KAPPA_B>,
    // Checks that C is a valid G1 curve point, leaking `secret` if not
    c_curve_check: curve_membership::Table<g1::Config>,
    // e(A,B)
    pairing: pairing::Table,
    // Leaks secret if evaluated on an input != 1 and its MAC
    pub conditional_reveal: conditional_reveal::Table<KAPPA_B>,
    // Groth16 verifier key
    verifier_key: VerifierKey,
}

impl Circuit<Input, InputMac, Output, IdentityMac<Output>> for Table {
    type Aux = VerifierKey;

    fn garble<P: PRG>(
        parent_nonce: &mut Nonce,
        prg_unseeded: &P::Unseeded,
        rng: &mut StdRng,
        aux: Option<&Self::Aux>,
        input_mac_key: &InputMacKey,
    ) -> (Self, <IdentityMac<Output> as Mac<Output>>::Key) {
        // The Groth16 verifier key
        let verifier_key = aux.unwrap();

        // Sample random MAC keys and auxiliary values for the pairing table:
        // - a_ec_mac_key: used as the output key for field_to_ec_mac of `A`
        // - b_ec_mac_key: used as the output key for field_to_ec_mac of `B`
        // - pairing_output_mac_key: used to derive the input key for conditional_reveal
        // - powers_of_tau: used in the pairing e(A,B) garbling table
        let ((a_ec_mac_key, b_ec_mac_key), pairing_output_mac_key, powers_of_tau) =
            pairing::Table::random_mac_keypair_and_tau(rng);

        // The secret that will be conditionally revealed.
        let secret = conditional_reveal::Table::derive_secret(&pairing_output_mac_key);

        // Construct c_ec_mac_key to use the same `phi`s as pairing_output_mac_key.
        // c_ec_mac_key is used as the output key for field_to_ec_mac of `C`
        let mut c_ec_mac_key = Vec::with_capacity(KAPPA_B);
        pairing_output_mac_key
            .iter()
            .for_each(|key| c_ec_mac_key.push(EndoMacKey::random_from_phi(rng, &key.phi)));
        let c_ec_mac_key: [EndoMacKey<g1::Config>; KAPPA_B] =
            c_ec_mac_key.try_into().unwrap_or_else(|_| panic!());

        // Garbling table for field_to_ec_mac of `A`
        let a_ec_mac = field_mac_to_ec_mac::Table::garble::<P>(
            parent_nonce,
            prg_unseeded,
            rng,
            // The a_ec_mac_key that we sampled earlier
            Some(&a_ec_mac_key),
            // Input MAC key for the (x, y) coordinates of `A`
            &input_mac_key.A,
        )
        .0;

        // Garbling table to check that `A` is a valid G1 curve point.
        // This table leaks `secret` if `A` is not a valid curve point.
        let a_curve_check = curve_membership::Table::garble::<P>(
            parent_nonce,
            prg_unseeded,
            rng,
            Some(&secret),
            // Input MAC key for the (x, y) coordinates of `A`
            &input_mac_key.A,
        )
        .0;

        // Garbling table for field_to_ec_mac of `B`
        let b_ec_mac = field_mac_to_ec_mac::Table::garble::<P>(
            parent_nonce,
            prg_unseeded,
            rng,
            // The b_ec_mac_key that we sampled earlier
            Some(&b_ec_mac_key),
            // Input MAC key for the (x, y) coordinates of `B`
            &input_mac_key.B,
        )
        .0;

        // Garbling table to check that `B` is a valid G2 curve point.
        // This table leaks `secret` if `B` is not a valid curve point.
        let b_curve_check = curve_membership::Table::garble::<P>(
            parent_nonce,
            prg_unseeded,
            rng,
            Some(&Fq2::from_base_prime_field(secret)),
            // Input MAC key for the (x, y) coordinates of `B`
            &input_mac_key.B,
        )
        .0;

        // Garbling table for field_to_ec_mac of `C`
        let c_ec_mac = field_mac_to_ec_mac::Table::garble::<P>(
            parent_nonce,
            prg_unseeded,
            rng,
            // The c_ec_mac_key that we sampled earlier
            Some(&c_ec_mac_key),
            // Input MAC key for the (x, y) coordinates of `C`
            &input_mac_key.C,
        )
        .0;

        // Garbling table to check that `C`` is a valid G1 curve point.
        // This table leaks `secret` if `C`` is not a valid curve point.
        let c_curve_check = curve_membership::Table::garble::<P>(
            parent_nonce,
            prg_unseeded,
            rng,
            Some(&secret),
            // Input MAC key for the (x, y) coordinates of `C`
            &input_mac_key.C,
        )
        .0;

        // Garbling table for the pairing e(A, B). It uses a_ec_mac_key,
        // b_ec_mac_key as input MAC keys, and pairing_output_mac_key as the output
        // MAC key.
        let pairing = pairing::Table::garble::<P>(
            parent_nonce,
            prg_unseeded,
            rng,
            Some(&(powers_of_tau, pairing_output_mac_key.clone())),
            &(a_ec_mac_key, b_ec_mac_key),
        )
        .0;

        // Garbling table for the conditional revelation of `secret`.
        let conditional_reveal = {
            // The input MAC key for e(A, B) e(C, H_1) * c. Recall that if
            // (A, B, C) is a valid Groth16 proof, it should satisfy
            //               e(A, B) e(C, H_1) = c^{-1}.
            // In other words, this should be the input MAC key for the value
            // Fq12::ONE.
            let input_mac_key = conditional_reveal::Table::groth16_input_mac_key(
                verifier_key,
                &pairing_output_mac_key,
                &c_ec_mac_key,
            );

            conditional_reveal::Table::garble::<P>(
                parent_nonce,
                prg_unseeded,
                rng,
                None,
                &input_mac_key,
            )
            .0
        };

        (
            Self {
                a_ec_mac,
                a_curve_check,
                b_ec_mac,
                b_curve_check,
                c_ec_mac,
                c_curve_check,
                pairing,
                conditional_reveal,
                verifier_key: *verifier_key,
            },
            IdentityMacKey,
        )
    }

    fn evaluate<P: PRG>(
        &self,
        parent_nonce: &mut Nonce,
        prg_unseeded: &P::Unseeded,
        input: &Input,
        input_mac: &InputMac,
    ) -> (Output, IdentityMac<Output>) {
        // The Groth16 verifier key.
        let VerifierKey { H1, .. } = self.verifier_key;

        // The Groth16 input (A, B, C).
        let Input { A, B, C } = *input;

        // Evaluate the field_to_ec_mac table for `A`.
        let a_ec_mac = self
            .a_ec_mac
            .evaluate::<P>(parent_nonce, prg_unseeded, &A.into(), &input_mac.A)
            .1;

        // Evaluate the curve validity check table for `A`. This evaluates to
        // `secret` if `A` is not a valid curve point.
        let a_curve_check = self
            .a_curve_check
            .evaluate::<P>(parent_nonce, prg_unseeded, &A.into(), &input_mac.A)
            .0;
        // If `A` is not a valid curve point, return `secret` early
        if !a_curve_check.is_zero() {
            return (a_curve_check, IdentityMac(PhantomData));
        }

        // Evaluate the field_to_ec_mac table for `B`.
        let b_ec_mac = self
            .b_ec_mac
            .evaluate::<P>(parent_nonce, prg_unseeded, &B.into(), &input_mac.B)
            .1;

        // Evaluate the curve validity check table for `B`. This evaluates to
        // `secret` if `B` is not a valid curve point.
        let b_curve_check = self
            .b_curve_check
            .evaluate::<P>(parent_nonce, prg_unseeded, &B.into(), &input_mac.B)
            .0;
        // If `B` is not a valid curve point, return `secret` early
        if !b_curve_check.is_zero() {
            return (
                b_curve_check.to_base_prime_field_elements().next().unwrap(),
                IdentityMac(PhantomData),
            );
        } else if let Some(s) = check_g2_membership(
            B,
            &b_ec_mac
                .iter()
                .map(|endo_mac| **endo_mac)
                .collect::<Vec<_>>(),
        ) {
            return (s, IdentityMac(PhantomData));
        }

        // Evaluate the field_to_ec_mac table for `C`.
        let c_ec_mac = self
            .c_ec_mac
            .evaluate::<P>(parent_nonce, prg_unseeded, &C.into(), &input_mac.C)
            .1;
        // Evaluate the curve validity check table for `C`. This evaluates to
        // `secret` if `C` is not a valid curve point.
        let c_curve_check = self
            .c_curve_check
            .evaluate::<P>(parent_nonce, prg_unseeded, &C.into(), &input_mac.C)
            .0;
        // If `C` is not a valid curve point, return `secret` early
        if !c_curve_check.is_zero() {
            return (c_curve_check, IdentityMac(PhantomData));
        }

        // Evaluate the pairing of (A, B)
        let AB_mac = self
            .pairing
            .evaluate::<P>(parent_nonce, prg_unseeded, &(A, B), &(a_ec_mac, b_ec_mac))
            .1;

        // Compute the pairing e(A,B) e(C, H_1) on the input A, B, C
        let conditional_reveal_input = self.verifier_key.pairing(input);

        // Compute the pairing e(A,B) e(C, H_1) on the MACs of A, B, C
        let conditional_reveal_mac: Box<[EndoMultGroupMac; KAPPA_B]> = {
            let mut conditional_reveal_mac = Vec::with_capacity(KAPPA_B);
            (0..KAPPA_B).for_each(|idx| {
                let AB = &AB_mac[idx];
                let C = &c_ec_mac[idx];
                let mac = (Bn254::multi_pairing(vec![**C], vec![H1]).0 * **AB).into();
                conditional_reveal_mac.push(mac);
            });
            conditional_reveal_mac
                .try_into()
                .unwrap_or_else(|_| panic!())
        };

        // Evaluate the conditional reveal circuit.
        self.conditional_reveal.evaluate::<P>(
            parent_nonce,
            prg_unseeded,
            &conditional_reveal_input,
            &conditional_reveal_mac.into(),
        )
    }
}

#[test]
/// Evaluating the circuit on a valid SNARK and its MAC should return Fq::ZERO,
/// i.e. it should not reveal the `secret`.
fn test_valid_snark() {
    use ark_ff::AdditiveGroup;
    use ark_ff::Field;

    const SEED: u64 = 0;

    let unseeded_prg = Aes128::init(&[]);

    let mut rng = StdRng::seed_from_u64(SEED);

    let input = InputMac::random_input(&mut rng);
    let G1 = G1Projective::rand(&mut rng);
    let H1 = G2Projective::rand(&mut rng);

    // Construct `c` for the verifier key s.t. our chosen inputs are valid,
    // i.e. e(A,B) e(C, H_1) = c^{-1}
    let c = Bn254::multi_pairing(
        vec![input.A.into_group(), input.C.into_group()],
        vec![input.B.into_group(), H1],
    )
    .0
    .inverse()
    .unwrap();

    let verifier_key = VerifierKey { G1, H1, c };

    let input_mac_key = InputMacKey::random_key(&mut rng);
    let table = Table::garble::<Aes128>(
        &mut Nonce::new(),
        &unseeded_prg,
        &mut rng,
        Some(&verifier_key),
        &input_mac_key,
    )
    .0;
    clear_set();

    let input_mac = input_mac_key.encode(&input);
    let output = table
        .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
        .0;
    clear_set();

    // Assert that the circuit returns Fq::ZERO when evaluated on a valid input
    assert_eq!(output, Fq::ZERO);
}

#[test]
/// Evaluating the circuit on an invalid SNARK and its MAC should reveal the
/// `secret`.
fn test_invalid_groth16() {
    use ark_ec::PrimeGroup;
    use ark_ff::{BigInteger, PrimeField};

    const SEED: u64 = 0;
    let unseeded_prg = Aes128::init(&[]);

    let mut rng = StdRng::seed_from_u64(SEED);

    let input = InputMac::random_input(&mut rng);

    let verifier_key = VerifierKey {
        G1: G1Projective::rand(&mut rng),
        H1: G2Projective::rand(&mut rng),
        c: Bn254::pairing(G1Projective::rand(&mut rng), G2Projective::generator()).0,
    };

    let input_mac_key = InputMacKey::random_key(&mut rng);

    let table = Table::garble::<Aes128>(
        &mut Nonce::new(),
        &unseeded_prg,
        &mut rng,
        Some(&verifier_key),
        &input_mac_key,
    )
    .0;
    clear_set();

    let input_mac = input_mac_key.encode(&input);

    let output = table
        .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
        .0;
    clear_set();

    let hash_s = table.conditional_reveal.hash_s;

    let bytes = output.into_bigint().to_bytes_le();
    let mut hasher = blake3::Hasher::new();
    hasher.update(&bytes);
    let hash_output = *hasher.finalize().as_bytes();

    assert_eq!(hash_output.to_vec(), hash_s);
}

#[test]
/// If any of the input (A, B, C) is an invalid curve point, the circuit should
/// leak the same `secret` as in the case where (A, B, C) are a well-formed but
/// invalid Groth16 input.
fn test_malformed_curve_points() {
    use crate::garble::g2_subgroup_check::get_curve_point_not_in_g2;

    const SEED: u64 = 0;

    let unseeded_prg = Aes128::init(&[]);

    let mut rng = StdRng::seed_from_u64(SEED);

    // Evaluate the same circuit on differently malformed inputs, check that
    // we recover the same secret
    let input = InputMac::random_input(&mut rng);
    let G1 = G1Projective::rand(&mut rng);
    let H1 = G2Projective::rand(&mut rng);

    // e(A,B) e(C, H_1) = c^{-1}
    let c = Bn254::multi_pairing(
        vec![input.A.into_group(), input.C.into_group()],
        vec![input.B.into_group(), H1],
    )
    .0
    .inverse()
    .unwrap();

    let verifier_key = VerifierKey { G1, H1, c };

    let input_mac_key = InputMacKey::random_key(&mut rng);

    let table = Table::garble::<Aes128>(
        &mut Nonce::new(),
        &unseeded_prg,
        &mut rng,
        Some(&verifier_key),
        &input_mac_key,
    )
    .0;
    clear_set();

    // Malformed A
    let malformed_a = {
        let invalid_A = G1Affine::new_unchecked(Fq::rand(&mut rng), Fq::rand(&mut rng));
        let mut input = input;
        input.A = invalid_A;
        let input_mac = input_mac_key.encode(&input);

        table
            .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
            .0
    };
    clear_set();

    // Invalid curve point B
    let invalid_curve_point_b = {
        let invalid_B = G2Affine::new_unchecked(Fq2::rand(&mut rng), Fq2::rand(&mut rng));
        let mut input = input;
        input.B = invalid_B;
        let input_mac = input_mac_key.encode(&input);

        table
            .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
            .0
    };
    clear_set();

    // Valid curve point B but not in G2 subgroup
    let valid_b_not_in_g2 = {
        let invalid_B = get_curve_point_not_in_g2(&mut rng);
        let mut input = input;
        input.B = invalid_B;
        let input_mac = input_mac_key.encode(&input);

        table
            .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
            .0
    };
    clear_set();

    assert_eq!(invalid_curve_point_b, valid_b_not_in_g2);
    let malformed_b = invalid_curve_point_b;

    // Malformed C
    let malformed_c = {
        let invalid_C = G1Affine::new_unchecked(Fq::rand(&mut rng), Fq::rand(&mut rng));
        let mut input = input;
        input.C = invalid_C;
        let input_mac = input_mac_key.encode(&input);

        table
            .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
            .0
    };
    clear_set();

    // Well-formed points but invalid proof
    let well_formed_invalid = {
        let input = InputMac::random_input(&mut rng);

        let input_mac = input_mac_key.encode(&input);

        table
            .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
            .0
    };
    clear_set();

    assert!(!malformed_a.is_zero());
    assert_eq!(malformed_a, malformed_b);
    assert_eq!(malformed_b, malformed_c);
    assert_eq!(malformed_c, well_formed_invalid);
}

#[test]
/// Roundtrip serialisation and deserialisation of the table.
fn test_groth16_serde_roundtrip() {
    use crate::garble::Circuit;
    use crate::garble::Nonce;
    use crate::mac::MacKey;
    use ark_bn254::{G1Projective, G2Projective};
    use ark_ec::PrimeGroup;
    use ark_ff::UniformRand;
    use ark_std::rand::rngs::StdRng;
    use std::fs::File;
    use std::io::{BufReader, BufWriter, Write};

    const SEED: u64 = 0;

    let unseeded_prg = Aes128::init(&[]);

    let mut rng = StdRng::seed_from_u64(SEED);

    let input = InputMac::random_input(&mut rng);

    let verifier_key = VerifierKey {
        G1: G1Projective::rand(&mut rng),
        H1: G2Projective::rand(&mut rng),
        c: Bn254::pairing(G1Projective::rand(&mut rng), G2Projective::generator()).0,
    };

    let (original_table, input_mac_key) = garble_groth16(SEED, verifier_key);

    let input_mac: InputMac = input_mac_key.encode(&input);

    let original_output = original_table
        .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
        .0;
    clear_set();

    std::fs::create_dir_all("bin").expect("Could not create directory /bin");
    let file_name = "bin/groth16_table.bin";
    let file = File::create(file_name).expect("Could not create file");
    let mut writer = BufWriter::new(file);
    original_table.serialize_compressed(&mut writer).unwrap();
    writer.flush().unwrap();

    drop(original_table);

    let file = File::open(file_name).expect("Could not open file");
    let mut reader = BufReader::new(file);
    let deserialised_table = Table::deserialize_compressed(&mut reader)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
        .unwrap();
    let deserialised_output = deserialised_table
        .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
        .0;
    clear_set();

    assert_eq!(original_output, deserialised_output);
}
