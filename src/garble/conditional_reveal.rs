//! The conditional_reveal table stores the MAC of 1 `one_mac`, and a hash
//! of the secret `hash_s`.
//!
//! When evaluated on an input != one and its MAC, it returns the `secret`.

use std::marker::PhantomData;

use ark_bn254::{Bn254, Fq, Fq12, g1::Config as G1Config};
use ark_ec::pairing::Pairing;
use ark_ff::{AdditiveGroup, BigInteger, Field, Fp12, PrimeField, batch_inversion};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::rngs::StdRng;

use crate::{
    garble::Nonce,
    garble::{Circuit, groth16_verifier},
    mac::{
        IdentityMac, IdentityMacKey, MacKey,
        endo::{EndoMacKey, EndoMultGroupMacKey, IteratedEndoMultGroupMac},
    },
    utils::{bits::pow_n_lincomb, prg::PRG},
};

/// Input is a single P::TargetField
type Input = Fq12;
pub(crate) type InputMac<const KAPPA: usize> = IteratedEndoMultGroupMac<KAPPA>;
type InputMacKey<const KAPPA: usize> = [EndoMultGroupMacKey; KAPPA];

pub(crate) type Output = Fq;

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct Table<const KAPPA: usize> {
    pub(crate) one_mac: InputMac<KAPPA>,
    pub hash_s: Vec<u8>,
}

impl<const KAPPA: usize> Table<KAPPA> {
    pub(crate) fn groth16_input_mac_key(
        verifier_key: &groth16_verifier::VerifierKey,
        pairing_output_mac_key: &[EndoMultGroupMacKey; KAPPA],
        c_ec_mac_key: &[EndoMacKey<G1Config>; KAPPA],
    ) -> InputMacKey<KAPPA> {
        let groth16_verifier::VerifierKey { H1, c, .. } = *verifier_key;

        let mut key = Vec::with_capacity(pairing_output_mac_key.len());

        // e(A,B) e(C, H_1) = c^{-1}
        pairing_output_mac_key
            .iter()
            .zip(c_ec_mac_key.iter())
            .for_each(|(AB, C)| {
                // All these MAC keys use the same phi
                let phi = C.phi.clone();

                let k = Bn254::multi_pairing(vec![C.k], vec![H1]).0;
                // c^{p^(2phi + 6)} = c^{-p^(2phi)}
                // This has the effect of multiplying the MAC'd value by c
                let c_term = c.frobenius_map(2 * phi.as_u8() as usize + 6);
                let k = k * AB.k * c_term;

                key.push(EndoMultGroupMacKey { phi, k });
            });

        key.try_into().unwrap_or_else(|_| panic!())
    }

    pub(crate) fn derive_secret(input_mac_key: &InputMacKey<KAPPA>) -> Fq {
        pow_n_lincomb::<_, Fq, Fq, 6>(input_mac_key.iter().map(|key| Fq::from(key.phi.as_u8())))
    }
}

impl<const KAPPA: usize> Circuit<Input, InputMac<KAPPA>, Output, IdentityMac<Output>>
    for Table<KAPPA>
{
    type Aux = ();

    fn garble<P: PRG>(
        _parent_nonce: &mut Nonce,
        _prg_unseeded: &P::Unseeded,
        _rng: &mut StdRng,
        _aux: Option<&Self::Aux>,
        input_mac_key: &InputMacKey<KAPPA>,
    ) -> (Self, IdentityMacKey) {
        let one_mac = input_mac_key.encode(&Fp12::ONE);

        let hash_s = {
            let s = Self::derive_secret(input_mac_key);
            let bytes = s.into_bigint().to_bytes_le();

            let mut hasher = blake3::Hasher::new();
            hasher.update(&bytes);
            *hasher.finalize().as_bytes()
        };

        (
            Self {
                one_mac,
                hash_s: hash_s.to_vec(),
            },
            IdentityMacKey,
        )
    }

    fn evaluate<P: PRG>(
        &self,
        _parent_nonce: &mut Nonce,
        _prg_unseeded: &P::Unseeded,
        input: &Input,
        input_mac: &InputMac<KAPPA>,
    ) -> (Output, IdentityMac<Output>) {
        if *input_mac == self.one_mac {
            return (Fq::ZERO, IdentityMac(PhantomData));
        }

        // Precompute the powers of input
        let mut input_exp = [*input; 6];
        input_exp
            .iter_mut()
            .enumerate()
            .for_each(|(exp, input_exp)| *input_exp = input_exp.frobenius_map(exp * 2));

        // If input_mac is not equal to one_mac, we can recover the `phi`s of
        // the input_mac_key in the following way:
        // - one_mac   = 1 ^ {phi_i} * k_i
        // - input_mac = input ^ {phi_i} * k_i
        // - input_mac / one_mac = input^{phi_i}
        let mut input_div_one: Vec<_> = self.one_mac.iter().map(|v| **v).collect();
        batch_inversion(&mut input_div_one);
        input_div_one
            .iter_mut()
            .zip(input_mac.iter())
            .for_each(|(one_inv, input)| *one_inv *= **input);

        // For each input^{phi_i}, brute force the `phi`. This is efficient
        // because there are only six possible values of `phi`.
        let phis = input_div_one
            .iter()
            .map(|input| Fq::from(input_exp.iter().position(|&x| x == *input).unwrap() as u8));

        // Recover the secret as the pow_6_lincomb of the phis.
        let s = pow_n_lincomb::<_, Fq, Fq, 6>(phis);
        (s, IdentityMac(PhantomData))
    }
}

#[test]
fn test_conditional_reveal() {
    use crate::garble::pairing::KAPPA_B;
    use aes::Aes128;
    use ark_std::{UniformRand, rand::SeedableRng};

    const SEED: u64 = 0;

    let unseeded_prg = Aes128::init(&[]);

    let mut rng = StdRng::seed_from_u64(SEED);

    let input = Input::rand(&mut rng);
    let input_mac_key = InputMacKey::random_key(&mut rng);
    let input_mac = input_mac_key.encode(&input);

    let table = Table::<KAPPA_B>::garble::<Aes128>(
        &mut Nonce::new(),
        &unseeded_prg,
        &mut rng,
        None,
        &input_mac_key,
    )
    .0;
    let expected_hash_s = &table.hash_s;

    let output = table
        .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
        .0;

    let hash_output = {
        let bytes = output.into_bigint().to_bytes_le();
        let mut hasher = blake3::Hasher::new();
        hasher.update(&bytes);
        *hasher.finalize().as_bytes()
    };

    assert_eq!(hash_output.as_slice(), expected_hash_s);
}
