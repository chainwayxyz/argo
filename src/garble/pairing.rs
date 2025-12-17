//! Given the input MACs of (A, B): (G1, G2), returns the MAC of e(A, B).

use core::panic;
use std::marker::PhantomData;

use ark_bn254::{Bn254, Fq12, Fr, G1Projective, G2Affine, G2Projective, g1, g2};
use ark_ec::{
    AffineRepr, CurveGroup, PrimeGroup, ScalarMul,
    pairing::{Pairing, PairingOutput},
    scalar_mul::glv::GLVConfig,
    short_weierstrass::Affine,
};
use ark_ff::{AdditiveGroup, Field, PrimeField, UniformRand};
use ark_std::rand::rngs::StdRng;
use lazy_static::lazy_static;
use num_bigint::BigInt;

use crate::{
    garble::{Circuit, Nonce},
    mac::{
        Mac, MacKey, Random,
        endo::{
            ENDO_SIZE, EndoMacKey, EndoMultGroupMac, EndoMultGroupMacKey, EndoScalar,
            IteratedEndoMac, IteratedEndoMultGroupMac,
        },
    },
    utils::{bits::pow_n_lincomb, prg::PRG},
};

lazy_static! {
    /// e(G, H)
    pub static ref TARGET_GENERATOR: PairingOutput<Bn254> = Bn254::multi_pairing(
            vec![<G1Projective as PrimeGroup>::generator()],
            vec![G2Projective::generator()],
        );
}

// KAPPA_B is log_7(2^128)
pub(crate) const KAPPA_B: usize = 46;
// KAPPA_A is ceil(log_7(r))
pub(crate) const KAPPA_A: usize = 92;

pub(crate) type Input = (Affine<g1::Config>, Affine<g2::Config>);
pub(crate) type InputMacKey = (
    [EndoMacKey<g1::Config>; KAPPA_A],
    [EndoMacKey<g2::Config>; KAPPA_B],
);
pub(crate) type InputMac = (
    IteratedEndoMac<g1::Config, KAPPA_A>,
    IteratedEndoMac<g2::Config, KAPPA_B>,
);

pub(crate) type Output = Fq12;
pub(crate) type OutputMacKey = [EndoMultGroupMacKey; KAPPA_B];
pub(crate) type OutputMac = IteratedEndoMultGroupMac<KAPPA_B>;

impl Mac<Input> for InputMac {
    type Key = InputMacKey;

    #[cfg(test)]
    fn random_input(rng: &mut StdRng) -> Input {
        use ark_bn254::{G1Affine, G2Affine};

        (G1Affine::rand(rng), G2Affine::rand(rng))
    }
}

impl MacKey<Input, InputMac> for InputMacKey {
    fn encode(&self, value: &Input) -> InputMac {
        (self.0.encode(&value.0), self.1.encode(&value.1))
    }

    #[cfg(test)]
    fn decode(&self, mac: &InputMac) -> Input {
        (self.0.decode(&mac.0), self.1.decode(&mac.1))
    }

    fn random_key(_rng: &mut StdRng) -> Self {
        unimplemented!()
    }
}

pub(crate) struct Table(pub(crate) Box<[G2Affine; KAPPA_A]>);

impl Circuit<Input, InputMac, Output, OutputMac> for Table {
    type Aux = (Box<[Fr; KAPPA_A]>, OutputMacKey);

    fn garble<P: PRG>(
        _parent_nonce: &mut Nonce,
        _prg_unseeded: &P::Unseeded,
        _rng: &mut StdRng,
        aux: Option<&Self::Aux>,
        _input_mac_key: &InputMacKey,
    ) -> (Self, OutputMacKey) {
        let (X, output_mac_key) = aux.unwrap();

        let mut commitment = Box::new([G2Affine::zero(); KAPPA_A]);

        let X = G2Projective::generator().batch_mul(X.as_slice());
        commitment.copy_from_slice(&X);

        (Self(commitment), output_mac_key.clone())
    }

    fn evaluate<P: PRG>(
        &self,
        _parent_nonce: &mut Nonce,
        _prg_unseeded: &P::Unseeded,
        input: &Input,
        input_mac: &InputMac,
    ) -> (Output, OutputMac) {
        let (A, B) = input;
        let (A_mac, B_mac) = input_mac;
        let output = Bn254::pairing(A, B);

        let S = pow_n_lincomb(
            A_mac.iter().map(|mac| mac.into_group()),
            EndoScalar::<g2::Config>::two_minus_omega(),
        );

        let mut output_mac = Vec::with_capacity(KAPPA_B);
        (0..KAPPA_B).for_each(|i| {
            let B_mac = *B_mac[i];
            let X = self.0[i];

            let pairing = Bn254::multi_pairing(vec![*A, S.into_affine()], vec![B_mac, X]);
            output_mac.push(EndoMultGroupMac::from(pairing.0));
        });

        let output_mac: Box<[EndoMultGroupMac; KAPPA_B]> =
            output_mac.try_into().unwrap_or_else(|_| panic!());

        (output.0, output_mac.into())
    }
}

impl Table {
    fn derive_phis_from_x(x: Fr) -> [EndoScalar<g1::Config>; KAPPA_A] {
        let phi_components = |z: u8| match z {
            0 => (0, 0),
            1 => (1, 0),
            2 => (0, 1),
            3 => (1, 1),
            4 => (-1, -1),
            5 => (0, -1),
            6 => (-1, 0),
            _ => unreachable!(),
        };

        let phi_mapping = |z: u8| match z {
            0 => EndoScalar::NegInf(PhantomData),
            1 => EndoScalar::Zero(PhantomData),
            2 => EndoScalar::Four(PhantomData),
            3 => EndoScalar::Five(PhantomData),
            4 => EndoScalar::Two(PhantomData),
            5 => EndoScalar::One(PhantomData),
            6 => EndoScalar::Three(PhantomData),
            _ => unreachable!(),
        };

        let mut endoscalars = Vec::with_capacity(KAPPA_A);

        // `a`, `b` are signed bigints
        let (mut a, mut b) = {
            let (a, b) = g1::Config::scalar_decomposition(x);

            let mut a_signed = BigInt::from(a.1.into_bigint());
            if !a.0 {
                a_signed = -a_signed
            };

            let mut b_signed = BigInt::from(b.1.into_bigint());
            if !b.0 {
                b_signed = -b_signed
            };

            let a = a_signed - &b_signed;
            let b = -b_signed;

            (a, b)
        };

        for _ in 0..KAPPA_A {
            let z: BigInt = (&a + 2 * &b) % ENDO_SIZE;
            let z = (BigInt::from(ENDO_SIZE) + z) % ENDO_SIZE;

            let z = z.to_bytes_le().1[0];
            let (d, e) = phi_components(z);

            let phi = phi_mapping(z);
            endoscalars.push(phi.clone());

            let a_prime: BigInt = a - d;
            let b_prime: BigInt = b - e;

            let u: BigInt = 3 * &a_prime - &b_prime;
            let v: BigInt = 2 * b_prime + a_prime;

            a = u / ENDO_SIZE;
            b = v / ENDO_SIZE;
        }

        endoscalars.try_into().unwrap_or_else(|_| panic!())
    }

    pub(crate) fn random_mac_keypair_and_vec(
        rng: &mut StdRng,
    ) -> (InputMacKey, OutputMacKey, Box<[Fr; KAPPA_A]>) {
        // Sample `KAPPA_A` random `r`s
        let mut rs = [Fr::ZERO; KAPPA_A];
        rs.iter_mut().for_each(|r| *r = Fr::rand(rng));

        // Sample x uniformly at random
        let x = Fr::rand(rng);
        let A_phis = Self::derive_phis_from_x(x);

        let mut A_ks_scalar = [Fr::ZERO; KAPPA_A];
        A_ks_scalar
            .iter_mut()
            .for_each(|A_k_scalar| *A_k_scalar = Fr::rand(rng));

        let A_ks = G1Projective::generator().batch_mul(&A_ks_scalar);

        #[cfg(test)]
        {
            let A_phis_scalar = A_phis.iter().map(|A_phi| A_phi.as_scalar());
            let expected_x = {
                pow_n_lincomb::<_, Fr, Fr>(
                    A_phis_scalar,
                    EndoScalar::<g1::Config>::two_minus_omega(),
                )
            };
            assert_eq!(x, expected_x);
        }

        let x_rs: Vec<_> = rs.iter().map(|r| *r * x).collect();

        let neg_B_ks = G2Projective::generator().batch_mul(&x_rs);

        let B_phis: Vec<_> = (0..KAPPA_B)
            .map(|_| EndoScalar::<g2::Config>::random(rng))
            .collect();

        let k =
            pow_n_lincomb::<_, Fr, Fr>(A_ks_scalar, EndoScalar::<g1::Config>::two_minus_omega());

        let output_ks = {
            let mut powers_of_2 = Vec::with_capacity(Fr::MODULUS_BIT_SIZE as usize);
            let mut p = TARGET_GENERATOR.0;
            powers_of_2.push(p);
            for _ in 1..Fr::MODULUS_BIT_SIZE {
                p.square_in_place();
                powers_of_2.push(p);
            }

            let mut output_ks = [Fq12::ZERO; KAPPA_B];
            output_ks
                .iter_mut()
                .zip(rs.iter())
                .for_each(|(output_k, r)| {
                    let exp = (*r * k).into_bigint();
                    *output_k = Fq12::pow_with_table(&powers_of_2, exp).unwrap();
                });

            output_ks
        };

        let A_mac_key = A_phis
            .iter()
            .zip(A_ks.iter())
            .map(|(phi, &k)| EndoMacKey {
                phi: phi.clone(),
                k,
            })
            .collect::<Vec<_>>()
            .try_into()
            .unwrap_or_else(|_| panic!());

        let B_mac_key: [EndoMacKey<g2::Config>; KAPPA_B] = B_phis
            .iter()
            .zip(neg_B_ks.iter())
            .map(|(phi, &neg_k)| EndoMacKey {
                phi: phi.clone(),
                k: -neg_k,
            })
            .collect::<Vec<_>>()
            .try_into()
            .unwrap_or_else(|_| panic!());

        let output_mac_key = B_phis
            .iter()
            .zip(output_ks.iter())
            .map(|(phi, &k)| EndoMultGroupMacKey {
                phi: phi.clone().into(),
                k,
            })
            .collect::<Vec<_>>()
            .try_into()
            .unwrap_or_else(|_| panic!());

        (
            (A_mac_key, B_mac_key),
            output_mac_key,
            rs.try_into().unwrap_or_else(|_| panic!()),
        )
    }
}

#[test]
fn test_pairing_circuit() {
    use crate::garble::Nonce;
    use aes::Aes128;
    use ark_std::rand::SeedableRng;

    for seed in 0..100 {
        let mut rng = StdRng::seed_from_u64(seed);
        let unseeded_prg = Aes128::init(&[]);

        let input = InputMac::random_input(&mut rng);
        let (input_mac_key, output_mac_key, X) = Table::random_mac_keypair_and_vec(&mut rng);
        let input_mac = input_mac_key.encode(&input);

        let output = Bn254::pairing(input.0, input.1).0;

        let expected_output_mac = output_mac_key.encode(&output);

        let circuit = Table::garble::<Aes128>(
            &mut Nonce::new(),
            &unseeded_prg,
            &mut rng,
            Some(&(X, output_mac_key)),
            &input_mac_key,
        )
        .0;

        let output_mac = circuit
            .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
            .1;

        for (expected, mac) in expected_output_mac.iter().zip(output_mac.iter()) {
            assert!(expected == mac);
        }
    }
}
