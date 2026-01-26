//! Given the input MACs of (A, B): (G1, G2), returns the MAC of e(A, B).

use core::panic;
use std::{marker::PhantomData, sync::LazyLock};

use ark_bn254::{Bn254, Fq2, Fq12, Fr, G1Projective, G2Affine, G2Projective, g1, g2};
use ark_ec::{
    AffineRepr, CurveGroup, PrimeGroup, ScalarMul,
    pairing::{Pairing, PairingOutput},
    scalar_mul::glv::GLVConfig,
    short_weierstrass::{Affine, Projective},
};
use ark_ff::{AdditiveGroup, Field, PrimeField, UniformRand, Zero};
use ark_std::rand::rngs::StdRng;
use num_bigint::BigInt;

use crate::{
    garble::groth16_verifier::{KAPPA_A, KAPPA_B},
    mac::{
        Mac, MacKey,
        endo::{
            ENDO_SIZE, EndoMacKey, EndoMultGroupMac, EndoMultGroupMacKey, EndoScalar,
            IteratedEndoMac, IteratedEndoMultGroupMac,
        },
    },
    utils::bits::pow_n_lincomb,
};

/// e(G, H)
pub static TARGET_GENERATOR: LazyLock<PairingOutput<Bn254>> = LazyLock::new(|| {
    Bn254::multi_pairing(
        vec![<G1Projective as PrimeGroup>::generator()],
        vec![G2Projective::generator()],
    )
});

/// [r]P, where P is a point on E'/Fp2 with x-coordinate = 1
/// This unwrap cannot fail since this is a valid curve point.
pub static R_MUL_COFACTOR_GENERATOR: LazyLock<Affine<g2::Config>> = LazyLock::new(|| {
    (Affine::get_point_from_x_unchecked(Fq2::ONE, true)
        .unwrap()
        .mul_bigint(Fr::MODULUS))
    .into()
});

// This works only up to KAPPA_B = 65;
pub(crate) const GENERATOR_EXPONENT: usize = 3;

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

pub(crate) struct Table(pub(crate) Box<[G2Affine; KAPPA_B]>);

impl Table {
    pub(crate) fn garble(aux: (Box<[Fr; KAPPA_B]>, OutputMacKey)) -> (Self, OutputMacKey) {
        let (X, output_mac_key) = aux;

        let mut commitment = Box::new([G2Affine::zero(); KAPPA_B]);

        let X = G2Projective::generator().batch_mul(X.as_slice());
        commitment.copy_from_slice(&X);

        (Self(commitment), output_mac_key.clone())
    }

    pub(crate) fn evaluate(&self, input: &Input, input_mac: &InputMac) -> (Output, OutputMac) {
        let (A, B) = input;
        let (A_mac, B_mac) = input_mac;

        // Subtract off GENERATOR_EXPONENT^i * R_MUL_COFACTOR_GENERATOR from B_mac
        let B_mac = B_mac.iter().enumerate().map(|(i, mac)| {
            **mac
                - R_MUL_COFACTOR_GENERATOR.mul_bigint(
                    BigInt::from(GENERATOR_EXPONENT)
                        .pow(i as u32)
                        .to_u64_digits()
                        .1,
                )
        });

        let output = Bn254::pairing(A, B);

        let S = pow_n_lincomb(
            A_mac.iter().map(|mac| mac.into_group()),
            EndoScalar::<g2::Config>::two_minus_omega(),
        );

        let mut output_mac = Vec::with_capacity(KAPPA_B);
        B_mac.zip(self.0.iter()).for_each(|(B_mac, X)| {
            let pairing =
                Bn254::multi_pairing(vec![*A, S.into_affine()], vec![B_mac.into_affine(), *X]);
            output_mac.push(EndoMultGroupMac::from(pairing.0));
        });

        // This unwrap cannot fail since both `B_mac` and `self.0` are fixed
        // arrays of size `KAPPA_B`.
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

        // This unwrap cannot fail since `endoscalars` was constructed in a for
        // loop ranging from 0..KAPPA_A.
        endoscalars.try_into().unwrap_or_else(|_| panic!())
    }

    // r (phi B + K) = phi (r B) + (r K) = r_m
    // -r (phi B + K) = phi (-r B) - (r K) = -r_m
    //
    // We will publish r_K = r K in "cofactor group"
    // We will CDS decode r_m and -r_m
    //    i.e. we will take r_m - r_K =? phi(r B) for any phi in { NegInf, ... }
    //         we will take -r_m - r_K =? phi(r B) for any phi in { NegInf, ... }
    // At least one of them will always decode
    //.   i.e. we will always find a phi in the first check or one in the second check
    // Sometimes, we will find they *both* decode
    //.   i.e. WLOG r_m is the correct one so
    //.        r_m - r_K = phi (r B)
    //.        BUT by coincidence on B
    //.        -r_m - r_K = -(phi (r B)) - 2 r_K = psi (r B)
    // This is because the cofactor group has small subgroups (technically could also occur in the Gt decoding)
    // To solve this, we need to ensure the r_K are chosen so that this can occur at most one time
    //    i.e. for any single r_K we can always cause a decoding error by maliciously choosing B
    //.        pick a psi != phi, and let r B = -2 / (phi + psi) r_K
    //.        but we can ensure that *at most* most one failure can occur
    // for all (phi, psi, phi', psi') where phi != -psi and phi' != -psi'
    // must have (phi + psi) r_(K,i) != (phi' + psi') r_(K,j)
    // Note: r_K does not have to be random
    // We should also choose the phi such that given n-1 of them we can always recover the last one
    // i.e. the sum mod 7 is always zero
    //
    // Let G' = r * G
    // Let K_i = (x * r_i * H) + GENERATOR_EXPONENT^i * G'
    fn derive_b_ks(g2_b_ks: &[Affine<g2::Config>; KAPPA_B]) -> [Affine<g2::Config>; KAPPA_B] {
        let mut b_ks = [Projective::<g2::Config>::zero(); KAPPA_B];
        b_ks.iter_mut()
            .zip(g2_b_ks.iter())
            .enumerate()
            .for_each(|(i, (b_k, g2_b_k))| {
                *b_k = *g2_b_k
                    + R_MUL_COFACTOR_GENERATOR.mul_bigint(
                        BigInt::from(GENERATOR_EXPONENT as u8)
                            .pow(i as u32)
                            .to_u64_digits()
                            .1,
                    )
            });

        let b_ks = Projective::<g2::Config>::normalize_batch(&b_ks);

        // This unwrap cannot fail since g2_b_ks is a fixed array of length
        // `KAPPA_B`.
        b_ks.try_into().unwrap()
    }

    pub(crate) fn random_mac_keypair_and_vec(
        rng: &mut StdRng,
    ) -> (InputMacKey, OutputMacKey, Box<[Fr; KAPPA_B]>) {
        // Sample `KAPPA_B` random `r`s
        let mut rs = Box::new([Fr::ZERO; KAPPA_B]);
        for r in rs.iter_mut() {
            *r = Fr::rand(rng);
        }

        // Sample x uniformly at random
        let x = Fr::rand(rng);
        let A_phis = Self::derive_phis_from_x(x);

        // Sample A_ks_scalar directly into array
        let mut A_ks_scalar = [Fr::ZERO; KAPPA_A];
        for A_k_scalar in A_ks_scalar.iter_mut() {
            *A_k_scalar = Fr::rand(rng);
        }

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

        // Compute g2_B_ks_scalar in-place without intermediate Vec
        let mut g2_B_ks_scalar = [Fr::ZERO; KAPPA_B];
        for (i, r) in rs.iter().enumerate() {
            g2_B_ks_scalar[i] = -*r * x;
        }

        let g2_B_ks = G2Projective::generator().batch_mul(&g2_B_ks_scalar);
        // This unwrap cannot fail since `g2_B_ks_scalar` is a fixed array of
        // length `KAPPA_B`.
        let B_ks = Self::derive_b_ks(&g2_B_ks.try_into().unwrap());

        let B_phis = EndoScalar::<g2::Config>::sample_endoscalars_with_checksum::<KAPPA_B>(rng);

        let k = pow_n_lincomb::<_, Fr, Fr>(
            A_ks_scalar.iter().copied(),
            EndoScalar::<g1::Config>::two_minus_omega(),
        );

        // Precompute powers of 2 once
        let powers_of_2 = {
            let mut powers = Vec::with_capacity(Fr::MODULUS_BIT_SIZE as usize);
            let mut p = TARGET_GENERATOR.0;
            powers.push(p);
            for _ in 1..Fr::MODULUS_BIT_SIZE {
                p.square_in_place();
                powers.push(p);
            }
            powers
        };

        // Compute output_ks
        let mut output_ks = [Fq12::ZERO; KAPPA_B];
        for (output_k, r) in output_ks.iter_mut().zip(rs.iter()) {
            let exp = (*r * k).into_bigint();
            // This unwrap cannot fail because we compute enough `powers_of_2`
            // for `Fr`.
            *output_k = Fq12::pow_with_table(&powers_of_2, exp).unwrap();
        }

        // Build A_mac_key efficiently with known size
        let mut A_mac_key = Vec::with_capacity(KAPPA_A);
        for ((phi, &k_scalar), &k) in A_phis.iter().zip(A_ks_scalar.iter()).zip(A_ks.iter()) {
            A_mac_key.push(EndoMacKey {
                phi: phi.clone(),
                k_scalar: Some(k_scalar),
                k,
            });
        }
        // This unwrap cannot fail because `A_phis`, `A_ks_scalar`, `A_ks` are
        // all fixed-size arrays of length `KAPPA_A`.
        let A_mac_key = A_mac_key.try_into().unwrap_or_else(|_| panic!());

        // Build B_mac_key efficiently with known size
        let mut B_mac_key = Vec::with_capacity(KAPPA_B);
        for ((phi, &k_scalar), &k) in B_phis.iter().zip(g2_B_ks_scalar.iter()).zip(B_ks.iter()) {
            B_mac_key.push(EndoMacKey {
                phi: phi.clone(),
                k_scalar: Some(k_scalar),
                k,
            });
        }
        // This unwrap cannot fail because `B_phis`, `g2_B_ks_scalar`, `B_ks`
        // are all of length `KAPPA_B`.
        let B_mac_key = B_mac_key.try_into().unwrap_or_else(|_| panic!());

        // Build output_mac_key efficiently with known size
        let mut output_mac_key = Vec::with_capacity(KAPPA_B);
        for (phi, k_val) in B_phis.iter().zip(output_ks.iter()) {
            output_mac_key.push(EndoMultGroupMacKey {
                phi: phi.clone().into(),
                k: *k_val,
            });
        }
        // This unwrap cannot fail because `B_phis`, `output_ks` are of length
        // `KAPPA_B`.
        let output_mac_key = output_mac_key.try_into().unwrap_or_else(|_| panic!());

        ((A_mac_key, B_mac_key), output_mac_key, rs)
    }
}

#[allow(unused)]
fn test_generator() {
    use ark_bn254::Fq2;
    use num_bigint::Sign;

    // # E'/Fp2 = n = h * r
    // Prime factors of h = 21888242871839275222246405745257275088844257914179612981679871602714643921549
    let prime_factors = [
        BigInt::from(10069),
        BigInt::from(5864401),
        // 1875725156269
        BigInt::new(Sign::Plus, vec![3119415213, 436]),
        // 197620364512881247228717050342013327560683201906968909
        BigInt::new(
            Sign::Plus,
            vec![
                2383251789, 2607680309, 252471992, 3923425572, 1462367600, 135217,
            ],
        ),
    ];

    let r = BigInt::from(Fr::MODULUS);

    // Check: (n / p_k) * G = G'_k != 0
    let g_ks = (0..4).map(|skip_idx| {
        let scalar = prime_factors
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != skip_idx)
            .map(|(_, factor)| factor)
            .product::<BigInt>();
        let g_k = R_MUL_COFACTOR_GENERATOR.mul_bigint(scalar.to_u64_digits().1);
        assert_ne!(g_k, Affine::<g2::Config>::identity());
        g_k
    });

    // Check: for all (i != j in KAPPA_B) that (phi + psi) i G'_k != (phi' + psi') j G'_k
    for g_k in g_ks {
        for i in 0..KAPPA_B {
            let i_g_k = g_k.mul_bigint(
                BigInt::from(GENERATOR_EXPONENT)
                    .pow(i as u32)
                    .to_u64_digits()
                    .1,
            );
            for j in (i + 1)..KAPPA_B {
                let j_g_k = g_k.mul_bigint(
                    BigInt::from(GENERATOR_EXPONENT)
                        .pow(j as u32)
                        .to_u64_digits()
                        .1,
                );
                for phi in EndoScalar::<g2::Config>::iter() {
                    for psi in EndoScalar::<g2::Config>::iter() {
                        for phi_prime in EndoScalar::<g2::Config>::iter() {
                            for psi_prime in EndoScalar::<g2::Config>::iter() {
                                let lhs = psi.endo(i_g_k.into()) + phi.endo(i_g_k.into());
                                let rhs =
                                    psi_prime.endo(j_g_k.into()) + phi_prime.endo(j_g_k.into());

                                let lhs: Affine<g2::Config> = lhs.into();
                                let rhs: Affine<g2::Config> = rhs.into();

                                if !lhs.is_zero() || !rhs.is_zero() {
                                    assert_ne!(lhs, rhs);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn test_pairing_circuit() {
    use ark_std::rand::SeedableRng;

    for seed in 0..100 {
        let mut rng = StdRng::seed_from_u64(seed);

        let input = InputMac::random_input(&mut rng);
        let (input_mac_key, output_mac_key, X) = Table::random_mac_keypair_and_vec(&mut rng);
        let input_mac = input_mac_key.encode(&input);

        let output = Bn254::pairing(input.0, input.1).0;

        let expected_output_mac = output_mac_key.encode(&output);

        let circuit = Table::garble((X, output_mac_key)).0;

        let output_mac = circuit.evaluate(&input, &input_mac).1;

        for (expected, mac) in expected_output_mac.iter().zip(output_mac.iter()) {
            assert!(expected == mac);
        }
    }
}
