//! Given the input MACs of (A, B): (G1, G2), returns the MAC of e(A, B).

use core::panic;

use ark_bn254::{Bn254, Fq12, Fr, G1Projective, G2Affine, G2Projective, g1, g2};
use ark_ec::{
    AffineRepr, CurveGroup, PrimeGroup, ScalarMul,
    pairing::{Pairing, PairingOutput},
    short_weierstrass::Projective,
};
use ark_ff::{Field, PrimeField, UniformRand};
use ark_std::rand::rngs::StdRng;
use lazy_static::lazy_static;

use crate::{
    garble::{Circuit, Nonce},
    mac::{
        Mac, MacKey, Random,
        endo::{
            EndoMacKey, EndoMultGroupMac, EndoMultGroupMacKey, EndoScalar, IteratedEndoMac,
            IteratedEndoMultGroupMac,
        },
    },
    utils::{prg::PRG, try_into_array},
};

lazy_static! {
    /// e(G, H)
    pub static ref TARGET_GENERATOR: PairingOutput<Bn254> = Bn254::multi_pairing(
            vec![<G1Projective as PrimeGroup>::generator()],
            vec![G2Projective::generator()],
        );
}

pub(crate) const KAPPA_B: usize = 50;
pub(crate) const KAPPA_A: usize = KAPPA_B * 3;

pub(crate) type Input = (Projective<g1::Config>, Projective<g2::Config>);
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
        use ark_bn254::G2Projective;

        (G1Projective::rand(rng), G2Projective::rand(rng))
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

/// Commitment to {KAPPA_A + KAPPA_B} powers of tau
pub(crate) struct Table(pub(crate) Box<[G2Affine; KAPPA_A + KAPPA_B]>);

impl Circuit<Input, InputMac, Output, OutputMac> for Table {
    type Aux = (Box<[Fr; KAPPA_A + KAPPA_B]>, OutputMacKey);

    fn garble<P: PRG>(
        _parent_nonce: &mut Nonce,
        _prg_unseeded: &P::Unseeded,
        _rng: &mut StdRng,
        aux: Option<&Self::Aux>,
        _input_mac_key: &InputMacKey,
    ) -> (Self, OutputMacKey) {
        let (powers_of_tau, output_mac_key) = aux.unwrap();

        let mut commitment = Box::new([G2Affine::zero(); KAPPA_A + KAPPA_B]);

        // Commit to powers of tau in G2
        let powers_of_tau = G2Projective::generator().batch_mul(powers_of_tau.as_slice());
        commitment.copy_from_slice(&powers_of_tau);

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

        let mut output_mac = Vec::with_capacity(KAPPA_B);
        (0..KAPPA_B).for_each(|i| {
            let B_mac = B_mac[i].into_affine();
            let X_row = (0..KAPPA_A).map(|j| self.0[i + j]);

            let pairing = Bn254::multi_pairing(
                std::iter::once(A.into_affine()).chain(A_mac.iter().map(|mac| mac.into_affine())),
                std::iter::once(B_mac).chain(X_row),
            );
            output_mac.push(EndoMultGroupMac::from(pairing.0));
        });

        let output_mac: Box<[EndoMultGroupMac; KAPPA_B]> =
            output_mac.try_into().unwrap_or_else(|_| panic!());

        (output.0, output_mac.into())
    }
}

impl Table {
    /// Given:
    ///   - an array of `KAPPA_A` endoscalars phi's,
    /// Generates:
    ///   - powers of a random scalar `tau`, and
    ///   - an array of `KAPPA_B` G2 elements K,
    /// s.t. K[i] = [-(x tau^i)]G, where x = ∑ phi_j tau^j.
    fn generate_random_tau_and_vec(
        rng: &mut StdRng,
        phis: &[EndoScalar<g1::Config>],
    ) -> (Box<[Fr; KAPPA_A + KAPPA_B]>, Box<[G2Projective; KAPPA_B]>) {
        assert_eq!(phis.len(), KAPPA_A);

        // Sample a random tau
        let tau = Fr::rand(rng);

        // Compute the powers of tau, [tau^0, ..., tau^{KAPPA_A + KAPPA_B - 1}]
        let powers_of_tau: Box<[Fr; KAPPA_A + KAPPA_B]> = (0..(KAPPA_A + KAPPA_B))
            .scan(Fr::ONE, |state, _| {
                let current = *state;
                *state *= tau;
                Some(current)
            })
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();

        // x = ∑ phi_j tau^j
        let x = phis
            .iter()
            .zip(powers_of_tau.iter())
            .map(|(phi, power)| phi.as_scalar() * power)
            .sum::<Fr>();

        // k[i] = -(x tau^i)
        let k_scalars = powers_of_tau
            .iter()
            .take(KAPPA_B)
            .map(|power| -x * power)
            .collect::<Vec<_>>();

        // K[i] = k[i]G
        let k_points = G2Projective::generator()
            .batch_mul(&k_scalars)
            .iter()
            .map(|point| point.into_group())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap_or_else(|_| panic!());

        (powers_of_tau, k_points)
    }

    pub(crate) fn random_mac_keypair_and_tau(
        rng: &mut StdRng,
    ) -> (InputMacKey, OutputMacKey, Box<[Fr; KAPPA_A + KAPPA_B]>) {
        // Sample `KAPPA_A` random endoscalars for `A_phis`
        let mut A_phis = Vec::with_capacity(KAPPA_A);
        for _ in 0..KAPPA_A {
            A_phis.push(EndoScalar::<g1::Config>::random(rng));
        }

        // Sample `KAPPA_A` random scalars for `A_ks`
        let A_ks_scalars = (0..KAPPA_A).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let A_ks = G1Projective::generator().batch_mul(&A_ks_scalars);

        // Sample `KAPPA_B` random endocscalars for `B_phis`
        let mut B_phis = Vec::with_capacity(KAPPA_B);
        for _ in 0..KAPPA_B {
            B_phis.push(EndoScalar::<g2::Config>::random(rng));
        }

        // Use `A_phis` to generate powers_of_tau, and B_ks to be used in
        // constructing B_mac_key.
        let (powers_of_tau, B_ks) = Self::generate_random_tau_and_vec(rng, &A_phis);

        // Compute output_ks as ∏ e(A_ks[j], [tau^{i + j}]G)
        let output_ks = (0..KAPPA_B).map(|i| {
            TARGET_GENERATOR.0.pow(
                (0..KAPPA_A)
                    .map(|j| A_ks_scalars[j] * powers_of_tau[i + j])
                    .sum::<Fr>()
                    .into_bigint(),
            )
        });

        // Construct A_mac_key
        let A_mac_key = try_into_array(A_phis.iter().zip(A_ks.iter()).map(|(phi, k)| EndoMacKey {
            phi: phi.clone(),
            k: k.into_group(),
        }));

        // Construct B_mac_key
        let B_mac_key =
            try_into_array(B_phis.iter().zip(B_ks.iter()).map(|(phi, &k)| EndoMacKey {
                phi: phi.clone(),
                k,
            }));

        // Construct output_mac_key
        let output_mac_key: [_; KAPPA_B] =
            try_into_array(B_phis.iter().zip(output_ks).map(|(phi, k)| {
                let phi: EndoScalar<g1::Config> = phi.clone().into();
                EndoMultGroupMacKey { phi, k }
            }));

        ((A_mac_key, B_mac_key), output_mac_key, powers_of_tau)
    }
}

#[test]
fn test_pairing_circuit() {
    use crate::garble::Nonce;
    use aes::Aes128;
    use ark_std::rand::SeedableRng;

    const SEED: u64 = 0;

    let mut rng = StdRng::seed_from_u64(SEED);
    let unseeded_prg = Aes128::init(&[]);

    let (input_mac_key, output_mac_key, X) = Table::random_mac_keypair_and_tau(&mut rng);

    let input = InputMac::random_input(&mut rng);
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
