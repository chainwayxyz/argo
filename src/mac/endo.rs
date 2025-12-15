//! Since the endoscalar assumes existence of sixth roots of unity, this works
//! only for curves with j-invariant = 0.
//!
//! Here, we fix the implementation to be for Bn254.

use std::marker::PhantomData;

use ark_bn254::{Fq12, g1, g2};
use ark_ec::{CurveGroup, scalar_mul::glv::GLVConfig, short_weierstrass::Affine};
use ark_ff::{BigInteger, Field, PrimeField, UniformRand};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::rngs::StdRng;
use derive_more::{Deref, From};

use crate::mac::{IteratedMac, Mac, MacKey, Random, RandomFromPhi};

#[derive(Clone, Copy)]
#[repr(u8)]
pub(crate) enum EndoScalar<G> {
    Zero(PhantomData<G>) = 0,
    One(PhantomData<G>),
    Two(PhantomData<G>),
    Three(PhantomData<G>),
    Four(PhantomData<G>),
    Five(PhantomData<G>),
}

impl<G> From<u8> for EndoScalar<G> {
    fn from(value: u8) -> Self {
        match value {
            0 => Self::Zero(PhantomData),
            1 => Self::One(PhantomData),
            2 => Self::Two(PhantomData),
            3 => Self::Three(PhantomData),
            4 => Self::Four(PhantomData),
            5 => Self::Five(PhantomData),
            _ => panic!(),
        }
    }
}

impl From<EndoScalar<g2::Config>> for EndoScalar<g1::Config> {
    fn from(value: EndoScalar<g2::Config>) -> Self {
        match value {
            EndoScalar::Zero(_) => EndoScalar::Zero(PhantomData),
            EndoScalar::One(_) => EndoScalar::One(PhantomData),
            EndoScalar::Two(_) => EndoScalar::Two(PhantomData),
            EndoScalar::Three(_) => EndoScalar::Three(PhantomData),
            EndoScalar::Four(_) => EndoScalar::Four(PhantomData),
            EndoScalar::Five(_) => EndoScalar::Five(PhantomData),
        }
    }
}

impl<G: GLVConfig> EndoScalar<G> {
    #[inline]
    fn p_sq_mod_r() -> G::ScalarField {
        let p = <<G::BaseField as Field>::BasePrimeField as PrimeField>::MODULUS;
        let p_mod_r = G::ScalarField::from_le_bytes_mod_order(&p.to_bytes_le());
        p_mod_r.square()
    }

    #[inline]
    fn scalar() -> G::ScalarField {
        // We need multiplication by this scalar to equal Frobenius squared
        // The same LAMBDA does not correspond to the same ENDO_COEFFS in G1
        // and G2
        // So, if you fix ENDO_COEFFS you get different LAMBDA (as BN254 does)
        // This reverses that and has the effect of fixing the endomorphism to
        // be multiplication by p^2.
        if -G::LAMBDA == Self::p_sq_mod_r() {
            -G::LAMBDA
        } else {
            -(G::LAMBDA.square())
        }
    }

    #[inline]
    fn base() -> G::BaseField {
        if -G::LAMBDA == Self::p_sq_mod_r() {
            -G::ENDO_COEFFS[0]
        } else {
            -(G::ENDO_COEFFS[0].square())
        }
    }

    pub(crate) fn as_u8(&self) -> u8 {
        match self {
            EndoScalar::Zero(_) => 0,
            EndoScalar::One(_) => 1,
            EndoScalar::Two(_) => 2,
            EndoScalar::Three(_) => 3,
            EndoScalar::Four(_) => 4,
            EndoScalar::Five(_) => 5,
        }
    }

    pub(crate) fn as_base(&self) -> G::BaseField {
        let exp = self.as_u8() as u64;
        Self::base().pow([exp])
    }

    pub(crate) fn as_scalar(&self) -> G::ScalarField {
        let exp = self.as_u8() as u64;
        Self::scalar().pow([exp])
    }

    /// phi^0 * P = P
    /// phi^1 * P = (-P).endo()
    /// phi^2 * P = P.endo().endo()
    /// phi^3 * P = -P
    /// phi^4 * P = P.endo()
    /// phi^5 * P = -P.endo().endo()
    pub(crate) fn endo(&self, point: Affine<G>) -> Affine<G> {
        let phi = self.as_base();
        let (a, b) = {
            let phi_sq = phi.square();
            (phi_sq.square(), phi_sq * phi)
        };

        Affine {
            x: a * point.x,
            y: b * point.y,
            infinity: point.infinity,
        }
    }

    #[cfg(test)]
    pub(crate) fn inverse(&self) -> Self {
        Self::from_u8(6 - self.as_u8())
    }

    #[cfg(test)]
    fn from_u8(variant: u8) -> Self {
        match variant {
            0 => EndoScalar::Zero(PhantomData),
            1 => EndoScalar::One(PhantomData),
            2 => EndoScalar::Two(PhantomData),
            3 => EndoScalar::Three(PhantomData),
            4 => EndoScalar::Four(PhantomData),
            5 => EndoScalar::Five(PhantomData),
            6 => EndoScalar::Zero(PhantomData),
            _ => unreachable!(),
        }
    }
}

#[derive(
    Clone, Copy, PartialEq, CanonicalSerialize, CanonicalDeserialize, Default, From, Deref,
)]
pub(crate) struct EndoMac<G: GLVConfig>(Affine<G>);
pub(crate) type IteratedEndoMac<G, const N: usize> = IteratedMac<Affine<G>, EndoMac<G>, N>;

impl<G> Mac<Affine<G>> for EndoMac<G>
where
    G: GLVConfig + Clone,
{
    type Key = EndoMacKey<G>;

    #[cfg(test)]
    fn random_input(rng: &mut StdRng) -> Affine<G> {
        use ark_std::rand::{distributions::Standard, prelude::Distribution as _};

        Standard.sample(rng)
    }
}

#[derive(Clone, Copy)]
pub(crate) struct EndoMacKey<G: GLVConfig> {
    pub(crate) phi: EndoScalar<G>,
    pub(crate) k: Affine<G>,
}

impl<G> RandomFromPhi<Affine<G>, EndoMac<G>> for EndoMacKey<G>
where
    G: GLVConfig + Clone,
{
    type Phi = EndoScalar<G>;

    fn random_from_phi(rng: &mut StdRng, phi: &EndoScalar<G>) -> Self {
        let k = Affine::<G>::rand(rng);

        Self {
            phi: phi.clone(),
            k,
        }
    }
}

impl<G> MacKey<Affine<G>, EndoMac<G>> for EndoMacKey<G>
where
    G: GLVConfig + Clone,
{
    fn encode(&self, value: &Affine<G>) -> EndoMac<G> {
        let phi_mul_value = self.phi.endo(*value);
        EndoMac((phi_mul_value + self.k).into_affine())
    }

    #[cfg(test)]
    fn decode(&self, mac: &EndoMac<G>) -> Affine<G> {
        let phi_mul_value = mac.0 - self.k;
        let phi_inv = self.phi.inverse();

        phi_inv.endo(phi_mul_value.into())
    }

    fn random_key(rng: &mut StdRng) -> Self {
        let phi = EndoScalar::random(rng);
        let k = Affine::<G>::rand(rng);

        Self { phi, k }
    }
}

/// G^phi * k, where phi is a bool and k is a curve point
#[derive(Clone, PartialEq, CanonicalDeserialize, CanonicalSerialize, Deref, From)]
pub(crate) struct EndoMultGroupMac(Fq12);

pub(crate) type IteratedEndoMultGroupMac<const N: usize> = IteratedMac<Fq12, EndoMultGroupMac, N>;

#[derive(Clone)]
pub(crate) struct EndoMultGroupMacKey {
    pub(crate) phi: EndoScalar<g1::Config>,
    pub(crate) k: Fq12,
}

impl RandomFromPhi<Fq12, EndoMultGroupMac> for EndoMultGroupMacKey {
    type Phi = EndoScalar<g1::Config>;

    fn random_from_phi(rng: &mut StdRng, phi: &Self::Phi) -> Self {
        Self {
            phi: phi.clone(),
            k: Fq12::rand(rng),
        }
    }
}

// Fr : u[0], u[1] are third roots
// Fp : v[0], v[1] are third
// G1
// G2
// Fp12
//
// Recall: p^12 = 1 mod r
// u[0] = p^4 OR p^8 mod r
// u[1] = p^8 OR p^4 mod r
//
// (x, y) in G1
// u[0] * (x, y) = (v[i] * x, y)
// u[1] * (x, y) = (v[1-i] * x, y)
//
// (x, y) in G2
// u[0] * (x, y) = (v[j] * x, y)
// u[1] * (x, y) = (v[1-j] * x, y)
//
// e(u[a] * P, Q) = e(P, Q)^u[a]
// e(P, u[a] * Q) = e(P, Q)^u[a]

// We store an integer phi from 0..5
// We have three functions: apply_g1, apply_g2, apply_gt
// We want to have the following properties:
//
// e(apply_g1(phi, P), Q) = apply_gt(phi, e(P, Q))
// e(P, apply_g2(phi, Q)) = apply_gt(phi, e(P, Q))
//
// apply_g1(phi, P) = (u[0]^a * P, (-1)^(phi mod 2) * P)
//.   a = phi mod 3 OR a = -phi mod 3
// apply_g2(phi, Q) = (u[0]^c * Q, (-1)^(phi mod 2) * Q)
//.   a = phi mod 3 OR a = -phi mod 3
// apply_gt(phi, g) = g^(p^phi)

impl MacKey<Fq12, EndoMultGroupMac> for EndoMultGroupMacKey {
    fn encode(&self, value: &Fq12) -> EndoMultGroupMac {
        (value.frobenius_map(self.phi.as_u8() as usize * 2) * self.k).into()
    }

    #[cfg(test)]
    fn decode(&self, _mac: &EndoMultGroupMac) -> Fq12 {
        // We cannot decode a EndoMultGroupMac
        unimplemented!()
    }

    fn random_key(rng: &mut StdRng) -> Self {
        Self {
            phi: EndoScalar::random(rng),
            k: Fq12::rand(rng),
        }
    }
}

impl Mac<Fq12> for EndoMultGroupMac {
    type Key = EndoMultGroupMacKey;

    #[cfg(test)]
    fn random_input(rng: &mut StdRng) -> Fq12 {
        Fq12::rand(rng)
    }
}

#[test]
fn test_endo() {
    use ark_bn254::g1::Config as G;
    use ark_ec::CurveConfig;
    use ark_ff::{BigInteger, PrimeField};

    {
        let endo = -(G::LAMBDA.square());

        let r = <G as CurveConfig>::ScalarField::MODULUS.to_bits_le();
        let r = <<G as CurveConfig>::BaseField as PrimeField>::BigInt::from_bits_le(&r);
        let r = <<G as CurveConfig>::BaseField as PrimeField>::from_bigint(r).unwrap();

        // This is p - r in Fp
        let neg_r = -r;
        let neg_r = neg_r.into_bigint().to_bits_le();
        let neg_r = <<G as CurveConfig>::ScalarField as PrimeField>::BigInt::from_bits_le(&neg_r);

        // This is p - r in Fr
        let neg_r = <<G as CurveConfig>::ScalarField as PrimeField>::from_bigint(neg_r);
        let p = neg_r.unwrap();
        let p_sq = p.square();

        assert_eq!(endo, p_sq);
    }
}

#[test]
fn test_gt_mac_commute() {
    use ark_bn254::g1::Config as G;
    use ark_bn254::{Bn254, G1Affine, G2Affine};
    use ark_ec::{AffineRepr, pairing::Pairing};
    use ark_std::rand::SeedableRng;

    const SEED: u64 = 0;

    let mut rng = StdRng::seed_from_u64(SEED);

    let point = G1Affine::rand(&mut rng);
    let gt_point = Bn254::pairing(point, G2Affine::generator()).0;

    let g1_mac_key = EndoMacKey::<G>::random_key(&mut rng);
    let g1_mac = g1_mac_key.encode(&point);

    let gt_mac_key = EndoMultGroupMacKey {
        phi: g1_mac_key.phi,
        k: Bn254::pairing(g1_mac_key.k, G2Affine::generator()).0,
    };
    let gt_mac = gt_mac_key.encode(&gt_point);

    assert_eq!(gt_mac.0, Bn254::pairing(g1_mac.0, G2Affine::generator()).0)
}
