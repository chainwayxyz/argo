//! The `field_to_ec_mac` on a G2 point produces `KAPPA_B` iterated MACs of the
//! form `m = phi * P + k`, where `k` is always in G2.
//!
//! This gadget leaks a secret if `P` is not in G2.
//!
//! If x is not in subgroup then r * m = y != 0, where `r` is group order
//! Solve dlog of y wrt z = r * P will return phi as long as gcd(h, 6) = 1,
//! where h is the cofactor

use ark_bn254::{Fq, Fr, g2};
use ark_ec::AffineRepr;
use ark_ec::{CurveGroup, short_weierstrass::Affine};
use ark_ff::PrimeField;
#[cfg(test)]
use ark_std::rand::rngs::StdRng;

use crate::{mac::endo::EndoScalar, utils::bits::pow_n_lincomb};

/// Leaks a secret if `curve_point` is not in G2.
pub(crate) fn check_g2_membership(
    curve_point: Affine<g2::Config>,
    g2_mac: &[Affine<g2::Config>],
) -> Option<Fq> {
    // For each possible value of `phi`, compute [r * phi] P
    let r_phi_P: Vec<_> = (0..6)
        .map(|phi| {
            let phi: EndoScalar<g2::Config> = phi.into();

            phi.endo(curve_point.mul_bigint(Fr::MODULUS).into_affine())
        })
        .collect();

    if r_phi_P.iter().all(|point| point.is_zero()) {
        return None;
    }

    let phis = g2_mac.iter().map(|point| {
        // m     =     phi * P +     k
        // r * m = r * phi * P + r * k = r * phi * P
        let y = point.mul_bigint(Fr::MODULUS);
        let endoscalar = r_phi_P.iter().position(|&x| x == y).unwrap();
        Fq::from(endoscalar as u8)
    });

    // Recover the secret as the pow_6_lincomb of the phis.
    let s = pow_n_lincomb::<_, Fq, Fq, 6>(phis);

    Some(s)
}

#[cfg(test)]
pub(crate) fn get_curve_point_not_in_g2(rng: &mut StdRng) -> Affine<g2::Config> {
    use ark_bn254::Fq2;
    use ark_ff::UniformRand;
    use ark_std::Zero;

    loop {
        let x = Fq2::rand(rng);
        if let Some(point) = Affine::get_point_from_x_unchecked(x, false)
            && !point.mul_bigint(Fr::MODULUS).is_zero()
        {
            return point;
        }
    }
}
