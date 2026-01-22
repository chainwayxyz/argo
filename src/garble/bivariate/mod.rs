//! Garbling tables for oblivious evaluation of bivariate polynomials

use ark_ec::{
    AffineRepr,
    short_weierstrass::{Affine, SWCurveConfig},
};
use ark_ff::Field;
use ark_std::rand::rngs::StdRng;

use crate::{
    mac::{Mac, MacKey, Random, digit::bytes},
    utils::bits::{FromBits, NumBits, ToBits},
};

/// Bivariate quadratic poly in coefficient form
/// Ordered according to weight:
/// [constant, x, y, x_sq, xy, y_sq]
pub(crate) struct Poly<F>(pub(crate) [F; 6]);

impl<F: Field> Poly<F> {
    pub(crate) fn new(coeffs: [F; 6]) -> Self {
        Self(coeffs)
    }

    pub(crate) fn leading_coeff(&self) -> F {
        self.0
            .iter()
            .rev()
            .find(|x| **x != F::ZERO)
            .copied()
            .unwrap_or(F::ZERO)
    }

    pub(crate) fn scale(&mut self, scalar: F) {
        self.0.iter_mut().for_each(|coeff| *coeff *= scalar);
    }
}

/// Garbling table for arbitrary bivariate quadratic polynomials
pub(crate) mod biquadratic;
/// Garbling table that takes as input the values (x, y), and leaks `secret`
/// if x^3 + B != y^2, i.e. if (x, y) are not affine coordinates of a valid
/// curve point.
pub(crate) mod curve_membership;

/// Bivariate input (x, y)
#[derive(Debug, Clone, Copy)]
pub struct Input<F> {
    pub(crate) x: F,
    pub(crate) y: F,
}

impl<C> From<C> for Input<C::BaseField>
where
    C: AffineRepr,
    C::BaseField: Field,
{
    /// # Panics
    ///
    /// Panics if `value` is the point at infinity.
    /// FIXME is this okay?
    fn from(value: C) -> Self {
        let (x, y) = value.xy().unwrap();
        Self { x, y }
    }
}

impl<C> From<Input<C::BaseField>> for Affine<C>
where
    C: SWCurveConfig,
    C::BaseField: Field,
{
    fn from(val: Input<C::BaseField>) -> Self {
        Affine::new(val.x, val.y)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct InputMac<F> {
    pub(crate) x: bytes::DigitMac<F>,
    pub(crate) y: bytes::DigitMac<F>,
}

#[derive(Clone)]
pub struct InputMacKey<F> {
    pub(crate) x: bytes::DigitMacKey<F>,
    pub(crate) y: bytes::DigitMacKey<F>,
}

impl<F> MacKey<Input<F>, InputMac<F>> for InputMacKey<F>
where
    F: Field + NumBits + ToBits + FromBits<F, F> + Random,
{
    fn encode(&self, value: &Input<F>) -> InputMac<F> {
        InputMac {
            x: self.x.encode(&value.x),
            y: self.y.encode(&value.y),
        }
    }

    #[cfg(test)]
    fn decode(&self, mac: &InputMac<F>) -> Input<F> {
        Input {
            x: self.x.decode(&mac.x),
            y: self.y.decode(&mac.y),
        }
    }

    fn random_key(rng: &mut StdRng) -> Self {
        Self {
            x: bytes::DigitMacKey::random_key(rng),
            y: bytes::DigitMacKey::random_key(rng),
        }
    }
}

impl<F> Mac<Input<F>> for InputMac<F>
where
    F: Field + NumBits + ToBits + FromBits<F, F> + Random,
{
    type Key = InputMacKey<F>;

    #[cfg(test)]
    fn random_input(rng: &mut StdRng) -> Input<F> {
        Input {
            x: F::rand(rng),
            y: F::rand(rng),
        }
    }
}
