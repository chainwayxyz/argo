use ark_bn254::{Fq, Fq2};
use ark_ff::{AdditiveGroup, BigInteger, Field, PrimeField};
use lazy_static::lazy_static;
use std::ops::{Add, Mul, MulAssign};

lazy_static! {
    static ref ALPHA: Fq2 = Fq2::from_base_prime_field_elems([Fq::ZERO, Fq::ONE]).unwrap();
}

pub(crate) trait FromBits<B, V> {
    /// Takes a bit representation where each bit is encoded as a value in `B`,
    /// and composes it to reconstruct a value in `V`.
    fn from_bits<I: DoubleEndedIterator<Item = B> + ExactSizeIterator>(bits: I) -> V;
}

pub trait NumBits {
    /// The number of bits in the bit representation of this value.
    fn num_bits() -> u32;
}

pub trait ToBits {
    /// Decomposes this value into its bit representation.
    fn to_bits(&self) -> Vec<bool>;
}

impl FromBits<bool, Fq> for Fq {
    fn from_bits<I: DoubleEndedIterator<Item = bool>>(bits: I) -> Fq {
        Fq::from_bigint(BigInteger::from_bits_le(&bits.collect::<Vec<_>>())).unwrap()
    }
}

impl FromBits<bool, Fq2> for Fq2 {
    fn from_bits<I: DoubleEndedIterator<Item = bool> + ExactSizeIterator>(mut bits: I) -> Fq2 {
        let fp_len = Fq::MODULUS_BIT_SIZE as usize;
        assert_eq!(bits.len(), fp_len * 2);

        let c0 = Fq::from_bits(bits.by_ref().take(fp_len));
        let c1 = Fq::from_bits(bits.skip(fp_len));

        Fq2::from_base_prime_field_elems([c0, c1]).unwrap()
    }
}

impl FromBits<Fq, Fq> for Fq {
    fn from_bits<I: DoubleEndedIterator<Item = Fq>>(bits: I) -> Fq {
        pow_two_lincomb::<_, Fq, Fq>(bits)
    }
}

impl FromBits<Fq, Fq2> for Fq2 {
    fn from_bits<I: DoubleEndedIterator<Item = Fq> + ExactSizeIterator>(mut bits: I) -> Fq2 {
        let fp_len = Fq::MODULUS_BIT_SIZE as usize;
        assert_eq!(bits.len(), fp_len * 2);

        let c0 = Fq::from_bits(bits.by_ref().take(fp_len));
        let c1 = Fq::from_bits(bits.skip(fp_len));

        Fq2::from_base_prime_field_elems([c0, c1]).unwrap()
    }
}

impl FromBits<Fq2, Fq2> for Fq2 {
    fn from_bits<I: DoubleEndedIterator<Item = Fq2> + ExactSizeIterator>(mut bits: I) -> Fq2 {
        let fp_len = Fq::MODULUS_BIT_SIZE as usize;
        assert_eq!(bits.len(), fp_len * 2);

        let c0_bits: Vec<_> = bits.by_ref().take(fp_len).collect();
        let c0 = pow_two_lincomb::<_, Fq2, Fq2>(c0_bits);

        let c1 = pow_two_lincomb::<_, Fq2, Fq2>(bits.map(|bit| bit * *ALPHA));

        c0 + c1
    }
}

impl NumBits for Fq {
    #[inline]
    fn num_bits() -> u32 {
        Self::MODULUS_BIT_SIZE
    }
}

impl NumBits for Fq2 {
    #[inline]
    fn num_bits() -> u32 {
        Fq::MODULUS_BIT_SIZE * 2
    }
}

impl ToBits for Fq {
    #[inline]
    fn to_bits(&self) -> Vec<bool> {
        let mut bits = Vec::with_capacity(Self::MODULUS_BIT_SIZE as usize);
        bits.extend_from_slice(&self.into_bigint().to_bits_le()[..Self::MODULUS_BIT_SIZE as usize]);
        bits
    }
}

impl ToBits for Fq2 {
    #[inline]
    fn to_bits(&self) -> Vec<bool> {
        let base_elems = self.to_base_prime_field_elements();
        let mut bits = Vec::with_capacity(
            base_elems.size_hint().0 * <Fq2 as Field>::BasePrimeField::MODULUS_BIT_SIZE as usize,
        );

        for v in base_elems {
            bits.extend_from_slice(&v.to_bits());
        }

        bits
    }
}

/// Computes ∑ v_i * 2^i using Horner's method
pub(crate) fn pow_two_lincomb<I, Value, Coeff>(values: I) -> Value
where
    I: IntoIterator<Item = Value>,
    I::IntoIter: DoubleEndedIterator,
    Value: Mul<Coeff, Output = Value>
        + MulAssign<Coeff>
        + Add<Value, Output = Value>
        + Copy
        + AdditiveGroup,
    Coeff: Copy,
{
    values
        .into_iter()
        .rev()
        .fold(Value::zero(), |acc, value| acc.double() + value)
}

/// Computes ∑ v_i * N^i using Horner's method
pub(crate) fn pow_n_lincomb<I, Value, Coeff>(values: I, n: Coeff) -> Value
where
    I: IntoIterator<Item = Value>,
    I::IntoIter: DoubleEndedIterator,
    Value: AdditiveGroup
        + Mul<Coeff, Output = Value>
        + MulAssign<Coeff>
        + Add<Value, Output = Value>
        + Copy,
    Coeff: Clone,
{
    values
        .into_iter()
        .rev()
        .fold(Value::zero(), |acc, value| acc * n.clone() + value)
}
