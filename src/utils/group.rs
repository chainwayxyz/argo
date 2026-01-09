use ark_bn254::{Fq, Fr};
use ark_ff::{BigInt, BigInteger, PrimeField};
use lazy_static::lazy_static;

use crate::utils::bits::ToBits;

lazy_static! {
    /// A 17-byte offset in Bn254::Fq, i.e. Fq(2^{8 * 17})
    /// Used in optimised from_le_bytes_mod_order for hash_to_field
    pub static ref FQ_LO_BYTES_OFFSET: Fq = Fq::from_bigint(BigInt::one() << (8 * 17)).unwrap();

    /// A 17-byte offset in Bn254::Fr, i.e. Fr(2^{8 * 17})
    /// Used in optimised from_le_bytes_mod_order for hash_to_field
    pub static ref FR_LO_BYTES_OFFSET: Fr = Fr::from_bigint(BigInt::one() << (8 * 17)).unwrap();
}

fn to_u64_array<const N: usize>(bytes: &[u8]) -> [u64; N] {
    let mut result = [0u64; N];

    for (i, chunk) in bytes.chunks(8).take(N).enumerate() {
        let mut u64_bytes = [0u8; 8];
        u64_bytes[..chunk.len()].copy_from_slice(chunk);
        result[i] = u64::from_le_bytes(u64_bytes);
    }

    result
}

pub trait FromLeBytesModOrder {
    fn from_le_bytes_mod_order(bytes: &[u8]) -> Self;
}

impl FromLeBytesModOrder for Fq {
    /// Optimized implementation of from_le_bytes_mod_order for Bn254::Fq
    fn from_le_bytes_mod_order(bytes: &[u8]) -> Fq {
        let num_modulus_bytes = Fq::MODULUS_BIT_SIZE.div_ceil(8) as usize;
        let num_bytes_to_directly_convert = ark_std::cmp::min(num_modulus_bytes - 1, bytes.len());

        // Copy the leading little-endian bytes directly into a field element.
        // The number of bytes directly converted must be less than the
        // number of bytes needed to represent the modulus, as we must begin
        // modular reduction once the data is of the same number of bytes as the
        // modulus.

        let (lo, hi) = bytes.split_at(bytes.len() - num_bytes_to_directly_convert);

        let lo: Fq = Fq::new_unchecked(BigInt::new(to_u64_array::<4>(lo)));
        let mut hi = Fq::new_unchecked(BigInt::new(to_u64_array::<4>(hi)));

        hi *= *FQ_LO_BYTES_OFFSET;
        hi += lo;

        // Sanity-check equivalence with arkworks unoptimised implementation
        #[cfg(test)]
        assert_eq!(
            hi * Fq::new_unchecked(Fq::R2),
            <Fq as PrimeField>::from_le_bytes_mod_order(bytes)
        );

        hi
    }
}

impl FromLeBytesModOrder for Fr {
    /// Optimized implementation of from_le_bytes_mod_order for Bn254::Fr
    fn from_le_bytes_mod_order(bytes: &[u8]) -> Fr {
        let num_modulus_bytes = Fr::MODULUS_BIT_SIZE.div_ceil(8) as usize;
        let num_bytes_to_directly_convert = ark_std::cmp::min(num_modulus_bytes - 1, bytes.len());

        // Copy the leading little-endian bytes directly into a field element.
        // The number of bytes directly converted must be less than the
        // number of bytes needed to represent the modulus, as we must begin
        // modular reduction once the data is of the same number of bytes as the
        // modulus.

        let (lo, hi) = bytes.split_at(bytes.len() - num_bytes_to_directly_convert);

        let lo: Fr = Fr::new_unchecked(BigInt::new(to_u64_array::<4>(lo)));
        let mut hi = Fr::new_unchecked(BigInt::new(to_u64_array::<4>(hi)));

        hi *= *FR_LO_BYTES_OFFSET;
        hi += lo;

        // Sanity-check equivalence with arkworks unoptimised implementation
        #[cfg(test)]
        assert_eq!(
            hi * Fr::new_unchecked(Fr::R2),
            <Fr as PrimeField>::from_le_bytes_mod_order(bytes)
        );

        hi
    }
}

pub(crate) fn is_valid<F>(value: F) -> bool
where
    F: PrimeField + ToBits,
{
    let value_bits = value.to_bits();
    let modulus_bits = F::MODULUS.to_bits_le();

    let N = F::MODULUS_BIT_SIZE as usize;
    for i in 0..N {
        if value_bits[N - i - 1] < modulus_bits[N - i - 1] {
            return true;
        } else if value_bits[N - i - 1] > modulus_bits[N - i - 1] {
            return false;
        }
    }
    false
}
