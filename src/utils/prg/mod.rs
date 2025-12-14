//! Pseudo-random generator used in the bit_adaptor table.
//! This implementation is tailored for AES128.

use ark_bn254::Fq;
use ark_ec::AdditiveGroup;
use ark_ff::{Field, PrimeField};

use crate::{garble::Nonce, utils::group::from_le_bytes_mod_order};

pub(crate) mod aes_ctr;

pub(crate) const AES_BLOCK_SIZE: usize = 16;

/// Number of random bytes needed for Fq hash_to_field at 128 bits of security
const FQ_LEN_PER_BASE_ELEM: usize = 48;

/// Pseudo-random generator with random access to output bytestream
pub trait PRG: Clone {
    type Unseeded;

    fn init(dst: &[u8]) -> Self::Unseeded;

    fn seed(unseeded: &Self::Unseeded, seed: &[u8; AES_BLOCK_SIZE]) -> Self;

    fn get_chunk(&self, counter: &[u8; AES_BLOCK_SIZE]) -> [u8; AES_BLOCK_SIZE];

    fn hash_to_field<F>(&self, parent_nonce: &Nonce) -> F
    where
        F: Field<BasePrimeField = Fq>,
    {
        let parent_nonce = &mut parent_nonce.clone();

        // The user requests an element of F_p^m, comprising `m` BasePrimeField elements.
        let m = F::extension_degree() as usize;
        let mut base_prime_field_elems = vec![F::BasePrimeField::ZERO; m];

        // The number of blocks required for each base prime field element
        let num_blocks_per_base_elem = FQ_LEN_PER_BASE_ELEM / AES_BLOCK_SIZE;

        base_prime_field_elems.iter_mut().for_each(|elem| {
            let mut bytes = Vec::with_capacity(FQ_LEN_PER_BASE_ELEM);
            (0..num_blocks_per_base_elem).for_each(|_| {
                // Fork the nonce to get a unique counter
                let nonce = parent_nonce.fork();

                // Get the chunk of the PRG corresponding to this counter
                let block = self.get_chunk(&nonce.bytes());

                bytes.extend_from_slice(&block);
            });

            *elem = from_le_bytes_mod_order(&bytes);
        });

        F::from_base_prime_field_elems(base_prime_field_elems).unwrap()
    }

    fn num_blocks<F: Field>() -> usize {
        let base_prime_field_num_bytes = F::BasePrimeField::MODULUS_BIT_SIZE.div_ceil(8) as usize;
        assert_eq!(base_prime_field_num_bytes & AES_BLOCK_SIZE, 0);
        let num_blocks = base_prime_field_num_bytes / AES_BLOCK_SIZE;

        F::extension_degree() as usize * num_blocks
    }
}
