#[cfg(test)]
use std::{cell::RefCell, collections::BTreeSet};
use std::{collections::BTreeMap, sync::Mutex};

use aes::{
    Aes128, Block,
    cipher::{BlockEncrypt, KeyInit},
};
use ark_ec::AdditiveGroup;
use ark_ff::{Field, PrimeField};
use once_cell::sync::Lazy;

use crate::{
    garble::Nonce,
    utils::{
        aes_hash::{AES_BLOCK_SIZE, Hash, LEN_PER_BASE_ELEM},
        group::FromLeBytesModOrder,
    },
};

static AES_CTRS: Lazy<Mutex<BTreeMap<[u8; AES_BLOCK_SIZE], Aes128>>> =
    Lazy::new(|| Mutex::new(BTreeMap::new()));

#[derive(Debug)]
pub(crate) struct AesCtr;

impl Hash for AesCtr {
    fn hash(key: &[u8; 16], msg: &[u8; 16]) -> [u8; 16] {
        let xof = AES_CTRS
            .lock()
            .unwrap()
            .entry(*key)
            .or_insert_with(|| Aes128::new(key.into()))
            .clone();

        let mut block: Block = (*msg).into();
        xof.encrypt_block(&mut block);

        #[cfg(test)]
        {
            // Check that each counter is used only once per key.
            // We cannot extract the key from the hasher, so use counter 0 value instead
            let mut zero_block = Block::default();
            xof.encrypt_block(&mut zero_block);
            assert!(COUNTERS_SET.with(|map| {
                map.borrow_mut()
                    .entry(zero_block.into())
                    .or_default()
                    .insert(*msg)
            }));
        }

        block.into()
    }
}

impl AesCtr {
    pub(crate) fn hash_to_field<F>(value: [u8; AES_BLOCK_SIZE], parent_nonce: &Nonce) -> F
    where
        F: Field,
        F::BasePrimeField: FromLeBytesModOrder,
    {
        let parent_nonce = &mut parent_nonce.clone();

        // The user requests an element of F_p^m, comprising `m` BasePrimeField elements.
        let m = F::extension_degree() as usize;
        let mut base_prime_field_elems = vec![F::BasePrimeField::ZERO; m];

        // The number of blocks required for each base prime field element
        let num_blocks_per_base_elem = LEN_PER_BASE_ELEM / AES_BLOCK_SIZE;

        base_prime_field_elems.iter_mut().for_each(|elem| {
            let mut bytes = Vec::with_capacity(LEN_PER_BASE_ELEM);
            (0..num_blocks_per_base_elem).for_each(|_| {
                // Fork the nonce to get a unique counter
                let nonce = parent_nonce.fork();

                // Get the chunk of the PRG corresponding to this counter
                let block = AesCtr::hash(&value, &nonce.bytes());

                bytes.extend_from_slice(&block);
            });

            *elem = <F::BasePrimeField as FromLeBytesModOrder>::from_le_bytes_mod_order(&bytes);
        });

        F::from_base_prime_field_elems(base_prime_field_elems).unwrap()
    }

    pub(crate) fn num_blocks<F: Field>() -> usize {
        let base_prime_field_num_bytes = F::BasePrimeField::MODULUS_BIT_SIZE.div_ceil(8) as usize;
        assert_eq!(base_prime_field_num_bytes & AES_BLOCK_SIZE, 0);
        let num_blocks = base_prime_field_num_bytes / AES_BLOCK_SIZE;

        F::extension_degree() as usize * num_blocks
    }
}

#[cfg(test)]
thread_local! {
    /// Mapping from seed to set of used counters; a counter should be used
    /// only once per seed
    pub(crate) static COUNTERS_SET: RefCell<BTreeMap<[u8; AES_BLOCK_SIZE], BTreeSet<[u8; AES_BLOCK_SIZE]>>> = const { RefCell::new(BTreeMap::new()) }
}

#[cfg(test)]
/// Clear the set
pub(crate) fn clear_set() {
    COUNTERS_SET.with(|map| {
        map.borrow_mut().clear();
    });
}
