#[cfg(test)]
use std::collections::{BTreeMap, BTreeSet};
use std::{cell::RefCell, rc::Rc};

use aes::{
    Aes128, Block,
    cipher::{BlockEncrypt, KeyInit},
};
use ahash::AHashMap;
use ark_ec::AdditiveGroup;
use ark_ff::{Field, PrimeField};

use crate::{
    garble::Nonce,
    utils::{
        aes_hash::{AES_BLOCK_SIZE, Hash, LEN_PER_BASE_ELEM},
        group::FromLeBytesModOrder,
    },
};

thread_local! {
    static AES_INSTANCES: RefCell<AHashMap<[u8; 16], Rc<Aes128>>> =
        RefCell::new(AHashMap::default());
}

#[derive(Debug)]
pub(crate) struct AesCtr;

impl AesCtr {
    pub(crate) fn get_instance(key: &[u8; 16]) -> Rc<Aes128> {
        AES_INSTANCES.with(|instances| {
            Rc::clone(
                instances
                    .borrow_mut()
                    .entry(*key)
                    .or_insert_with(|| Rc::new(Aes128::new(&(*key).into()))),
            )
        })
    }
}

impl Hash for AesCtr {
    #[inline(always)]
    fn hash(key: &[u8; 16], msg: &[u8; 16]) -> [u8; 16] {
        let mut block = Block::from(*msg);
        Self::get_instance(key).encrypt_block(&mut block);

        #[cfg(test)]
        {
            // Check that each counter is used only once per key.
            // We cannot extract the key from the hasher, so use counter 0 value instead
            let mut zero_block = Block::default();
            Aes128::new(key.into()).encrypt_block(&mut zero_block);
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

        let m = F::extension_degree() as usize;
        let mut base_prime_field_elems = vec![F::BasePrimeField::ZERO; m];
        let num_blocks_per_base_elem = LEN_PER_BASE_ELEM / AES_BLOCK_SIZE;

        let aes = AesCtr::get_instance(&value);

        // Reuse a single buffer instead of allocating per element
        let mut bytes = [0u8; LEN_PER_BASE_ELEM];

        for elem in base_prime_field_elems.iter_mut() {
            // Generate blocks directly into the buffer
            for i in 0..num_blocks_per_base_elem {
                let nonce = parent_nonce.fork();
                let mut block = Block::from(nonce.bytes());
                aes.encrypt_block(&mut block);

                let start = i * AES_BLOCK_SIZE;
                bytes[start..start + AES_BLOCK_SIZE].copy_from_slice(&block);
            }

            *elem = <F::BasePrimeField as FromLeBytesModOrder>::from_le_bytes_mod_order(&bytes);
        }

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
