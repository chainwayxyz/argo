#[cfg(test)]
use std::{cell::RefCell, collections::BTreeSet};
use std::{collections::BTreeMap, sync::Mutex};

use aes::{
    Aes128, Block,
    cipher::{BlockEncrypt, KeyInit},
};
use once_cell::sync::Lazy;

use crate::utils::{PRG, prg::AES_BLOCK_SIZE};

static AES_CTRS: Lazy<Mutex<BTreeMap<[u8; AES_BLOCK_SIZE], Aes128>>> =
    Lazy::new(|| Mutex::new(BTreeMap::new()));

impl PRG for Aes128 {
    type Unseeded = &'static Lazy<Mutex<BTreeMap<[u8; AES_BLOCK_SIZE], Aes128>>>;

    fn init(_dst: &[u8]) -> Self::Unseeded {
        AES_CTRS.lock().unwrap().clear();
        &AES_CTRS
    }

    fn seed(unseeded: &Self::Unseeded, seed: &[u8; AES_BLOCK_SIZE]) -> Self {
        unseeded
            .lock()
            .unwrap()
            .entry(*seed)
            .or_insert_with(|| Aes128::new(&(*seed).into()))
            .clone()
    }

    fn get_chunk(&self, counter: &[u8; AES_BLOCK_SIZE]) -> [u8; AES_BLOCK_SIZE] {
        let mut block: Block = (*counter).into();
        self.encrypt_block(&mut block);

        #[cfg(test)]
        {
            // Check that each counter is used only once per seed.
            // We cannot extract the seed from the PRG, so use counter 0 value instead
            let mut block = Block::default();
            self.encrypt_block(&mut block);
            assert!(COUNTERS_SET.with(|map| {
                map.borrow_mut()
                    .entry(block.into())
                    .or_default()
                    .insert(*counter)
            }));
        }

        block.into()
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
