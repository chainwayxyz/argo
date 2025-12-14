use crate::mac::{Mac, MacKey};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::{RngCore as _, rngs::StdRng};
use derive_more::{Deref, From};

/// MAC from a bit to a 16-byte word
#[derive(Debug, Clone, Copy, CanonicalSerialize, CanonicalDeserialize, PartialEq, From, Deref)]
pub(crate) struct Bit16ByteMac([u8; 16]);

#[derive(Debug, Clone, Copy, CanonicalSerialize, CanonicalDeserialize)]
pub(crate) struct Bit16ByteMacKey {
    phi: [u8; 16],
    k: [u8; 16],
}

impl Bit16ByteMacKey {
    pub(crate) fn phi_xor_k(&self) -> [u8; 16] {
        let mut bytes = self.phi;
        bytes
            .iter_mut()
            .zip(self.k.iter())
            .for_each(|(phi, k)| *phi ^= k);

        bytes
    }
}

impl MacKey<bool, Bit16ByteMac> for Bit16ByteMacKey {
    fn encode(&self, value: &bool) -> Bit16ByteMac {
        if *value {
            self.phi_xor_k().into()
        } else {
            self.k.into()
        }
    }

    #[cfg(test)]
    fn decode(&self, mac: &Bit16ByteMac) -> bool {
        if mac.0 == self.k {
            false
        } else {
            assert_eq!(mac.0, self.phi_xor_k());
            true
        }
    }

    fn random_key(rng: &mut StdRng) -> Self {
        let mut phi = [0; 16];
        rng.fill_bytes(&mut phi);

        let mut k = [0; 16];
        rng.fill_bytes(&mut k);

        Self { phi, k }
    }
}

impl Mac<bool> for Bit16ByteMac {
    type Key = Bit16ByteMacKey;

    #[cfg(test)]
    fn random_input(rng: &mut StdRng) -> bool {
        use ark_std::rand::Rng;
        rng.gen_bool(0.5)
    }
}
