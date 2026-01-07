use ark_std::rand::rngs::StdRng;

use crate::{
    garble::{
        Circuit, Nonce,
        bit_adaptor::{DeriveOutputMacKey, Table},
    },
    mac::{
        Encrypt, MacKey,
        bit::bytes::{Bit16ByteMac, Bit16ByteMacKey},
    },
    utils::PRG,
};

/// A bit adaptor table mapping from a Bit16ByteMac to a Bit16ByteMac.
impl Circuit<bool, Bit16ByteMac, bool, Bit16ByteMac> for Table<Bit16ByteMac, Bit16ByteMac> {
    type Aux = [u8; 16];

    fn garble<P: PRG>(
        parent_nonce: &mut Nonce,
        prg_unseeded: &P::Unseeded,
        _rng: &mut StdRng,
        aux: Option<&Self::Aux>,
        input_mac_key: &Bit16ByteMacKey,
    ) -> (Self, Bit16ByteMacKey) {
        // The chosen `phi` for the output_mac_key.
        let phi = aux.unwrap();

        // The output_mac_key is derived s.t.
        // - output_mac_key.phi = phi, and
        // - output_mac_key.k = hash(input_mac_key.encode(false)).
        let output_mac_key = <Bit16ByteMacKey as DeriveOutputMacKey<
            bool,
            bool,
            Bit16ByteMac,
            Bit16ByteMac,
        >>::derive_output_mac_key::<P>(
            &input_mac_key, prg_unseeded, &parent_nonce, phi
        );

        // We encrypt the encoding of output_mac_key.encode(true), using
        // input_mac_key.encode(true) as the key. This ciphertext is stored
        // in the bit adaptor garbling table.
        let key = input_mac_key.encode(&true);
        let msg = output_mac_key.encode(&true);

        let cyph = key.encrypt::<P>(prg_unseeded, parent_nonce, &msg);

        (Self(cyph), output_mac_key)
    }

    fn evaluate<P: PRG>(
        &self,
        parent_nonce: &mut Nonce,
        prg_unseeded: &P::Unseeded,
        input: &bool,
        input_mac: &Bit16ByteMac,
    ) -> (bool, Bit16ByteMac) {
        let output_mac = if !(*input) {
            // If the input is `false`, the `output_mac` will be
            //                  phi * 0 + k = k;
            // but since output_mac_key.k = hash(input_mac_key.encode(false)),
            // we can directly derive `k` by hashing the `input_mac`.
            P::seed(prg_unseeded, input_mac)
                .get_chunk(&parent_nonce.bytes())
                .into()
        } else {
            // If the input is `true`, we have to decrypt the table entry
            // using the `input_mac` as a key, to recover the plaintext
            // `output_mac`.
            input_mac.decrypt::<P>(prg_unseeded, parent_nonce, &self.0)
        };

        (*input, output_mac)
    }
}

impl DeriveOutputMacKey<bool, bool, Bit16ByteMac, Bit16ByteMac> for Bit16ByteMacKey {
    type Phi = [u8; 16];

    fn derive_output_mac_key<P: crate::utils::PRG>(
        &self,
        prg_unseeded: &P::Unseeded,
        parent_nonce: &Nonce,
        phi: &Self::Phi,
    ) -> Bit16ByteMacKey {
        let mac = self.encode(&false);
        let k = P::seed(prg_unseeded, &mac).get_chunk(&parent_nonce.bytes());

        Self { phi: *phi, k }
    }
}

impl<Mac1, Mac2> Encrypt<Mac1> for Mac2
where
    Mac1: Copy + Into<Bit16ByteMac> + From<Bit16ByteMac>,
    Mac2: Copy + Into<Bit16ByteMac> + From<Bit16ByteMac>,
{
    type Ciphertext = [u8; 16];

    fn encrypt<P: PRG>(
        &self,
        prg_unseeded: &P::Unseeded,
        parent_nonce: &mut Nonce,
        message: &Mac1,
    ) -> Self::Ciphertext {
        let seed: Bit16ByteMac = (*self).into();
        let message: Bit16ByteMac = (*message).into();

        // Seed the PRG with the input_mac.
        let xof = P::seed(prg_unseeded, &seed);

        let mut ciphertext = xof.get_chunk(&parent_nonce.bytes());
        ciphertext
            .iter_mut()
            .zip(message.iter())
            .for_each(|(pad, msg)| *pad ^= msg);

        ciphertext
    }

    fn decrypt<P: PRG>(
        &self,
        prg_unseeded: &P::Unseeded,
        parent_nonce: &mut Nonce,
        ciphertext: &Self::Ciphertext,
    ) -> Mac1 {
        let seed: Bit16ByteMac = (*self).into();

        // Seed the PRG with the input_mac.
        let xof = P::seed(prg_unseeded, &seed);

        let mut plaintext = xof.get_chunk(&parent_nonce.bytes());
        plaintext
            .iter_mut()
            .zip(ciphertext.iter())
            .for_each(|(pad, ciphertext)| *pad ^= ciphertext);

        Bit16ByteMac(plaintext).into()
    }
}
