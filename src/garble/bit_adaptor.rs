//! Given an `input_mac_key` and a value `phi`, derives an `output_mac_key`.
//! Maps the MAC of a single bit under input_mac_key to its MAC under the
//! derived output_mac_key.

use ark_bn254::Fq;
use ark_ff::{BigInteger, Field, PrimeField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::rngs::StdRng;

use crate::{
    garble::{Circuit, Nonce},
    mac::{
        Encrypt, Mac, MacKey, Random,
        bit::{
            bytes::{Bit16ByteMac, Bit16ByteMacKey},
            group::BitGroupMac,
        },
        group::GroupMacKey,
    },
    utils::prg::{AES_BLOCK_SIZE, PRG},
};

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub(crate) struct Table<InputMac, OutputMac>(
    pub(crate) <InputMac as Encrypt<OutputMac>>::Ciphertext,
)
where
    InputMac: Mac<bool> + Encrypt<OutputMac>,
    OutputMac: Mac<bool>,
    InputMac::Key: DeriveOutputMacKey<bool, bool, InputMac, OutputMac>;

/// A bit adaptor table mapping from a Bit16ByteMac to a BitGroupMac.
impl<G> Circuit<bool, Bit16ByteMac, bool, BitGroupMac<G>> for Table<Bit16ByteMac, BitGroupMac<G>>
where
    G: Field<BasePrimeField = Fq> + Random,
{
    type Aux = G;

    fn garble<P: PRG>(
        parent_nonce: &mut Nonce,
        prg_unseeded: &P::Unseeded,
        _rng: &mut StdRng,
        aux: Option<&Self::Aux>,
        input_mac_key: &Bit16ByteMacKey,
    ) -> (Self, GroupMacKey<G>) {
        let mut nonce = parent_nonce.fork();
        // The chosen `phi` for the output_mac_key.
        let phi = aux.unwrap();

        // The output_mac_key is derived s.t.
        // - output_mac_key.phi = phi, and
        // - output_mac_key.k = hash(input_mac_key.encode(false)).
        let output_mac_key =
            input_mac_key.derive_output_mac_key::<P>(prg_unseeded, &mut nonce, phi);

        // We encrypt the encoding of output_mac_key.encode(true), using
        // input_mac_key.encode(true) as the key. This ciphertext is stored
        // in the bit adaptor garbling table.
        let key = input_mac_key.encode(&true);
        let msg = output_mac_key.encode(&true);
        let cyph = key.encrypt::<P>(prg_unseeded, &mut nonce, &msg);

        (Self(cyph), output_mac_key)
    }

    fn evaluate<P: PRG>(
        &self,
        parent_nonce: &mut Nonce,
        prg_unseeded: &P::Unseeded,
        input: &bool,
        input_mac: &Bit16ByteMac,
    ) -> (bool, BitGroupMac<G>) {
        let mut nonce = parent_nonce.fork();

        let output_mac = if !(*input) {
            // If the input is `false`, the `output_mac` will be
            //                  phi * 0 + k = k;
            // but since output_mac_key.k = hash(input_mac_key.encode(false)),
            // we can directly derive `k` by hashing the `input_mac`.
            P::seed(prg_unseeded, input_mac)
                .hash_to_field::<G>(&nonce)
                .into()
        } else {
            // If the input is `true`, we have to decrypt the table entry
            // using the `input_mac` as a key, to recover the plaintext
            // `output_mac`.
            input_mac.decrypt::<P>(prg_unseeded, &mut nonce, &self.0)
        };

        (*input, output_mac)
    }
}

pub(crate) trait DeriveOutputMacKey<Input, Output, InputMac, OutputMac>
where
    InputMac: Mac<Input>,
    OutputMac: Mac<Output>,
{
    type Phi;

    fn derive_output_mac_key<P: PRG>(
        &self,
        prg_unseeded: &P::Unseeded,
        parent_nonce: &mut Nonce,
        phi: &Self::Phi,
    ) -> OutputMac::Key;
}

/// Given an input_mac_key: Bit16ByteMacKey and a group element phi: G,
/// derive an output_mac_key: GroupMacKey with s.t.
/// - output_mac_key.phi = phi, and
/// - output_mac_key.k = hash(input_mac_key.encode(false)).
impl<G> DeriveOutputMacKey<bool, bool, Bit16ByteMac, BitGroupMac<G>> for Bit16ByteMacKey
where
    G: Field<BasePrimeField = Fq> + Random,
{
    type Phi = G;

    fn derive_output_mac_key<P: PRG>(
        &self,
        prg_unseeded: &P::Unseeded,
        parent_nonce: &mut Nonce,
        phi: &G,
    ) -> GroupMacKey<G> {
        let mac = self.encode(&false);
        let k = P::seed(prg_unseeded, &mac).hash_to_field(parent_nonce);

        GroupMacKey { phi: *phi, k }
    }
}

/// Encrypt a BitGroupMac(G) using a Bit16ByteMac([u8; 16]) as the key.
/// Produces a ciphertext Vec<u8> of the same length as G::NUM_BYTES.
impl<G> Encrypt<BitGroupMac<G>> for Bit16ByteMac
where
    G: Field<BasePrimeField = Fq>,
{
    type Ciphertext = Vec<u8>;

    fn encrypt<P: PRG>(
        &self,
        prg_unseeded: &P::Unseeded,
        parent_nonce: &mut Nonce,
        message: &BitGroupMac<G>,
    ) -> Self::Ciphertext {
        // Seed the PRG with the input_mac.
        let xof = P::seed(prg_unseeded, self);

        // Initialise the ciphertext as an empty vector with the desired
        // number of AES blocks.
        let mut ciphertext = Vec::with_capacity(P::num_blocks::<G>() * AES_BLOCK_SIZE);
        (0..P::num_blocks::<G>()).for_each(|_| {
            // Use a new block counter for each block.
            let nonce = parent_nonce.fork();

            // Get the XOF output for this block counter.
            let block = xof.get_chunk(&nonce.bytes());

            // Copy the output block to the ciphertext.
            ciphertext.extend_from_slice(&block);
        });

        // Get the byte representation of the output_mac group element.
        let message_bytes = message
            .to_base_prime_field_elements()
            .flat_map(|elem| elem.into_bigint().to_bytes_le());

        // Encrypt the message using the ciphertext as a one-time pad.
        ciphertext
            .iter_mut()
            .zip(message_bytes)
            .for_each(|(pad, msg)| *pad ^= msg);

        ciphertext
    }

    fn decrypt<P: PRG>(
        &self,
        prg_unseeded: &P::Unseeded,
        parent_nonce: &mut Nonce,
        ciphertext: &Self::Ciphertext,
    ) -> BitGroupMac<G> {
        // Seed the XOF with the input_mac.
        let xof = P::seed(prg_unseeded, self);

        // Initialise the plaintext as an empty vector with the desired
        // number of AES blocks.
        let mut plaintext = Vec::with_capacity(P::num_blocks::<G>() * AES_BLOCK_SIZE);

        // Recover the one-time pad used to encrypt the ciphertext.
        (0..P::num_blocks::<G>()).for_each(|_| {
            let nonce = parent_nonce.fork();

            let block = xof.get_chunk(&nonce.bytes());
            plaintext.extend_from_slice(&block);
        });

        // XOR the ciphertext with the one-time pad to recover the plaintext
        // byte representation of the group element.
        plaintext
            .iter_mut()
            .zip(ciphertext.iter())
            .for_each(|(pad, ciph)| *pad ^= ciph);

        // Cast the plaintext bytes as a group element to recover output_mac.
        G::from_random_bytes(&plaintext).unwrap().into()
    }
}

#[cfg(test)]
fn test_bit_field_adaptor<G>(rng: &mut StdRng)
where
    G: Field<BasePrimeField = Fq> + Random,
{
    use crate::utils::prg::aes_ctr::clear_set;
    use aes::Aes128;

    let unseeded_prg = Aes128::init(&[]);

    let input_mac_key = Bit16ByteMacKey::random_key(rng);
    let phi = G::rand(rng);

    let (garbling_table, output_mac_key) = Table::<Bit16ByteMac, BitGroupMac<G>>::garble::<Aes128>(
        &mut Nonce::new(),
        &unseeded_prg,
        rng,
        Some(&phi),
        &input_mac_key,
    );
    clear_set();

    let input = Bit16ByteMac::random_input(rng);
    let input_mac = input_mac_key.encode(&input);

    let output = input;
    let output_mac: BitGroupMac<G> = output_mac_key.encode(&output);

    let evaluated_output_mac: BitGroupMac<G> = garbling_table
        .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
        .1;
    clear_set();

    assert!(output_mac == evaluated_output_mac);
}

#[test]
fn test_field_2_field_adaptor() {
    use ark_bn254::Fq as F;
    use ark_std::rand::SeedableRng;

    const N: usize = 128;
    const SEED: u64 = 0;

    let mut rng = StdRng::seed_from_u64(SEED);

    for _ in 0..N {
        test_bit_field_adaptor::<F>(&mut rng);
    }
}

#[test]
fn test_fp2_2_fp2_adaptor() {
    use ark_bn254::Fq2;
    use ark_std::rand::SeedableRng;

    const N: usize = 128;
    const SEED: u64 = 0;

    let mut rng = StdRng::seed_from_u64(SEED);

    for _ in 0..N {
        test_bit_field_adaptor::<Fq2>(&mut rng);
    }
}
