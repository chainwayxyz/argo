//! Given an `input_mac_key` and a value `phi`, derives an `output_mac_key`.
//! Maps the MAC of a single bit under input_mac_key to its MAC under the
//! derived output_mac_key.

pub(crate) mod bytes;
pub(crate) mod group;

use std::fmt::Debug;

use ark_bn254::Fq;
use ark_ff::Field;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

use crate::{
    garble::Nonce,
    mac::{
        Encrypt, Mac, MacKey, Random,
        bit::{
            bytes::{Bit16ByteMac, Bit16ByteMacKey},
            group::BitGroupMac,
        },
        group::GroupMacKey,
    },
    utils::prg::PRG,
};

#[derive(CanonicalSerialize, CanonicalDeserialize, Debug)]
pub(crate) struct Table<InputMac, OutputMac>(
    pub(crate) <InputMac as Encrypt<OutputMac>>::Ciphertext,
)
where
    InputMac: Mac<bool> + Encrypt<OutputMac>,
    OutputMac: Mac<bool>,
    InputMac::Key: DeriveOutputMacKey<bool, bool, InputMac, OutputMac>,
    <InputMac as Encrypt<OutputMac>>::Ciphertext: Debug;

pub(crate) trait DeriveOutputMacKey<Input, Output, InputMac, OutputMac>
where
    InputMac: Mac<Input>,
    OutputMac: Mac<Output>,
{
    type Phi;

    fn derive_output_mac_key<P: PRG>(
        &self,
        prg_unseeded: &P::Unseeded,
        parent_nonce: &Nonce,
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
        parent_nonce: &Nonce,
        phi: &G,
    ) -> GroupMacKey<G> {
        let mac = self.encode(&false);
        let k = P::seed(prg_unseeded, &mac).hash_to_field(parent_nonce);

        GroupMacKey { phi: *phi, k }
    }
}
