//! XOR garbling table

use std::marker::PhantomData;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::rngs::StdRng;
use serde::{Deserialize, Serialize};

use crate::{
    garble::{Circuit, Nonce},
    mac::{
        Mac,
        binary::{BinaryInput, BinaryInputMac, BinaryInputMacKey, BitMac, BitMacKey, Xor},
        bit::bytes::Bit16ByteMacKey,
    },
    utils::PRG,
};

pub(crate) type InputMacKey = BinaryInputMacKey<Xor>;

pub(crate) type Input = BinaryInput<Xor>;

pub(crate) type InputMac = BinaryInputMac<Xor>;

type OutputMac = BitMac<Xor>;
// phi has to be the same as InputMacKey
// k is the XOR of InputMacKey.left.k, InputMacKey.right.k
type OutputMacKey = BitMacKey<Xor>;

impl InputMacKey {
    pub(crate) fn output_mac_key(&self) -> OutputMacKey {
        let phi = self.left.phi;
        assert_eq!(phi, self.right.phi);

        let mut k = [0; 16];
        k.iter_mut()
            .zip(self.left.k.iter())
            .zip(self.right.k.iter())
            .for_each(|((k, left), right)| *k = left ^ right);

        // OutputMacKey phi has to be the same as InputMacKey
        // OutputMacKey k is the XOR of InputMacKey.left.k, InputMacKey.right.k
        let output_mac_key = BitMacKey {
            key: Bit16ByteMacKey { phi, k },
            _marker: PhantomData,
        };

        output_mac_key
    }
}

#[derive(Debug, Clone, CanonicalSerialize, CanonicalDeserialize, Serialize, Deserialize)]
pub(crate) struct Table;

impl Circuit<Input, InputMac, bool, OutputMac> for Table {
    type Aux = ();

    fn garble<P: PRG>(
        _parent_nonce: &mut Nonce,
        _prg_unseeded: &P::Unseeded,
        _rng: &mut StdRng,
        _aux: Option<&Self::Aux>,
        input_mac_key: &<BinaryInputMac<Xor> as Mac<BinaryInput<Xor>>>::Key,
    ) -> (Self, <OutputMac as Mac<bool>>::Key) {
        (Self, input_mac_key.output_mac_key())
    }

    fn evaluate<P: PRG>(
        &self,
        _parent_nonce: &mut Nonce,
        _prg_unseeded: &P::Unseeded,
        input: &BinaryInput<Xor>,
        input_mac: &BinaryInputMac<Xor>,
    ) -> (bool, OutputMac) {
        (input.xor(), input_mac.xor())
    }
}

/// FIXME test more cases
#[test]
fn test_xor_circuit() {
    use crate::mac::MacKey;
    use aes::Aes128;
    use ark_std::rand::SeedableRng;

    const SEED: u64 = 0;

    let mut rng = StdRng::seed_from_u64(SEED);
    let unseeded_prg = Aes128::init(&[]);

    let input = InputMac::random_input(&mut rng);
    let input_mac_key = InputMacKey::random_key(&mut rng);
    let output_mac_key = input_mac_key.output_mac_key();
    let input_mac = input_mac_key.encode(&input);

    let circuit = Table::garble::<Aes128>(
        &mut Nonce::new(),
        &unseeded_prg,
        &mut rng,
        None,
        &input_mac_key,
    )
    .0;
    let output_mac = circuit
        .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
        .1;

    let output = input.xor();
    let expected_output_mac = output_mac_key.encode(&output);

    assert_eq!(output_mac, expected_output_mac);
}
