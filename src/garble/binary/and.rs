//! AND garbling table

use std::marker::PhantomData;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::rngs::StdRng;

#[cfg(test)]
use crate::utils::prg::aes_ctr::clear_set;
use crate::{
    garble::{Circuit, Nonce, bit_adaptor},
    mac::{
        binary::{And, BinaryInput, BinaryInputMac, BinaryInputMacKey, BitMac, BitMacKey},
        bit::bytes::{Bit16ByteMac, Bit16ByteMacKey},
    },
    utils::PRG,
};

pub(crate) type InputMac = BinaryInputMac<And>;

pub(crate) type Input = BinaryInput<And>;
pub(crate) type InputMacKey = BinaryInputMacKey<And>;

type OutputMac = BitMac<And>;
// phi has to be the same as InputMacKey
// k is the And of InputMacKey.left.k, InputMacKey.right.k
type OutputMacKey = BitMacKey<And>;

#[derive(CanonicalSerialize, CanonicalDeserialize, Debug)]
pub(crate) struct Table(bit_adaptor::Table<Bit16ByteMac, Bit16ByteMac>);

impl Circuit<Input, InputMac, bool, OutputMac> for Table {
    type Aux = ();

    fn garble<P: PRG>(
        parent_nonce: &mut Nonce,
        prg_unseeded: &P::Unseeded,
        rng: &mut StdRng,
        _aux: Option<&Self::Aux>,
        input_mac_key: &InputMacKey,
    ) -> (Self, OutputMacKey) {
        let mut nonce = parent_nonce.fork();

        // Bit adaptor table mapping from a MAC under input_mac_key.right to a
        // MAC under a mac_key s.t.
        // - mac_key.phi = input_mac_key.left.k,
        // - mac_key.k = hash(input_mac_key.right.encode(0))
        // This stores the encryption of mac_key.encode(1),
        // using input_mac_key.right.encode(1) as the key
        let (adaptor_table, mut output_mac_key): (
            bit_adaptor::Table<Bit16ByteMac, Bit16ByteMac>,
            Bit16ByteMacKey,
        ) = bit_adaptor::Table::garble::<P>(
            &mut nonce,
            prg_unseeded,
            rng,
            // `k` from the input_mac_key of the left wire
            Some(&input_mac_key.left.k),
            &input_mac_key.right.key,
        );

        output_mac_key.phi = input_mac_key.left.phi;

        (
            Self(adaptor_table),
            OutputMacKey {
                key: output_mac_key,
                _marker: PhantomData,
            },
        )
    }

    fn evaluate<P: PRG>(
        &self,
        parent_nonce: &mut Nonce,
        prg_unseeded: &P::Unseeded,
        input: &Input,
        input_mac: &InputMac,
    ) -> (bool, OutputMac) {
        let mut nonce = parent_nonce.fork();

        let mac_left = if input.right() {
            // If input.right = 1, use input_mac.left as mac_left, i.e.
            input_mac.left()
        } else {
            // If input.right = 0, use the all-zeros string as mac_left
            BitMac {
                mac: Bit16ByteMac::from([0; 16]),
                _marker: PhantomData,
            }
        };

        let mac_right = if input.right() {
            // If input.right = 1, use the input_mac.right as the key
            // to decrypt the mac_right, to get
            // input_mac_key.left.k ^ hash(input_mac_key.right.encode(0))
            self.0
                .evaluate::<P>(&mut nonce, prg_unseeded, &input.right(), &input_mac.right())
                .1
        } else {
            // If input.right = 0, derive mac_right by hashing input_mac.right
            // to get hash(input_mac_key.right.encode(0))
            let xof = P::seed(prg_unseeded, &input_mac.right());
            xof.get_chunk(&nonce.bytes()).into()
        };

        // Recall that our output_mac_key is such that:
        // - output_mac_key.phi = input_mac_key.left.phi,
        // - output_mac_key.k = hash(input_mac_key.right.encode(0))
        //
        // Enumerating through the cases:
        // 1. left = 0, right = 0, output = 0
        //    mac_left = [0; 16]
        //    mac_right = hash(input_mac_key.right.encode(0))
        //    mac_left ^ mac_right = hash(input_mac_key.right.encode(0)),
        //    which recovers the MAC of 0 under output_mac_key.
        // 2. left = 0, right = 1, output = 0
        //    mac_left = input_mac.left = input_mac_key.left.k
        //    mac_right = input_mac_key.left.k ^ hash(input_mac_key.right.encode(0))
        //    mac_left ^ mac_right = hash(input_mac_key.right.encode(0)),
        //    which recovers the MAC of 0 under output_mac_key.
        // 3. left = 1, right = 0, output = 0
        //    mac_left = [0; 16]
        //    mac_right = hash(input_mac_key.right.encode(0))
        //    mac_left ^ mac_right = hash(input_mac_key.right.encode(0)),
        //    which recovers the MAC of 0 under output_mac_key.
        // 4. left = 1, right = 1, output = 1
        //    mac_left = input_mac.left = input_mac_key.left.phi ^ input_mac_key.left.k
        //    mac_right = input_mac_key.left.k ^ hash(input_mac_key.right.encode(0))
        //    mac_left ^ mac_right = input_mac_key.left.phi ^ hash(input_mac_key.right.encode(0)),
        //    which recovers the MAC of 1 under output_mac_key.
        (input.and(), mac_left ^ mac_right.into())
    }
}

#[test]
fn test_and_circuit() {
    use crate::mac::{Mac, MacKey};
    use aes::Aes128;
    use ark_std::rand::SeedableRng;

    const N: usize = 100;
    const SEED: u64 = 0;

    let mut rng = StdRng::seed_from_u64(SEED);
    // FIXME use TCCR for this
    let unseeded_prg = Aes128::init(&[]);

    for _ in 0..N {
        let input = InputMac::random_input(&mut rng);

        let input_mac_key = InputMacKey::random_key(&mut rng);

        let input_mac = input_mac_key.encode(&input);

        let (circuit, output_mac_key) = Table::garble::<Aes128>(
            &mut Nonce::new(),
            &unseeded_prg,
            &mut rng,
            None,
            &input_mac_key,
        );
        clear_set();

        let output_mac = circuit
            .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
            .1;
        clear_set();

        let output = input.and();
        let expected_output_mac = output_mac_key.encode(&output);

        assert_eq!(output_mac, expected_output_mac);
    }
}
