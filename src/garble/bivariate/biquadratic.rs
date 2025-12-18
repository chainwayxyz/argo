//! Garbling table for the oblivious evaluation of a quadratic equation in
//! two variables.

use ark_bn254::Fq;
use ark_ff::Field;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::rngs::StdRng;
use std::{collections::BTreeSet, marker::PhantomData};

use crate::{
    garble::{
        Circuit, Nonce,
        bivariate::{Input, InputMac, InputMacKey},
        digit_adaptor,
    },
    mac::Random,
    mac::{IdentityMac, IdentityMacKey},
    utils::{
        bits::{FromBits, NumBits, ToBits},
        prg::PRG,
    },
};

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub(crate) struct Table<F>
where
    F: Field<BasePrimeField = Fq> + NumBits + FromBits<bool, F> + Random,
{
    c_primes: [Option<F>; 6],
    // Same x input to both adaptors
    // c'_7, c'_9
    x_adaptors: [Option<digit_adaptor::Table<F>>; 2],
    // Same y input to all three adaptors
    // c'_6, c'_8, c'_10
    y_adaptors: [Option<digit_adaptor::Table<F>>; 3],
}

impl<F> Table<F>
where
    F: Field<BasePrimeField = Fq> + NumBits + FromBits<bool, F> + Random,
{
    pub(crate) fn c_prime_indices(cs: &[bool; 6]) -> BTreeSet<usize> {
        //   c'0 + c'1 X + c'2 Y + c'3 X^2 + c'4 XY + c'5 Y^2 + c'6(Y) X
        // + c'7(X) X + c'8(Y) Y + c'9(X) + c'10(Y)
        let mut c_prime_indices = BTreeSet::new();

        let [constant, x, y, x_sq, xy, y_sq] = cs;
        if *constant {
            c_prime_indices.insert(0);
        }

        // # x   -> (c_1, c_9)
        if *x {
            c_prime_indices.insert(1);
            c_prime_indices.insert(9);
        }

        // # y   -> (c_2, c_10)
        if *y {
            c_prime_indices.insert(2);
            c_prime_indices.insert(10);
        }

        // # x^2 -> (c_3, c_7, c_9)
        if *x_sq {
            c_prime_indices.insert(3);
            c_prime_indices.insert(7);
            c_prime_indices.insert(9);
        }

        // # x y -> (c_4, c_6, c_9)
        if *xy {
            c_prime_indices.insert(4);
            c_prime_indices.insert(6);
            c_prime_indices.insert(9);
        }

        // # y^2 -> (c_5, c_8, c_10)
        if *y_sq {
            c_prime_indices.insert(5);
            c_prime_indices.insert(8);
            c_prime_indices.insert(10);
        }

        c_prime_indices
    }

    //   c'0 + c'1 X + c'2 Y + c'3 X^2 + c'4 XY + c'5 Y^2 + c'6(Y) X
    // + c'7(X) X + c'8(Y) Y + c'9(X) + c'10(Y)
    fn coord(input: &Input<F>, c_primes: &[Option<F>; 11]) -> F {
        let Input { x, y } = input;

        let [
            c_prime_0,
            c_prime_1,
            c_prime_2,
            c_prime_3,
            c_prime_4,
            c_prime_5,
            c_prime_6,
            c_prime_7,
            c_prime_8,
            c_prime_9,
            c_prime_10,
        ] = c_primes.map(|val| val.map_or(F::ZERO, |v| v));

        c_prime_0
            + c_prime_1 * x
            + c_prime_2 * y
            + c_prime_3 * x.square()
            + c_prime_4 * x * y
            + c_prime_5 * y.square()
            + c_prime_6 * x
            + c_prime_7 * x
            + c_prime_8 * y
            + c_prime_9
            + c_prime_10
    }

    // c_0 + c_1 X + c_2 Y + c_3 X^2 + c_4 XY + c_5 Y^2
    #[cfg(test)]
    fn coord_plaintext(input: &Input<F>, c: &[F; 6]) -> F {
        let Input { x, y } = input;

        c[0] + c[1] * x + c[2] * y + c[3] * x.square() + c[4] * x * y + c[5] * y.square()
    }
}

impl<F> Circuit<Input<F>, InputMac<F>, F, IdentityMac<F>> for Table<F>
where
    F: Field<BasePrimeField = Fq> + NumBits + ToBits + FromBits<F, F> + FromBits<bool, F> + Random,
{
    // c_i's
    type Aux = ([F; 6], [bool; 6]);

    fn garble<P: PRG>(
        parent_nonce: &mut Nonce,
        prg_unseeded: &P::Unseeded,
        rng: &mut StdRng,
        aux: Option<&Self::Aux>,
        input_mac_key: &InputMacKey<F>,
    ) -> (Self, IdentityMacKey) {
        let mut nonce = parent_nonce.fork();

        let (cs, cs_flags) = aux.unwrap();
        let c_prime_indices = Self::c_prime_indices(cs_flags);

        let mut rs = [F::ZERO; 11];

        for i in 1..=5 {
            if c_prime_indices.contains(&i) {
                rs[i] = F::rand(rng);
            }
        }

        // Construct an adaptor table from input_mac_key.y, to
        // a triple of DigitMacKey's using:
        // - { phi: -r_4, k: r_6 }
        // - { phi: -r_5, k: r_8 }
        // - { phi: -(r_2 + r_8), k: r_10 }
        let y_adaptors = {
            // c'_6
            let table_0 = if c_prime_indices.contains(&6) {
                let (table, output_mac_key) = digit_adaptor::Table::garble::<P>(
                    &mut nonce,
                    prg_unseeded,
                    rng,
                    Some(&-rs[4]),
                    &input_mac_key.y,
                );
                rs[6] = output_mac_key.from_bits_k();

                Some(table)
            } else {
                assert_eq!(rs[4], F::ZERO);
                assert_eq!(rs[6], F::ZERO);

                None
            };

            // c'_8
            let table_1 = if c_prime_indices.contains(&8) {
                let (table, output_mac_key) = digit_adaptor::Table::garble::<P>(
                    &mut nonce,
                    prg_unseeded,
                    rng,
                    Some(&-rs[5]),
                    &input_mac_key.y,
                );
                rs[8] = output_mac_key.from_bits_k();

                Some(table)
            } else {
                assert_eq!(rs[5], F::ZERO);
                assert_eq!(rs[8], F::ZERO);

                None
            };

            // c'_10
            let table_2 = if c_prime_indices.contains(&10) {
                let (table, output_mac_key) = digit_adaptor::Table::garble::<P>(
                    &mut nonce,
                    prg_unseeded,
                    rng,
                    Some(&-(rs[2] + rs[8])),
                    &input_mac_key.y,
                );
                rs[10] = output_mac_key.from_bits_k();

                Some(table)
            } else {
                assert_eq!(rs[10], F::ZERO);
                assert_eq!(-(rs[2] + rs[8]), F::ZERO);

                None
            };

            [table_0, table_1, table_2]
        };

        // Construct an adaptor table from input_mac_key.x, to
        // a pair of DigitMacKey's using:
        // - { phi: -r_3, k: r_7 }
        // - { phi: -(r_1 + r_6 + r_7), k: r_9 }
        let x_adaptors = {
            // c'_7
            let table_0 = if c_prime_indices.contains(&7) {
                let (table, output_mac_key) = digit_adaptor::Table::garble::<P>(
                    &mut nonce,
                    prg_unseeded,
                    rng,
                    Some(&-rs[3]),
                    &input_mac_key.x,
                );
                rs[7] = output_mac_key.from_bits_k();

                Some(table)
            } else {
                assert_eq!(rs[3], F::ZERO);
                assert_eq!(rs[7], F::ZERO);
                None
            };

            // c'_9
            let table_1 = if c_prime_indices.contains(&9) {
                let (table, output_mac_key) = digit_adaptor::Table::garble::<P>(
                    &mut nonce,
                    prg_unseeded,
                    rng,
                    Some(&-(rs[1] + rs[6] + rs[7])),
                    &input_mac_key.x,
                );
                rs[9] = output_mac_key.from_bits_k();

                Some(table)
            } else {
                assert_eq!(rs[9], F::ZERO);
                assert_eq!(-(rs[1] + rs[6] + rs[7]), F::ZERO);
                None
            };
            [table_0, table_1]
        };

        rs[0] = -rs[10] - rs[9];

        let mut c_primes = [None; 6];
        cs.iter()
            .zip(rs.iter())
            .enumerate()
            .for_each(|(idx, (c, r))| {
                let val = if c_prime_indices.contains(&idx) {
                    Some(*c + r)
                } else {
                    assert_eq!(*c + *r, F::ZERO);
                    None
                };
                c_primes[idx] = val;
            });

        (
            Self {
                c_primes,
                x_adaptors,
                y_adaptors,
            },
            IdentityMacKey,
        )
    }

    fn evaluate<P: PRG>(
        &self,
        parent_nonce: &mut Nonce,
        prg_unseeded: &P::Unseeded,
        input: &Input<F>,
        input_mac: &InputMac<F>,
    ) -> (F, IdentityMac<F>) {
        let mut nonce = parent_nonce.fork();

        let [
            c_prime_0,
            c_prime_1,
            c_prime_2,
            c_prime_3,
            c_prime_4,
            c_prime_5,
        ] = self.c_primes;

        let evaluate_x = |table: &digit_adaptor::Table<F>, nonce: &mut Nonce| {
            let digit_mac = table
                .evaluate::<P>(nonce, prg_unseeded, &input.x, &input_mac.x)
                .1;
            let field_elems = digit_mac.iter().map(|mac| **mac);
            F::from_bits(field_elems)
        };

        let evaluate_y = |table: &digit_adaptor::Table<F>, nonce: &mut Nonce| {
            let digit_mac = table
                .evaluate::<P>(nonce, prg_unseeded, &input.y, &input_mac.y)
                .1;
            let field_elems = digit_mac.iter().map(|mac| **mac);
            F::from_bits(field_elems)
        };

        let c_prime_6 = self.y_adaptors[0]
            .as_ref()
            .map(|table| evaluate_y(table, &mut nonce));

        let c_prime_8 = self.y_adaptors[1]
            .as_ref()
            .map(|table| evaluate_y(table, &mut nonce));

        let c_prime_10 = self.y_adaptors[2]
            .as_ref()
            .map(|table| evaluate_y(table, &mut nonce));

        let c_prime_7 = self.x_adaptors[0]
            .as_ref()
            .map(|table| evaluate_x(table, &mut nonce));

        let c_prime_9 = self.x_adaptors[1]
            .as_ref()
            .map(|table| evaluate_x(table, &mut nonce));

        let c_primes = [
            c_prime_0, c_prime_1, c_prime_2, c_prime_3, c_prime_4, c_prime_5, c_prime_6, c_prime_7,
            c_prime_8, c_prime_9, c_prime_10,
        ];

        let coord: F = Self::coord(input, &c_primes);

        (coord, IdentityMac(PhantomData))
    }
}

#[cfg(test)]
fn test_random_biquadratic<F>()
where
    F: Field<BasePrimeField = Fq> + NumBits + FromBits<bool, F> + FromBits<F, F> + ToBits + Random,
{
    use ark_std::rand::SeedableRng;

    use crate::mac::{Mac, MacKey};
    use crate::utils::prg::aes_ctr::clear_set;
    use aes::Aes128;

    const SEED: u64 = 0;

    let mut rng = StdRng::seed_from_u64(SEED);
    let unseeded_prg = Aes128::init(&[]);

    let mut c_s = [F::default(); 6];
    (0..6).for_each(|i| c_s[i] = F::rand(&mut rng));
    let cs_flags = [true; 6];

    let input_mac_key = InputMacKey::random_key(&mut rng);

    let garbling_table = Table::garble::<Aes128>(
        &mut Nonce::new(),
        &unseeded_prg,
        &mut rng,
        Some(&(c_s, cs_flags)),
        &input_mac_key,
    )
    .0;
    clear_set();

    let input = <InputMac<F> as Mac<Input<F>>>::random_input(&mut rng);
    let input_mac = input_mac_key.encode(&input);

    let output = Table::coord_plaintext(&input, &c_s);

    let evaluated_output = garbling_table
        .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
        .0;
    clear_set();

    assert_eq!(output, evaluated_output);
}

#[test]
fn test_random_biquadratic_fq() {
    use ark_bn254::Fq;
    test_random_biquadratic::<Fq>();
}

#[test]
fn test_random_biquadratic_fp2() {
    use ark_bn254::Fq2;
    test_random_biquadratic::<Fq2>();
}
