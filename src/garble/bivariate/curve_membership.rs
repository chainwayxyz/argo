//! This table takes in (P.x, P.y) affine coordinates and obliviously evaluates
//! the polynomial
//!                               (B + x^3 - y^2) s
//! on them. This evaluates to 0 when B + x^3 - y^2 = 0 i.e. when (x,y) is a
//! valid curve point, and `c * s` otherwise, where B + x^3 - y^2 = c != 0.

use ark_ec::{AdditiveGroup, short_weierstrass::SWCurveConfig};
use ark_ff::{Field, Zero};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{UniformRand, rand::rngs::StdRng};
use std::marker::PhantomData;

use crate::{
    garble::{
        Circuit, Nonce,
        bivariate::{Input, InputMac, InputMacKey},
        digit_adaptor,
    },
    mac::{IdentityMac, IdentityMacKey, Random},
    utils::{
        bits::{FromBits, NumBits, ToBits},
        group::FromLeBytesModOrder,
        prg::PRG,
    },
};

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub(crate) struct Table<C>
where
    C: SWCurveConfig,
    C::BaseField: Field + NumBits + FromBits<bool, C::BaseField> + Random,
    <C::BaseField as Field>::BasePrimeField: FromLeBytesModOrder,
{
    // c'0 + c'1 X^3 + c'1 Y^2 + c'2(X)X^2 + c'3(Y)Y + c'4(X)X + c'5(Y) + c'6(X)
    c_primes: [C::BaseField; 3],
    // c'_2, c'_4, c'_6
    x_adaptors: [digit_adaptor::Table<C::BaseField>; 3],
    // c'_3, c'_5
    y_adaptors: [digit_adaptor::Table<C::BaseField>; 2],
}

impl<C> Table<C>
where
    C: SWCurveConfig,
    C::BaseField: Field + NumBits + FromBits<bool, C::BaseField> + Random,
    <C::BaseField as Field>::BasePrimeField: FromLeBytesModOrder,
{
    // c'0 + c'1 X^3 + c'2 Y^2 + c'3(X)X^2 + c'4(Y)Y + c'5(X)X + c'6(Y) + c'7(X)
    fn coord(input: &Input<C::BaseField>, c_primes: &[C::BaseField; 8]) -> C::BaseField {
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
        ] = c_primes;

        let x_sq = x.square();
        *c_prime_0
            + *c_prime_1 * x_sq * x
            + *c_prime_2 * y.square()
            + *c_prime_3 * x_sq
            + *c_prime_4 * y
            + *c_prime_5 * x
            + *c_prime_6
            + *c_prime_7
    }

    // B + X^3 -Y^2
    fn coord_plaintext(input: &Input<C::BaseField>) -> C::BaseField {
        let Input { x, y } = input;

        C::COEFF_B + x.pow([3]) - y.square()
    }
}

impl<C>
    Circuit<Input<C::BaseField>, InputMac<C::BaseField>, C::BaseField, IdentityMac<C::BaseField>>
    for Table<C>
where
    C: SWCurveConfig,
    C::BaseField: Field
        + NumBits
        + ToBits
        + FromBits<C::BaseField, C::BaseField>
        + FromBits<bool, C::BaseField>
        + Random,
    <C::BaseField as Field>::BasePrimeField: FromLeBytesModOrder,
{
    // s
    type Aux = C::BaseField;

    fn garble<P: PRG>(
        parent_nonce: &mut Nonce,
        prg_unseeded: &P::Unseeded,
        rng: &mut StdRng,
        aux: Option<&Self::Aux>,
        input_mac_key: &InputMacKey<C::BaseField>,
    ) -> (Self, IdentityMacKey) {
        let mut nonce = parent_nonce.fork();

        // Bs + s x^3 + (-s) y^2
        // c'0 + c'1 X^3 + c'2 Y^2 + c'3(X)X^2 + c'4(Y)Y + c'5(X)X + c'6(Y) + c'7(X)
        let s = aux.unwrap();

        let mut rs = [C::BaseField::ZERO; 8];
        rs[1] = C::BaseField::rand(rng);
        rs[2] = C::BaseField::rand(rng);

        // Construct an adaptor table from input_mac_key.y, to
        // two DigitMacKeys using:
        // - { phi: -r_2, k: r_4 }
        // - { phi: -r_4, k: r_6 }
        let y_adaptors = {
            // c'_4
            let table_0 = {
                let (table, output_mac_key) = digit_adaptor::Table::garble::<P>(
                    &mut nonce,
                    prg_unseeded,
                    rng,
                    Some(&-rs[2]),
                    &input_mac_key.y,
                );
                rs[4] = output_mac_key.from_bits_k();

                table
            };

            // c'_6
            let table_1 = {
                let (table, output_mac_key) = digit_adaptor::Table::garble::<P>(
                    &mut nonce,
                    prg_unseeded,
                    rng,
                    Some(&-rs[4]),
                    &input_mac_key.y,
                );

                rs[6] = output_mac_key.from_bits_k();

                table
            };

            [table_0, table_1]
        };

        // Construct an adaptor table from input_mac_key.x, to
        // a triple of DigitMacKey's using:
        // - { phi: -r_1, k: r_3 }
        // - { phi: -r_3, k: r_5 }
        // - { phi: -r_5, k: r_7 }
        let x_adaptors = {
            // c'_3
            let table_0 = {
                let (table, output_mac_key) = digit_adaptor::Table::garble::<P>(
                    &mut nonce,
                    prg_unseeded,
                    rng,
                    Some(&-rs[1]),
                    &input_mac_key.x,
                );

                rs[3] = output_mac_key.from_bits_k();

                table
            };

            // c'_5
            let table_1 = {
                let (table, output_mac_key) = digit_adaptor::Table::garble::<P>(
                    &mut nonce,
                    prg_unseeded,
                    rng,
                    Some(&-rs[3]),
                    &input_mac_key.x,
                );
                rs[5] = output_mac_key.from_bits_k();

                table
            };

            // c'_7
            let table_2 = {
                let (table, output_mac_key) = digit_adaptor::Table::garble::<P>(
                    &mut nonce,
                    prg_unseeded,
                    rng,
                    Some(&-rs[5]),
                    &input_mac_key.x,
                );
                rs[7] = output_mac_key.from_bits_k();

                table
            };

            [table_0, table_1, table_2]
        };

        rs[0] = -rs[6] - rs[7];

        let c_primes = {
            // c'0 = Bs + r[0]
            let c_prime_0 = C::COEFF_B * s + rs[0];

            // c'1 = s + r[1]
            let c_prime_1 = *s + rs[1];

            // c'2 = -s + r[2]
            let c_prime_2 = -*s + rs[2];

            [c_prime_0, c_prime_1, c_prime_2]
        };

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
        input: &Input<C::BaseField>,
        input_mac: &InputMac<C::BaseField>,
    ) -> (C::BaseField, IdentityMac<C::BaseField>) {
        let mut nonce = parent_nonce.fork();

        let [c_prime_0, c_prime_1, c_prime_2] = self.c_primes;

        let evaluate_x = |table: &digit_adaptor::Table<C::BaseField>, nonce: &mut Nonce| {
            let digit_mac = table
                .evaluate::<P>(nonce, prg_unseeded, &input.x, &input_mac.x)
                .1;
            let field_elems = digit_mac.iter().map(|mac| **mac);
            C::BaseField::from_bits(field_elems)
        };

        let evaluate_y = |table: &digit_adaptor::Table<C::BaseField>, nonce: &mut Nonce| {
            let digit_mac = table
                .evaluate::<P>(nonce, prg_unseeded, &input.y, &input_mac.y)
                .1;
            let field_elems = digit_mac.iter().map(|mac| **mac);
            C::BaseField::from_bits(field_elems)
        };

        let c_prime_4 = evaluate_y(&self.y_adaptors[0], &mut nonce);

        let c_prime_6 = evaluate_y(&self.y_adaptors[1], &mut nonce);

        let c_prime_3 = evaluate_x(&self.x_adaptors[0], &mut nonce);

        let c_prime_5 = evaluate_x(&self.x_adaptors[1], &mut nonce);

        let c_prime_7 = evaluate_x(&self.x_adaptors[2], &mut nonce);

        let c_primes = [
            c_prime_0, c_prime_1, c_prime_2, c_prime_3, c_prime_4, c_prime_5, c_prime_6, c_prime_7,
        ];

        let coord: C::BaseField = Self::coord(input, &c_primes);
        let curve_eqn = Self::coord_plaintext(input);

        let secret = if curve_eqn.is_zero() {
            assert!(coord.is_zero());
            C::BaseField::ZERO
        } else {
            coord / curve_eqn
        };

        (secret, IdentityMac(PhantomData))
    }
}

#[cfg(test)]
fn test_curve_membership<C>(input: Input<C::BaseField>, rng: &mut StdRng)
where
    C: SWCurveConfig,
    C::BaseField: Field
        + NumBits
        + ToBits
        + FromBits<C::BaseField, C::BaseField>
        + FromBits<bool, C::BaseField>
        + Random,
    <C::BaseField as Field>::BasePrimeField: FromLeBytesModOrder,
{
    use crate::{mac::MacKey, utils::prg::aes_ctr::clear_set};
    use aes::Aes128;

    let s = C::BaseField::rand(rng);
    let unseeded_prg = Aes128::init(&[]);

    let input_mac_key = InputMacKey::random_key(rng);

    let garbling_table = Table::<C>::garble::<Aes128>(
        &mut Nonce::new(),
        &unseeded_prg,
        rng,
        Some(&s),
        &input_mac_key,
    )
    .0;
    clear_set();

    let input_mac = input_mac_key.encode(&input);

    let curve_eqn = Table::<C>::coord_plaintext(&input);
    let evaluated_output = garbling_table
        .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
        .0;
    clear_set();

    if curve_eqn == C::BaseField::ZERO {
        use ark_ff::Zero;
        assert!(evaluated_output.is_zero());
    } else {
        assert_eq!(s, evaluated_output);
    }
}

#[test]
fn test_valid_curve_points() {
    use crate::mac::Mac;
    use ark_bn254::{Fq, Fq2};
    use ark_std::rand::SeedableRng;

    const SEED: u64 = 0;

    let mut rng = StdRng::seed_from_u64(SEED);

    let g1_input = InputMac::<Fq>::random_input(&mut rng);
    test_curve_membership::<ark_bn254::g1::Config>(g1_input, &mut rng);

    let g2_input = InputMac::<Fq2>::random_input(&mut rng);
    test_curve_membership::<ark_bn254::g2::Config>(g2_input, &mut rng);
}

#[test]
fn test_invalid_curve_points() {
    use ark_bn254::{Fq, Fq2};
    use ark_std::rand::SeedableRng;

    const SEED: u64 = 0;

    let mut rng = StdRng::seed_from_u64(SEED);

    let g1_input = Input {
        x: Fq::rand(&mut rng),
        y: Fq::rand(&mut rng),
    };
    test_curve_membership::<ark_bn254::g1::Config>(g1_input, &mut rng);

    let g2_input = Input {
        x: Fq2::rand(&mut rng),
        y: Fq2::rand(&mut rng),
    };
    test_curve_membership::<ark_bn254::g2::Config>(g2_input, &mut rng);
}
