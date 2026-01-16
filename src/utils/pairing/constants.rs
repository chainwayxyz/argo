use ark_bn254::{Fq, Fq2};
use ark_ff::{BigInteger, Field, PrimeField};
use lazy_static::lazy_static;
use num_bigint::BigUint;

/// 6x+2 for in NAF form
/// x = 4965661367192848881 for Bn254
pub const SIX_X_PLUS_2_NAF: [i8; 65] = [
    0, 0, 0, 1, 0, 1, 0, -1, 0, 0, 1, -1, 0, 0, 1, 0, 0, 1, 1, 0, -1, 0, 0, 1, 0, -1, 0, 0, 0, 0,
    1, 1, 1, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 0, -1, 0,
    0, 1, 0, 1, 1,
];

lazy_static! {
    pub static ref U_PLUS_NINE: Fq2 = Fq2::new(Fq::from(9), Fq::ONE);

    // Fq2(u + 9)**(((q^1) - 1) / 2)
    pub static ref XI_TO_Q_MINUS_1_OVER_2: Fq2 = U_PLUS_NINE.pow(Fq::MODULUS_MINUS_ONE_DIV_TWO);

    pub static ref FROBENIUS_COEFF_FQ6_C1: [Fq2; 6] = {
        let q = BigUint::from_bytes_le(&Fq::MODULUS.to_bytes_le());
        let mut q_pow = q.clone();

        // (q_pow - 1) / 3
        let minus_one_div_three = |q_pow: &BigUint| {
            let one = BigUint::from(1u8);
            let three = BigUint::from(3u8);
            ((q_pow.clone() - one) / three).to_u64_digits()
        };

        [
        // Fq2(u + 9)**(((q^0) - 1) / 3)
        Fq2::ONE,
        // Fq2(u + 9)**(((q^1) - 1) / 3)
        {
            U_PLUS_NINE.pow(minus_one_div_three(&q_pow))
        },
        // Fq2(u + 9)**(((q^2) - 1) / 3)
        {
            q_pow *= q.clone();
            U_PLUS_NINE.pow(minus_one_div_three(&q_pow))
        },
        // Fq2(u + 9)**(((q^3) - 1) / 3)
        {
            q_pow *= q.clone();
            U_PLUS_NINE.pow(minus_one_div_three(&q_pow))
        },
        // Fq2(u + 9)**(((q^4) - 1) / 3)
        {
            q_pow *= q.clone();
            U_PLUS_NINE.pow(minus_one_div_three(&q_pow))
        },
        // Fq2(u + 9)**(((q^5) - 1) / 3)
        {
            q_pow *= q.clone();
            U_PLUS_NINE.pow(minus_one_div_three(&q_pow))
        },
    ]};
}
