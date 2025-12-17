//! Maps the MAC of an affine curve point (P.x, P.y) to its Jacobian form.
//! Takes in MAC of bits (P.x, P.y); evaluates to iterated MAC of Jacobian
//! coordinates of P, in the form MAC(P) = phi * P + k, where phi is an
//! endoscalar.
//!
//! This is generic over inputs in Fq and Fq2.

use ark_bn254::Fq;
use ark_ec::{AdditiveGroup, AffineRepr, scalar_mul::glv::GLVConfig, short_weierstrass::Affine};
use ark_ff::{Field, batch_inversion};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::rngs::StdRng;

use crate::{
    garble::{
        Circuit, Nonce,
        bivariate::{Input, InputMac, InputMacKey, biquadratic},
    },
    mac::{
        Random,
        endo::{EndoMac, EndoMacKey, EndoScalar, IteratedEndoMac},
    },
    utils::{
        bits::{FromBits, NumBits, ToBits},
        prg::PRG,
    },
};

pub(crate) type Output<Config> = Affine<Config>;
pub(crate) type OutputMac<Config, const KAPPA: usize> = IteratedEndoMac<Config, KAPPA>;
pub(crate) type OutputMacKey<Config, const KAPPA: usize> = [EndoMacKey<Config>; KAPPA];

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub(crate) struct Table<C, const KAPPA: usize>
where
    C: GLVConfig,
    C::BaseField: Field<BasePrimeField = Fq> + NumBits + FromBits<bool, C::BaseField> + Random,
{
    x: Vec<biquadratic::Table<C::BaseField>>,
    y: Vec<biquadratic::Table<C::BaseField>>,
    z: Vec<biquadratic::Table<C::BaseField>>,
}

impl<Config, F, const KAPPA: usize>
    Circuit<Input<F>, InputMac<F>, Output<Config>, OutputMac<Config, KAPPA>>
    for Table<Config, KAPPA>
where
    Config: GLVConfig<BaseField = F> + Clone,
    F: Field<BasePrimeField = Fq> + NumBits + ToBits + FromBits<F, F> + FromBits<bool, F> + Random,
{
    type Aux = OutputMacKey<Config, KAPPA>;

    fn garble<P: PRG>(
        parent_nonce: &mut Nonce,
        prg_unseeded: &P::Unseeded,
        rng: &mut StdRng,
        aux: Option<&Self::Aux>,
        input_mac_key: &InputMacKey<F>,
    ) -> (Self, OutputMacKey<Config, KAPPA>) {
        assert_eq!(Config::COEFF_A, F::ZERO);

        let mut nonce = parent_nonce.fork();
        let output_mac_key = aux.unwrap();

        // c_0, c_1 * a, c_2 * b, c_3 * a * b, c_4 * a^2, c_5 * b^2
        // The `c`s are coeffs of a polynomial to evaluate P -> P + k
        // The endo_coeffs converts it into a polynomial to evaluate P -> phi * P + k
        fn endo_coeffs<Config: GLVConfig>(
            phi: &EndoScalar<Config>,
            cs: [Config::BaseField; 6],
            r: Config::BaseField,
            k: Config::BaseField,
        ) -> [Config::BaseField; 6] {
            // If phi = NegInf, just return k
            match phi {
                EndoScalar::NegInf(_) => [
                    k * r,
                    Config::BaseField::ZERO,
                    Config::BaseField::ZERO,
                    Config::BaseField::ZERO,
                    Config::BaseField::ZERO,
                    Config::BaseField::ZERO,
                ],
                _ => {
                    let phi = EndoScalar::<Config>::base().pow([phi.as_i8() as u64]);
                    let (a, b) = {
                        let phi_sq = phi.square();
                        (phi_sq.square(), phi_sq * phi)
                        // phi * (X, Y) = (a * X, b * Y)
                        // c_0, c_1 * a, c_2 * b, c_3 * a * b, c_4 * a^2, c_5 * b^2
                    };
                    let mut coeffs = [Config::BaseField::ZERO; 6];
                    coeffs
                        .iter_mut()
                        .zip([Config::BaseField::ONE, a, b, a * b, a.square(), b.square()].iter())
                        .zip(cs.iter())
                        .for_each(|((coeff, v), c)| *coeff = *v * c * r);
                    coeffs
                }
            }
            // For x, k = P.x
            // For y, k = P.y
            // For z, k = P.z
        }

        fn x_aux<C, F>(mac_key: &EndoMacKey<C>, a2: F) -> ([F; 6], [bool; 6])
        where
            C: GLVConfig<BaseField = F>,
            F: Field,
        {
            let B = C::COEFF_B;
            let cs_flags = [true, true, true, false, true, false];
            let (x, y) = mac_key.k.xy().unwrap();
            // f(X, Y) = phi * (X, Y) + K = (a X, b Y) + K
            let cs = {
                //  (6, 1)                 # c_0
                let c_0 = B.double();
                //  (xQ_^2, xP_)           # c_1
                let c_1 = x.square();
                //  (-2*yQ_, yP_)          # c_2
                let c_2 = -y.double();
                //  (0, xP_ * yP_)         # c_3 = 0
                let c_3 = F::ZERO;
                //  (xQ_, xP_^2)           # c_4
                let c_4 = x;
                //  (0, yP_^2)             # c_5 = 0
                let c_5 = F::ZERO;

                // scale by random a2 = a^2
                [c_0 * a2, c_1 * a2, c_2 * a2, c_3 * a2, c_4 * a2, c_5 * a2]
            };

            let cs = endo_coeffs(&mac_key.phi, cs, a2, x);

            (cs, cs_flags)
        }

        fn y_aux<C, F>(mac_key: &EndoMacKey<C>, a3: F) -> ([F; 6], [bool; 6])
        where
            C: GLVConfig<BaseField = F>,
            F: Field,
        {
            let B = C::COEFF_B;
            let (x, y) = mac_key.k.xy().unwrap();

            let cs_flags = [true, false, true, true, true, true];
            let cs = {
                //  (9*yQ_, 1)                  # c_0
                let c_0 = F::from(3) * B * y;
                //  (0,   xP_),                 # c_1 = 0
                let c_1 = F::ZERO;
                //  (-(yQ_^2 + 9), yP_)         # c_2
                let c_2 = -(y.square() + F::from(3) * B);
                //  (-3 * xQ_^2, xP_ * yP_)     # c_3
                let c_3 = -F::from(3) * x.square();
                //  (3 * xQ_ * yQ_, xP_^2),     # c_4
                let c_4 = F::from(3) * x * y;
                //  (yQ_^2, yP_),               # c_5
                let c_5 = y;

                // scale by random a3 = a^3
                [c_0 * a3, c_1 * a3, c_2 * a3, c_3 * a3, c_4 * a3, c_5 * a3]
            };
            let cs = endo_coeffs(&mac_key.phi, cs, a3, y);

            (cs, cs_flags)
        }

        fn z_aux<C, F>(mac_key: &EndoMacKey<C>, a: F) -> ([F; 6], [bool; 6])
        where
            C: GLVConfig<BaseField = F>,
            F: Field,
        {
            let (x, _) = mac_key.k.xy().unwrap();
            let cs_flags = [true, true, false, false, false, false];

            let cs = {
                //  (-xQ_,      1)]       # c_0
                let c_0 = -x;
                //  (1, xP_)              # c_1
                let c_1 = F::ONE;
                let c_2 = F::ZERO;
                let c_3 = F::ZERO;
                let c_4 = F::ZERO;
                let c_5 = F::ZERO;

                // scale by random `a`
                [c_0 * a, c_1 * a, c_2 * a, c_3 * a, c_4 * a, c_5 * a]
            };

            let cs = endo_coeffs(&mac_key.phi, cs, a, F::ONE);

            (cs, cs_flags)
        }

        let mut x = Vec::with_capacity(output_mac_key.len());
        let mut y = Vec::with_capacity(output_mac_key.len());
        let mut z = Vec::with_capacity(output_mac_key.len());

        output_mac_key.iter().for_each(|key| {
            // Generate random scalar `a`
            let a = F::rand(rng);

            // Scale z coordinate by `a`
            let z_table = biquadratic::Table::garble::<P>(
                &mut nonce,
                prg_unseeded,
                rng,
                Some(&z_aux(key, a)),
                input_mac_key,
            )
            .0;
            z.push(z_table);

            // Scale x coordinate by a^2
            let a2 = a.square();
            let x_table = biquadratic::Table::garble::<P>(
                &mut nonce,
                prg_unseeded,
                rng,
                Some(&x_aux(key, a2)),
                input_mac_key,
            )
            .0;
            x.push(x_table);

            // Scale y coordinate by a^3
            let a3 = a2 * a;
            let y_table = biquadratic::Table::garble::<P>(
                &mut nonce,
                prg_unseeded,
                rng,
                Some(&y_aux(key, a3)),
                input_mac_key,
            )
            .0;
            y.push(y_table);
        });

        (Self { x, y, z }, output_mac_key.clone())
    }

    fn evaluate<P: PRG>(
        &self,
        parent_nonce: &mut Nonce,
        prg_unseeded: &P::Unseeded,
        input: &Input<Config::BaseField>,
        input_mac: &InputMac<Config::BaseField>,
    ) -> (Affine<Config>, OutputMac<Config, KAPPA>) {
        let mut nonce = parent_nonce.fork();

        let mut xs = Vec::with_capacity(self.x.len());
        let mut ys = Vec::with_capacity(self.y.len());
        let mut zs = Vec::with_capacity(self.z.len());

        self.x
            .iter()
            .zip(self.y.iter())
            .zip(self.z.iter())
            .for_each(|((x_table, y_table), z_table)| {
                let z_table = z_table
                    .evaluate::<P>(&mut nonce, prg_unseeded, input, input_mac)
                    .0;
                zs.push(z_table);

                let x_table = x_table
                    .evaluate::<P>(&mut nonce, prg_unseeded, input, input_mac)
                    .0;
                xs.push(x_table);

                let y_table = y_table
                    .evaluate::<P>(&mut nonce, prg_unseeded, input, input_mac)
                    .0;
                ys.push(y_table);
            });

        let mut zs_inv = zs;
        batch_inversion(&mut zs_inv);

        let mut macs = Vec::with_capacity(KAPPA);
        xs.into_iter()
            .zip(ys)
            .zip(zs_inv)
            .for_each(|((x, y), z_inv)| {
                // Use Jacobian coordinates, so x / z^2, y / z^3
                let z2 = z_inv.square();
                let point = Affine::new_unchecked(x * z2, y * z_inv * z2);
                let mac = EndoMac::from(point);
                macs.push(mac);
            });

        let macs: Box<[EndoMac<Config>; KAPPA]> = macs.try_into().unwrap_or_else(|_| panic!());

        (Affine::new_unchecked(input.x, input.y), macs.into())
    }
}

#[cfg(test)]
fn test_field_mac_to_ec_mac<C>()
where
    C: GLVConfig + Clone,
    C::BaseField: Field<BasePrimeField = Fq>
        + NumBits
        + FromBits<bool, C::BaseField>
        + Random
        + ToBits
        + FromBits<C::BaseField, C::BaseField>,
{
    use crate::{garble::pairing::KAPPA_B, mac::MacKey, utils::prg::aes_ctr::clear_set};
    use aes::Aes128;
    use ark_std::{UniformRand, rand::SeedableRng};

    const SEED: u64 = 0;

    let mut rng = StdRng::seed_from_u64(SEED);
    let unseeded_prg = Aes128::init(&[]);

    let input_mac_key = InputMacKey::random_key(&mut rng);
    let output_mac_key = OutputMacKey::<_, KAPPA_B>::random_key(&mut rng);

    let input_point = Affine::<C>::rand(&mut rng);

    let input = input_point.into();
    let input_mac = input_mac_key.encode(&input);

    let output = input_point;
    let output_mac = output_mac_key.encode(&output);

    let garbling_table = Table::garble::<Aes128>(
        &mut Nonce::new(),
        &unseeded_prg,
        &mut rng,
        Some(&output_mac_key),
        &input_mac_key,
    )
    .0;
    clear_set();

    let evaluated_output = garbling_table
        .evaluate::<Aes128>(&mut Nonce::new(), &unseeded_prg, &input, &input_mac)
        .1;
    clear_set();

    for (expected, evaluated) in output_mac.iter().zip(evaluated_output.iter()) {
        assert_eq!(**expected, **evaluated);
    }
}

#[test]
fn test_field_mac_to_ec_mac_g1() {
    test_field_mac_to_ec_mac::<ark_bn254::g1::Config>();
}

#[test]
fn test_field_mac_to_ec_mac_g2() {
    test_field_mac_to_ec_mac::<ark_bn254::g2::Config>();
}
