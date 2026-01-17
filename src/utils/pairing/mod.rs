//! Adapted from https://github.com/matter-labs/zksync-crypto/blob/main/crates/pairing/src/bn254/pairing_certificate.rs
//! to work with arkworks types

use ark_bn254::{Bn254, Fq2, Fq12, G1Affine, G2Affine};
use ark_ec::pairing::{MillerLoopOutput, Pairing, PairingOutput};
use ark_ff::{AdditiveGroup, Field, One, batch_inversion};

mod constants;
use constants::*;

pub type Line = Vec<(Fq2, Fq2)>;

pub fn multi_pairing(g1_points: &[G1Affine], g2_lines: &[&Line]) -> PairingOutput<Bn254> {
    let miller_loop_f = miller_loop_with_prepared_lines(g1_points, g2_lines);
    Bn254::final_exponentiation(MillerLoopOutput(miller_loop_f)).unwrap()
}

fn double_with_alpha(value: &mut G2Affine, alpha: Fq2) {
    //x_r = alpha^2 - 2 * x
    let mut x_r = alpha;
    x_r.square_in_place();
    x_r -= value.x;
    x_r -= value.x;

    // y_r = alpha * (x - x_r) - y
    let mut y_r = value.x;
    y_r -= x_r;
    y_r *= alpha;
    y_r -= value.y;

    value.x = x_r;
    value.y = y_r;
}

fn add_assign_mixed_with_alpha(value: &mut G2Affine, other: &G2Affine, alpha: Fq2) {
    let x1 = value.x;
    let y1 = value.y;
    let x2 = other.x;

    // x_r = alpha^2 - x1 - x2
    let mut x_r = alpha;
    x_r.square_in_place();
    x_r -= x1;
    x_r -= x2;

    // y_r = alpha * (x1 - x_r) - y1
    let mut y_r = x1;
    y_r -= x_r;
    y_r *= alpha;
    y_r -= y1;

    value.x = x_r;
    value.y = y_r;
}

pub(crate) fn prepare_all_line_functions(q: G2Affine) -> Line {
    fn line_double(t: G2Affine) -> (Fq2, Fq2) {
        let mut alpha = t.x;
        let mut mu = t.y;

        // alpha = 3 * x^2 / 2 * y
        alpha.square_in_place();
        let mut tmp = alpha;
        tmp.double_in_place();
        alpha += tmp;

        let mut tmp = t.y;
        tmp.double_in_place();
        let tmp_inv = tmp.inverse().unwrap();
        alpha *= tmp_inv;

        let mut tmp = t.x;
        tmp *= alpha;
        mu -= tmp;

        (alpha, mu)
    }

    fn line_add(t: G2Affine, p: G2Affine) -> (Fq2, Fq2) {
        let x1 = t.x;
        let y1 = t.y;
        let x2 = p.x;
        let y2 = p.y;

        let mut alpha = y2;
        let mut mu = y1;

        // alpha = (y2 - y1) / (x2 - x1)
        alpha -= y1;
        let mut tmp = x2;
        tmp -= x1;
        let tmp_inv = tmp.inverse().unwrap();
        alpha *= tmp_inv;

        let mut tmp = x1;
        tmp *= alpha;
        mu -= tmp;

        (alpha, mu)
    }

    let mut l = vec![];
    let mut t = q;
    let q_negated = -q;

    for i in (1..SIX_X_PLUS_2_NAF.len()).rev() {
        let (alpha, mu) = line_double(t);
        double_with_alpha(&mut t, alpha);

        l.push((alpha, mu));

        let bit = SIX_X_PLUS_2_NAF[i - 1];

        if bit != 0 {
            let q_t = if bit == 1 { q } else { q_negated };
            let (alpha, mu) = line_add(t, q_t);
            add_assign_mixed_with_alpha(&mut t, &q_t, alpha);
            l.push((alpha, mu));
        }
    }

    // Frobenius map calculations for BN254
    let mut pi_1_q_x = q.x;
    let mut pi_1_q_y = q.y;
    pi_1_q_x.conjugate_in_place();
    pi_1_q_x *= FROBENIUS_COEFF_FQ6_C1[1];
    pi_1_q_y.conjugate_in_place();
    pi_1_q_y *= &*XI_TO_Q_MINUS_1_OVER_2;
    let pi_1_q = G2Affine::new(pi_1_q_x, pi_1_q_y);

    let mut pi_2_q_x = q.x;
    pi_2_q_x *= FROBENIUS_COEFF_FQ6_C1[2];

    let pi_2_q = G2Affine::new(pi_2_q_x, q.y);

    let (alpha, mu) = line_add(t, pi_1_q);
    add_assign_mixed_with_alpha(&mut t, &pi_1_q, alpha);
    l.push((alpha, mu));

    let (alpha, mu) = line_add(t, pi_2_q);
    add_assign_mixed_with_alpha(&mut t, &pi_2_q, alpha);
    l.push((alpha, mu));

    l
}

pub(crate) fn prepare_g1_points(points: &[G1Affine]) -> Vec<G1Affine> {
    let mut y_invs = points.iter().map(|p| p.y).collect::<Vec<_>>();
    batch_inversion(&mut y_invs);

    let ys = y_invs.iter().map(|&y| -y);

    points
        .iter()
        .zip(ys)
        .map(|(p, y)| {
            let x = p.x * y;
            G1Affine::new_unchecked(x, y)
        })
        .collect()
}

fn prepare_g1_point(p: G1Affine) -> G1Affine {
    // we "prepare" p by recomputing (x, y) -> (x', y') where: x' = - p.x / p.y; y' = -1 /p.y
    // it is required to enforce that in c0c3c4, c0 = 1
    let y_new = -p.y.inverse().unwrap();
    let mut x_new = p.x;
    x_new *= y_new;

    G1Affine::new_unchecked(x_new, y_new)
}

fn miller_loop_with_prepared_lines(eval_points: &[G1Affine], lines: &[&Line]) -> Fq12 {
    assert_eq!(eval_points.len(), lines.len());

    fn line_evaluation(alpha: &Fq2, mu: &Fq2, p: &G1Affine) -> (Fq2, Fq2, Fq2) {
        let mut c3 = *alpha;
        c3.c0 *= p.x;
        c3.c1 *= p.x;
        let mut c4 = *mu;
        c4.c0 *= p.y;
        c4.c1 *= p.y;
        let c0 = Fq2::one();

        (c0, c3, c4)
    }

    let mut f = Fq12::one();

    let mut lc = 0;
    for i in (1..SIX_X_PLUS_2_NAF.len()).rev() {
        if i != SIX_X_PLUS_2_NAF.len() - 1 {
            f.square_in_place();
        }
        let x = SIX_X_PLUS_2_NAF[i - 1];
        for (p, l) in eval_points.iter().zip(lines.iter()) {
            let (alpha, mu) = l[lc];
            let (c0, c1, c2) = line_evaluation(&alpha, &mu, p);
            f.mul_by_034(&c0, &c1, &c2);

            if x * x == 1 {
                let (alpha, bias) = l[lc + 1];
                let (c0, c1, c2) = line_evaluation(&alpha, &bias, p);
                f.mul_by_034(&c0, &c1, &c2);
            }
        }

        if x == 0 {
            lc += 1;
        } else {
            lc += 2;
        }
    }

    // Frobenius map part: p - p^2 + p^3 steps
    // this part runs through each eval point and applies
    // three additional line evaluations with special conditions.
    for (p, l) in eval_points.iter().zip(lines.iter()) {
        for k in 0..2 {
            let (alpha, mu) = l[lc + k];

            let (c0, c1, c2) = line_evaluation(&alpha, &mu, p);
            f.mul_by_034(&c0, &c1, &c2);
        }
    }

    lc += 2;

    // Check we consumed all lines
    if !lines.is_empty() {
        assert_eq!(lc, lines[0].len());
    }

    f
}

#[test]
fn test_precomputed_lines() {
    use ark_ff::UniformRand;
    use ark_std::rand::{SeedableRng, rngs::StdRng};

    const SEED: u64 = 0;
    let mut rng = StdRng::seed_from_u64(SEED);

    let g1_points = [G1Affine::rand(&mut rng), G1Affine::rand(&mut rng)];
    let g2_points = [G2Affine::rand(&mut rng), G2Affine::rand(&mut rng)];

    let lhs = Bn254::multi_pairing(&g1_points, &g2_points);

    let g1_points: Vec<_> = g1_points.into_iter().map(prepare_g1_point).collect();

    let g2_lines: Vec<_> = g2_points
        .into_iter()
        .map(prepare_all_line_functions)
        .collect();
    let g2_lines: Vec<&Vec<_>> = g2_lines.iter().collect();

    let rhs = multi_pairing(&g1_points, &g2_lines);

    assert_eq!(lhs, rhs);
}
