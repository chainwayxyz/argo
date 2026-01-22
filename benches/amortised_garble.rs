use argo::{
    garble::groth16_verifier::{GID, InputMacKey, VerifierKey, garble_groth16},
    mac::{MacKey, Random},
};
use ark_bn254::{Bn254, Fq, G1Projective, G2Projective};
use ark_ec::{PrimeGroup, pairing::Pairing};
use ark_ff::UniformRand;
use ark_std::rand::{SeedableRng, rngs::StdRng};
use criterion::{Criterion, criterion_group, criterion_main};
use rayon::{iter::ParallelIterator, slice::ParallelSlice};

fn amortised_garble(c: &mut Criterion) {
    const N: u64 = 400;

    let mut rng = StdRng::seed_from_u64(0);

    let verifier_key = VerifierKey {
        G1: G1Projective::rand(&mut rng),
        H1: G2Projective::rand(&mut rng),
        H2: G2Projective::rand(&mut rng),
        c: Bn254::pairing(G1Projective::rand(&mut rng), G2Projective::generator()).0,
    };

    let mut garble_group = c.benchmark_group("amortised-garble");
    garble_group.sample_size(10);
    garble_group.bench_function("amortised_garble", |b| {
        b.iter(|| {
            let chunk_size = (N as usize).div_ceil(rayon::current_num_threads());

            (0..N)
                .collect::<Vec<_>>()
                .par_chunks(chunk_size)
                .for_each(|seeds| {
                    for &seed in seeds {
                        let mut rng = StdRng::seed_from_u64(seed);
                        let input_mac_key = InputMacKey::random_key(&mut rng);
                        let secret = Fq::rand(&mut rng);
                        let gid = GID::random(&mut rng);

                        garble_groth16(seed, secret, verifier_key, &input_mac_key, gid).unwrap();
                    }
                });
        })
    });
    garble_group.finish();
}

criterion_group!(benches, amortised_garble);
criterion_main!(benches);
