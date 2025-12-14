use std::marker::PhantomData;

use ark_bn254::{Fq, Fq2};
use ark_ec::scalar_mul::glv::GLVConfig;
use ark_ff::UniformRand;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::{Rng, RngCore, rngs::StdRng};
use derive_more::Deref;

use crate::{garble::Nonce, mac::endo::EndoScalar, utils::prg::PRG};

pub(crate) mod bit;
pub(crate) mod digit;
pub(crate) mod endo;
pub(crate) mod group;

pub trait Mac<Value>: Sized {
    // phi, k
    type Key: MacKey<Value, Self>;

    #[cfg(test)]
    fn random_input(rng: &mut StdRng) -> Value;
}

pub trait MacKey<V, M> {
    // phi * value + k
    fn encode(&self, value: &V) -> M;
    #[cfg(test)]
    fn decode(&self, mac: &M) -> V;
    fn random_key(rng: &mut StdRng) -> Self;
}

pub trait RandomFromPhi<V, M>: MacKey<V, M> {
    type Phi: Random;

    fn random_from_phi(rng: &mut StdRng, phi: &Self::Phi) -> Self;
}

pub trait Encrypt<Message> {
    type Ciphertext: CanonicalDeserialize + CanonicalSerialize;

    fn encrypt<P: PRG>(
        &self,
        prg_unseeded: &P::Unseeded,
        parent_nonce: &mut Nonce,
        message: &Message,
    ) -> Self::Ciphertext;

    fn decrypt<P: PRG>(
        &self,
        prg_unseeded: &P::Unseeded,
        parent_nonce: &mut Nonce,
        ciphertext: &Self::Ciphertext,
    ) -> Message;
}

#[derive(Deref)]
pub(crate) struct IteratedMac<Value, M: Mac<Value>, const N: usize> {
    #[deref]
    pub(crate) value: Box<[M; N]>,
    _marker: PhantomData<Value>,
}

impl<Value, M: Mac<Value>, const N: usize> From<Box<[M; N]>> for IteratedMac<Value, M, N> {
    fn from(value: Box<[M; N]>) -> Self {
        Self {
            value,
            _marker: PhantomData,
        }
    }
}

impl<Value, M, const N: usize> PartialEq for IteratedMac<Value, M, N>
where
    M: PartialEq + Mac<Value>,
{
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

/// Single value, multiple MAC keys
impl<BaseValue, BaseMac, BaseMacKey, const KAPPA: usize> Mac<BaseValue>
    for IteratedMac<BaseValue, BaseMac, KAPPA>
where
    BaseValue: PartialEq + Clone,
    BaseMac: Mac<BaseValue, Key = BaseMacKey> + Clone,
    BaseMacKey: MacKey<BaseValue, BaseMac> + Clone,
{
    type Key = [BaseMacKey; KAPPA];

    #[cfg(test)]
    fn random_input(rng: &mut StdRng) -> BaseValue {
        BaseMac::random_input(rng)
    }
}

impl<BaseValue, BaseMac, BaseMacKey, const KAPPA: usize>
    MacKey<BaseValue, IteratedMac<BaseValue, BaseMac, KAPPA>> for [BaseMacKey; KAPPA]
where
    BaseMacKey: MacKey<BaseValue, BaseMac>,
    BaseMac: Mac<BaseValue> + Clone,
    BaseValue: PartialEq + Clone,
{
    fn encode(&self, value: &BaseValue) -> IteratedMac<BaseValue, BaseMac, KAPPA> {
        let mut macs = Vec::with_capacity(KAPPA);
        self.iter()
            .for_each(|base_key| macs.push(base_key.encode(value)));

        let macs: Box<[BaseMac; KAPPA]> = macs.try_into().unwrap_or_else(|_| panic!());

        IteratedMac::from(macs)
    }

    #[cfg(test)]
    fn decode(&self, mac: &IteratedMac<BaseValue, BaseMac, KAPPA>) -> BaseValue {
        let mut values = self
            .iter()
            .zip(mac.iter())
            .map(|(key, mac)| key.decode(mac));

        let first_value = values.next().unwrap();
        assert!(values.all(|value| value == first_value));

        first_value
    }

    fn random_key(rng: &mut StdRng) -> Self {
        let mut keys = Vec::with_capacity(KAPPA);
        for _ in 0..KAPPA {
            keys.push(BaseMacKey::random_key(rng));
        }
        keys.try_into().unwrap_or_else(|_| panic!())
    }
}

#[derive(Debug, Clone, Copy)]
pub struct IdentityMac<V>(pub PhantomData<V>);

#[derive(Debug, Clone)]
pub struct IdentityMacKey;

impl<V, M> MacKey<V, M> for IdentityMacKey {
    fn encode(&self, _value: &V) -> M {
        unimplemented!()
    }

    #[cfg(test)]
    fn decode(&self, _mac: &M) -> V {
        unimplemented!()
    }

    fn random_key(_rng: &mut StdRng) -> Self {
        unimplemented!()
    }
}

impl<V> Mac<V> for IdentityMac<V> {
    type Key = IdentityMacKey;

    #[cfg(test)]
    fn random_input(_rng: &mut StdRng) -> V {
        unimplemented!()
    }
}

pub trait Random {
    fn random(rng: &mut StdRng) -> Self;
}

impl Random for () {
    fn random(_rng: &mut StdRng) -> Self {}
}

impl Random for Fq {
    fn random(rng: &mut StdRng) -> Self {
        Fq::rand(rng)
    }
}

impl Random for Fq2 {
    fn random(rng: &mut StdRng) -> Self {
        Fq2::rand(rng)
    }
}

impl Random for [u8; 16] {
    fn random(rng: &mut StdRng) -> Self {
        let mut bytes = [0u8; 16];
        rng.fill_bytes(&mut bytes);
        bytes
    }
}

impl<G: GLVConfig> Random for EndoScalar<G> {
    fn random(rng: &mut StdRng) -> EndoScalar<G> {
        let exp: u8 = rng.gen_range(0..6);
        match exp {
            0 => Self::Zero(PhantomData),
            1 => Self::One(PhantomData),
            2 => Self::Two(PhantomData),
            3 => Self::Three(PhantomData),
            4 => Self::Four(PhantomData),
            5 => Self::Five(PhantomData),
            _ => unreachable!(),
        }
    }
}
