use ark_bn254::G2Affine;
use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate,
};

use crate::{
    garble::pairing::{self, KAPPA_A},
    mac::{IteratedMac, Mac},
};
use std::fmt::Debug;

impl<Value, M, const N: usize> CanonicalSerialize for IteratedMac<Value, M, N>
where
    M: Mac<Value> + CanonicalSerialize,
{
    fn serialize_with_mode<W: std::io::Write>(
        &self,
        mut writer: W,
        compress: Compress,
    ) -> Result<(), SerializationError> {
        (*self.value).serialize_with_mode(&mut writer, compress)?;
        Ok(())
    }

    fn serialized_size(&self, compress: Compress) -> usize {
        (*self.value).serialized_size(compress)
    }
}

impl<Value, M, const N: usize> CanonicalDeserialize for IteratedMac<Value, M, N>
where
    M: Mac<Value> + CanonicalDeserialize + Sync,
    Value: Sync,
{
    fn deserialize_with_mode<R: std::io::Read>(
        mut reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let mut macs = Vec::with_capacity(N);
        for _ in 0..N {
            macs.push(M::deserialize_with_mode(&mut reader, compress, validate)?);
        }
        let macs: Box<[_; N]> = macs
            .into_boxed_slice()
            .try_into()
            .unwrap_or_else(|_| panic!());

        Ok(IteratedMac::from(macs))
    }
}

impl<Value, M, const N: usize> Valid for IteratedMac<Value, M, N>
where
    Value: Sync,
    M: Mac<Value> + Sync,
{
    fn check(&self) -> Result<(), SerializationError> {
        Ok(())
    }
}

impl CanonicalSerialize for pairing::Table {
    #[inline(never)]
    fn serialize_with_mode<W: std::io::Write>(
        &self,
        mut writer: W,
        compress: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        self.0.serialize_with_mode(&mut writer, compress)?;
        Ok(())
    }

    fn serialized_size(&self, compress: ark_serialize::Compress) -> usize {
        self.0.serialized_size(compress)
    }
}

impl CanonicalDeserialize for pairing::Table {
    #[inline(never)]
    fn deserialize_with_mode<R: std::io::Read>(
        mut reader: R,
        compress: ark_serialize::Compress,
        validate: ark_serialize::Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        fn deserialize_box_vec<T, const N: usize, R>(
            mut reader: R,
            compress: Compress,
            validate: Validate,
        ) -> Result<Box<[T; N]>, SerializationError>
        where
            T: CanonicalDeserialize + Debug,
            R: std::io::Read,
        {
            let mut vec = Vec::with_capacity(N);
            for _ in 0..N {
                vec.push(T::deserialize_with_mode(&mut reader, compress, validate)?);
            }

            Ok(vec.into_boxed_slice().try_into().expect("N elements"))
        }

        let X = deserialize_box_vec::<G2Affine, { KAPPA_A }, _>(&mut reader, compress, validate)?;

        Ok(Self(X))
    }
}

impl Valid for pairing::Table {
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        Ok(())
    }
}
