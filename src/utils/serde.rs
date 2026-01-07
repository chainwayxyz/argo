use ark_bn254::{Fq, G2Affine, g2};
use ark_ec::short_weierstrass::Affine;
use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate,
};

use crate::{
    garble::{
        g2_subgroup_check,
        pairing::{self, KAPPA_B},
    },
    mac::{IteratedMac, Mac},
};
use std::fmt::Debug;

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

impl CanonicalSerialize for g2_subgroup_check::Table {
    fn serialize_with_mode<W: std::io::Write>(
        &self,
        mut writer: W,
        compress: Compress,
    ) -> Result<(), SerializationError> {
        self.enc_s.serialize_with_mode(&mut writer, compress)?;
        self.r_ks.serialize_with_mode(&mut writer, compress)?;

        Ok(())
    }

    fn serialized_size(&self, compress: Compress) -> usize {
        self.enc_s.serialized_size(compress) + self.r_ks.serialized_size(compress)
    }
}

impl CanonicalDeserialize for g2_subgroup_check::Table {
    fn deserialize_with_mode<R: std::io::Read>(
        mut reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let enc_s = Fq::deserialize_with_mode(&mut reader, compress, validate)?;

        // These are not g2 points, so don't validate that they are in the
        // subgroup
        let r_ks = deserialize_box_vec::<Affine<g2::Config>, { KAPPA_B }, _>(
            &mut reader,
            compress,
            Validate::No,
        )?;

        Ok(Self { enc_s, r_ks })
    }
}

impl Valid for g2_subgroup_check::Table {
    fn check(&self) -> Result<(), SerializationError> {
        Ok(())
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
        let X = deserialize_box_vec::<G2Affine, { KAPPA_B }, _>(&mut reader, compress, validate)?;

        Ok(Self(X))
    }
}

impl Valid for pairing::Table {
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        Ok(())
    }
}

/// Serializers for elements that are Ark-Works serializable but not serde
/// serializable.
///
/// Many cryptographic objects (e.g. finite field elements) are foreign types
/// that we cannot apply [tagged] or `#[derive(Deserialize, Serialize)]` to.
/// Instead, use `#[serde(with = "canonical")]` at the point where the object is
/// used inside a struct or enum definition.
///
/// [tagged]: tagged_base64::tagged
pub mod canonical {
    use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
    use ark_std::format;
    use serde::{
        Deserialize, Serialize,
        de::{Deserializer, Error as DeError},
        ser::{Error as SerError, Serializer},
    };
    use tagged_base64::TaggedBase64;

    pub fn serialize<S: Serializer, T: CanonicalSerialize>(
        elem: &T,
        serializer: S,
    ) -> Result<S::Ok, S::Error> {
        let mut bytes = Vec::new();
        T::serialize_compressed(elem, &mut bytes)
            .map_err(|e| S::Error::custom(format!("{e:?}")))?;
        Serialize::serialize(&TaggedBase64::new("FIELD", &bytes).unwrap(), serializer)
    }

    pub fn deserialize<'de, D: Deserializer<'de>, T: CanonicalDeserialize>(
        deserializer: D,
    ) -> Result<T, D::Error> {
        let tb64 = <TaggedBase64 as Deserialize>::deserialize(deserializer)?;
        if tb64.tag() == "FIELD" {
            T::deserialize_compressed_unchecked(tb64.as_ref())
                .map_err(|e| D::Error::custom(format!("{e:?}")))
        } else {
            Err(D::Error::custom(format!(
                "incorrect tag (expected FIELD, got {})",
                tb64.tag()
            )))
        }
    }
}
