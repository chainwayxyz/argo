//! Utility methods

/// Bit composition and bitwise operations
pub(crate) mod bits;

/// Group operations
pub(crate) mod group;

/// Pseudo-random generator
pub(crate) mod prg;
pub use prg::PRG;

/// Serialization and deserialization
pub(crate) mod serde;

pub(crate) fn try_into_array<I: Iterator<Item = T>, T, const N: usize>(iterator: I) -> [T; N] {
    iterator
        .collect::<Vec<_>>()
        .try_into()
        .unwrap_or_else(|_| panic!())
}
