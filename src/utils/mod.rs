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

pub(crate) fn unwrap_single<T>(iter: impl Iterator<Item = Option<T>>) -> Option<T> {
    let mut flattened = iter.flatten();

    match (flattened.next(), flattened.next()) {
        (Some(value), None) => Some(value),
        (None, _) => None,
        (Some(_), Some(_)) => None,
    }
}
