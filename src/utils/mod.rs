//! Utility methods

/// Bit composition and bitwise operations
pub mod bits;

/// Error types
pub mod error;

/// Group operations
pub mod group;

/// Hash function used in half-gates
pub(crate) mod aes_hash;

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
