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
