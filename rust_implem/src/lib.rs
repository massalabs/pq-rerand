//! Post-quantum re-randomizable encryption (BFV/LPR + CRT noise-limb trick).
//!
//! # ⚠️ WARNING: NOT PRODUCTION READY ⚠️
//!
//! This is a research prototype. NOT audited, NOT constant-time,
//! NOT safe against side-channel attacks.

pub mod params;
pub mod ntt;
pub mod poly;
pub mod sampling;
pub mod encoding;
pub mod keygen;
pub mod encrypt;
pub mod rerandomize;
pub mod decrypt;
pub mod serialize;
