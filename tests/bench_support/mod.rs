//! In-repo support for consuming the cross-language filter benchmark standard.
//!
//! The authoritative standard lives in the sibling repository `filter-bench`
//! (`data/`: cases, golden anchors, harness policy, gate definitions). mfs is a
//! **published** crate whose CI checks out only mfs, so a
//! `path = "../filter-bench/rust"` dev-dependency would break `cargo test` for
//! every solo clone, crates.io consumer and downstream CI. The crate switch only
//! becomes viable once `filter-bench` has a fetchable remote.
//!
//! Until then the standard's `data/` tree is **vendored** into `bench-data/` by
//! the standard's own tool (byte-for-byte copy, plus `MANIFEST.sha256`):
//!
//! ```text
//! python filter-bench/tools/vendor.py --into mfs/bench-data
//! ```
//!
//! `tests/benchmark_standard.rs` guards it in two tiers: byte-for-byte against
//! the authoritative source when the sibling directory exists, otherwise against
//! `MANIFEST.sha256`.
//!
//! The public shapes here deliberately mirror the frozen `filterbench` Rust crate
//! API, so migrating to the crate later is a change of `use` statements only.

#![allow(dead_code)]

pub mod data;
pub mod gates;
pub mod sha256;
