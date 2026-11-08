//! Python bindings for `xtalforge`.
//!
//! The pyo3 module, the `xtalforge` package layout and the type stubs are
//! introduced in Task 0.5.4 and grow with each phase. Keeping this crate in
//! the workspace from the start makes `cargo test --workspace` cover the
//! binding build and keeps the version lockstep automatic.

#![cfg_attr(
    test,
    allow(
        clippy::unwrap_used,
        clippy::expect_used,
        clippy::panic,
        clippy::indexing_slicing
    )
)]
