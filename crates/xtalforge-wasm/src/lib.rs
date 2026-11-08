//! WebAssembly bindings for `xtalforge`.
//!
//! The wasm-bindgen exports, the generated TypeScript definitions and the
//! size budget are introduced in Task 0.5.5. The crate stays in the workspace
//! so `cargo check --target wasm32-unknown-unknown` covers it from day one.

#![cfg_attr(
    test,
    allow(
        clippy::unwrap_used,
        clippy::expect_used,
        clippy::panic,
        clippy::indexing_slicing
    )
)]
