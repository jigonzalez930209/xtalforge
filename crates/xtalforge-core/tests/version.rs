//! Integration test for the workspace version accessor.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

#[test]
fn version_matches_the_package_metadata() {
    assert_eq!(xtalforge_core::version(), env!("CARGO_PKG_VERSION"));
    assert!(!xtalforge_core::version().is_empty());
    assert!(xtalforge_core::version().contains('.'));
}
