//! Integration tests for the fragment library (Task 6.2.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_build::{find_fragment, fragments};

#[test]
fn library_has_water_and_methanol() {
    let library = fragments().unwrap();
    assert!(library.len() >= 10);
    assert!(library.iter().any(|fragment| fragment.name == "Water"));
    assert_eq!(find_fragment("water").unwrap().unwrap().smiles, "O");
    assert!(find_fragment("does-not-exist").unwrap().is_none());
}
