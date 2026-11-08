//! Integration tests for the C60 builder (Task 6.7).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_build::fullerene;

#[test]
fn c60_has_sixty_carbons_and_bonds() {
    let ball = fullerene(60).unwrap();
    assert_eq!(ball.sites().len(), 60);
    assert_eq!(ball.formula(), "C60");
    assert!(ball.is_p1());
    let graph = ball.bond_graph().unwrap();
    assert_eq!(graph.bonds.len(), 90);
    assert!(fullerene(70).is_err());
}
