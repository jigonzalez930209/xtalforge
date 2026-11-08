//! Integration tests for the Pareto front (Task 7.9.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::pareto_front;

#[test]
fn keeps_only_non_dominated_points() {
    let points = [(1.0, 4.0), (2.0, 2.0), (3.0, 1.0), (4.0, 4.0), (2.0, 3.0)];
    assert_eq!(pareto_front(&points), vec![0, 1, 2]);
    assert_eq!(pareto_front(&[]), Vec::<usize>::new());
    assert_eq!(pareto_front(&[(1.0, 1.0)]), vec![0]);
}
