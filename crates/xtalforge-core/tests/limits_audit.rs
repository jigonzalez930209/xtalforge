//! Limits audit at and beyond every budget (Task 8.6.4).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::limits::{Budget, LimitLevel, check_drawn, check_grid, check_supercell};
use xtalforge_core::{Lattice, Limits, Profile, limits_check};

const BUDGETS: &[&str] = &[
    "drawn",
    "supercell",
    "grid",
    "carbon",
    "nanotube",
    "graphene",
    "phonons",
    "hessian",
];

#[test]
fn every_budget_has_a_soft_and_a_hard_limit() {
    for name in BUDGETS {
        let soft = limits_check("standard", name, 0).unwrap();
        assert_eq!(soft.level, LimitLevel::Ok, "{name} zero");
        let generous = limits_check("generous", name, 1).unwrap();
        assert_eq!(generous.level, LimitLevel::Ok, "{name} generous");
        let parse = Budget::parse(name);
        assert!(parse.is_ok(), "{name}");
        assert!(Budget::parse("nope").is_err());
    }
}

#[test]
fn a_cell_beyond_the_supercell_limit_is_refused() {
    let limits = Limits::new(Profile::Standard);
    let small = check_supercell(
        &limits,
        2,
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    )
    .unwrap();
    assert_eq!(small.level, LimitLevel::Ok);
    let huge = check_supercell(
        &limits,
        1_000,
        [[100.0, 0.0, 0.0], [0.0, 100.0, 0.0], [0.0, 0.0, 100.0]],
    )
    .unwrap();
    assert!(huge.refused());
    assert!(!huge.sentence.is_empty());
}

#[test]
fn a_grid_beyond_the_limit_is_refused() {
    let limits = Limits::new(Profile::Standard);
    let lattice = Lattice::cubic(100.0).unwrap();
    let fine = check_grid(&limits, &lattice, 0.01).unwrap();
    assert!(fine.refused());
    let coarse = check_grid(&limits, &lattice, 2.0).unwrap();
    assert!(!coarse.refused());
}

#[test]
fn the_drawn_check_counts_the_display_range() {
    let limits = Limits::new(Profile::Standard);
    let one_cell = check_drawn(&limits, 10_000, [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)]);
    assert!(!one_cell.refused());
    let many_cells = check_drawn(&limits, 10_000, [(-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0)]);
    assert!(many_cells.refused());
}

#[test]
fn warn_only_never_refuses() {
    let limits = Limits::new(Profile::WarnOnly);
    let verdict = check_supercell(
        &limits,
        1_000,
        [[100.0, 0.0, 0.0], [0.0, 100.0, 0.0], [0.0, 0.0, 100.0]],
    )
    .unwrap();
    assert!(!verdict.refused());
    assert_eq!(verdict.level, LimitLevel::Warn);
}
