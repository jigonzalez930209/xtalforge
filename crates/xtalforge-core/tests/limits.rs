//! Integration tests for the size-limit model (Task 1.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::Lattice;
use xtalforge_core::limits::{
    Budget, LimitLevel, Limits, Profile, check_drawn, check_grid, check_hessian, check_phonons,
    check_supercell, count, drawn_atoms, grid_points, multiplier, supercell_atoms,
};

#[test]
fn profile_bounds() {
    let standard = Limits::new(Profile::Standard);
    assert_eq!(standard.bounds(Budget::Supercell), (50_000, Some(200_000)));
    let generous = Limits::new(Profile::Generous);
    assert_eq!(generous.bounds(Budget::Supercell), (100_000, Some(400_000)));
    let warn_only = Limits::new(Profile::WarnOnly);
    assert_eq!(warn_only.bounds(Budget::Supercell), (50_000, None));
}

#[test]
fn counts_are_readable() {
    assert_eq!(count(999), "999");
    assert_eq!(count(52_000), "52,000");
    assert_eq!(count(1_400_000), "1.4 M");
}

#[test]
fn estimates_are_arithmetic() {
    assert_eq!(
        drawn_atoms(1_000, [(0.0, 2.0), (0.0, 2.0), (0.0, 2.0)]),
        8_000
    );
    assert_eq!(
        multiplier([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]).unwrap(),
        8
    );
    assert_eq!(
        supercell_atoms(1_000, [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]).unwrap(),
        8_000
    );
    let lattice = Lattice::cubic(10.0).unwrap();
    assert_eq!(grid_points(&lattice, 0.5).unwrap(), 8_000);
    assert!(multiplier([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]).is_err());
}

#[test]
fn verdict_levels_and_sentences() {
    let limits = Limits::new(Profile::Standard);
    let p = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];
    let ok = check_supercell(&limits, 1_000, p).unwrap();
    assert_eq!(ok.level, LimitLevel::Ok);
    assert!(!ok.warned());

    let warn = check_supercell(&limits, 10_000, p).unwrap();
    assert_eq!(warn.level, LimitLevel::Warn);
    assert!(warn.sentence.contains("80,000"));
    assert!(!warn.refused());

    let refuse = check_supercell(&limits, 30_000, p).unwrap();
    assert_eq!(refuse.level, LimitLevel::Refuse);
    assert!(refuse.refused());
    assert!(refuse.sentence.contains("200,000"));
}

#[test]
fn warn_only_never_refuses() {
    let limits = Limits::new(Profile::WarnOnly);
    let p = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];
    let verdict = check_supercell(&limits, 1_000_000, p).unwrap();
    assert_eq!(verdict.level, LimitLevel::Warn);
    assert!(!verdict.refused());
}

#[test]
fn other_checks_use_their_budgets() {
    let limits = Limits::new(Profile::Standard);
    let lattice = Lattice::cubic(4.0).unwrap();
    assert_eq!(
        check_drawn(&limits, 1_000, [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)]).level,
        LimitLevel::Ok
    );
    assert_eq!(
        check_grid(&limits, &lattice, 0.4).unwrap().level,
        LimitLevel::Ok
    );
    assert_eq!(
        check_phonons(
            &limits,
            100,
            [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]]
        )
        .unwrap()
        .level,
        LimitLevel::Warn
    );
    assert_eq!(check_hessian(&limits, 301).level, LimitLevel::Refuse);
}

#[test]
fn limit_error_has_a_stable_code() {
    let limits = Limits::new(Profile::Standard);
    let verdict = check_hessian(&limits, 301);
    let error = xtalforge_core::Error::LimitExceeded(verdict);
    assert_eq!(error.code(), "limit_exceeded");
    assert!(error.message().contains("Hessian"));
}
