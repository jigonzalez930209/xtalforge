//! Integration tests for the `huge` limit preset (Task 15.2.5).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::limits::{Budget, LimitLevel, check_hessian, check_supercell};
use xtalforge_core::{Limits, Profile, limits_bounds, limits_check};

#[test]
fn the_huge_preset_is_named_and_parsed() {
    assert_eq!(Profile::parse("huge"), Profile::Huge);
    assert_eq!(Profile::Huge.name(), "huge");
    assert_eq!(
        limits_bounds("huge", "supercell").unwrap(),
        (2_000_000, Some(20_000_000))
    );
}

#[test]
fn million_atom_supercells_are_inside_the_huge_preset() {
    let huge = Limits::new(Profile::Huge);
    let standard = Limits::new(Profile::Standard);
    let p = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];
    let verdict = check_supercell(&huge, 250_000, p).unwrap();
    assert_eq!(verdict.level, LimitLevel::Ok);
    assert!(!verdict.warned());
    assert!(check_supercell(&standard, 250_000, p).unwrap().refused());

    let warning = check_supercell(&huge, 1_500_000, p).unwrap();
    assert_eq!(warning.level, LimitLevel::Warn);
    let refusal = check_supercell(&huge, 3_000_000, p).unwrap();
    assert!(refusal.refused());
    assert!(refusal.sentence.contains("huge limits preset"));
}

#[test]
fn every_budget_has_a_huge_bound() {
    for name in [
        "drawn",
        "supercell",
        "grid",
        "carbon",
        "nanotube",
        "graphene",
        "phonons",
        "hessian",
    ] {
        let (soft, hard) = limits_bounds("huge", name).unwrap();
        let hard = hard.expect("huge refuses above the hard limit");
        assert!(hard > soft, "{name}");
        assert_eq!(
            limits_check("huge", name, soft).unwrap().level,
            LimitLevel::Ok
        );
        assert_eq!(Budget::parse(name).unwrap().name(), name);
    }
}

#[test]
fn the_hessian_stays_capped_even_in_huge() {
    let huge = Limits::new(Profile::Huge);
    assert_eq!(check_hessian(&huge, 1_000).level, LimitLevel::Ok);
    assert!(check_hessian(&huge, 3_001).refused());
    assert!(huge.bounds(Budget::Hessian).0 >= 1_000);
}
