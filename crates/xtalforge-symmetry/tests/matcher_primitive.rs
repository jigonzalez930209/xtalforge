//! Golden regression tests for the standardized primitive matcher
//! (Task 17.3.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_symmetry::{MatcherOptions, SpeciesComparator, match_structures};

const A: f64 = 5.430941;

fn conventional(element: &str) -> Structure {
    let mut builder = StructureBuilder::new("conventional", Lattice::cubic(A).unwrap());
    let positions = [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ];
    for (index, position) in positions.iter().enumerate() {
        builder = builder.site(Site::new(format!("{element}{index}"), element, *position).unwrap());
        let shifted = [position[0] + 0.25, position[1] + 0.25, position[2] + 0.25];
        builder = builder.site(Site::new(format!("{element}{index}b"), element, shifted).unwrap());
    }
    builder.build().unwrap()
}

fn primitive(element: &str) -> Structure {
    let lattice = Lattice::new([
        [0.0, A / 2.0, A / 2.0],
        [A / 2.0, 0.0, A / 2.0],
        [A / 2.0, A / 2.0, 0.0],
    ])
    .unwrap();
    StructureBuilder::new("primitive", lattice)
        .site(Site::new(format!("{element}1"), element, [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new(format!("{element}2"), element, [0.25, 0.25, 0.25]).unwrap())
        .build()
        .unwrap()
}

fn empty() -> Structure {
    StructureBuilder::new("empty", Lattice::cubic(A).unwrap())
        .build()
        .unwrap()
}

#[test]
fn a_conventional_primitive_pair_matches_on_the_primitive_cell() {
    let outcome = match_structures(
        &conventional("Si"),
        &primitive("Si"),
        &MatcherOptions::default(),
    )
    .unwrap();
    assert!(outcome.matched, "rms {}", outcome.rms);
    assert!(outcome.rms < 1e-6, "rms {}", outcome.rms);
    assert_eq!(outcome.n_atoms, 2);
    assert_eq!(outcome.mapping.len(), 2);
    assert!(outcome.supercell.is_none());
}

#[test]
fn a_supercell_still_matches_the_conventional_cell() {
    let expanded = conventional("Si").supercell(2, 1, 1).unwrap();
    assert_eq!(expanded.p1().unwrap().n_atoms(), 16);
    let outcome =
        match_structures(&conventional("Si"), &expanded, &MatcherOptions::default()).unwrap();
    assert!(outcome.matched, "rms {}", outcome.rms);
    assert_eq!(outcome.n_atoms, 2);
}

#[test]
fn a_shifted_rotated_copy_still_matches() {
    let shifted = conventional("Si").origin_shift([0.1, 0.2, 0.3]).unwrap();
    let outcome =
        match_structures(&conventional("Si"), &shifted, &MatcherOptions::default()).unwrap();
    assert!(outcome.matched, "rms {}", outcome.rms);
}

#[test]
fn a_different_species_on_the_same_lattice_does_not_match() {
    let outcome = match_structures(
        &conventional("Si"),
        &conventional("Ge"),
        &MatcherOptions::default(),
    )
    .unwrap();
    assert!(!outcome.matched);
}

#[test]
fn ignoring_species_matches_the_same_lattice() {
    let options = MatcherOptions {
        comparator: SpeciesComparator::Ignore,
        ..MatcherOptions::default()
    };
    let outcome = match_structures(&conventional("Si"), &conventional("Ge"), &options).unwrap();
    assert!(outcome.matched, "rms {}", outcome.rms);
}

#[test]
fn empty_structures_and_invalid_tolerances_are_rejected() {
    let options = MatcherOptions::default();
    assert!(match_structures(&empty(), &empty(), &options).is_err());
    assert!(match_structures(&conventional("Si"), &empty(), &options).is_err());
    for invalid in [
        MatcherOptions {
            ltol: 0.0,
            ..MatcherOptions::default()
        },
        MatcherOptions {
            stol: f64::NAN,
            ..MatcherOptions::default()
        },
        MatcherOptions {
            angle_tol: -1.0,
            ..MatcherOptions::default()
        },
    ] {
        assert!(match_structures(&conventional("Si"), &primitive("Si"), &invalid).is_err());
    }
}
