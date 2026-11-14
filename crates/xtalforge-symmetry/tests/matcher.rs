//! Integration tests for the structure matcher and equality classes
//! (Task 9.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_symmetry::{
    MatcherOptions, SpeciesComparator, group_structures, interpolate_matched, match_structures,
};

const A: f64 = 5.430941;

fn conventional_silicon() -> Structure {
    let mut builder = StructureBuilder::new("Si-conventional", Lattice::cubic(A).unwrap());
    let positions = [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ];
    for (index, position) in positions.iter().enumerate() {
        builder = builder.site(Site::new(format!("Si{index}"), "Si", *position).unwrap());
        let shifted = [position[0] + 0.25, position[1] + 0.25, position[2] + 0.25];
        builder = builder.site(Site::new(format!("Si{index}b"), "Si", shifted).unwrap());
    }
    builder.build().unwrap()
}

fn primitive_silicon() -> Structure {
    let lattice = Lattice::new([
        [0.0, A / 2.0, A / 2.0],
        [A / 2.0, 0.0, A / 2.0],
        [A / 2.0, A / 2.0, 0.0],
    ])
    .unwrap();
    StructureBuilder::new("Si-primitive", lattice)
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Si2", "Si", [0.25, 0.25, 0.25]).unwrap())
        .build()
        .unwrap()
}

fn p1_three_atoms() -> Structure {
    StructureBuilder::new(
        "random",
        Lattice::new([[4.1, 0.1, 0.0], [0.2, 4.3, 0.1], [0.1, 0.3, 4.7]]).unwrap(),
    )
    .site(Site::new("Na1", "Na", [0.1, 0.2, 0.3]).unwrap())
    .site(Site::new("Cl1", "Cl", [0.6, 0.15, 0.75]).unwrap())
    .site(Site::new("O1", "O", [0.35, 0.8, 0.45]).unwrap())
    .build()
    .unwrap()
}

#[test]
fn conventional_and_primitive_silicon_match() {
    let options = MatcherOptions::default();
    let outcome =
        match_structures(&conventional_silicon(), &primitive_silicon(), &options).unwrap();
    assert!(outcome.matched, "rms {}", outcome.rms);
    assert!(outcome.rms < 0.1);
    assert_eq!(outcome.n_atoms, 2);
    assert_eq!(outcome.mapping.len(), 2);
}

#[test]
fn a_shifted_permuted_copy_matches() {
    let first = p1_three_atoms();
    let second = StructureBuilder::new("copy", *first.lattice())
        .site(Site::new("O1", "O", [0.35 + 0.2, 0.8 - 0.1, 0.45]).unwrap())
        .site(Site::new("Na1", "Na", [0.1 + 0.2, 0.2 - 0.1, 0.3]).unwrap())
        .site(Site::new("Cl1", "Cl", [0.6 + 0.2, 0.15 - 0.1, 0.75]).unwrap())
        .build()
        .unwrap();
    let options = MatcherOptions {
        ..MatcherOptions::default()
    };
    let outcome = match_structures(&first, &second, &options).unwrap();
    assert!(outcome.matched, "rms {}", outcome.rms);
    assert!(outcome.rms < 0.1);
}

#[test]
fn a_different_species_does_not_match() {
    let first = p1_three_atoms();
    let second = StructureBuilder::new("copy", *first.lattice())
        .site(Site::new("C1", "C", [0.1, 0.2, 0.3]).unwrap())
        .site(Site::new("Cl1", "Cl", [0.6, 0.15, 0.75]).unwrap())
        .site(Site::new("O1", "O", [0.35, 0.8, 0.45]).unwrap())
        .build()
        .unwrap();
    let options = MatcherOptions {
        ..MatcherOptions::default()
    };
    assert!(!match_structures(&first, &second, &options).unwrap().matched);
    let ignored = MatcherOptions {
        comparator: SpeciesComparator::Ignore,
        ..MatcherOptions::default()
    };
    assert!(match_structures(&first, &second, &ignored).unwrap().matched);
}

#[test]
fn a_one_percent_strain_still_matches_as_a_supercell() {
    let first = conventional_silicon();
    let strained = first
        .clone()
        .with_lattice(Lattice::cubic(A * 1.01).unwrap())
        .as_p1()
        .unwrap();
    let options = MatcherOptions {
        ltol: 0.05,
        stol: 0.2,
        ..MatcherOptions::default()
    };
    let outcome = match_structures(&first, &strained, &options).unwrap();
    assert!(outcome.matched, "rms {}", outcome.rms);
}

#[test]
fn grouping_is_deterministic() {
    let first = p1_three_atoms();
    let second = first.origin_shift([0.1, 0.0, 0.0]).unwrap();
    let third = primitive_silicon();
    let options = MatcherOptions {
        ..MatcherOptions::default()
    };
    let groups = group_structures(&[first, second, third], &options).unwrap();
    assert_eq!(groups.len(), 2);
    assert_eq!(groups[0], vec![0, 1]);
    assert_eq!(groups[1], vec![2]);
}

#[test]
fn matched_interpolation_returns_images() {
    let first = p1_three_atoms();
    let second = first.origin_shift([0.1, 0.0, 0.0]).unwrap();
    let options = MatcherOptions {
        ..MatcherOptions::default()
    };
    let frames = interpolate_matched(&first, &second, 3, &options).unwrap();
    assert_eq!(frames.len(), 3);
    let sodium = frames[1]
        .sites()
        .iter()
        .find(|site| site.element() == "Na")
        .unwrap();
    assert!((0.1..=0.2).contains(&sodium.frac()[0]));
    assert_ne!(frames[0].sites(), frames[2].sites());
}
