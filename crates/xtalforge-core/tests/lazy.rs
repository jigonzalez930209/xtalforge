//! Integration tests for the lazy symmetry expansion (Task 15.2.4).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, LazyExpansion, Site, StructureBuilder, SymOp};

fn structure() -> xtalforge_core::Structure {
    StructureBuilder::new("lazy", Lattice::cubic(6.0).unwrap())
        .hall("P -1")
        .operation(SymOp::identity())
        .operation(SymOp::parse("-x,-y,-z").unwrap())
        .site(Site::new("C1", "C", [0.1, 0.2, 0.3]).unwrap())
        .site(
            Site::new("O1", "O", [0.4, 0.5, 0.6])
                .unwrap()
                .with_occupancy(0.5)
                .unwrap()
                .with_charge(-1.0)
                .unwrap(),
        )
        .site(Site::new("C2", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn counts_match_the_materialised_cell() {
    let structure = structure();
    let p1 = structure.p1().unwrap();
    let lazy = LazyExpansion::new(&structure, 0.05).unwrap();
    assert_eq!(lazy.n_sites(), 3);
    assert_eq!(lazy.n_atoms().unwrap(), p1.n_atoms());
    let carbon = p1.atoms.iter().filter(|atom| atom.element == "C").count();
    assert_eq!(lazy.composition().unwrap().get("C"), Some(&carbon));
    assert_eq!(lazy.composition().unwrap().get("O"), Some(&2));
    let expected: f64 = p1.atoms.iter().map(|atom| atom.occupancy).sum();
    assert!((lazy.occupancy_sum().unwrap() - expected).abs() < 1e-12);
    assert!((lazy.charge_sum().unwrap() - structure.total_charge()).abs() < 1e-12);
    assert!((lazy.charge_sum().unwrap() + 1.0).abs() < 1e-12);
}

#[test]
fn orbits_report_multiplicity_and_disorder() {
    let structure = structure();
    let lazy = LazyExpansion::new(&structure, 0.05).unwrap();
    let orbits = lazy.orbits().unwrap();
    assert_eq!(orbits.len(), 3);
    assert_eq!(orbits[0].multiplicity, 2);
    assert_eq!(orbits[1].multiplicity, 2);
    assert_eq!(orbits[1].occupancy, 0.5);
    assert_eq!(
        orbits[2].multiplicity, 1,
        "the special position is deduplicated"
    );
    assert_eq!(orbits[2].site_idx, 2);
    assert_eq!(orbits[0].element, "C");
}

#[test]
fn single_site_images_and_the_callback_agree() {
    let structure = structure();
    let p1 = structure.p1().unwrap();
    let lazy = LazyExpansion::new(&structure, 0.05).unwrap();
    assert_eq!(lazy.site_images(0).unwrap().len(), 2);
    let mut fracs = Vec::new();
    assert!(
        lazy.for_each_atom(|atom| {
            fracs.push(atom.frac);
            true
        })
        .unwrap()
    );
    assert_eq!(fracs.len(), p1.n_atoms());
    let mut expected: Vec<u64> = p1
        .atoms
        .iter()
        .flat_map(|atom| atom.frac.iter().map(|value| value.to_bits()))
        .collect();
    let mut actual: Vec<u64> = fracs
        .iter()
        .flat_map(|frac| frac.iter().map(|value| value.to_bits()))
        .collect();
    expected.sort_unstable();
    actual.sort_unstable();
    assert_eq!(actual, expected);
    let mut visited = 0;
    assert!(
        !lazy
            .for_each_atom(|_| {
                visited += 1;
                false
            })
            .unwrap()
    );
    assert_eq!(visited, 1);
}

#[test]
fn a_p1_structure_uses_the_identity_only() {
    let p1 = StructureBuilder::new("p1", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [0.25, 0.25, 0.25]).unwrap())
        .site(Site::new("C2", "C", [0.75, 0.75, 0.75]).unwrap())
        .build()
        .unwrap();
    let lazy = LazyExpansion::new(&p1, 0.05).unwrap();
    assert_eq!(lazy.n_atoms().unwrap(), 2);
    assert_eq!(lazy.composition().unwrap().get("C"), Some(&2));
}

#[test]
fn invalid_arguments_are_reported() {
    let structure = structure();
    assert!(LazyExpansion::new(&structure, -1.0).is_err());
    assert!(LazyExpansion::new(&structure, f64::INFINITY).is_err());
    let lazy = LazyExpansion::new(&structure, 0.05).unwrap();
    assert!(lazy.orbit(99).is_err());
    assert!(lazy.site_images(99).is_err());
}
