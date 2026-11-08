//! Integration tests for scalars and atom descriptions (Task 1.14.2/1.14.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder, describe_atom};
use xtalforge_core::{atom_charges, atom_masses, bond_lengths, coordination_numbers};

fn diamond() -> xtalforge_core::Structure {
    StructureBuilder::new("diamond", Lattice::cubic(3.566_986).unwrap())
        .site(
            Site::new("C1", "C", [0.0, 0.0, 0.0])
                .unwrap()
                .with_charge(-1.0)
                .unwrap(),
        )
        .site(Site::new("C2", "C", [0.25, 0.25, 0.25]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn scalars_follow_the_bond_graph() {
    let structure = diamond();
    let lengths = bond_lengths(&structure).unwrap();
    assert_eq!(lengths.len(), 1);
    assert!((lengths[0] - 1.544_550_245_471_719_7).abs() < 1e-9);
    assert_eq!(coordination_numbers(&structure).unwrap(), vec![1, 1]);
    assert_eq!(atom_charges(&structure), vec![Some(-1.0), None]);
    let masses = atom_masses(&structure).unwrap();
    assert_eq!(masses.len(), 2);
    assert!((masses[0] - 12.0107).abs() < 1e-6);
}

#[test]
fn atom_description_contains_the_site() {
    let structure = diamond();
    let text = describe_atom(&structure, 0).unwrap();
    assert!(text.contains("C1"));
    assert!(text.contains("coordination=1"));
    assert!(text.contains("charge=-1.00"));
    assert!(describe_atom(&structure, 9).is_err());
}
