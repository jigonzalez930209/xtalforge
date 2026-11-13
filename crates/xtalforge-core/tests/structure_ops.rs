//! Integration tests for merge-duplicates and inversion (Task 14.3.3 support).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{
    ExplicitBond, Lattice, Site, StructureBuilder, SymOp, TopologyBond, invert, is_centrosymmetric,
    merge_duplicates,
};

fn nearly_coincident() -> xtalforge_core::Structure {
    StructureBuilder::new("duplicates", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("C2", "C", [0.000_025, 0.0, 0.0]).unwrap())
        .site(Site::new("O1", "O", [0.5, 0.5, 0.5]).unwrap())
        .explicit_bond(ExplicitBond {
            i: 1,
            j: 2,
            image: [0, 0, 0],
            order: 2.0,
            remove: false,
        })
        .topology_bond(TopologyBond {
            i: 1,
            j: 2,
            image: [0, 0, 0],
        })
        .build()
        .unwrap()
}

#[test]
fn merge_duplicates_keeps_the_first_site_and_remaps_bonds() {
    let structure = nearly_coincident();
    let merged = merge_duplicates(&structure, 0.05).unwrap();
    assert_eq!(merged.sites().len(), 2);
    assert_eq!(merged.sites()[0].element(), "C");
    assert_eq!(merged.sites()[1].element(), "O");
    assert_eq!(merged.explicit_bonds().len(), 1);
    assert_eq!(
        (merged.explicit_bonds()[0].i, merged.explicit_bonds()[0].j),
        (0, 1)
    );
    assert!((merged.explicit_bonds()[0].order - 2.0).abs() < 1e-12);
    assert_eq!(
        (merged.topology_bonds()[0].i, merged.topology_bonds()[0].j),
        (0, 1)
    );
}

#[test]
fn merge_duplicates_leaves_disorder_and_clean_cells_alone() {
    let disorder = StructureBuilder::new("disorder", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("O1", "O", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    assert_eq!(merge_duplicates(&disorder, 0.05).unwrap().sites().len(), 2);

    let clean = StructureBuilder::new("clean", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("O1", "O", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap();
    assert_eq!(merge_duplicates(&clean, 0.05).unwrap(), clean);
    assert_eq!(
        merge_duplicates(&nearly_coincident(), 1e-9)
            .unwrap()
            .sites()
            .len(),
        3
    );
    assert!(merge_duplicates(&clean, 0.0).is_err());
    assert!(merge_duplicates(&clean, -1.0).is_err());
}

#[test]
fn invert_negates_p1_coordinates_and_keeps_p1() {
    let structure = StructureBuilder::new("p1", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [0.1, 0.2, 0.3]).unwrap())
        .explicit_bond(ExplicitBond {
            i: 0,
            j: 0,
            image: [1, 0, 0],
            order: 1.0,
            remove: false,
        })
        .build()
        .unwrap();
    assert!(!is_centrosymmetric(&structure));
    let inverted = invert(&structure).unwrap();
    assert!(inverted.is_p1());
    assert_eq!(inverted.sites().len(), 1);
    let frac = inverted.sites()[0].frac();
    for (actual, expected) in frac.iter().zip([0.9, 0.8, 0.7].iter()) {
        assert!((actual - expected).abs() < 1e-12);
    }
    assert_eq!(inverted.explicit_bonds()[0].image, [-1, 0, 0]);
}

#[test]
fn invert_leaves_a_centrosymmetric_crystal_unchanged() {
    let structure = StructureBuilder::new("centred", Lattice::cubic(4.0).unwrap())
        .hall("P -1")
        .operation(SymOp::identity())
        .operation(SymOp::parse("-x,-y,-z").unwrap())
        .site(Site::new("C1", "C", [0.1, 0.2, 0.3]).unwrap())
        .build()
        .unwrap();
    assert!(is_centrosymmetric(&structure));
    let inverted = invert(&structure).unwrap();
    assert_eq!(inverted, structure);
}

#[test]
fn invert_drops_a_setting_that_cannot_describe_the_other_hand() {
    let structure = StructureBuilder::new("chiral", Lattice::cubic(4.0).unwrap())
        .hall("P 21")
        .operation(SymOp::identity())
        .site(Site::new("C1", "C", [0.1, 0.2, 0.3]).unwrap())
        .build()
        .unwrap();
    let inverted = invert(&structure).unwrap();
    assert!(inverted.is_p1());
    assert_eq!(inverted.hall(), "P 1");
    assert!(inverted.operations().is_empty());
    let frac = inverted.sites()[0].frac();
    for (actual, expected) in frac.iter().zip([0.9, 0.8, 0.7].iter()) {
        assert!((actual - expected).abs() < 1e-12);
    }
}

#[test]
fn invert_rejects_an_empty_structure() {
    let empty = StructureBuilder::new("empty", Lattice::cubic(4.0).unwrap())
        .build()
        .unwrap();
    assert!(invert(&empty).is_err());
}
