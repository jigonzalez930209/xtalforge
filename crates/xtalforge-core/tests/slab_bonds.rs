//! Integration tests for slab bond carrying and the severed-bond count
//! (Task 14.3.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{
    BondRules, ExplicitBond, Lattice, Site, Structure, StructureBuilder, make_slab,
};

/// A one-dimensional C-C chain: a 1.4 A bond and a one-cell image at 1.0 A.
fn chain() -> Structure {
    StructureBuilder::new("chain", Lattice::orthorhombic(2.0, 2.0, 2.4).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("C2", "C", [0.0, 0.0, 1.4 / 2.4]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn slab_reports_the_bonds_the_surfaces_sever() {
    let structure = chain();
    let graph = structure.bond_graph().unwrap();
    assert_eq!(graph.bonds.len(), 2, "two images of the 0.7 A contact");

    let slab = make_slab(&structure, [0, 0, 1], 2, 5.0, 0.0).unwrap();
    assert_eq!(slab.structure.sites().len(), 4);
    assert_eq!(slab.cut, 2, "one severed bond at each surface");
    assert_eq!(
        slab.structure.explicit_bonds().len(),
        3,
        "the three internal bonds are carried"
    );
    for bond in slab.structure.explicit_bonds() {
        assert!(!bond.remove);
        assert_eq!(bond.image[2], 0);
    }
    let rebuilt = slab.structure.bond_graph().unwrap();
    assert_eq!(rebuilt.bonds.len(), 3);
    assert!(rebuilt.bonds.iter().all(|bond| bond.explicit));
}

#[test]
fn slab_keeps_the_bond_order_the_user_stated() {
    let structure = chain();
    let mut builder = StructureBuilder::new(structure.name(), *structure.lattice());
    for site in structure.sites() {
        builder = builder.site(site.clone());
    }
    let stated = builder
        .explicit_bond(ExplicitBond {
            i: 0,
            j: 1,
            image: [0, 0, 0],
            order: 2.0,
            remove: false,
        })
        .build()
        .unwrap();
    let slab = make_slab(&stated, [0, 0, 1], 1, 5.0, 0.0).unwrap();
    assert_eq!(slab.structure.sites().len(), 2);
    assert_eq!(slab.cut, 2);
    let orders: Vec<f64> = slab
        .structure
        .explicit_bonds()
        .iter()
        .map(|b| b.order)
        .collect();
    assert!(orders.iter().any(|order| (order - 2.0).abs() < 1e-12));
}

#[test]
fn slab_carries_a_suppression_instead_of_perceiving_it_again() {
    let structure = chain();
    let mut builder = StructureBuilder::new(structure.name(), *structure.lattice());
    for site in structure.sites() {
        builder = builder.site(site.clone());
    }
    let suppressed = builder
        .explicit_bond(ExplicitBond {
            i: 0,
            j: 1,
            image: [0, 0, 0],
            order: 1.0,
            remove: true,
        })
        .build()
        .unwrap();
    assert_eq!(suppressed.bond_graph().unwrap().bonds.len(), 1);
    let slab = make_slab(&suppressed, [0, 0, 1], 1, 5.0, 0.0).unwrap();
    let rebuilt = slab.structure.bond_graph().unwrap();
    assert_eq!(
        rebuilt.bonds.len(),
        0,
        "the carried suppression covers the only bond the slab geometry has"
    );
}

#[test]
fn slab_with_a_termination_shift_still_carries_bonds() {
    let structure = chain();
    let slab = make_slab(&structure, [0, 0, 1], 3, 8.0, 0.25).unwrap();
    assert_eq!(slab.structure.sites().len(), 6);
    assert!(slab.cut >= 2);
    assert_eq!(slab.structure.explicit_bonds().len(), 5);
}

#[test]
fn slab_copies_the_bond_rules_metadata() {
    let rules = BondRules::default().with_metal_metal(true);
    let mut builder = StructureBuilder::new("rules", Lattice::orthorhombic(2.0, 2.0, 2.4).unwrap())
        .metadata("bond_rules", rules.to_json());
    for site in chain().sites() {
        builder = builder.site(site.clone());
    }
    let structure = builder.build().unwrap();
    let slab = make_slab(&structure, [0, 0, 1], 1, 5.0, 0.0).unwrap();
    assert_eq!(BondRules::from_structure(&slab.structure).unwrap(), rules);
}
