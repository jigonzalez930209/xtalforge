//! Integration tests for the sketch data model and SMILES writer (Tasks 6.1.5
//! and 6.1.6).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_build::{RDKIT_BEHAVIOUR_INVENTORY, Sketch, SketchAtom, SketchBond, from_smiles};

#[test]
fn ethanol_round_trips_through_smiles() {
    let sketch = from_smiles("CCO").unwrap();
    assert_eq!(sketch.atoms.len(), 3);
    assert_eq!(sketch.bonds.len(), 2);
    let text = sketch.to_smiles().unwrap();
    let again = from_smiles(&text).unwrap();
    assert_eq!(again.atoms.len(), 3);
    assert_eq!(again.bonds.len(), 2);
    assert_eq!(
        again
            .atoms
            .iter()
            .map(|atom| atom.element.as_str())
            .collect::<Vec<_>>(),
        ["C", "C", "O"]
    );
}

#[test]
fn benzene_round_trips_with_aromatic_bond_orders() {
    let sketch = from_smiles("c1ccccc1").unwrap();
    assert_eq!(sketch.atoms.len(), 6);
    assert_eq!(sketch.bonds.len(), 6);
    assert!(sketch.bonds.iter().all(|bond| bond.order == 1.5));
    let text = sketch.to_smiles().unwrap();
    let again = from_smiles(&text).unwrap();
    assert_eq!(again.atoms.len(), 6);
    assert_eq!(again.bonds.len(), 6);
    assert!(again.bonds.iter().all(|bond| bond.order == 1.5));
}

#[test]
fn atom_removal_renumbers_bonds() {
    let mut sketch = Sketch::new();
    for (index, element) in ["C", "N", "C", "O"].iter().enumerate() {
        sketch.add_atom(SketchAtom {
            element: (*element).to_owned(),
            position: [index as f64, 0.0],
            charge: 0,
        });
    }
    sketch
        .add_bond(SketchBond {
            a: 0,
            b: 1,
            order: 1.0,
        })
        .unwrap();
    sketch
        .add_bond(SketchBond {
            a: 1,
            b: 2,
            order: 1.0,
        })
        .unwrap();
    sketch
        .add_bond(SketchBond {
            a: 2,
            b: 3,
            order: 2.0,
        })
        .unwrap();
    sketch.remove_atom(1).unwrap();
    assert_eq!(sketch.atoms.len(), 3);
    assert_eq!(sketch.bonds.len(), 1);
    assert_eq!(sketch.bonds[0].a, 1);
    assert_eq!(sketch.bonds[0].b, 2);
    assert_eq!(sketch.bonds[0].order, 2.0);
    assert_eq!(sketch.to_smiles().unwrap(), "C.C=O");
}

#[test]
fn charges_survive_the_round_trip() {
    let sketch = from_smiles("[NH4+]").unwrap();
    assert_eq!(sketch.atoms.len(), 1);
    assert_eq!(sketch.atoms[0].charge, 1);
    let text = sketch.to_smiles().unwrap();
    assert!(text.contains('+'));
    let again = from_smiles(&text).unwrap();
    assert_eq!(again.atoms[0].element, "N");
    assert_eq!(again.atoms[0].charge, 1);
}

#[test]
fn disconnected_sketches_keep_both_components() {
    let sketch = from_smiles("C.C").unwrap();
    assert_eq!(sketch.atoms.len(), 2);
    assert!(sketch.bonds.is_empty());
    let again = from_smiles(&sketch.to_smiles().unwrap()).unwrap();
    assert_eq!(again.atoms.len(), 2);
}

#[test]
fn invalid_operations_are_rejected() {
    let mut sketch = Sketch::new();
    sketch.add_atom(SketchAtom {
        element: "C".to_owned(),
        position: [0.0, 0.0],
        charge: 0,
    });
    assert!(
        sketch
            .add_bond(SketchBond {
                a: 0,
                b: 1,
                order: 1.0,
            })
            .is_err()
    );
    assert!(
        sketch
            .add_bond(SketchBond {
                a: 0,
                b: 0,
                order: 1.0,
            })
            .is_err()
    );
    assert!(
        sketch
            .add_bond(SketchBond {
                a: 0,
                b: 0,
                order: 1.0,
            })
            .is_err()
    );
    assert!(sketch.remove_atom(7).is_err());
    assert!(Sketch::new().to_smiles().is_err());
}

#[test]
fn the_rdkit_inventory_is_data() {
    assert!(RDKIT_BEHAVIOUR_INVENTORY.len() >= 5);
    assert!(
        RDKIT_BEHAVIOUR_INVENTORY
            .iter()
            .any(|(name, _)| *name == "aromaticity perception")
    );
}
