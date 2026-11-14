//! Integration tests for the SMILES graph parser (Task 6.1.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_build::parse_smiles;

#[test]
fn parses_ethanol_and_benzene() {
    let ethanol = parse_smiles("CCO").unwrap();
    assert_eq!(ethanol.atoms.len(), 3);
    assert_eq!(ethanol.bonds.len(), 2);
    assert_eq!(ethanol.atoms[0].element, "C");
    assert_eq!(ethanol.atoms[2].element, "O");
    let benzene = parse_smiles("c1ccccc1").unwrap();
    assert_eq!(benzene.atoms.len(), 6);
    assert_eq!(benzene.bonds.len(), 6);
    assert!(benzene.atoms.iter().all(|atom| atom.aromatic));
    assert!(benzene.bonds.iter().all(|bond| bond.order == 1.5));
}

#[test]
fn aromatic_bonds_default_to_one_and_a_half() {
    let benzene = parse_smiles("c1ccccc1").unwrap();
    assert_eq!(benzene.bonds.len(), 6);
    assert!(benzene.bonds.iter().all(|bond| bond.order == 1.5));
    let pyridine = parse_smiles("c1ccncc1").unwrap();
    assert_eq!(pyridine.atoms[3].element, "N");
    assert_eq!(pyridine.bonds.len(), 6);
    assert!(pyridine.bonds.iter().all(|bond| bond.order == 1.5));
    assert!(pyridine.bonds.iter().any(|bond| {
        let first = &pyridine.atoms[bond.a];
        let second = &pyridine.atoms[bond.b];
        (first.element == "N" && second.element == "C")
            || (first.element == "C" && second.element == "N")
    }));
}

#[test]
fn explicit_bond_symbols_override_the_aromatic_order() {
    let single = parse_smiles("c1ccccc-1").unwrap();
    let closure = single
        .bonds
        .iter()
        .find(|bond| bond.a == 0 && bond.b == 5)
        .unwrap();
    assert_eq!(closure.order, 1.0);
    let double = parse_smiles("c1ccccc=1").unwrap();
    let closure = double
        .bonds
        .iter()
        .find(|bond| bond.a == 0 && bond.b == 5)
        .unwrap();
    assert_eq!(closure.order, 2.0);
}

#[test]
fn branches_charges_and_errors() {
    let graph = parse_smiles("C(C)(N)O").unwrap();
    assert_eq!(graph.atoms.len(), 4);
    assert_eq!(graph.bonds.len(), 3);
    let ammonium = parse_smiles("[NH4+]").unwrap();
    assert_eq!(ammonium.atoms[0].element, "N");
    assert_eq!(ammonium.atoms[0].charge, 1);
    assert!(parse_smiles("").is_err());
    assert!(parse_smiles("C(").is_err());
    assert!(parse_smiles("C1CC").is_err());
    assert!(parse_smiles("Z").is_err());
}
