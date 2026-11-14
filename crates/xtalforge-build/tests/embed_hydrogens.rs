//! Integration tests for hydrogen inference from SMILES bond orders.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_build::{EmbedOptions, embed_smiles, embed_smiles_with_hydrogens, parse_smiles};
use xtalforge_core::Structure;

fn embed(smiles: &str, infer_hydrogens: bool) -> Structure {
    let graph = parse_smiles(smiles).unwrap();
    let options = EmbedOptions::default();
    if infer_hydrogens {
        embed_smiles_with_hydrogens(&graph, &options).unwrap()
    } else {
        embed_smiles(&graph, &options).unwrap()
    }
}

fn distance(first: [f64; 3], second: [f64; 3]) -> f64 {
    ((first[0] - second[0]).powi(2)
        + (first[1] - second[1]).powi(2)
        + (first[2] - second[2]).powi(2))
    .sqrt()
}

/// Index of the closest heavy atom for every hydrogen site.
fn hydrogen_owners(structure: &Structure) -> Vec<(usize, f64)> {
    let points = structure.cart_sites();
    let heavy: Vec<usize> = structure
        .sites()
        .iter()
        .enumerate()
        .filter_map(|(index, site)| (site.element() != "H").then_some(index))
        .collect();
    let mut owners = Vec::new();
    for (index, site) in structure.sites().iter().enumerate() {
        if site.element() != "H" {
            continue;
        }
        let mut best = heavy[0];
        let mut best_distance = f64::INFINITY;
        for &candidate in &heavy {
            let length = distance(points[index], points[candidate]);
            if length < best_distance {
                best_distance = length;
                best = candidate;
            }
        }
        owners.push((best, best_distance));
    }
    owners
}

#[test]
fn benzene_gets_one_hydrogen_per_carbon() {
    let structure = embed("c1ccccc1", true);
    assert_eq!(structure.sites().len(), 12);
    assert_eq!(structure.formula(), "C6 H6");
    let owners = hydrogen_owners(&structure);
    assert_eq!(owners.len(), 6);
    let mut per_carbon = [0_usize; 6];
    for (owner, length) in owners {
        assert!((length - 1.09).abs() < 1e-6, "C-H length {length}");
        per_carbon[owner] += 1;
    }
    assert_eq!(per_carbon, [1; 6]);
}

#[test]
fn pyridine_keeps_the_nitrogen_bare() {
    let structure = embed("c1ccncc1", true);
    assert_eq!(structure.sites().len(), 11);
    assert_eq!(
        structure
            .sites()
            .iter()
            .filter(|site| site.element() == "N")
            .count(),
        1
    );
    let owners = hydrogen_owners(&structure);
    assert_eq!(owners.len(), 5);
    assert!(
        owners
            .iter()
            .all(|(owner, length)| { *owner != 3 && (length - 1.09).abs() < 1e-6 })
    );
}

#[test]
fn naphthalene_keeps_the_fusion_carbons_bare() {
    let structure = embed("c1ccc2ccccc2c1", true);
    assert_eq!(structure.sites().len(), 18);
    let owners = hydrogen_owners(&structure);
    assert_eq!(owners.len(), 8);
    let mut per_carbon = [0_usize; 10];
    for (owner, length) in owners {
        assert!((length - 1.09).abs() < 1e-6, "C-H length {length}");
        per_carbon[owner] += 1;
    }
    assert_eq!(per_carbon[3], 0);
    assert_eq!(per_carbon[8], 0);
    assert_eq!(per_carbon.iter().sum::<usize>(), 8);
}

#[test]
fn explicit_ammonium_keeps_its_four_hydrogens() {
    let structure = embed("[NH4+]", true);
    assert_eq!(structure.sites().len(), 5);
    let owners = hydrogen_owners(&structure);
    assert_eq!(owners.len(), 4);
    for (owner, length) in owners {
        assert_eq!(owner, 0);
        assert!((length - 1.01).abs() < 1e-6, "N-H length {length}");
    }
}

#[test]
fn inference_is_opt_in() {
    let structure = embed("c1ccccc1", false);
    assert_eq!(structure.sites().len(), 6);
    assert!(structure.sites().iter().all(|site| site.element() == "C"));
}

#[test]
fn impossible_hydrogen_count_is_an_error() {
    let graph = parse_smiles("[O+3]").unwrap();
    let options = EmbedOptions::default();
    assert!(embed_smiles_with_hydrogens(&graph, &options).is_err());
}
