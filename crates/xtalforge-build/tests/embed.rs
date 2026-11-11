//! Integration tests for the reduced distance-geometry embedder and metal
//! coordination (Tasks 6.1.2 and 6.1.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_build::{
    EmbedOptions, apply_coordination, coordination_geometry, embed_smiles, parse_smiles,
};

fn distance(points: &[[f64; 3]], first: usize, second: usize) -> f64 {
    let a = points[first];
    let b = points[second];
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

#[test]
fn ethanol_embeds_three_atoms_with_reasonable_bonds() {
    let graph = parse_smiles("CCO").unwrap();
    let structure = embed_smiles(&graph, &EmbedOptions::default()).unwrap();
    assert_eq!(structure.sites().len(), 3);
    assert_eq!(structure.formula(), "C2 O");
    let points = structure.cart_sites();
    assert!(points.iter().flatten().all(|value| value.is_finite()));
    let carbon_carbon = distance(&points, 0, 1);
    let carbon_oxygen = distance(&points, 1, 2);
    assert!((carbon_carbon - 1.54).abs() < 0.15, "C-C {carbon_carbon}");
    assert!((carbon_oxygen - 1.43).abs() < 0.15, "C-O {carbon_oxygen}");
}

#[test]
fn benzene_ring_closes() {
    let graph = parse_smiles("c1ccccc1").unwrap();
    let structure = embed_smiles(&graph, &EmbedOptions::default()).unwrap();
    assert_eq!(structure.sites().len(), 6);
    let points = structure.cart_sites();
    for index in 0..6 {
        let other = (index + 1) % 6;
        let length = distance(&points, index, other);
        assert!(
            (length - 1.43).abs() < 0.3,
            "ring bond {index}-{other} is {length}"
        );
    }
}

#[test]
fn embedding_is_deterministic_for_a_seed() {
    let graph = parse_smiles("CC(C)O").unwrap();
    let options = EmbedOptions {
        seed: 42,
        ..EmbedOptions::default()
    };
    let first = embed_smiles(&graph, &options).unwrap();
    let second = embed_smiles(&graph, &options).unwrap();
    assert_eq!(first.cart_sites(), second.cart_sites());
}

#[test]
fn invalid_embedding_options_are_rejected() {
    let graph = parse_smiles("CCO").unwrap();
    let bad_length = EmbedOptions {
        bond_length: f64::NAN,
        ..EmbedOptions::default()
    };
    assert!(embed_smiles(&graph, &bad_length).is_err());
    let bad_iterations = EmbedOptions {
        iterations: 0,
        ..EmbedOptions::default()
    };
    assert!(embed_smiles(&graph, &bad_iterations).is_err());
}

#[test]
fn metal_geometries_are_tabulated() {
    let iron = coordination_geometry("Fe").unwrap();
    assert_eq!(iron.coordination, 6);
    assert_eq!(iron.shape, "octahedral");
    let zinc = coordination_geometry("Zn").unwrap();
    assert_eq!(zinc.coordination, 4);
    assert_eq!(zinc.shape, "tetrahedral");
    let copper = coordination_geometry("Cu").unwrap();
    assert_eq!(copper.coordination, 4);
    assert_eq!(copper.shape, "square planar");
    let magnesium = coordination_geometry("Mg").unwrap();
    assert_eq!(magnesium.coordination, 6);
    assert!(coordination_geometry("C").is_none());
    assert!(coordination_geometry("Xx").is_none());
}

#[test]
fn coordination_places_neighbours_on_the_shape() {
    let graph = parse_smiles("N[Pt](N)(Cl)Cl").unwrap();
    let geometry = coordination_geometry("Pt").unwrap();
    assert_eq!(geometry.shape, "square planar");
    let (structure, warnings) =
        apply_coordination(&graph, &geometry, &EmbedOptions::default()).unwrap();
    assert!(warnings.is_empty());
    assert_eq!(structure.sites().len(), 5);
    let points = structure.cart_sites();
    let metal = points[1];
    for neighbour in [0, 2, 3, 4] {
        let point = points[neighbour];
        let length = ((point[0] - metal[0]).powi(2)
            + (point[1] - metal[1]).powi(2)
            + (point[2] - metal[2]).powi(2))
        .sqrt();
        assert!((length - 1.5).abs() < 1e-6, "metal bond {length}");
    }
}

#[test]
fn coordination_warns_when_the_degree_differs() {
    let graph = parse_smiles("N[Pt]N").unwrap();
    let geometry = coordination_geometry("Pt").unwrap();
    let (structure, warnings) =
        apply_coordination(&graph, &geometry, &EmbedOptions::default()).unwrap();
    assert_eq!(structure.sites().len(), 3);
    assert_eq!(warnings.len(), 1);
    assert!(warnings[0].contains("2 neighbours"));
    assert!(warnings[0].contains("square planar"));
}

/// Bond lengths and bond angles of the embedded heavy-atom skeleton.
fn skeleton_geometry(smiles: &str) -> (Vec<f64>, Vec<f64>) {
    let graph = parse_smiles(smiles).unwrap();
    let structure = embed_smiles(&graph, &EmbedOptions::default()).unwrap();
    let points = structure.cart_sites();
    let distance = |a: usize, b: usize| {
        let (p, q) = (points[a], points[b]);
        ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
    };
    let lengths = graph
        .bonds
        .iter()
        .map(|bond| distance(bond.a, bond.b))
        .collect();
    let mut angles = Vec::new();
    for center in 0..graph.atoms.len() {
        let arms: Vec<usize> = graph
            .bonds
            .iter()
            .filter_map(|bond| match (bond.a == center, bond.b == center) {
                (true, _) => Some(bond.b),
                (_, true) => Some(bond.a),
                _ => None,
            })
            .collect();
        for (index, &first) in arms.iter().enumerate() {
            for &second in &arms[index + 1..] {
                let (a, b, c) = (
                    distance(center, first),
                    distance(center, second),
                    distance(first, second),
                );
                angles.push(
                    ((a * a + b * b - c * c) / (2.0 * a * b))
                        .acos()
                        .to_degrees(),
                );
            }
        }
    }
    (lengths, angles)
}

#[test]
fn embedded_rings_keep_their_bond_angles() {
    for (smiles, angle, tolerance) in [
        ("c1ccccc1", 120.0, 3.0),
        ("c1ccc2ccccc2c1", 120.0, 6.0),
        ("C1CCCCC1", 109.47, 6.0),
        ("CC(C)(C)C", 109.47, 2.0),
        ("C#CC", 180.0, 6.0),
    ] {
        let (lengths, angles) = skeleton_geometry(smiles);
        assert!(
            lengths.iter().all(|length| (1.15..1.62).contains(length)),
            "{smiles}: {lengths:?}"
        );
        assert!(
            angles.iter().all(|value| (value - angle).abs() < tolerance),
            "{smiles}: {angles:?}"
        );
    }
}
