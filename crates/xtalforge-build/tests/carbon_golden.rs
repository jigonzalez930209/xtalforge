//! Statistical golden test for the disordered-carbon recipe (Task 6.8.5).
//!
//! The Crystal-Builder reference uses a different RNG, remesh and ribbon
//! path, so the exact atom list cannot match. The test compares what the
//! recipe promises: density, atom count, edge chemistry, one percolating
//! piece and a six-membered ring census.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::Value;
use xtalforge_build::carbon_recipe::{CarbonRecipe, build_carbon};

const GOLDEN: &str = include_str!("../../../tests/golden/phase-6/carbon/dia.json");

#[test]
fn dia_recipe_matches_the_reference_statistics() {
    let golden: Value = serde_json::from_str(GOLDEN).unwrap();
    let recipe = CarbonRecipe {
        repeat: [1, 1, 1],
        seed: golden["seed"].as_u64().unwrap_or(0),
        ..CarbonRecipe::default()
    };
    let build = build_carbon(&recipe).unwrap();
    let atoms = build.structure.sites().len() as f64;
    let expected_atoms = golden["atoms"].as_f64().unwrap();
    assert!(
        (atoms - expected_atoms).abs() <= 0.1 * expected_atoms,
        "atoms {atoms} vs {expected_atoms}"
    );
    let expected_density = golden["built_density"].as_f64().unwrap();
    assert!(
        (build.density - expected_density).abs() <= 0.03,
        "density {} vs {expected_density}",
        build.density
    );
    assert_eq!(build.pieces, 1);
    assert_eq!(build.periodicity, 3);
    let edge = golden["edge_fraction"].as_f64().unwrap();
    assert!(
        (build.edge_fraction - edge).abs() <= 0.08,
        "edge fraction {} vs {edge}",
        build.edge_fraction
    );
    let expected = &golden["ratios"];
    for (index, key) in ["H/C", "F/C", "O/C"].iter().enumerate() {
        let want = expected[*key].as_f64().unwrap();
        let got = build.ratios[index];
        assert!((got - want).abs() <= 0.03, "{key}: {got} vs {want}");
    }
    assert!(
        (1.25..=1.56).contains(&build.mean_bond),
        "mean C-C bond {}",
        build.mean_bond
    );
    let expected_closed = golden["closed_carbons"].as_f64().unwrap();
    assert!(
        (build.closed_carbons as f64 - expected_closed).abs() <= 0.15 * expected_closed,
        "closed carbons {} vs {expected_closed}",
        build.closed_carbons
    );
    let total: usize = build.census.iter().map(|(_, count)| *count).sum();
    let reference: usize = golden["census"]
        .as_object()
        .unwrap()
        .values()
        .map(|value| value.as_u64().unwrap() as usize)
        .sum();
    assert!(
        total * 2 >= reference,
        "ring census {total} is far below the reference {reference}"
    );
    let sixes = build
        .census
        .iter()
        .find(|(size, _)| *size == 6)
        .map_or(0, |(_, count)| *count);
    let largest = build
        .census
        .iter()
        .map(|(_, count)| *count)
        .max()
        .unwrap_or(0);
    assert_eq!(sixes, largest, "six-membered rings are not the mode");
}
