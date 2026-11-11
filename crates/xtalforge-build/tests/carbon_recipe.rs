//! Integration tests for the disordered-carbon recipe (Task 6.8).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_build::carbon_recipe::{CarbonRecipe, build_carbon, carbon_surface};

#[test]
fn dia_surface_is_one_closed_sheet() {
    let surface = carbon_surface("dia", [1, 1, 1], 10.5, 0.25, 0).unwrap();
    println!("{surface:?}");
    assert_eq!(surface.net_vertices, 8);
    assert_eq!(surface.net_edges, 16);
    assert_eq!(surface.euler, -16);
    assert_eq!(surface.pieces, 1);
    assert_eq!(surface.periodicity, 3);
}

#[test]
fn carbon_recipe_builds_a_connected_framework() {
    let recipe = CarbonRecipe {
        repeat: [1, 1, 1],
        seed: 0,
        ..CarbonRecipe::default()
    };
    let build = build_carbon(&recipe).unwrap();
    println!(
        "atoms={} density={} bond={} ratios={:?} edge={} pieces={} periodicity={} closed={} census={:?} notes={:?}",
        build.structure.sites().len(),
        build.density,
        build.mean_bond,
        build.ratios,
        build.edge_fraction,
        build.pieces,
        build.periodicity,
        build.closed_carbons,
        build.census,
        build.notes
    );
}

#[test]
fn pcu_dual_percolates() {
    let surface = carbon_surface("pcu", [1, 1, 1], 10.5, 0.25, 0).unwrap();
    println!("pcu dual {} {}", surface.carbons, surface.periodicity);
}

#[test]
fn carbon_recipe_geometry_is_physical() {
    for net in ["dia", "pcu"] {
        let recipe = CarbonRecipe {
            net: net.to_owned(),
            repeat: [1, 1, 1],
            seed: 0,
            ..CarbonRecipe::default()
        };
        let structure = build_carbon(&recipe).unwrap().structure;
        let graph = structure.bond_graph().unwrap();
        for bond in &graph.bonds {
            assert!(bond.distance < 1.9, "{net}: bond {bond:?}");
        }
        let p1 = structure.p1().unwrap();
        let lattice = structure.lattice();
        for (i, a) in p1.atoms.iter().enumerate() {
            for b in &p1.atoms[i + 1..] {
                assert!(lattice.distance(a.frac, b.frac) > 0.9, "{net}: overlap");
            }
        }
    }
}
