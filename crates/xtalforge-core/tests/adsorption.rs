//! Integration tests for adsorption sites below the top layer (Task 17.3.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{
    AdsorptionKind, Lattice, Site, StructureBuilder, TerminationOptions, adsorption_sites,
    adsorption_sites_with_options,
};

fn two_layer_slab() -> xtalforge_core::Structure {
    StructureBuilder::new("two-layer", Lattice::cubic(10.0).unwrap())
        .site(Site::new("Pt1", "Pt", [0.1, 0.1, 0.6]).unwrap())
        .site(Site::new("Pt2", "Pt", [0.6, 0.1, 0.6]).unwrap())
        .site(Site::new("O1", "O", [0.1, 0.1, 0.3]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn the_default_search_is_the_top_layer_only() {
    let slab = two_layer_slab();
    let default = adsorption_sites(&slab, 0.05).unwrap();
    let explicit =
        adsorption_sites_with_options(&slab, 0.05, &TerminationOptions::default()).unwrap();
    assert_eq!(default, explicit);
    assert!(
        default
            .iter()
            .all(|site| (site.frac[2] - 0.6).abs() < 1e-12)
    );
    assert!(default.iter().all(|site| site.kind == AdsorptionKind::Top));
    assert_eq!(default.len(), 2);
}

#[test]
fn site_layers_reaches_the_layer_below_the_top() {
    let slab = two_layer_slab();
    let options = TerminationOptions {
        site_layers: 2,
        ..TerminationOptions::default()
    };
    let sites = adsorption_sites_with_options(&slab, 0.05, &options).unwrap();
    let top = sites
        .iter()
        .filter(|site| (site.frac[2] - 0.6).abs() < 1e-12)
        .count();
    let lower = sites
        .iter()
        .filter(|site| (site.frac[2] - 0.3).abs() < 1e-12)
        .count();
    assert_eq!(top, 2);
    assert_eq!(lower, 1);
    assert_eq!(sites.len(), 3);
}

#[test]
fn asking_for_more_layers_than_exist_stops_at_the_bottom() {
    let slab = two_layer_slab();
    let options = TerminationOptions {
        site_layers: 8,
        ..TerminationOptions::default()
    };
    let sites = adsorption_sites_with_options(&slab, 0.05, &options).unwrap();
    assert_eq!(sites.len(), 3);
    let heights: Vec<f64> = sites.iter().map(|site| site.frac[2]).collect();
    assert!(
        heights
            .iter()
            .all(|z| (z - 0.6).abs() < 1e-12 || (z - 0.3).abs() < 1e-12)
    );
}

#[test]
fn invalid_adsorption_inputs_are_rejected() {
    let slab = two_layer_slab();
    assert!(adsorption_sites(&slab, 0.0).is_err());
    assert!(adsorption_sites(&slab, f64::NAN).is_err());
    let no_layers = TerminationOptions {
        site_layers: 0,
        ..TerminationOptions::default()
    };
    assert!(adsorption_sites_with_options(&slab, 0.05, &no_layers).is_err());
    let empty = StructureBuilder::new("empty", Lattice::cubic(10.0).unwrap())
        .build()
        .unwrap();
    assert!(adsorption_sites(&empty, 0.05).is_err());
}
