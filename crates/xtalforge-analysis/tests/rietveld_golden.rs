//! Rietveld refinement against the Phase 7 golden rutile pattern (Task 7.8).
//!
//! The pattern comes from the Crystal-Builder `_rutile_xy` fixture: a rutile
//! cell, Cu K alpha 1 + K alpha 2, 20-80 degrees, Poisson noise with seed 0.
//! The structure is started 0.03 fractional units away from the generating
//! oxygen position, so the test proves that the McCusker stages move the site
//! back along its free Wyckoff direction.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{RefineOptions, RietveldOptions, free_directions, refine_phase};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_io::xy::PatternPoint;

const RUTILE: &str = include_str!("../../../tests/corpus/simple/TiO2.cif");
const GOLDEN: &str = include_str!("../../../tests/golden/phase-7/rietveld/rutile.xy");

/// Generating oxygen position in the golden pattern.
const X_TRUE: f64 = 0.29;
/// Generating cell.
const A_TRUE: f64 = 4.5940;
const C_TRUE: f64 = 2.9590;

fn rutile(x: f64) -> Structure {
    let base = xtalforge_io::cif::read_cif(RUTILE).unwrap();
    let mut builder = StructureBuilder::new(
        "rutile",
        Lattice::from_parameters(A_TRUE, A_TRUE, C_TRUE, 90.0, 90.0, 90.0).unwrap(),
    )
    .hall(base.hall());
    builder = builder.site(Site::new("Ti", "Ti", [0.0, 0.0, 0.0]).unwrap());
    builder = builder.site(Site::new("O", "O", [x, x, 0.0]).unwrap());
    for operation in base.operations() {
        builder = builder.operation(*operation);
    }
    builder.build().unwrap()
}

fn observed() -> Vec<PatternPoint> {
    xtalforge_io::read_xy(GOLDEN).unwrap()
}

fn value(result: &xtalforge_analysis::RefinedPhase, name: &str) -> f64 {
    result
        .parameters
        .iter()
        .find(|(key, _)| key == name)
        .map(|(_, value)| *value)
        .unwrap_or_else(|| panic!("missing parameter {name}"))
}

fn esd(result: &xtalforge_analysis::RefinedPhase, name: &str) -> f64 {
    result
        .parameter_esds
        .iter()
        .find(|(key, _)| key == name)
        .map(|(_, value)| *value)
        .unwrap_or_else(|| panic!("missing esd {name}"))
}

#[test]
fn refinement_returns_the_oxygen_to_its_golden_position() {
    let start = rutile(0.26);
    let options = RefineOptions {
        max_iterations: 30,
        tolerance: 1e-6,
        refine_positions: true,
        refine_cell: true,
        refine_zero: true,
        refine_width: true,
        doublet: true,
        ..RefineOptions::default()
    };
    let result = refine_phase(&start, &observed(), &RietveldOptions::default(), &options).unwrap();
    assert_eq!(
        result.stages,
        vec!["background_scale", "zero_cell", "width", "structure"]
    );
    assert!(result.converged, "the refinement did not converge");
    assert!(result.rwp < 0.25, "Rwp {}", result.rwp);

    let oxygen = result
        .structure
        .sites()
        .iter()
        .find(|site| site.label() == "O")
        .expect("the oxygen site must survive");
    let frac = oxygen.frac();
    assert!(
        (frac[0] - X_TRUE).abs() < 3e-3,
        "x = {} is not {X_TRUE}",
        frac[0]
    );
    // The site stays on (x, x, 0): only its free Wyckoff direction moved.
    assert!((frac[1] - frac[0]).abs() < 1e-10);
    assert!(frac[2].abs() < 1e-10);

    let parameters = result.structure.lattice().parameters();
    assert!(
        (parameters[0] - A_TRUE).abs() < 2e-3,
        "a = {} is not {A_TRUE}",
        parameters[0]
    );
    // The tetragonal tie survives the cell refinement.
    assert!((parameters[1] - parameters[0]).abs() < 1e-12);
    assert!(
        (parameters[2] - C_TRUE).abs() < 2e-3,
        "c = {} is not {C_TRUE}",
        parameters[2]
    );
    assert!(
        value(&result, "zero").abs() < 0.02,
        "zero shift {}",
        value(&result, "zero")
    );
    assert_eq!(result.structure.sites().len(), 2);

    for name in ["zero", "cell.0", "cell.2", "width_scale", "O.d0"] {
        let uncertainty = esd(&result, name);
        assert!(
            uncertainty.is_finite() && uncertainty > 0.0,
            "{name} esd {uncertainty}"
        );
    }
}

#[test]
fn the_oxygen_free_direction_is_the_wyckoff_one() {
    let base = xtalforge_io::cif::read_cif(RUTILE).unwrap();
    let directions = free_directions(base.operations(), [X_TRUE, X_TRUE, 0.0]);
    assert_eq!(directions.len(), 1);
    let direction = directions[0];
    let norm = direction
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    assert!((norm - 1.0).abs() < 1e-12);
    assert!((direction[0] - direction[1]).abs() < 1e-12);
    assert!((direction[2]).abs() < 1e-12);
}
