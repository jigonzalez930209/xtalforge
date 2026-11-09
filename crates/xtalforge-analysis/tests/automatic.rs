//! Integration tests for the automatic peaks-index-Le Bail-Rietveld run
//! (Task 7.10).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{AutomaticOptions, PxrdOptions, Radiation, calculate_pxrd, run_automatic};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_io::xy::PatternPoint;

fn silicon() -> Structure {
    StructureBuilder::new("Si", Lattice::cubic(5.430941).unwrap())
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

fn pattern() -> Vec<PatternPoint> {
    calculate_pxrd(
        &silicon(),
        &PxrdOptions {
            radiation: Radiation::CuKa,
            two_theta_range: [10.0, 80.0],
            step: 0.01,
            fwhm: 0.05,
        },
    )
    .unwrap()
    .points
    .iter()
    .map(|point| PatternPoint {
        two_theta: point.two_theta,
        intensity: point.intensity,
        uncertainty: None,
    })
    .collect()
}

#[test]
fn automatic_run_records_steps_in_order() {
    let observed = pattern();
    let result = run_automatic(&silicon(), &observed, &AutomaticOptions::default()).unwrap();
    assert_eq!(
        result.steps,
        vec![
            "peaks",
            "index",
            "leading_cells",
            "le_bail",
            "pawley",
            "rietveld"
        ]
    );
    assert!(!result.peaks.is_empty());
    assert!(result.lebail.is_some());
    assert!(result.pawley.is_some());
    assert!(result.phase.is_some());
    let phase = result.phase.unwrap();
    assert!(phase.rwp.is_finite());
    assert!(phase.rwp >= 0.0);
    if let Some(index) = result.index {
        assert!(index.cell.a.is_finite());
        assert!(index.completeness >= AutomaticOptions::default().min_completeness);
        assert_eq!(index.system, xtalforge_analysis::LatticeSystem::Cubic);
    }
    let pawley = result.pawley.unwrap();
    // The default profile is a pseudo-Voigt while the synthetic pattern is
    // Gaussian, so only finiteness and the reflection list are asserted here.
    assert!(pawley.rwp.is_finite() && pawley.rwp >= 0.0);
    assert!(!pawley.reflections.is_empty());
    assert!(pawley.cell[0].is_finite());
    assert!(!result.leading_cells.is_empty());
    let table = &result.leading_cells;
    assert!(
        table
            .iter()
            .any(|row| row.index.system == xtalforge_analysis::LatticeSystem::Cubic)
    );
    assert!(table.iter().any(|row| row.pawley_rwp.is_some()));
    assert!(table.iter().any(|row| row.pawley_cell_esds.is_some()));
    let rwps: Vec<f64> = table.iter().filter_map(|row| row.pawley_rwp).collect();
    assert!(rwps.windows(2).all(|pair| pair[0] <= pair[1]));
}

#[test]
fn empty_and_peakless_inputs_error() {
    let structure = silicon();
    assert!(run_automatic(&structure, &[], &AutomaticOptions::default()).is_err());
    let flat: Vec<PatternPoint> = (0..100)
        .map(|index| PatternPoint {
            two_theta: 10.0 + index as f64 * 0.02,
            intensity: 0.0,
            uncertainty: None,
        })
        .collect();
    assert!(run_automatic(&structure, &flat, &AutomaticOptions::default()).is_err());
    assert!(
        run_automatic(
            &structure,
            &flat,
            &AutomaticOptions {
                min_completeness: 2.0,
                ..AutomaticOptions::default()
            }
        )
        .is_err()
    );
}

#[test]
fn automatic_run_indexes_a_tetragonal_cell() {
    let structure = StructureBuilder::new(
        "TiO2",
        Lattice::new([[4.5937, 0.0, 0.0], [0.0, 4.5937, 0.0], [0.0, 0.0, 2.9587]]).unwrap(),
    )
    .site(Site::new("Ti1", "Ti", [0.0, 0.0, 0.0]).unwrap())
    .site(Site::new("O1", "O", [0.3, 0.3, 0.0]).unwrap())
    .build()
    .unwrap();
    let observed: Vec<PatternPoint> = calculate_pxrd(
        &structure,
        &PxrdOptions {
            radiation: Radiation::CuKa,
            two_theta_range: [10.0, 70.0],
            step: 0.01,
            fwhm: 0.05,
        },
    )
    .unwrap()
    .points
    .iter()
    .map(|point| PatternPoint {
        two_theta: point.two_theta,
        intensity: point.intensity,
        uncertainty: None,
    })
    .collect();
    let result = run_automatic(&structure, &observed, &AutomaticOptions::default()).unwrap();
    let index = result.index.expect("an indexed cell");
    assert_eq!(index.system, xtalforge_analysis::LatticeSystem::Tetragonal);
    assert!((index.cell.a - 4.5937).abs() < 0.05);
    assert!((index.cell.c - 2.9587).abs() < 0.05);
    assert!(
        result
            .leading_cells
            .iter()
            .any(|row| row.index.system == xtalforge_analysis::LatticeSystem::Tetragonal)
    );
}
