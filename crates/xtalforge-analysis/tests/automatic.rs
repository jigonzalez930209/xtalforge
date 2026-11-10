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
    assert!(
        rwps.windows(2).all(|pair| pair[0] <= pair[1] * 1.01),
        "Rwp ties must stay within one percent"
    );
    assert!(table.iter().all(|row| row.class_rank >= 1));
    assert!(table.iter().all(|row| !row.class_symbol.is_empty()));
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

#[test]
fn automatic_run_screens_space_groups_from_absences() {
    let cif = include_str!("../../../tests/corpus/simple/TiO2.cif");
    let base = xtalforge_io::cif::read_cif(cif).unwrap();
    let mut builder = StructureBuilder::new(
        "rutile",
        Lattice::from_parameters(4.5940, 4.5940, 2.9590, 90.0, 90.0, 90.0).unwrap(),
    )
    .hall(base.hall());
    builder = builder.site(Site::new("Ti", "Ti", [0.0, 0.0, 0.0]).unwrap());
    builder = builder.site(Site::new("O", "O", [0.29, 0.29, 0.0]).unwrap());
    for operation in base.operations() {
        builder = builder.operation(*operation);
    }
    let structure = builder.build().unwrap();
    let observed = xtalforge_io::read_xy(include_str!(
        "../../../tests/golden/phase-7/rietveld/rutile.xy"
    ))
    .unwrap();
    let result = run_automatic(&structure, &observed, &AutomaticOptions::default()).unwrap();
    for row in &result.leading_cells {
        println!(
            "{} class {} {} delta {:.2} rwp {:?} groups {:?}",
            row.index.system.name(),
            row.class_rank,
            row.class_symbol,
            row.delta_bic,
            row.pawley_rwp,
            row.space_groups
                .iter()
                .map(|group| group.number)
                .collect::<Vec<_>>()
        );
    }
    assert!(
        result.leading_cells[0]
            .space_groups
            .iter()
            .any(|group| group.number == 136),
        "the rutile class is not the first row: {:?}",
        result.leading_cells[0]
    );
    assert_eq!(result.leading_cells[0].class_rank, 1);
}

/// A backgrounded, noisy copy of the silicon pattern with a deterministic
/// pseudo-random generator (linear congruential), so the test never changes.
fn noisy_pattern() -> Vec<PatternPoint> {
    let clean = calculate_pxrd(
        &silicon(),
        &PxrdOptions {
            radiation: Radiation::CuKa,
            two_theta_range: [20.0, 80.0],
            step: 0.04,
            fwhm: 0.1,
        },
    )
    .unwrap();
    let maximum = clean
        .points
        .iter()
        .map(|point| point.intensity)
        .fold(0.0_f64, f64::max);
    let mut state: u64 = 12_345;
    clean
        .points
        .iter()
        .map(|point| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let uniform = (state >> 33) as f64 / f64::from(1_u32 << 31) - 0.5;
            let signal = 8000.0 * point.intensity / maximum;
            let background = 150.0 + 0.5 * (80.0 - point.two_theta);
            PatternPoint {
                two_theta: point.two_theta,
                intensity: signal + background + 40.0 * uniform,
                uncertainty: None,
            }
        })
        .collect()
}

#[test]
fn a_noisy_pattern_on_a_raised_baseline_is_not_indexed_on_its_noise() {
    let result = run_automatic(&silicon(), &noisy_pattern(), &AutomaticOptions::default()).unwrap();
    // The one-atom cubic cell has 17 reflections between 20 and 80 degrees;
    // noise maxima must not be taken for peaks (before the background fix the
    // search kept the 30 strongest noise maxima), and the cell is the true one.
    assert!(result.peaks.len() <= 20, "{} peaks", result.peaks.len());
    let a = result.index.as_ref().unwrap().cell.a;
    assert!((a - 5.43).abs() < 0.05, "indexed a = {a}");
}
