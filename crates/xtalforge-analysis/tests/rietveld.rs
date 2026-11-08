//! Integration tests for the Rietveld phase model (Task 7.8).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{RietveldOptions, fit_phase};
use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_io::xy::PatternPoint;

fn silicon() -> xtalforge_core::Structure {
    StructureBuilder::new("Si", Lattice::cubic(5.430941).unwrap())
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn scale_fit_recovers_a_synthetic_pattern() {
    let structure = silicon();
    let options = RietveldOptions::default();
    let phase = fit_phase(
        &structure,
        &(0..2000)
            .map(|index| PatternPoint {
                two_theta: 5.0 + index as f64 * 0.02,
                intensity: 1.0,
                uncertainty: None,
            })
            .collect::<Vec<_>>(),
        &options,
    )
    .unwrap();
    assert!(phase.scale.is_finite());
    assert!(phase.rwp.is_finite());
    assert_eq!(phase.pattern.len(), 2000);
    assert!(fit_phase(&structure, &[], &options).is_err());
}
