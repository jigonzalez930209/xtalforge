//! Integration tests for the persistent powder parameter sets (Task 14.1.5).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{ParameterSet, PawleyOptions, PawleyResult, RefinedPhase};
use xtalforge_core::{Lattice, Site, StructureBuilder};

fn cubic() -> xtalforge_core::Structure {
    StructureBuilder::new("CsCl", Lattice::cubic(4.0).unwrap())
        .site(Site::new("Cs1", "Cs", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

fn pawley_fit() -> PawleyResult {
    PawleyResult {
        reflections: Vec::new(),
        cell: [4.0, 4.0, 4.0, 90.0, 90.0, 90.0],
        cell_esds: [1e-4; 6],
        zero_shift: 0.02,
        zero_shift_esd: 1e-3,
        size: 0.12,
        strain: 0.03,
        rp: 0.02,
        rwp: 0.03,
        rexp: 0.01,
        chi_squared: 1.2,
        pattern: Vec::new(),
        converged: true,
        iterations: 3,
    }
}

fn refined() -> RefinedPhase {
    RefinedPhase {
        structure: cubic(),
        rwp: 0.1,
        parameters: vec![("cell.0".to_owned(), 4.0), ("zero".to_owned(), 0.02)],
        parameter_esds: Vec::new(),
        stages: vec!["zero_cell".to_owned()],
        converged: true,
        iterations: 2,
    }
}

#[test]
fn refined_and_pawley_results_become_persistent_sets() {
    let set = ParameterSet::from_refined_phase("phase", &refined());
    assert_eq!(set.get("cell.0"), Some(4.0));
    assert_eq!(set.get("zero"), Some(0.02));
    assert_eq!(set.get("rwp"), Some(0.1));
    let text = set.to_json().unwrap();
    assert_eq!(ParameterSet::from_json(&text).unwrap(), set);

    let fit = pawley_fit();
    let pawley = ParameterSet::from_pawley("phase", &fit, 0.7);
    assert_eq!(pawley.get("eta"), Some(0.7));
    assert_eq!(pawley.get("fwhm"), Some(0.12));
    assert_eq!(pawley.get("strain"), Some(0.03));
    assert_eq!(pawley.get("cell.a"), Some(4.0));
    assert_eq!(pawley.get("cell.gamma"), Some(90.0));

    let options = pawley.pawley_options(&PawleyOptions::default());
    assert_eq!(options.eta, 0.7);
    assert_eq!(options.fwhm, 0.12);

    let mut broken = ParameterSet::new("broken");
    broken.set("eta", 2.0).unwrap();
    broken.set("fwhm", -1.0).unwrap();
    let unchanged = broken.pawley_options(&PawleyOptions::default());
    assert_eq!(unchanged.eta, PawleyOptions::default().eta);
    assert_eq!(unchanged.fwhm, PawleyOptions::default().fwhm);
}
