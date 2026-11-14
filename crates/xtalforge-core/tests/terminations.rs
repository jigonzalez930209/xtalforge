//! Integration tests for slab termination options (Task 17.3.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{
    Lattice, Site, StructureBuilder, TerminationOptions, slab_cut_shifts, slab_terminations,
    slab_terminations_with_options,
};

fn cscl() -> xtalforge_core::Structure {
    StructureBuilder::new("CsCl", Lattice::cubic(4.2).unwrap())
        .site(Site::new("Cs1", "Cs", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Cl1", "Cl", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn default_cuts_are_the_layer_midpoints() {
    let shifts = slab_cut_shifts(&cscl(), [1, 0, 0], &TerminationOptions::default()).unwrap();
    assert_eq!(shifts.len(), 2);
    assert!((shifts[0] - 0.25).abs() < 1e-9, "{shifts:?}");
    assert!((shifts[1] - 0.75).abs() < 1e-9, "{shifts:?}");
}

#[test]
fn default_terminations_match_the_task_9_4_report() {
    let report = slab_terminations(&cscl(), [1, 0, 0], 1, 10.0).unwrap();
    assert_eq!(report.hkl, [1, 0, 0]);
    assert_eq!(report.terminations.len(), 2);
    assert_eq!(report.terminations[0].slab.structure.sites().len(), 2);
    assert_eq!(report.terminations[1].slab.structure.sites().len(), 2);
    assert!(report.terminations[0].shift < report.terminations[1].shift);
    for termination in &report.terminations {
        assert!(termination.slab.structure.is_p1());
    }
}

#[test]
fn intra_layer_cuts_add_the_gap_quarters() {
    let options = TerminationOptions {
        include_intra_layer: true,
        ..TerminationOptions::default()
    };
    let shifts = slab_cut_shifts(&cscl(), [1, 0, 0], &options).unwrap();
    let expected = [0.125, 0.25, 0.375, 0.625, 0.75, 0.875];
    assert_eq!(shifts.len(), expected.len());
    for (shift, expected) in shifts.iter().zip(expected) {
        assert!((shift - expected).abs() < 1e-9, "{shifts:?}");
    }
    let report = slab_terminations_with_options(&cscl(), [1, 0, 0], 1, 10.0, &options).unwrap();
    assert_eq!(report.terminations.len(), 6);
}

#[test]
fn explicit_cuts_are_wrapped_and_sorted() {
    let options = TerminationOptions {
        cuts: vec![0.3, 1.25, -0.25],
        ..TerminationOptions::default()
    };
    let shifts = slab_cut_shifts(&cscl(), [1, 0, 0], &options).unwrap();
    // 1.25 wraps to 0.25 and -0.25 wraps to 0.75; both duplicates are merged.
    let expected = [0.25, 0.3, 0.75];
    assert_eq!(shifts.len(), expected.len(), "{shifts:?}");
    for (shift, expected) in shifts.iter().zip(expected) {
        assert!((shift - expected).abs() < 1e-9, "{shifts:?}");
    }
}

#[test]
fn a_single_layer_still_has_one_midpoint_cut() {
    let structure = StructureBuilder::new("C", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let shifts = slab_cut_shifts(&structure, [1, 0, 0], &TerminationOptions::default()).unwrap();
    assert_eq!(shifts.len(), 1);
    assert!((shifts[0] - 0.5).abs() < 1e-9);
}

#[test]
fn invalid_termination_options_are_rejected() {
    let structure = cscl();
    for options in [
        TerminationOptions {
            cuts: vec![f64::NAN],
            ..TerminationOptions::default()
        },
        TerminationOptions {
            layer_tolerance: 0.0,
            ..TerminationOptions::default()
        },
        TerminationOptions {
            site_layers: 0,
            ..TerminationOptions::default()
        },
    ] {
        assert!(slab_cut_shifts(&structure, [1, 0, 0], &options).is_err());
    }
    assert!(slab_terminations(&structure, [0, 0, 0], 1, 10.0).is_err());
    assert!(slab_terminations(&structure, [1, 0, 0], 0, 10.0).is_err());
    assert!(slab_terminations(&structure, [1, 0, 0], 1, -1.0).is_err());
    let empty = StructureBuilder::new("empty", Lattice::cubic(4.0).unwrap())
        .build()
        .unwrap();
    assert!(slab_cut_shifts(&empty, [1, 0, 0], &TerminationOptions::default()).is_err());
}
