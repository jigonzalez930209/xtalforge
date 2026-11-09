//! Integration tests for the preparation pipeline (Task 4.7).
//!
//! The water case exercises the default pipeline: FIRE, QEq charges and the
//! report fields.  Ewald and vibrations are checked as separate opt-in steps so
//! a failure names the step that broke rather than the whole workflow.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Site, Structure};
use xtalforge_ff::{OptimizeOptions, PrepareOptions, PrepareStep, UffCalculator, prepare};

fn molecule(name: &str, atoms: &[(&str, [f64; 3])]) -> Structure {
    let points: Vec<(String, [f64; 3])> = atoms
        .iter()
        .map(|(element, point)| ((*element).to_owned(), *point))
        .collect();
    Structure::molecule(name, &points, 8.0).expect("molecule")
}

fn distorted_water() -> Structure {
    molecule(
        "distorted-water",
        &[
            ("O", [0.0, 0.0, 0.0]),
            ("H", [1.08, 0.0, 0.0]),
            ("H", [-0.27, 1.04, 0.0]),
        ],
    )
}

fn options() -> PrepareOptions {
    PrepareOptions {
        minimize: true,
        charges: true,
        ewald: false,
        vibrations: false,
        optimize: OptimizeOptions {
            max_steps: 300,
            force_tolerance: 0.05,
            ..OptimizeOptions::default()
        },
        ..PrepareOptions::default()
    }
}

fn names(steps: &[PrepareStep]) -> Vec<&str> {
    steps.iter().map(|step| step.name.as_str()).collect()
}

#[test]
fn prepare_relaxes_water_and_records_every_step() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let report = prepare(&calculator, &distorted_water(), &options()).unwrap();
    assert_eq!(
        names(&report.steps),
        [
            "duplicates",
            "disorder",
            "solvent",
            "hydrogens",
            "minimize",
            "charges",
            "ewald",
            "vibrations",
            "close_contacts",
            "charge_balance",
        ]
    );
    assert!(report.steps[4].done && report.steps[5].done);
    assert!(!report.steps[6].done && !report.steps[7].done);
    assert!(!report.steps[8].done && !report.steps[9].done);
    assert!(report.energy.is_finite(), "energy {}", report.energy);
    assert!(report.total_charge.abs() < 1e-8, "{}", report.total_charge);
    assert!(report.ewald_energy.is_none());
    assert!(report.frequencies.is_empty());
    let optimized = report.optimized.expect("minimizer result");
    assert!(optimized.steps > 0);
}

#[test]
fn prepare_can_evaluate_ewald_on_the_assigned_charges() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let mut settings = options();
    settings.minimize = false;
    settings.ewald = true;
    let report = prepare(&calculator, &distorted_water(), &settings).unwrap();
    assert!(report.steps[6].done);
    let energy = report.ewald_energy.expect("Ewald energy");
    assert!(energy.is_finite(), "Ewald energy {energy}");
    assert!(report.total_charge.abs() < 1e-8);
    assert!(report.optimized.is_none());
}

#[test]
fn prepare_can_skip_optional_steps() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let settings = PrepareOptions {
        minimize: false,
        charges: false,
        ..PrepareOptions::default()
    };
    let report = prepare(&calculator, &distorted_water(), &settings).unwrap();
    assert!(report.steps.iter().all(|step| !step.done));
    assert!(report.energy.is_finite());
    assert!(report.optimized.is_none());
    assert!(report.ewald_energy.is_none());
    assert!(report.frequencies.is_empty());
}

#[test]
fn prepare_can_build_gamma_frequencies() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let settings = PrepareOptions {
        minimize: false,
        charges: false,
        vibrations: true,
        ..PrepareOptions::default()
    };
    let report = prepare(&calculator, &distorted_water(), &settings).unwrap();
    assert!(report.steps[7].done);
    assert_eq!(report.frequencies.len(), 9);
    assert!(report.frequencies.iter().all(|value| value.is_finite()));
}

#[test]
fn prepare_rejects_invalid_options() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let invalid_minimizer = PrepareOptions {
        optimize: OptimizeOptions {
            max_steps: 0,
            ..OptimizeOptions::default()
        },
        ..PrepareOptions::default()
    };
    assert!(prepare(&calculator, &distorted_water(), &invalid_minimizer).is_err());
    let invalid_charges = PrepareOptions {
        charge_options: xtalforge_ff::ChargeOptions {
            total_charge: f64::INFINITY,
            ..xtalforge_ff::ChargeOptions::default()
        },
        ..PrepareOptions::default()
    };
    assert!(prepare(&calculator, &distorted_water(), &invalid_charges).is_err());
}

fn cleanup_options() -> PrepareOptions {
    PrepareOptions {
        merge_duplicates: true,
        order_disorder: true,
        remove_solvent: true,
        check_charge_balance: true,
        minimize: false,
        charges: false,
        ..PrepareOptions::default()
    }
}

fn cell_with(sites: Vec<Site>) -> Structure {
    let mut builder = xtalforge_core::StructureBuilder::new(
        "cleanup",
        xtalforge_core::Lattice::cubic(10.0).unwrap(),
    );
    for site in sites {
        builder = builder.site(site);
    }
    builder.build().unwrap()
}

#[test]
fn cleanup_merges_duplicate_sites() {
    let structure = cell_with(vec![
        Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap(),
        Site::new("C2", "C", [0.0, 0.0, 0.0]).unwrap(),
    ]);
    let calculator = UffCalculator::with_defaults().unwrap();
    let report = prepare(&calculator, &structure, &cleanup_options()).unwrap();
    assert_eq!(report.removed_sites, 1);
    assert_eq!(report.charge_balance.expect("balance").total_charge, 0.0);
}

#[test]
fn cleanup_chooses_one_disorder_component() {
    let structure = cell_with(vec![
        Site::new("O1", "O", [0.0, 0.0, 0.0])
            .unwrap()
            .with_occupancy(0.5)
            .unwrap()
            .with_disorder_group(1),
        Site::new("O2", "O", [0.05, 0.0, 0.0])
            .unwrap()
            .with_occupancy(0.5)
            .unwrap()
            .with_disorder_group(2),
        Site::new("C1", "C", [0.5, 0.5, 0.5]).unwrap(),
    ]);
    let calculator = UffCalculator::with_defaults().unwrap();
    let mut options = cleanup_options();
    options.remove_solvent = false;
    let report = prepare(&calculator, &structure, &options).unwrap();
    assert_eq!(report.removed_sites, 1);
}

#[test]
fn cleanup_removes_a_solvent_molecule() {
    let structure = cell_with(vec![
        Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap(),
        Site::new("O1", "O", [0.5, 0.5, 0.5]).unwrap(),
    ]);
    let calculator = UffCalculator::with_defaults().unwrap();
    let report = prepare(&calculator, &structure, &cleanup_options()).unwrap();
    assert_eq!(report.removed_sites, 1);
    assert!(
        report
            .steps
            .iter()
            .any(|step| step.name == "solvent" && step.done)
    );
}

#[test]
fn cleanup_reports_a_non_neutral_cell() {
    let structure = cell_with(vec![
        Site::new("Na1", "Na", [0.0, 0.0, 0.0])
            .unwrap()
            .with_charge(1.0)
            .unwrap(),
    ]);
    let calculator = UffCalculator::with_defaults().unwrap();
    let report = prepare(&calculator, &structure, &cleanup_options()).unwrap();
    let balance = report.charge_balance.expect("balance");
    assert!(!balance.neutral);
    assert!((balance.total_charge - 1.0).abs() < 1e-12);
}

#[test]
fn prepared_golden_samples_are_ordered_and_neutral() {
    use std::fs;
    use std::path::Path;

    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/corpus/prepared");
    let calculator = UffCalculator::with_defaults().unwrap();
    let mut checked = 0;
    for entry in fs::read_dir(root).unwrap() {
        let path = entry.unwrap().path();
        if path.extension().and_then(|value| value.to_str()) != Some("cif") {
            continue;
        }
        let structure = xtalforge_io::read_cif(&fs::read_to_string(&path).unwrap()).unwrap();
        if structure.p1().unwrap().n_atoms() > 600 {
            continue;
        }
        let report = prepare(&calculator, &structure, &cleanup_options()).unwrap();
        assert!(
            report
                .structure
                .sites()
                .iter()
                .all(|site| (site.occupancy() - 1.0).abs() < 1e-6),
            "{}: partial occupancy left",
            path.display()
        );
        assert!(
            report.charge_balance.is_some(),
            "{}: no charge balance",
            path.display()
        );
        checked += 1;
    }
    assert!(checked >= 5, "prepared golden samples checked: {checked}");
}

#[test]
fn cleanup_adds_missing_hydrogens() {
    let structure = cell_with(vec![Site::new("C1", "C", [0.5, 0.5, 0.5]).unwrap()]);
    let calculator = UffCalculator::with_defaults().unwrap();
    let options = PrepareOptions {
        hydrogens: true,
        minimize: false,
        charges: false,
        ..PrepareOptions::default()
    };
    let report = prepare(&calculator, &structure, &options).unwrap();
    assert_eq!(report.structure.sites().len(), 5);
    assert!(
        report
            .steps
            .iter()
            .any(|step| step.name == "hydrogens" && step.done)
    );
}

#[test]
fn cleanup_lists_close_contacts() {
    let structure = cell_with(vec![
        Site::new("H1", "H", [0.0, 0.0, 0.0]).unwrap(),
        Site::new("H2", "H", [0.08, 0.0, 0.0]).unwrap(),
    ]);
    let calculator = UffCalculator::with_defaults().unwrap();
    let options = PrepareOptions {
        check_close_contacts: true,
        minimize: false,
        charges: false,
        ..PrepareOptions::default()
    };
    let report = prepare(&calculator, &structure, &options).unwrap();
    assert!(!report.close_contacts.is_empty());
}

#[test]
fn graph_planner_keeps_a_connected_disorder_unit() {
    let structure = cell_with(vec![
        Site::new("O1", "O", [0.0, 0.0, 0.0])
            .unwrap()
            .with_occupancy(0.6)
            .unwrap(),
        Site::new("O2", "O", [0.18, 0.0, 0.0])
            .unwrap()
            .with_occupancy(0.4)
            .unwrap(),
    ]);
    let calculator = UffCalculator::with_defaults().unwrap();
    let mut options = cleanup_options();
    options.remove_solvent = false;
    let report = prepare(&calculator, &structure, &options).unwrap();
    assert_eq!(report.removed_sites, 1);
    assert_eq!(report.structure.sites().len(), 1);
}

#[test]
fn disorder_choice_follows_the_declared_formula() {
    let structure = xtalforge_core::StructureBuilder::new(
        "formula",
        xtalforge_core::Lattice::cubic(10.0).unwrap(),
    )
    .site(Site::new("C1", "C", [0.5, 0.5, 0.5]).unwrap())
    .site(
        Site::new("O1", "O", [0.0, 0.0, 0.0])
            .unwrap()
            .with_occupancy(0.6)
            .unwrap()
            .with_disorder_group(1),
    )
    .site(
        Site::new("N1", "N", [0.18, 0.0, 0.0])
            .unwrap()
            .with_occupancy(0.4)
            .unwrap()
            .with_disorder_group(2),
    )
    .metadata("chemical_formula", serde_json::json!("C1 N1"))
    .build()
    .unwrap();
    let calculator = UffCalculator::with_defaults().unwrap();
    let mut options = cleanup_options();
    options.remove_solvent = false;
    let report = prepare(&calculator, &structure, &options).unwrap();
    assert_eq!(report.structure.sites().len(), 2);
    assert!(
        report
            .structure
            .sites()
            .iter()
            .any(|site| site.element() == "N")
    );
    assert!(
        report
            .structure
            .sites()
            .iter()
            .all(|site| site.element() != "O")
    );
}

#[test]
fn disorder_choice_prefers_a_neutral_component() {
    let structure = cell_with(vec![
        Site::new("C1", "C", [0.5, 0.5, 0.5]).unwrap(),
        Site::new("O1", "O", [0.0, 0.0, 0.0])
            .unwrap()
            .with_occupancy(0.6)
            .unwrap()
            .with_charge(1.0)
            .unwrap()
            .with_disorder_group(1),
        Site::new("N1", "N", [0.18, 0.0, 0.0])
            .unwrap()
            .with_occupancy(0.4)
            .unwrap()
            .with_disorder_group(2),
    ]);
    let calculator = UffCalculator::with_defaults().unwrap();
    let mut options = cleanup_options();
    options.remove_solvent = false;
    let report = prepare(&calculator, &structure, &options).unwrap();
    assert_eq!(report.structure.sites().len(), 2);
    assert!(
        report
            .structure
            .sites()
            .iter()
            .any(|site| site.element() == "N")
    );
}
