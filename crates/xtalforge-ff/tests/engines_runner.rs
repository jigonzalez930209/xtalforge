//! DFTB+ and xTB calculators through the process runner (Task 4.9).

#![cfg(unix)]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::PathBuf;
use std::time::Duration;

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_ff::{Calculator, DftbCalculator, DftbOptions, XtbCalculator, XtbOptions};

fn script(name: &str, body: &str) -> PathBuf {
    let path =
        std::env::temp_dir().join(format!("xtalforge-engine-{}-{name}.sh", std::process::id()));
    fs::write(&path, format!("#!/bin/sh\n{body}\n")).unwrap();
    let mut permissions = fs::metadata(&path).unwrap().permissions();
    permissions.set_mode(0o755);
    fs::set_permissions(&path, permissions).unwrap();
    path
}

fn dimer() -> Structure {
    StructureBuilder::new("dimer", Lattice::cubic(10.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("C2", "C", [0.1, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn dftb_calculator_runs_and_parses_detailed_out() {
    let program = script(
        "dftb",
        "test -f dftb_in.hsd || exit 3\n\
         cat > detailed.out <<'OUT'\n\
         Total energy: -4.0 H\n\
         Total Forces\n\
         1 0.1 0.2 0.0\n\
         2 -0.1 0.0 0.3\n\
         OUT\n",
    );
    let calculator = DftbCalculator::new(DftbOptions::default())
        .with_program(program.to_string_lossy())
        .with_timeout(Duration::from_secs(5));
    let (energy, forces) = calculator.energy_forces(&dimer()).unwrap();
    assert!((energy - (-4.0 * xtalforge_ff::HARTREE_TO_EV)).abs() < 1e-9);
    assert_eq!(forces.len(), 2);
    assert!((forces[0][0] - 0.1 * 51.422_067_476_3).abs() < 1e-6);
    fs::remove_file(program).unwrap();
}

#[test]
fn xtb_calculator_runs_and_parses_stdout() {
    let program = script(
        "xtb",
        "test -f coord.xyz || exit 3\n\
         cat <<'OUT'\n\
         TOTAL ENERGY -5.0 Eh\n\
         TOTAL CHARGE 0.0 e\n\
         :: forces (Eh/a0) ::\n\
         1 0.1 0.0 0.0\n\
         2 0.0 -0.2 0.0\n\
         OUT\n",
    );
    let calculator = XtbCalculator::new(XtbOptions::default())
        .with_program(program.to_string_lossy())
        .with_timeout(Duration::from_secs(5));
    let (energy, forces) = calculator.energy_forces(&dimer()).unwrap();
    assert!((energy - (-5.0 * xtalforge_ff::HARTREE_TO_EV)).abs() < 1e-9);
    assert_eq!(forces.len(), 2);
    assert!((forces[1][1] + 0.2 * 51.422_067_476_3).abs() < 1e-6);
    fs::remove_file(program).unwrap();
}

#[test]
fn a_nonzero_exit_is_an_error() {
    let program = script("fail", "echo broken >&2\nexit 7");
    let calculator =
        DftbCalculator::new(DftbOptions::default()).with_program(program.to_string_lossy());
    let error = calculator.energy(&dimer()).unwrap_err();
    assert!(error.message().contains("Failed"));
    fs::remove_file(program).unwrap();
}

#[test]
fn an_unsupported_xtb_method_is_refused() {
    let calculator = XtbCalculator::new(XtbOptions {
        method: "not-a-method".to_owned(),
        ..XtbOptions::default()
    });
    assert!(calculator.energy(&dimer()).is_err());
}
