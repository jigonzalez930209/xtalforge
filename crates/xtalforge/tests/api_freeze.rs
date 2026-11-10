//! Representative public signatures of the 1.0 surface (Task 8.8.1).
//!
//! These are compile-time assertions. Changing one of these signatures
//! without a major-version decision fails this test, which is the point: the
//! facade is the stable API, and a rename or a new argument is a breaking
//! change on every surface.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use xtalforge::agent::{Answer, DiagnosticReport, Session};
use xtalforge::core::{Error, Lattice, Structure, StructureSummary};
use xtalforge::view::ViewSettings;

#[test]
fn the_core_hello_world_surface_is_frozen() {
    let version: fn() -> &'static str = xtalforge::core::version;
    let volume: fn(f64, f64, f64, f64, f64, f64) -> Result<Lattice, Error> =
        xtalforge::core::lattice_info;
    let summary: fn(&Structure) -> Result<StructureSummary, Error> = Structure::summary;
    assert!(!version().is_empty());
    let lattice = volume(5.0, 5.0, 5.0, 90.0, 90.0, 90.0).unwrap();
    assert_eq!(lattice.volume(), 125.0);
    let _ = summary;
}

#[test]
fn the_io_and_view_surface_is_frozen() {
    let read: fn(&str) -> Result<Structure, Error> = xtalforge::io::read_cif;
    let write: fn(&Structure) -> Result<String, Error> = xtalforge::io::write_cif;
    let scene: fn(&Structure, &ViewSettings) -> Result<xtalforge::view::SceneData, Error> =
        xtalforge::view::build_scene;
    let text = "data_x\n_cell_length_a 5\n_cell_length_b 5\n_cell_length_c 5\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                _symmetry_space_group_name_H-M 'P 1'\nloop_\n_atom_site_label\n\
                _atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n\
                _atom_site_fract_z\nC1 C 0 0 0\n";
    let structure = read(text).unwrap();
    let written = write(&structure).unwrap();
    assert!(written.contains("_cell_length_a"));
    let _ = scene(&structure, &ViewSettings::default()).unwrap();
}

#[test]
fn the_agent_surface_is_frozen() {
    let diagnose: fn(&Structure) -> Result<DiagnosticReport, xtalforge::agent::Error> =
        xtalforge::agent::diagnose;
    let batch: fn(&str) -> Result<Vec<Answer>, xtalforge::agent::Error> =
        xtalforge::agent::run_verbs;
    let catalogue: fn() -> &'static [xtalforge::agent::VerbInfo] = xtalforge::agent::catalogue;
    let text = "data_x\n_cell_length_a 5\n_cell_length_b 5\n_cell_length_c 5\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                _symmetry_space_group_name_H-M 'P 1'\nloop_\n_atom_site_label\n\
                _atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n\
                _atom_site_fract_z\nC1 C 0 0 0\n";
    let structure = xtalforge::io::read_cif(text).unwrap();
    let report = diagnose(&structure).unwrap();
    assert!(report.counts.contains_key("error"));
    assert!(!catalogue().is_empty());
    let script = "[]";
    assert!(batch(script).unwrap().is_empty());

    let mut session = Session::new();
    let answer = session.run("capabilities", &serde_json::json!({})).unwrap();
    assert_eq!(answer.kind(), "capabilities");
    assert_eq!(
        answer.to_json().unwrap().get("kind"),
        Some(&serde_json::json!("capabilities"))
    );
}

#[test]
fn the_error_codes_are_frozen() {
    assert_eq!(Error::InvalidInput(String::new()).code(), "invalid_input");
    assert_eq!(xtalforge::agent::Error::NoDocument.code(), "no_document");
    assert_eq!(
        xtalforge::ff::Error::UnknownCalculator(String::new()).code(),
        "unknown_calculator"
    );
}
