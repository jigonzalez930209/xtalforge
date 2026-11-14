//! Phase 19 CLI surfaces.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::process::{Command, Output};

const NACL: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tests/corpus/simple/NaCl.cif"
);
const MCIF: &str = r#"data_mn_test
_cell_length_a 5.0
_cell_length_b 4.0
_cell_length_c 8.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M 'P 1'
_space_group_name_Hall 'P 1'
loop_
_space_group_symop_magn_operation_xyz
'x,y,z,+1'
'-x,-y,z,-1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Mn1 Mn 0.0 0.0 0.0 1.0
Mn2 Mn 0.5 0.5 0.0 0.5
Fe2 Fe 0.5 0.5 0.0 0.5
O1 O 0.25 0.25 0.25 1.0
loop_
_atom_site_moment_label
_atom_site_moment_crystalaxis_x
_atom_site_moment_crystalaxis_y
_atom_site_moment_crystalaxis_z
Mn1 0.0 0.0 4.0
Mn2 0.0 0.0 -4.0
"#;

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args(args)
        .output()
        .expect("the CLI binary must run")
}

fn json(output: Output) -> serde_json::Value {
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).expect("stdout must be JSON")
}

fn write_mcif() -> String {
    let path = std::env::temp_dir().join(format!("xtalforge-phase19-{}.mcif", std::process::id()));
    std::fs::write(&path, MCIF).unwrap();
    path.to_string_lossy().into_owned()
}

#[test]
fn magnetism_disorder_and_derivatives_run() {
    let path = write_mcif();
    let orderings = json(run(&[
        "--json",
        "features",
        "magnetic-orderings",
        &path,
        "--supercell",
        "1",
        "--max",
        "2",
    ]));
    assert!(orderings.as_array().is_some_and(|items| !items.is_empty()));
    let disorder = json(run(&["--json", "features", "disorder", &path]));
    assert!(disorder["groups"].as_array().is_some());
    let derivatives = json(run(&[
        "--json",
        "features",
        "derivatives",
        &path,
        "--max-index",
        "1",
        "--max",
        "2",
    ]));
    assert!(derivatives.as_array().is_some());
    let _ = std::fs::remove_file(path);
}

#[test]
fn neutron_and_electron_patterns_run() {
    let neutron = json(run(&[
        "--json",
        "analysis",
        "neutron",
        NACL,
        "--wavelength",
        "1.54",
        "--step",
        "0.1",
    ]));
    assert!(neutron["reflections"].as_array().is_some());
    let image = std::env::temp_dir().join(format!("xtalforge-saed-{}.pgm", std::process::id()));
    let electron = json(run(&[
        "--json",
        "analysis",
        "electron",
        NACL,
        "--zone",
        "0,0,1",
        "--image",
        &image.to_string_lossy(),
    ]));
    assert!(electron["spots"].as_array().is_some());
    assert!(image.is_file());
    let _ = std::fs::remove_file(image);
}

#[test]
fn spectroscopy_commands_are_registered() {
    for command in ["raman", "ir-born"] {
        let output = run(&["ff", command, "--help"]);
        assert!(output.status.success(), "ff {command} --help");
    }
}
