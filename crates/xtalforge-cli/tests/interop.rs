//! CLI interoperability surface (Task 13.2/13.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::process::{Command, Output};

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args(args)
        .output()
        .expect("the CLI binary must run")
}

#[test]
fn formats_lists_the_phase_13_entries() {
    let output = run(&["formats"]);
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).unwrap();
    for name in [
        "mmcif", "extxyz", "pwscf", "castep", "fhi-aims", "cp2k", "abinit", "crystal", "mol2",
        "sdf", "outcar", "vasprun", "optimade",
    ] {
        assert!(text.contains(name), "missing format {name}");
    }
}

#[test]
fn convert_reads_a_castep_cell() {
    let directory =
        std::env::temp_dir().join(format!("xtalforge-cli-interop-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let input = directory.join("Si.cell");
    std::fs::write(
        &input,
        "%BLOCK LATTICE_CART\n5.43 0 0\n0 5.43 0\n0 0 5.43\n%ENDBLOCK LATTICE_CART\n\
         %BLOCK POSITIONS_FRAC\nSi 0 0 0\n%ENDBLOCK POSITIONS_FRAC\n",
    )
    .unwrap();
    let output = run(&["convert", input.to_str().unwrap(), "--format", "cif"]);
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).unwrap();
    assert!(text.contains("_cell_length_a"));
    assert!(text.contains("Si1"));
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn convert_writes_an_optimade_entry() {
    let directory =
        std::env::temp_dir().join(format!("xtalforge-cli-optimade-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let input = directory.join("Si.cif");
    std::fs::write(
        &input,
        "data_Si\n_cell_length_a 5.43\n_cell_length_b 5.43\n_cell_length_c 5.43\n\
         _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
         _symmetry_space_group_name_H-M 'F d -3 m'\nloop_\n_atom_site_label\n\
         _atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n\
         _atom_site_fract_z\nSi1 Si 0 0 0\n",
    )
    .unwrap();
    let output = run(&["convert", input.to_str().unwrap(), "--format", "optimade"]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["type"], "structures");
    assert!(value["attributes"]["lattice_vectors"].is_array());
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn mesh_writes_obj_and_gltf() {
    let cif = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/corpus/simple/Si.cif"
    );
    let obj = run(&["mesh", cif, "--mesh-format", "obj"]);
    assert!(obj.status.success());
    let text = String::from_utf8(obj.stdout).unwrap();
    assert!(text.contains("v "));
    assert!(text.contains("vn "));
    assert!(text.contains("f "));

    let gltf = run(&["mesh", cif, "--mesh-format", "gltf"]);
    assert!(gltf.status.success());
    let value: serde_json::Value = serde_json::from_slice(&gltf.stdout).unwrap();
    assert_eq!(value["asset"]["version"], "2.0");
    assert!(value["meshes"].is_array());
}

#[test]
fn modules_lists_the_builtin_registry() {
    let output = run(&["--json", "modules"]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let items = value.as_array().unwrap();
    assert!(items.iter().any(|item| item["name"] == "version"));
}
