//! Parser robustness (Task 8.6.3).
//!
//! A deterministic mutation loop over every text parser and a set of malformed
//! inputs. A parser may accept or reject a mutated file, but it must not
//! panic: the test fails on a panic. This is the reproducible stand-in for a
//! coverage-guided fuzzer, which needs nightly and a separate runner.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_io::{
    read_cif, read_cube, read_gen, read_lammps, read_pmg_json, read_poscar, read_project,
    read_shelx, read_xyz, write_gen, write_lammps, write_pmg_json, write_poscar, write_shelx,
    write_xyz,
};

const NACL: &str = include_str!("../../../tests/corpus/simple/NaCl.cif");

fn structure() -> xtalforge_core::Structure {
    read_cif(NACL).expect("NaCl")
}

fn next(seed: &mut u64) -> u64 {
    *seed = seed
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    *seed
}

fn mutate(text: &str, seed: &mut u64) -> String {
    let mut bytes = text.as_bytes().to_vec();
    if bytes.is_empty() {
        return String::new();
    }
    for _ in 0..1 + next(seed) % 4 {
        let index = (next(seed) as usize) % bytes.len();
        bytes[index] = (next(seed) % 128) as u8;
    }
    String::from_utf8_lossy(&bytes).into_owned()
}

#[test]
fn malformed_inputs_are_rejected_without_panicking() {
    let cif = [
        "",
        "data_",
        "data_x\n_cell_length_a nope\n",
        "data_x\nloop_\n_atom_site_label\n",
        "data_x\n_cell_length_a 5\n_cell_length_b -1\n",
        "data_x\n_symmetry_space_group_name_H-M 'not a group'\n",
    ];
    for text in cif {
        let _ = read_cif(text);
    }
    let poscar = [
        "",
        "NaCl\n1\n",
        "NaCl\n1.0\n",
        "NaCl\n1.0\n1 1 1\n",
        "x\n1\n5 5 5\n",
    ];
    for text in poscar {
        let _ = read_poscar(text);
    }
    let xyz = [
        "",
        "3\nwater\nO 0 0 0\n",
        "2\nx\nXx 0 0 0\nH 1 0 0\n",
        "abc\n",
    ];
    for text in xyz {
        let _ = read_xyz(text);
    }
    for text in [
        "",
        "cube",
        "cube 1 0 0",
        "cube\n1 0 0 0 0\n0 0 0\n1 0 0 0\n",
    ] {
        let _ = read_cube(text);
    }
    for text in ["", "ATOM", "HETATM", "ATOM 1\n"] {
        let _ = read_shelx(text);
    }
    for text in ["", "[", "{\"sites\": 1}", "null"] {
        let _ = read_pmg_json(text);
    }
    for text in ["", "LAMMPS", "LAMMPS data\n\n0 atoms\n"] {
        let _ = read_lammps(text);
    }
    for text in ["", "gen", "2 C\n"] {
        let _ = read_gen(text);
    }
    for bytes in [&b""[..], b"not a zip", &[0x50, 0x4b, 0x03, 0x04][..]] {
        let _ = read_project(bytes);
    }
}

#[test]
fn mutated_cif_poscar_and_xyz_never_panic() {
    let structure = structure();
    let poscar = write_poscar(&structure).expect("POSCAR");
    let xyz = write_xyz(&structure).expect("XYZ");
    let gen_text = write_gen(&structure).expect("GEN");
    let lammps = write_lammps(&structure).expect("LAMMPS");
    let shelx = write_shelx(&structure).expect("SHELX");
    let pmg = write_pmg_json(&structure).expect("pmsg JSON");
    let mut seed = 0x5eed_1234_5678_9abc_u64;
    for _ in 0..200 {
        let _ = read_cif(&mutate(NACL, &mut seed));
        let _ = read_poscar(&mutate(&poscar, &mut seed));
        let _ = read_xyz(&mutate(&xyz, &mut seed));
        let _ = read_gen(&mutate(&gen_text, &mut seed));
        let _ = read_lammps(&mutate(&lammps, &mut seed));
        let _ = read_shelx(&mutate(&shelx, &mut seed));
        let _ = read_pmg_json(&mutate(&pmg, &mut seed));
    }
}

#[test]
fn a_mutated_project_archive_never_panics() {
    let structure = structure();
    let parts =
        xtalforge_io::ProjectParts::from_structure(&structure, None, None, None, None, None)
            .unwrap();
    let bytes = xtalforge_io::write_project(&parts).expect("project");
    let mut seed = 0x1234_5678_9abc_def0_u64;
    for _ in 0..100 {
        let mut copy = bytes.clone();
        if !copy.is_empty() {
            let index = (next(&mut seed) as usize) % copy.len();
            copy[index] ^= (next(&mut seed) & 0xff) as u8;
        }
        let _ = read_project(&copy);
    }
}
