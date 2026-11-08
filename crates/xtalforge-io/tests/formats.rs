//! Integration tests for POSCAR, XYZ and format dispatch (Task 1.16).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_io::{
    read_poscar, read_structure, read_xyz, write_poscar, write_structure, write_xyz,
};

const POSCAR: &str =
    "Si\n1.0\n4.0 0.0 0.0\n0.0 4.0 0.0\n0.0 0.0 4.0\nSi\n2\nDirect\n0.0 0.0 0.0\n0.25 0.25 0.25\n";

#[test]
fn poscar_round_trip() {
    let structure = read_poscar(POSCAR).unwrap();
    assert_eq!(structure.formula(), "Si2");
    assert!((structure.lattice().volume() - 64.0).abs() < 1e-9);
    let text = write_poscar(&structure).unwrap();
    let back = read_poscar(&text).unwrap();
    assert_eq!(back.formula(), "Si2");
    assert!((back.lattice().volume() - 64.0).abs() < 1e-9);
}

#[test]
fn negative_scale_means_target_volume() {
    let text = POSCAR.replace("1.0\n", "-64.0\n");
    let structure = read_poscar(&text).unwrap();
    assert!((structure.lattice().volume() - 64.0).abs() < 1e-9);
}

#[test]
fn poscar_without_species_is_refused() {
    let text = POSCAR.replace("Si\n2\n", "2\n");
    assert!(read_poscar(&text).is_err());
}

#[test]
fn extended_xyz_round_trip() {
    let structure = StructureBuilder::new("C", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let text = write_xyz(&structure).unwrap();
    let back = read_xyz(&text).unwrap();
    assert_eq!(back.formula(), "C");
    assert!((back.lattice().volume() - 64.0).abs() < 1e-6);
}

#[test]
fn plain_xyz_gets_a_padded_box() {
    let text = "2\nplain\nH 0.0 0.0 0.0\nH 1.0 0.0 0.0\n";
    let structure = read_xyz(text).unwrap();
    assert_eq!(structure.formula(), "H2");
    assert!((structure.lattice().lengths()[0] - 11.0).abs() < 1e-9);
}

#[test]
fn dispatch_uses_names_and_extensions() {
    let structure = read_structure(POSCAR, ".vasp").unwrap();
    assert_eq!(structure.formula(), "Si2");
    let xyz = write_structure(&structure, "xyz").unwrap();
    assert!(xyz.contains("Lattice="));
    assert!(read_structure(&xyz, "xyz").is_ok());
    assert!(read_structure(POSCAR, "unknown").is_err());
}
