//! Integration tests for trajectory indexing and XDATCAR (Task 3.6).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_io::{
    XyzIndex, read_geo_end, read_lammps_dump, read_lammps_dump_with_types, read_trajectory,
    read_xdatcar,
};

const XYZ: &str = "2\nLattice=\"4 0 0 0 4 0 0 0 4\" energy=-1.5 step=1\n\
                   C 0 0 0\nH 1 0 0\n\
                   2\nLattice=\"4 0 0 0 4 0 0 0 4\" energy=-1.25 step=2\n\
                   C 0.1 0 0\nH 1.1 0 0\n";

const XDATCAR: &str = "Si C\n1.0\n4 0 0\n0 4 0\n0 0 4\nSi C\n1 1\n\
                       Direct configuration= 1\n0 0 0\n0.5 0.5 0.5\n\
                       Direct configuration= 2\n0.1 0 0\n0.6 0.5 0.5\n";

#[test]
fn xyz_index_reads_frames_on_demand() {
    let index = XyzIndex::new(XYZ).unwrap();
    assert_eq!(index.len(), 2);
    assert_eq!(index.elements(), ["C", "H"]);
    assert_eq!(index.energies(), vec![Some(-1.5), Some(-1.25)]);
    assert_eq!(index.steps(), vec![Some(1), Some(2)]);
    let first = index.frame(0).unwrap();
    assert_eq!(first.frac[0], [0.0, 0.0, 0.0]);
    let second = index.structure(1).unwrap();
    assert_eq!(second.sites().len(), 2);
    assert!((second.sites()[0].frac()[0] - 0.025).abs() < 1e-12);
    assert!(index.frame(2).is_err());
}

#[test]
fn xdatcar_reads_direct_frames() {
    let trajectory = read_xdatcar(XDATCAR).unwrap();
    assert_eq!(trajectory.elements, ["Si", "C"]);
    assert_eq!(trajectory.frames.len(), 2);
    assert_eq!(trajectory.frames[0].frac[1], [0.5, 0.5, 0.5]);
    assert_eq!(trajectory.frames[1].frac[0], [0.1, 0.0, 0.0]);
    let structure = trajectory.structure(0).unwrap();
    assert_eq!(structure.sites().len(), 2);
}

#[test]
fn read_trajectory_detects_the_format() {
    assert_eq!(read_trajectory(XYZ).unwrap().frames.len(), 2);
    assert_eq!(read_trajectory(XDATCAR).unwrap().frames.len(), 2);
}

#[test]
fn a_bad_frame_is_an_error() {
    assert!(XyzIndex::new("2\ncomment\nC 0 0 0\n").is_err());
    assert!(read_xdatcar("only a header\n").is_err());
}

const LAMMPS: &str = "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n\
                      ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n\
                      ITEM: ATOMS id type element x y z\n\
                      1 1 C 0 0 0\n2 1 C 1 0 0\n\
                      ITEM: TIMESTEP\n1\nITEM: NUMBER OF ATOMS\n2\n\
                      ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n\
                      ITEM: ATOMS id type element x y z\n\
                      1 1 C 0.5 0 0\n2 1 C 1.5 0 0\n";

const GEO_END: &str = "2\nenergy=-1.0\nC 0 0 0\nH 1 0 0\n\
                       2\nenergy=-0.9\nC 0.1 0 0\nH 1.1 0 0\n";

#[test]
fn lammps_dump_reads_cartesian_frames() {
    let trajectory = read_lammps_dump(LAMMPS).unwrap();
    assert_eq!(trajectory.elements, ["C", "C"]);
    assert_eq!(trajectory.frames.len(), 2);
    assert_eq!(trajectory.frames[0].step, Some(0));
    assert_eq!(trajectory.frames[1].step, Some(1));
    assert!((trajectory.frames[0].frac[1][0] - 0.1).abs() < 1e-12);
    assert!((trajectory.frames[1].frac[0][0] - 0.05).abs() < 1e-12);
    let structure = trajectory.structure(0).unwrap();
    assert_eq!(structure.sites().len(), 2);
}

#[test]
fn lammps_dump_maps_numeric_types() {
    let text = LAMMPS
        .replace("type element", "type")
        .replace("1 1 C", "1 1")
        .replace("2 1 C", "2 2");
    let trajectory = read_lammps_dump_with_types(&text, &["C", "H"]).unwrap();
    assert_eq!(trajectory.elements, ["C", "H"]);
    assert!(read_lammps_dump(&text).is_err());
}

#[test]
fn geo_end_boxes_plain_xyz() {
    let trajectory = read_geo_end(GEO_END).unwrap();
    assert_eq!(trajectory.frames.len(), 2);
    assert!(trajectory.lattice.is_some());
    let structure = trajectory.structure(0).unwrap();
    assert_eq!(structure.sites().len(), 2);
}

#[test]
fn read_trajectory_detects_lammps_and_geo_end() {
    assert_eq!(read_trajectory(LAMMPS).unwrap().frames.len(), 2);
    assert_eq!(read_trajectory(GEO_END).unwrap().frames.len(), 2);
}
