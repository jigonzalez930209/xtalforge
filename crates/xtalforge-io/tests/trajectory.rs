//! Integration tests for trajectory indexing and XDATCAR (Task 3.6).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_io::{XyzIndex, read_trajectory, read_xdatcar};

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
