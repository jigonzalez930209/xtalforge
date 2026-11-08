//! Tests for twists, spirals, ZTC and polymers (Tasks 6.5.3, 6.7.1, 6.8, 6.9).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_build::{
    PolymerOptions, ZtcOptions, build_polymer, build_schwarzite, build_ztc, commensurate_twist,
    face_spiral, nanoribbon, pack_chains, polymer_protocol, schwarzite_surface_mesh,
    spiral_sequence, ztc_surface_mesh,
};

#[test]
fn twists_and_ribbons_build() {
    let flat = commensurate_twist(2, 2, 0.0, 2).unwrap();
    assert_eq!(flat.sites().len(), 16);
    let twisted = commensurate_twist(2, 2, 58.0, 2).unwrap();
    assert_eq!(twisted.sites().len(), 16);
    assert_eq!(twisted.metadata()["twist_approximation"], 60.0);
    let ribbon = nanoribbon(2, 3, 12.0).unwrap();
    assert!(!ribbon.sites().is_empty());
    assert!(ribbon.lattice().lengths()[1] > 12.0);
    assert!(commensurate_twist(0, 2, 10.0, 1).is_err());
    assert!(nanoribbon(0, 2, 10.0).is_err());
}

#[test]
fn spirals_cover_c20_and_c60() {
    let c20 = face_spiral(20).unwrap();
    assert_eq!(c20.sites().len(), 20);
    assert_eq!(spiral_sequence(20).unwrap().len(), 12);
    let c60 = face_spiral(60).unwrap();
    assert_eq!(c60.sites().len(), 60);
    assert!(face_spiral(70).is_err());
}

#[test]
fn ztc_surfaces_build_and_mesh() {
    let options = ZtcOptions {
        samples: 8,
        cell: 12.0,
        isovalue: 0.0,
    };
    let ztc = build_ztc(&options).unwrap();
    assert!(!ztc.sites().is_empty());
    let schwarz = build_schwarzite(&options).unwrap();
    assert!(!schwarz.sites().is_empty());
    assert!(!ztc_surface_mesh(&options).unwrap().is_empty());
    assert!(!schwarzite_surface_mesh(&options).unwrap().is_empty());
    assert!(
        ztc_surface_mesh(&ZtcOptions {
            samples: 1,
            ..options
        })
        .is_err()
    );
    assert!(
        ztc_surface_mesh(&ZtcOptions {
            cell: 0.0,
            ..options
        })
        .is_err()
    );
}

#[test]
fn polymer_chains_pack_and_record_steps() {
    let options = PolymerOptions { density: 1.0 };
    let build = build_polymer("C=C", 3, &options).unwrap();
    assert!(!build.structure.sites().is_empty());
    assert!(!build.steps.is_empty());
    let chain = build_polymer("C=C", 2, &options).unwrap().structure;
    let packed = pack_chains(&[chain], 0.8).unwrap();
    assert!(!packed.sites().is_empty());
    assert!(pack_chains(&[], 0.8).is_err());
    assert!(!polymer_protocol().unwrap().is_empty());
}
