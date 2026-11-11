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
    // (i, j) = (2, 1): 21.787 degrees and 14 atoms per layer.
    let twisted = commensurate_twist(1, 1, 21.79, 2).unwrap();
    assert_eq!(twisted.sites().len(), 28);
    assert_eq!(
        twisted.metadata()["twist_indices"],
        serde_json::json!([2, 1])
    );
    let angle = twisted.metadata()["twist_angle"].as_f64().unwrap();
    assert!((angle - 21.786_789).abs() < 1e-5, "{angle}");
    // (3, 2): 13.174 degrees, 38 atoms per layer, repeated 2 x 1.
    let small = commensurate_twist(2, 1, 13.17, 2).unwrap();
    assert_eq!(small.sites().len(), 2 * 2 * 38);
    // Every atom keeps three in-plane neighbours at the graphene bond length.
    let p1 = twisted.p1().unwrap();
    let lattice = twisted.lattice();
    for (index, atom) in p1.atoms.iter().enumerate() {
        let close = p1
            .atoms
            .iter()
            .enumerate()
            .filter(|(other, b)| {
                *other != index && (b.frac[2] - atom.frac[2]).abs() < 1e-9 && {
                    let d = lattice.distance(atom.frac, b.frac);
                    (d - 1.42).abs() < 0.01
                }
            })
            .count();
        assert_eq!(close, 3, "atom {index}");
    }
    let sixty = commensurate_twist(1, 1, 58.0, 2).unwrap();
    let sixty_angle = sixty.metadata()["twist_angle"].as_f64().unwrap();
    assert!((sixty_angle - 60.0).abs() < 1e-9);
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
