//! Integration tests for schwarzite level-set meshes (Task 6.8).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_build::{ZtcOptions, schwarzite_surface_mesh, ztc_surface_mesh};

fn small() -> ZtcOptions {
    ZtcOptions {
        samples: 8,
        cell: 12.0,
        isovalue: 0.0,
    }
}

#[test]
fn gyroid_and_schwarz_p_meshes_are_not_empty() {
    let options = small();
    let gyroid = ztc_surface_mesh(&options).unwrap();
    assert!(!gyroid.is_empty());
    assert!(!gyroid.vertices.is_empty());
    assert_eq!(gyroid.triangles.len() * 3, gyroid.vertices.len());
    let schwarz = schwarzite_surface_mesh(&options).unwrap();
    assert!(!schwarz.is_empty());
    assert!(!schwarz.vertices.is_empty());
}

#[test]
fn meshes_are_deterministic() {
    let options = small();
    assert_eq!(
        ztc_surface_mesh(&options).unwrap(),
        ztc_surface_mesh(&options).unwrap()
    );
    assert_eq!(
        schwarzite_surface_mesh(&options).unwrap(),
        schwarzite_surface_mesh(&options).unwrap()
    );
}

#[test]
fn invalid_samples_and_cell_error() {
    assert!(
        ztc_surface_mesh(&ZtcOptions {
            samples: 1,
            ..small()
        })
        .is_err()
    );
    assert!(
        schwarzite_surface_mesh(&ZtcOptions {
            cell: 0.0,
            ..small()
        })
        .is_err()
    );
    assert!(
        ztc_surface_mesh(&ZtcOptions {
            isovalue: f64::NAN,
            ..small()
        })
        .is_err()
    );
}
