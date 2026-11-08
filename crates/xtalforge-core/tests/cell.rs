//! Integration tests for cell volumes, including the shared vectors.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde::Deserialize;
use xtalforge_core::{Error, cell_volume};

#[derive(Debug, Deserialize)]
struct Vectors {
    cell_volume: Vec<Case>,
}

#[derive(Debug, Deserialize)]
struct Case {
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
    expected: f64,
    tolerance: f64,
}

fn vectors() -> Vectors {
    serde_json::from_str(include_str!("../../../tests/vectors/hello_world.json"))
        .expect("shared hello-world vectors must be valid JSON")
}

#[test]
fn shared_cell_volume_vectors() {
    let vectors = vectors();
    assert!(!vectors.cell_volume.is_empty());
    for case in vectors.cell_volume {
        let got = cell_volume(case.a, case.b, case.c, case.alpha, case.beta, case.gamma)
            .expect("valid vector case");
        assert!(
            (got - case.expected).abs() <= case.tolerance,
            "cell_volume({}, {}, {}, {}, {}, {}) = {got}, expected {} +/- {}",
            case.a,
            case.b,
            case.c,
            case.alpha,
            case.beta,
            case.gamma,
            case.expected,
            case.tolerance
        );
    }
}

#[test]
fn cubic_and_hexagonal_cells() {
    assert!((cell_volume(1.0, 1.0, 1.0, 90.0, 90.0, 90.0).unwrap() - 1.0).abs() < 1e-12);
    assert!(
        (cell_volume(1.0, 1.0, 2.0, 90.0, 90.0, 120.0).unwrap() - 1.732_050_807_568_877_6).abs()
            < 1e-12
    );
}

#[test]
fn invalid_cells_are_rejected() {
    let cases = [
        (0.0, 1.0, 1.0, 90.0, 90.0, 90.0),
        (-1.0, 1.0, 1.0, 90.0, 90.0, 90.0),
        (1.0, 1.0, 1.0, 0.0, 90.0, 90.0),
        (1.0, 1.0, 1.0, 180.0, 90.0, 90.0),
        (1.0, 1.0, 1.0, f64::NAN, 90.0, 90.0),
        (1.0, 1.0, 1.0, f64::INFINITY, 90.0, 90.0),
        (1.0, 1.0, 1.0, 180.0, 180.0, 0.0),
    ];
    for (a, b, c, alpha, beta, gamma) in cases {
        let error = cell_volume(a, b, c, alpha, beta, gamma)
            .expect_err("invalid cell must not produce a volume");
        assert_eq!(error.code(), "invalid_input");
    }
}

#[test]
fn degenerate_cell_is_rejected() {
    let error = cell_volume(1.0, 1.0, 1.0, 170.0, 170.0, 170.0)
        .expect_err("degenerate cell must not produce a volume");
    assert_eq!(error.code(), "invalid_input");
    assert!(error.message().contains("degenerate"));
}

#[test]
fn error_type_is_cloneable_and_comparable() {
    let error = cell_volume(0.0, 1.0, 1.0, 90.0, 90.0, 90.0).unwrap_err();
    assert_eq!(error.clone(), error);
    assert!(matches!(error, Error::InvalidInput(_)));
}
