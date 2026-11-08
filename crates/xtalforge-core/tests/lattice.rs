//! Integration tests for the lattice conventions and geometry (Task 1.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::Lattice;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual}, expected={expected}, tolerance={tolerance}"
    );
}

fn close_point(actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    for axis in 0..3 {
        close(actual[axis], expected[axis], tolerance);
    }
}

#[test]
fn cubic_cell_has_the_standard_orientation() {
    let lattice = Lattice::cubic(4.0).unwrap();
    let expected = Lattice::new([[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]]).unwrap();
    assert!(lattice.almost_equal(&expected, 1e-12));
    close(lattice.volume(), 64.0, 1e-12);
    let parameters = lattice.parameters();
    for (actual, expected) in parameters.iter().zip([4.0, 4.0, 4.0, 90.0, 90.0, 90.0]) {
        close(*actual, expected, 1e-12);
    }
    assert!(lattice.is_right_handed());
}

#[test]
fn triclinic_cell_round_trips_its_parameters() {
    let lattice = Lattice::from_parameters(3.0, 4.0, 5.0, 60.0, 70.0, 80.0).unwrap();
    close(lattice.volume(), 48.827_604_320_084_98, 1e-9);
    let parameters = lattice.parameters();
    for (actual, expected) in parameters.iter().zip([3.0, 4.0, 5.0, 60.0, 70.0, 80.0]) {
        close(*actual, expected, 1e-9);
    }
    close(
        lattice.d_spacing([1.0, 1.0, 1.0]).unwrap(),
        2.589_161_046_877_841_4,
        1e-9,
    );
    close_point(
        lattice.plane_normal([1.0, 1.0, 1.0]).unwrap(),
        [
            0.863_053_682_292_613_7,
            0.505_096_107_376_449_4,
            -0.004_032_839_485_793_732,
        ],
        1e-9,
    );
}

#[test]
fn hexagonal_cell_uses_the_120_degree_setting() {
    let lattice = Lattice::from_parameters(3.0, 3.0, 5.0, 90.0, 90.0, 120.0).unwrap();
    close(lattice.volume(), 38.971_143_170_299_754, 1e-9);
    let [a, b, c] = lattice.lengths();
    close(a, 3.0, 1e-12);
    close(b, 3.0, 1e-12);
    close(c, 5.0, 1e-12);
    let angles = lattice.angles();
    close(angles[0], 90.0, 1e-12);
    close(angles[1], 90.0, 1e-12);
    close(angles[2], 120.0, 1e-12);
    close(
        lattice.d_spacing([1.0, 0.0, 0.0]).unwrap(),
        2.598_076_211_353_316,
        1e-9,
    );
}

#[test]
fn coordinate_conversions_are_inverse() {
    let lattice = Lattice::from_parameters(3.0, 4.0, 5.0, 60.0, 70.0, 80.0).unwrap();
    let frac = [0.125, 0.375, 0.625];
    let cart = lattice.to_cart(frac);
    close_point(lattice.to_frac(cart), frac, 1e-12);
    let block = [frac, [0.1, 0.2, 0.3]];
    let carts = lattice.to_cart_many(&block);
    let back = lattice.to_frac_many(&carts);
    close_point(back[0], block[0], 1e-12);
    close_point(back[1], block[1], 1e-12);
}

#[test]
fn wrapping_maps_near_one_to_zero() {
    close_point(Lattice::wrap([1.2, -0.2, 0.5]), [0.2, 0.8, 0.5], 1e-12);
    assert_eq!(
        Lattice::wrap_with_tolerance([0.999_999_99, -1e-9, 0.5], 1e-6),
        [0.0, 0.0, 0.5]
    );
}

#[test]
fn minimum_image_matches_the_reference_values() {
    let lattice = Lattice::from_parameters(3.0, 4.0, 5.0, 60.0, 70.0, 80.0).unwrap();
    let (delta, distance) = lattice.minimum_image([0.9, 0.4, -0.7]);
    close_point(delta, [-0.1, -0.6, 0.3], 1e-12);
    close(distance, 2.107_661_084_436_318, 1e-9);
}

#[test]
fn minimum_image_matches_brute_force_for_a_skewed_cell() {
    let lattice = Lattice::new([[3.0, 0.0, 0.0], [2.8, 1.0, 0.0], [2.8, 0.0, 1.0]]).unwrap();
    let deltas = [
        [0.5, 0.5, 0.5],
        [0.9, 0.4, -0.7],
        [0.1, 0.9, 0.9],
        [-0.9, 0.8, 0.2],
        [0.49, 0.51, 0.5],
    ];
    for delta in deltas {
        let mut best = f64::INFINITY;
        for i in -4..=4 {
            for j in -4..=4 {
                for k in -4..=4 {
                    let candidate = [
                        delta[0] - f64::from(i),
                        delta[1] - f64::from(j),
                        delta[2] - f64::from(k),
                    ];
                    let distance = norm(lattice.to_cart(candidate));
                    if distance < best {
                        best = distance;
                    }
                }
            }
        }
        let (_, distance) = lattice.minimum_image(delta);
        close(distance, best, 1e-12);
    }
}

#[test]
fn reciprocal_lattice_uses_the_crystallographic_convention() {
    let lattice = Lattice::from_parameters(3.0, 4.0, 5.0, 60.0, 70.0, 80.0).unwrap();
    let reciprocal = lattice.reciprocal(false).unwrap();
    let row_a = lattice.matrix()[0];
    let [ra, rb, rc] = reciprocal.matrix();
    close(dot(row_a, ra), 1.0, 1e-12);
    close(dot(row_a, rb), 0.0, 1e-12);
    close(dot(row_a, rc), 0.0, 1e-12);
    let reciprocal_two_pi = lattice.reciprocal(true).unwrap();
    close(
        dot(row_a, reciprocal_two_pi.matrix()[0]),
        2.0 * std::f64::consts::PI,
        1e-12,
    );
}

#[test]
fn transforms_scale_and_reject_bad_matrices() {
    let lattice = Lattice::cubic(2.0).unwrap();
    let supercell = lattice
        .transform([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]])
        .unwrap();
    close(supercell.volume(), 8.0 * lattice.volume(), 1e-12);
    assert!(
        lattice
            .transform([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]])
            .is_err()
    );
    assert!(
        lattice
            .transform([[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
            .is_err()
    );
}

#[test]
fn strain_and_volume_scaling() {
    let lattice = Lattice::cubic(2.0).unwrap();
    let strained = lattice.strained_isotropic(0.01).unwrap();
    close(
        strained.volume(),
        lattice.volume() * 1.01_f64.powi(3),
        1e-12,
    );
    let scaled = lattice.scaled_to_volume(27.0).unwrap();
    close(scaled.volume(), 27.0, 1e-12);
    close(scaled.lengths()[0], 3.0, 1e-12);
    assert!(lattice.scaled_to_volume(0.0).is_err());
}

#[test]
fn validation_reports_the_offending_values() {
    let cases = [
        (0.0, 1.0, 1.0, 90.0, 90.0, 90.0),
        (1.0, 1.0, 1.0, 0.0, 90.0, 90.0),
        (1.0, 1.0, 1.0, 180.0, 90.0, 90.0),
        (1.0, 1.0, 1.0, 170.0, 170.0, 170.0),
        (f64::NAN, 1.0, 1.0, 90.0, 90.0, 90.0),
    ];
    for (a, b, c, alpha, beta, gamma) in cases {
        let error = Lattice::from_parameters(a, b, c, alpha, beta, gamma)
            .expect_err("invalid cell must fail");
        assert_eq!(error.code(), "invalid_input");
    }
    let singular = Lattice::new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]])
        .expect_err("singular matrix must fail");
    assert!(singular.message().contains("singular"));
    assert!(Lattice::from_parameters(f64::INFINITY, 1.0, 1.0, 90.0, 90.0, 90.0).is_err());
}

#[test]
fn serde_round_trip_keeps_the_matrix() {
    let lattice = Lattice::from_parameters(3.0, 4.0, 5.0, 60.0, 70.0, 80.0).unwrap();
    let json = serde_json::to_string(&lattice).unwrap();
    let back: Lattice = serde_json::from_str(&json).unwrap();
    assert!(lattice.almost_equal(&back, 1e-12));
}

#[test]
fn shared_lattice_vectors() {
    let vectors: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/vectors/lattice.json")).unwrap();
    for case in vectors["cases"].as_array().unwrap() {
        let value = |name: &str| case[name].as_f64().unwrap();
        let lattice = Lattice::from_parameters(
            value("a"),
            value("b"),
            value("c"),
            value("alpha"),
            value("beta"),
            value("gamma"),
        )
        .unwrap();
        close(lattice.volume(), value("volume"), 1e-9);
        let delta: Vec<f64> = case["minimum_image"]["delta"]
            .as_array()
            .unwrap()
            .iter()
            .map(|value| value.as_f64().unwrap())
            .collect();
        let (_, distance) = lattice.minimum_image([delta[0], delta[1], delta[2]]);
        close(
            distance,
            case["minimum_image"]["distance"].as_f64().unwrap(),
            1e-9,
        );
    }
}

fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    first[0] * second[0] + first[1] * second[1] + first[2] * second[2]
}

fn norm(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}
