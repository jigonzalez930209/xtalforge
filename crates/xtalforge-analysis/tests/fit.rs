//! Integration tests for the least-squares engine (Task 7.4.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{FitOptions, least_squares};

#[test]
fn linear_fit_recovers_parameters() {
    let result = least_squares(&[0.0, 0.0], &FitOptions::default(), |parameters| {
        let [slope, intercept] = parameters else {
            return Err(xtalforge_analysis::Error::InvalidInput(
                "two parameters".to_owned(),
            ));
        };
        Ok((0..10)
            .map(|index| {
                let x = index as f64;
                2.5 * x + 1.0 - (slope * x + intercept)
            })
            .collect())
    })
    .unwrap();
    assert!((result.parameters[0] - 2.5).abs() < 1e-5);
    assert!((result.parameters[1] - 1.0).abs() < 1e-5);
    assert!(result.chi_squared < 1e-8);
    assert!(result.converged);
}

#[test]
fn invalid_inputs_are_rejected() {
    assert!(least_squares(&[], &FitOptions::default(), |_| Ok(vec![1.0])).is_err());
    assert!(least_squares(&[1.0], &FitOptions::default(), |_| Ok(Vec::new())).is_err());
    assert!(
        least_squares(
            &[1.0],
            &FitOptions {
                tolerance: 0.0,
                ..Default::default()
            },
            |_| Ok(vec![1.0])
        )
        .is_err()
    );
}
