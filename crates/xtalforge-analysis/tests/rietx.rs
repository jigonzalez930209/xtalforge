//! Integration tests for the rietx bridge inventory (Task 7.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{
    ACCEPTANCE_MAX_POSITION_ERROR, ACCEPTANCE_MAX_RWP, AcceptanceMetrics, RIETX_LICENSE_NOTE,
    check_acceptance, rietx_inventory,
};

#[test]
fn inventory_lists_the_bridge_capabilities() {
    let bridge = rietx_inventory();
    assert_eq!(bridge.name, "rietx");
    assert!(!bridge.capabilities.is_empty());
    for capability in [
        "pattern_io",
        "background",
        "profile",
        "le_bail",
        "pawley",
        "rietveld",
        "indexing",
    ] {
        assert!(
            bridge.capabilities.contains(&capability),
            "missing capability {capability}"
        );
    }
    assert!(bridge.license.contains("MIT"));
    assert!(!bridge.source.is_empty());
    assert!(RIETX_LICENSE_NOTE.contains("MIT"));
}

#[test]
fn acceptance_accepts_a_good_fit_and_rejects_bad_metrics() {
    assert!(
        check_acceptance(&AcceptanceMetrics {
            max_rwp: 0.05,
            max_position_error: 0.01,
        })
        .is_ok()
    );
    assert!(
        check_acceptance(&AcceptanceMetrics {
            max_rwp: ACCEPTANCE_MAX_RWP,
            max_position_error: ACCEPTANCE_MAX_POSITION_ERROR,
        })
        .is_ok()
    );
    assert!(
        check_acceptance(&AcceptanceMetrics {
            max_rwp: ACCEPTANCE_MAX_RWP + 0.01,
            max_position_error: 0.01,
        })
        .is_err()
    );
    assert!(
        check_acceptance(&AcceptanceMetrics {
            max_rwp: 0.05,
            max_position_error: ACCEPTANCE_MAX_POSITION_ERROR + 0.01,
        })
        .is_err()
    );
    assert!(
        check_acceptance(&AcceptanceMetrics {
            max_rwp: -0.1,
            max_position_error: 0.01,
        })
        .is_err()
    );
    assert!(
        check_acceptance(&AcceptanceMetrics {
            max_rwp: 0.05,
            max_position_error: f64::NAN,
        })
        .is_err()
    );
}
