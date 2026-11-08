//! Integration tests for the optional in-process quantum adapter (Task 4.10).
//!
//! The backend is not vendored in this workspace, so the adapter must fail with
//! a typed `invalid_input` error while the capability table stays truthful.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use xtalforge_ff::{Error, QuantumCapabilities, QuantumMethod, adapter, capabilities};

#[test]
fn capability_table_matches_the_documented_limits() {
    for method in [
        QuantumMethod::Gfn0,
        QuantumMethod::Gfn1,
        QuantumMethod::Gfn2,
    ] {
        let caps = capabilities(method);
        assert!(caps.periodic, "{} should be periodic", method.name());
        assert!(caps.forces, "{} should return forces", method.name());
        assert!(!caps.stress, "{} has no analytic stress", method.name());
        assert_eq!(
            caps,
            QuantumCapabilities {
                periodic: true,
                forces: true,
                stress: false,
            }
        );
    }
}

#[test]
fn method_names_are_stable() {
    assert_eq!(QuantumMethod::Gfn0.name(), "gfn0");
    assert_eq!(QuantumMethod::Gfn1.name(), "gfn1");
    assert_eq!(QuantumMethod::Gfn2.name(), "gfn2");
}

#[test]
fn adapter_returns_a_typed_error_without_the_backend() {
    for method in [
        QuantumMethod::Gfn0,
        QuantumMethod::Gfn1,
        QuantumMethod::Gfn2,
    ] {
        let error = adapter(method).expect_err("backend is not linked");
        assert_eq!(error.code(), "invalid_input");
        match error {
            Error::InvalidInput(message) => {
                assert!(
                    message.contains("sci-form"),
                    "unexpected message: {message}"
                );
                assert!(
                    message.contains("rebuild") || message.contains("not linked"),
                    "unexpected message: {message}"
                );
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }
}
