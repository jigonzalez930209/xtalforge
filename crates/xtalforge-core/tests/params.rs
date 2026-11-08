//! Integration tests for typed parameter descriptions (Task 3.8.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::json;
use xtalforge_core::{ParamSet, ParamSpec, ParamType};

#[test]
fn parameters_validate_and_expose_a_schema() {
    let set = ParamSet::new()
        .with_param(ParamSpec {
            name: "symprec".to_owned(),
            kind: ParamType::Float {
                min: Some(0.0),
                max: Some(1.0),
                unit: Some("Angstrom".to_owned()),
            },
            default: json!(0.001),
            help: "Position tolerance".to_owned(),
            required: false,
            group: Some("Symmetry".to_owned()),
        })
        .unwrap()
        .with_param(ParamSpec {
            name: "mode".to_owned(),
            kind: ParamType::Choice {
                choices: vec!["detect".to_owned(), "impose".to_owned()],
            },
            default: json!("detect"),
            help: "Operation mode".to_owned(),
            required: true,
            group: None,
        })
        .unwrap();
    let schema = set.schema();
    assert_eq!(schema["properties"]["symprec"]["type"], "number");
    assert_eq!(schema["required"][0], "mode");
    let mut set = set;
    assert!(set.set("symprec", json!(0.5)).is_ok());
    assert!(set.set("symprec", json!(2.0)).is_err());
    assert!(set.set("mode", json!("bad")).is_err());
    assert!(set.set("unknown", json!(1)).is_err());
    assert!(set.check_required().is_ok());
}

#[test]
fn duplicate_and_invalid_defaults_are_rejected() {
    let duplicate = ParamSet::new()
        .with_param(ParamSpec {
            name: "x".to_owned(),
            kind: ParamType::Boolean,
            default: json!(true),
            help: String::new(),
            required: false,
            group: None,
        })
        .unwrap()
        .with_param(ParamSpec {
            name: "x".to_owned(),
            kind: ParamType::Boolean,
            default: json!(false),
            help: String::new(),
            required: false,
            group: None,
        });
    assert!(duplicate.is_err());
    let invalid = ParamSet::new().with_param(ParamSpec {
        name: "x".to_owned(),
        kind: ParamType::Integer {
            min: Some(1),
            max: None,
        },
        default: json!(0),
        help: String::new(),
        required: false,
        group: None,
    });
    assert!(invalid.is_err());
}
