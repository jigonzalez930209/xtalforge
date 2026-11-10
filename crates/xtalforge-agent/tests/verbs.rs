//! The verb catalogue (Task 8.1.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeSet;

use xtalforge_agent::verbs::{ParameterKind, VERBS, catalogue, lookup, names};

#[test]
fn names_are_unique_and_findable() {
    let mut seen = BTreeSet::new();
    for verb in catalogue() {
        assert!(seen.insert(verb.name), "duplicate verb {}", verb.name);
        assert!(lookup(verb.name).is_some());
    }
    assert_eq!(names().len(), VERBS.len());
}

#[test]
fn examples_are_json_objects_that_only_name_declared_parameters() {
    for verb in catalogue() {
        let value: serde_json::Value = serde_json::from_str(verb.example).expect("example");
        let object = value.as_object().expect("example must be an object");
        for key in object.keys() {
            assert!(
                verb.parameters
                    .iter()
                    .any(|parameter| parameter.name == key),
                "{}.{} is in the example but not in the catalogue",
                verb.name,
                key
            );
        }
    }
}

#[test]
fn required_parameters_have_no_default() {
    for verb in catalogue() {
        for parameter in verb.parameters {
            if parameter.required {
                assert!(
                    parameter.default.is_none(),
                    "{}.{} is required and has a default",
                    verb.name,
                    parameter.name
                );
            }
        }
    }
}

#[test]
fn every_kind_is_used_by_at_least_one_verb() {
    let kinds: BTreeSet<ParameterKind> = catalogue()
        .iter()
        .flat_map(|verb| verb.parameters.iter().map(|parameter| parameter.kind))
        .collect();
    for kind in [
        ParameterKind::String,
        ParameterKind::Integer,
        ParameterKind::Number,
        ParameterKind::Array,
    ] {
        assert!(kinds.contains(&kind), "{kind:?} is declared but unused");
    }
}

#[test]
fn mutations_are_marked() {
    for name in ["supercell", "wrap", "hydrogens", "optimize", "undo", "redo"] {
        assert!(lookup(name).expect("verb").mutation, "{name} must mutate");
    }
    for name in ["info", "diagnostics", "symmetry", "bonds", "pxrd"] {
        assert!(
            !lookup(name).expect("verb").mutation,
            "{name} must not mutate"
        );
    }
}
