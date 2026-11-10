//! Runtime format registry (Task 8.3.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_io::registry::{Capabilities, FORMATS, FormatRegistry, FormatSpec};
use xtalforge_io::{read_xyz, write_xyz};

fn spec(name: &'static str, extensions: &'static [&'static str]) -> FormatSpec {
    FormatSpec {
        name,
        extensions,
        read: Some(read_xyz),
        write: Some(write_xyz),
        capabilities: Capabilities {
            cell: true,
            labels: true,
            ..Capabilities::default()
        },
    }
}

#[test]
fn a_registry_starts_with_the_builtins() {
    let registry = FormatRegistry::new();
    assert_eq!(registry.formats().len(), FORMATS.len());
    assert!(registry.find("cif").is_some());
    assert!(registry.find(".cif").is_some());
    assert!(registry.find("nope").is_none());
    assert!(FormatRegistry::empty().formats().is_empty());
}

#[test]
fn a_plugin_format_reads_and_writes() {
    let mut registry = FormatRegistry::new();
    registry.register(spec("xyz-copy", &[".xyzcopy"])).unwrap();
    let text = "1\nNa\nNa 0 0 0\n";
    let structure = registry.read_with(text, ".xyzcopy").unwrap();
    assert_eq!(structure.sites().len(), 1);
    let written = registry.write_with(&structure, "xyz-copy").unwrap();
    assert!(written.contains("Na"));
    assert!(registry.find("xyzcopy").is_some());
}

#[test]
fn duplicates_and_invalid_specs_are_refused() {
    let mut registry = FormatRegistry::new();
    assert!(registry.register(spec("cif", &[".cif2"])).is_err());
    assert!(registry.register(spec("new", &[".cif"])).is_err());
    assert!(registry.register(spec(" ", &[".x"])).is_err());
    assert!(registry.register(spec("no-extension", &[])).is_err());
    let read_only = FormatSpec {
        name: "read-only",
        extensions: &[".ro"],
        read: Some(read_xyz),
        write: None,
        capabilities: Capabilities::default(),
    };
    registry.register(read_only).unwrap();
    assert!(
        registry
            .write_with(
                &xtalforge_io::read_xyz("1\nNa\nNa 0 0 0\n").unwrap(),
                "read-only"
            )
            .is_err()
    );
}
