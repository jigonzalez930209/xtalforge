//! Model cache management (Task 8.4.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeSet;

use xtalforge_ff::models::{CATALOGUE, cache_directory, fetch, list, read_record, remove, verify};

fn cache() -> std::path::PathBuf {
    let mut path = std::env::temp_dir();
    path.push(format!("xtalforge-models-{}", std::process::id()));
    path
}

#[test]
fn the_catalogue_is_unique_and_licensed() {
    let mut names = BTreeSet::new();
    for spec in CATALOGUE {
        assert!(names.insert(spec.name), "duplicate {}", spec.name);
        assert!(!spec.license.is_empty(), "{} has no license", spec.name);
        assert!(spec.homepage.starts_with("https://"));
    }
    assert_eq!(CATALOGUE.len(), 3);
    assert!(!cache_directory().as_os_str().is_empty());
}

#[test]
fn fetch_verify_corrupt_and_remove() {
    let cache = cache();
    let _ = std::fs::remove_dir_all(&cache);
    std::fs::create_dir_all(&cache).unwrap();
    let source = cache.join("source.bin");
    std::fs::write(&source, b"not a real model, but a deterministic payload").unwrap();
    let url = format!("file://{}", source.display());

    let statuses = list(&cache).unwrap();
    assert!(statuses.iter().all(|status| status.record.is_none()));

    let record = fetch(&cache, "mace-mp-0", &url).unwrap();
    assert_eq!(record.name, "mace-mp-0");
    assert_eq!(record.sha256.len(), 64);
    assert_eq!(record.license, "MIT");
    assert!(record.bytes > 0);
    assert!(read_record(&cache, "mace-mp-0").unwrap().is_some());
    assert_eq!(verify(&cache, "mace-mp-0").unwrap(), record);

    let statuses = list(&cache).unwrap();
    let mace = statuses
        .iter()
        .find(|status| status.spec.name == "mace-mp-0")
        .unwrap();
    assert_eq!(mace.verified, Some(true));
    assert!(mace.record.is_some());

    let artifact = cache.join("mace-mp-0.model");
    let mut bytes = std::fs::read(&artifact).unwrap();
    bytes.push(0);
    std::fs::write(&artifact, bytes).unwrap();
    assert!(verify(&cache, "mace-mp-0").is_err());

    assert!(remove(&cache, "mace-mp-0").unwrap());
    assert!(!remove(&cache, "mace-mp-0").unwrap());
    let _ = std::fs::remove_dir_all(&cache);
}

#[test]
fn unknown_names_and_empty_urls_are_refused() {
    let cache = cache();
    assert!(fetch(&cache, "nope", "file:///tmp/x").is_err());
    assert!(fetch(&cache, "orb-v3", " ").is_err());
    assert!(verify(&cache, "mattersim").is_err());
    assert!(remove(&cache, "nope").is_err());
}
