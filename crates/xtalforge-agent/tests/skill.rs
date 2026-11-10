//! Skill documents and their installation (Task 8.2.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_agent::skill::{DOCUMENTS, default_directory, find, install};

#[test]
fn the_catalogue_is_complete_and_unique() {
    let mut names = std::collections::BTreeSet::new();
    for document in DOCUMENTS {
        assert!(names.insert(document.name), "duplicate {}", document.name);
        assert!(
            !document.body.trim().is_empty(),
            "{} is empty",
            document.name
        );
        assert!(
            document.body.starts_with('#'),
            "{} must start with a markdown heading",
            document.name
        );
    }
    assert!(DOCUMENTS.len() >= 12);
}

#[test]
fn documents_are_found_by_path_or_stem() {
    assert_eq!(
        find("SKILL.md").map(|document| document.name),
        Some("SKILL.md")
    );
    assert_eq!(
        find("SKILL").map(|document| document.name),
        Some("SKILL.md")
    );
    assert_eq!(find("/mcp").map(|document| document.name), Some("mcp.md"));
    assert_eq!(
        find("references/porosity.md").map(|document| document.name),
        Some("references/porosity.md")
    );
    assert!(find("nope").is_none());
}

#[test]
fn install_writes_every_document() {
    let mut directory = std::env::temp_dir();
    directory.push(format!("xtalforge-skill-{}", std::process::id()));
    let written = install(&directory).expect("install");
    assert_eq!(written.len(), DOCUMENTS.len());
    for path in &written {
        assert!(path.is_file(), "{}", path.display());
    }
    let skill = std::fs::read_to_string(directory.join("SKILL.md")).expect("SKILL.md");
    assert!(skill.contains("capabilities"));
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn the_default_directory_is_absolute_or_relative() {
    let path = default_directory();
    assert!(!path.as_os_str().is_empty());
}
