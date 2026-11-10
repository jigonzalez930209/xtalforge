//! Atomic writes to destinations that are not regular files.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::path::Path;

#[cfg(unix)]
#[test]
fn writing_to_a_device_does_not_replace_it() {
    let null = Path::new("/dev/null");
    xtalforge_io::atomic::write_text(null, "discarded").unwrap();
    // The node is still a character device and not a regular file.
    let metadata = std::fs::metadata(null).unwrap();
    assert!(!metadata.is_file());
    assert!(!metadata.is_dir());
}

#[test]
fn writing_to_a_directory_is_an_error() {
    let directory = std::env::temp_dir();
    let error = xtalforge_io::atomic::write_text(&directory, "x").unwrap_err();
    assert_eq!(error.code(), "invalid_input");
}
