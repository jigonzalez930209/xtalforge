//! Integration tests for text decoding (Task 1.15.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_io::decode_text;

#[test]
fn utf8_and_bom_are_decoded() {
    assert_eq!(decode_text(b"hello"), "hello");
    assert_eq!(decode_text(b"\xEF\xBB\xBFhello"), "hello");
}

#[test]
fn latin1_fallback_never_fails() {
    let bytes = b"caf\xE9";
    assert_eq!(decode_text(bytes), "café");
}
