//! Integration tests for embedded and external data bundles.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::io::Write;

use flate2::Compression;
use flate2::write::GzEncoder;
use serde::Deserialize;
use xtalforge_core::data::{EmbeddedBundle, JsonBundle, load_external_gzip, load_external_json};

const GREETING_SHA256: &str = "318290e761461f3e217f6b62a5c4c8b5412dbd1f28f06ca3d11b410227f10c10";
const ZERO_SHA256: &str = "0000000000000000000000000000000000000000000000000000000000000000";
const GREETING_BYTES: &[u8] = include_bytes!("fixtures/greeting.json.gz");
const GREETING_JSON_BYTES: &[u8] = br#"{"greeting":"hello","answer":42}"#;

static GREETING: EmbeddedBundle = EmbeddedBundle::new("greeting", GREETING_BYTES, GREETING_SHA256);
static GREETING_JSON: JsonBundle<Greeting> = JsonBundle::new(EmbeddedBundle::new(
    "greeting",
    GREETING_BYTES,
    GREETING_SHA256,
));
static BAD_CHECKSUM: EmbeddedBundle =
    EmbeddedBundle::new("bad-checksum", GREETING_BYTES, ZERO_SHA256);
static BAD_DIGEST: EmbeddedBundle =
    EmbeddedBundle::new("bad-digest", GREETING_BYTES, "not-a-digest");
static NOT_GZIP: EmbeddedBundle = EmbeddedBundle::new("not-gzip", b"plain text", ZERO_SHA256);

#[derive(Debug, Deserialize, PartialEq, Eq)]
struct Greeting {
    greeting: String,
    answer: u32,
}

fn gzip(payload: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder
        .write_all(payload)
        .expect("writing to an in-memory encoder cannot fail");
    encoder
        .finish()
        .expect("finishing an in-memory encoder cannot fail")
}

fn hex_digest(payload: &[u8]) -> String {
    use sha2::{Digest, Sha256};

    let digest = Sha256::digest(payload);
    let mut text = String::with_capacity(64);
    for byte in digest {
        use std::fmt::Write as _;
        write!(&mut text, "{byte:02x}").expect("writing to a String cannot fail");
    }
    text
}

#[test]
fn embedded_bundle_decompresses_once_and_verifies_checksum() {
    assert_eq!(GREETING.name(), "greeting");
    let first = GREETING.bytes().expect("valid embedded bundle");
    let second = GREETING.bytes().expect("cached embedded bundle");
    assert_eq!(first, GREETING_JSON_BYTES);
    assert!(std::ptr::eq(first.as_ptr(), second.as_ptr()));
}

#[test]
fn json_bundle_decodes_once() {
    assert_eq!(GREETING_JSON.name(), "greeting");
    let first = GREETING_JSON.get().expect("valid embedded JSON");
    let second = GREETING_JSON.get().expect("cached embedded JSON");
    assert_eq!(
        first,
        &Greeting {
            greeting: "hello".to_owned(),
            answer: 42,
        }
    );
    assert!(std::ptr::eq(first, second));
}

#[test]
fn checksum_mismatch_is_reported() {
    let error = BAD_CHECKSUM.bytes().expect_err("checksum must fail");
    assert_eq!(error.code(), "data_checksum");
    assert!(error.to_string().contains("bad-checksum"));
}

#[test]
fn malformed_digest_is_reported() {
    let error = BAD_DIGEST.bytes().expect_err("digest must fail");
    assert_eq!(error.code(), "data_digest");
}

#[test]
fn corrupt_gzip_is_reported() {
    let error = NOT_GZIP.bytes().expect_err("gzip must fail");
    assert_eq!(error.code(), "data_compression");
}

#[test]
fn external_bundle_is_verified() {
    let payload = b"external payload";
    let compressed = gzip(payload);
    let digest = hex_digest(payload);
    assert_eq!(
        load_external_gzip("external", &compressed, &digest).expect("valid external bundle"),
        payload
    );

    let error =
        load_external_gzip("external", &compressed, ZERO_SHA256).expect_err("checksum must fail");
    assert_eq!(error.code(), "data_checksum");
}

#[test]
fn external_json_is_decoded_and_validated() {
    let compressed = gzip(GREETING_JSON_BYTES);
    let digest = hex_digest(GREETING_JSON_BYTES);
    let greeting: Greeting =
        load_external_json("external-json", &compressed, &digest).expect("valid external JSON");
    assert_eq!(
        greeting,
        Greeting {
            greeting: "hello".to_owned(),
            answer: 42,
        }
    );

    let invalid = gzip(b"not json");
    let error = load_external_json::<Greeting>("bad-json", &invalid, &hex_digest(b"not json"))
        .expect_err("invalid JSON must fail");
    assert_eq!(error.code(), "data_json");
}
