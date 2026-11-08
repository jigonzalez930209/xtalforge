//! Integration tests for the stable error codes and messages.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::Error;

#[test]
fn invalid_input_code_and_message() {
    let error = Error::InvalidInput("negative volume".to_owned());
    assert_eq!(error.code(), "invalid_input");
    assert_eq!(error.message(), "invalid input: negative volume");
    assert_eq!(error.to_string(), error.message());
}

#[test]
fn data_error_codes_are_stable() {
    let compression = Error::DataCompression {
        name: "elements".to_owned(),
        reason: "bad gzip".to_owned(),
    };
    let checksum = Error::DataChecksum {
        name: "elements".to_owned(),
    };
    let digest = Error::DataDigest {
        name: "elements".to_owned(),
    };
    let json = Error::DataJson {
        name: "elements".to_owned(),
        reason: "bad json".to_owned(),
    };

    assert_eq!(compression.code(), "data_compression");
    assert_eq!(checksum.code(), "data_checksum");
    assert_eq!(digest.code(), "data_digest");
    assert_eq!(json.code(), "data_json");
    assert!(compression.message().contains("elements"));
    assert!(checksum.message().contains("checksum"));
    assert!(digest.message().contains("digest"));
    assert!(json.message().contains("JSON"));
}
