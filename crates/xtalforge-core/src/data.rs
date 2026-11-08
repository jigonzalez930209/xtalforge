//! Embedded and lazily loaded data bundles.
//!
//! Large tables live in compressed files instead of Rust literals (AGENTS.md
//! rule 1). This module decodes a gzip bundle, verifies its SHA-256 checksum
//! against the manifest value and keeps the result for the lifetime of the
//! process, so repeated lookups never decompress or parse twice.
//!
//! The checksum covers the **decompressed** payload, not the compressed file.
//! Compute it with:
//!
//! ```text
//! gzip -dc data/example.json.gz | sha256sum
//! ```
//!
//! Embedded bundles use [`EmbeddedBundle`] and [`JsonBundle`] with
//! `include_bytes!`; WASM callers and large external databases use
//! [`load_external_gzip`] or [`load_external_json`] with bytes fetched by the
//! host.

use std::io::Read;
use std::sync::OnceLock;

use serde::de::DeserializeOwned;
use sha2::{Digest, Sha256};

use crate::Error;

/// A gzip-compressed blob embedded in the binary at compile time.
///
/// The first call to [`EmbeddedBundle::bytes`] decompresses and verifies the
/// payload; later calls return the cached bytes. A failure is cached as well.
#[derive(Debug)]
pub struct EmbeddedBundle {
    name: &'static str,
    gzip: &'static [u8],
    sha256: &'static str,
    decoded: OnceLock<Result<Vec<u8>, Error>>,
}

impl EmbeddedBundle {
    /// Creates a bundle description.
    ///
    /// `sha256` is the lowercase or uppercase hexadecimal SHA-256 digest of
    /// the decompressed payload, as printed by `sha256sum`.
    #[must_use]
    pub const fn new(name: &'static str, gzip: &'static [u8], sha256: &'static str) -> Self {
        Self {
            name,
            gzip,
            sha256,
            decoded: OnceLock::new(),
        }
    }

    /// Name of the bundle, used in error messages and manifests.
    #[must_use]
    pub const fn name(&self) -> &'static str {
        self.name
    }

    /// Returns the decompressed and verified payload.
    ///
    /// # Errors
    ///
    /// Returns [`Error::DataCompression`] for invalid gzip data,
    /// [`Error::DataDigest`] for a malformed manifest digest and
    /// [`Error::DataChecksum`] when the payload does not match the digest.
    pub fn bytes(&self) -> Result<&[u8], Error> {
        let cached = self.decoded.get_or_init(|| self.decode());
        match cached {
            Ok(bytes) => Ok(bytes),
            Err(error) => Err(error.clone()),
        }
    }

    fn decode(&self) -> Result<Vec<u8>, Error> {
        decode_gzip(self.name, self.gzip, self.sha256)
    }
}

/// A JSON bundle decoded once into a typed value.
///
/// `T` must be `Send + Sync` so the bundle can live in a `static` and be
/// shared between threads without locking.
#[derive(Debug)]
pub struct JsonBundle<T> {
    bundle: EmbeddedBundle,
    parsed: OnceLock<Result<T, Error>>,
}

impl<T> JsonBundle<T> {
    /// Wraps an embedded bundle.
    #[must_use]
    pub const fn new(bundle: EmbeddedBundle) -> Self {
        Self {
            bundle,
            parsed: OnceLock::new(),
        }
    }
}

impl<T> JsonBundle<T>
where
    T: DeserializeOwned + Send + Sync,
{
    /// Returns the decoded value, decoding it on the first call only.
    ///
    /// # Errors
    ///
    /// Returns any error from [`EmbeddedBundle::bytes`] or
    /// [`Error::DataJson`] when the payload is not valid JSON for `T`.
    pub fn get(&self) -> Result<&T, Error> {
        let cached = self.parsed.get_or_init(|| {
            let payload = self.bundle.bytes()?;
            serde_json::from_slice(payload).map_err(|error| Error::DataJson {
                name: self.bundle.name.to_owned(),
                reason: error.to_string(),
            })
        });
        match cached {
            Ok(value) => Ok(value),
            Err(error) => Err(error.clone()),
        }
    }

    /// Name of the underlying bundle.
    #[must_use]
    pub fn name(&self) -> &'static str {
        self.bundle.name()
    }
}

/// Decompresses and verifies a gzip bundle supplied at runtime.
///
/// Hosts that cannot embed data (WASM) fetch the bytes and call this function;
/// native callers use it for databases loaded from disk. The result is not
/// cached: keep it in the caller's own structure when it is reused.
///
/// # Errors
///
/// Returns [`Error::DataCompression`] for invalid gzip data,
/// [`Error::DataDigest`] for a malformed digest and [`Error::DataChecksum`]
/// when the payload does not match the digest.
pub fn load_external_gzip(name: &str, gzip: &[u8], sha256: &str) -> Result<Vec<u8>, Error> {
    decode_gzip(name, gzip, sha256)
}

/// Decompresses, verifies and decodes a JSON bundle supplied at runtime.
///
/// # Errors
///
/// Returns the errors of [`load_external_gzip`] and [`Error::DataJson`].
pub fn load_external_json<T>(name: &str, gzip: &[u8], sha256: &str) -> Result<T, Error>
where
    T: DeserializeOwned,
{
    let payload = decode_gzip(name, gzip, sha256)?;
    serde_json::from_slice(&payload).map_err(|error| Error::DataJson {
        name: name.to_owned(),
        reason: error.to_string(),
    })
}

fn decode_gzip(name: &str, gzip: &[u8], sha256: &str) -> Result<Vec<u8>, Error> {
    let mut decoder = flate2::read::GzDecoder::new(gzip);
    let mut payload = Vec::new();
    decoder
        .read_to_end(&mut payload)
        .map_err(|error| Error::DataCompression {
            name: name.to_owned(),
            reason: error.to_string(),
        })?;
    verify_sha256(name, &payload, sha256)?;
    Ok(payload)
}

fn verify_sha256(name: &str, payload: &[u8], expected: &str) -> Result<(), Error> {
    let expected = decode_digest(name, expected)?;
    let actual = Sha256::digest(payload);
    if actual.as_slice() != expected.as_slice() {
        return Err(Error::DataChecksum {
            name: name.to_owned(),
        });
    }
    Ok(())
}

fn decode_digest(name: &str, hex: &str) -> Result<[u8; 32], Error> {
    let bytes = hex.as_bytes();
    if bytes.len() != 64 {
        return Err(Error::DataDigest {
            name: name.to_owned(),
        });
    }

    let mut digest = [0u8; 32];
    for (slot, pair) in digest.iter_mut().zip(bytes.chunks_exact(2)) {
        let invalid = || Error::DataDigest {
            name: name.to_owned(),
        };
        let (high, rest) = pair.split_first().ok_or_else(invalid)?;
        let low = rest.first().ok_or_else(invalid)?;
        let high = hex_nibble(*high).ok_or_else(invalid)?;
        let low = hex_nibble(*low).ok_or_else(invalid)?;
        *slot = (high << 4) | low;
    }
    Ok(digest)
}

fn hex_nibble(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'a'..=b'f' => Some(byte - b'a' + 10),
        b'A'..=b'F' => Some(byte - b'A' + 10),
        _ => None,
    }
}
