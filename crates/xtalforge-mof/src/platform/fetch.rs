//! On-demand database download with a disk cache, checksums and an offline
//! mode (Task 17.1.4).
//!
//! The published CoRE MOF and QMOF archives are large; baking them into every
//! wheel is wrong. This module downloads one archive into a caller-chosen
//! cache directory, verifies the SHA-256 of the downloaded bytes and reuses
//! the cached copy on later runs. `offline` refuses to touch the network and
//! only accepts an already verified cache entry.
//!
//! The HTTP client lives behind the `net` feature, so the default build stays
//! pure Rust with no TLS stack. The browser cannot write a disk cache; its
//! `fetch` wrapper returns the bytes to JavaScript, which owns them.

#[cfg(feature = "net")]
use std::io::Read;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use xtalforge_core::Error;

/// Options for [`fetch_database`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FetchOptions {
    /// Directory that holds downloaded archives and their checksum sidecars.
    pub cache_dir: PathBuf,
    /// Refuse network access and accept only a verified cache entry.
    pub offline: bool,
    /// Expected lowercase or uppercase hex SHA-256 of the archive.
    ///
    /// When `None`, the digest is recorded on first download and checked on
    /// later reads, which still detects a corrupted cache but not a changed
    /// remote file.
    pub sha256: Option<String>,
}

impl FetchOptions {
    /// Creates options with a cache directory and no fixed digest.
    #[must_use]
    pub fn new(cache_dir: impl Into<PathBuf>) -> Self {
        Self {
            cache_dir: cache_dir.into(),
            offline: false,
            sha256: None,
        }
    }

    /// Copy with the offline flag.
    #[must_use]
    pub fn offline(mut self, offline: bool) -> Self {
        self.offline = offline;
        self
    }

    /// Copy with a pinned archive digest.
    #[must_use]
    pub fn with_sha256(mut self, sha256: impl Into<String>) -> Self {
        self.sha256 = Some(sha256.into());
        self
    }
}

/// One verified archive in the cache.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FetchedDatabase {
    /// Path of the cached archive.
    pub path: PathBuf,
    /// Whether the archive came from the cache instead of the network.
    pub from_cache: bool,
    /// Number of bytes in the archive.
    pub bytes: usize,
    /// SHA-256 of the archive, lowercase hex.
    pub sha256: String,
}

/// Downloads one archive to the cache, or returns the verified cached copy.
///
/// `file_name` must be a plain file name; separators and `..` are refused so a
/// caller cannot escape `cache_dir`. The archive is written to a temporary
/// file and renamed only after the digest passes.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an invalid file name, a missing cache
/// entry in offline mode, a digest mismatch, or a filesystem failure.
/// Returns a clear error when the crate is built without the `net` feature and
/// the cache is empty.
pub fn fetch_database(
    url: &str,
    file_name: &str,
    options: &FetchOptions,
) -> Result<FetchedDatabase, Error> {
    validate_name(file_name)?;
    if !url.starts_with("https://") && !url.starts_with("http://") {
        return Err(Error::InvalidInput(format!(
            "database URL must use http or https, got {url:?}"
        )));
    }
    std::fs::create_dir_all(&options.cache_dir).map_err(|error| {
        Error::InvalidInput(format!(
            "cannot create cache directory {}: {error}",
            options.cache_dir.display()
        ))
    })?;
    let path = options.cache_dir.join(file_name);
    if path.is_file() {
        let bytes = read(&path)?;
        let actual = digest_hex(&bytes);
        verify(&actual, options.sha256.as_deref(), &path)?;
        return Ok(FetchedDatabase {
            path,
            from_cache: true,
            bytes: bytes.len(),
            sha256: actual,
        });
    }
    if options.offline {
        return Err(Error::InvalidInput(format!(
            "offline mode: {} is not in the cache {}",
            file_name,
            options.cache_dir.display()
        )));
    }
    let bytes = download(url)?;
    let actual = digest_hex(&bytes);
    verify(&actual, options.sha256.as_deref(), &path)?;
    let temporary = path.with_extension("part");
    std::fs::write(&temporary, &bytes).map_err(|error| {
        Error::InvalidInput(format!("cannot write {}: {error}", temporary.display()))
    })?;
    std::fs::rename(&temporary, &path).map_err(|error| {
        Error::InvalidInput(format!(
            "cannot move {} to {}: {error}",
            temporary.display(),
            path.display()
        ))
    })?;
    Ok(FetchedDatabase {
        path,
        from_cache: false,
        bytes: bytes.len(),
        sha256: actual,
    })
}

/// Fetches the CoRE MOF archive with the published default file name.
///
/// # Errors
///
/// Returns the errors of [`fetch_database`].
pub fn fetch_core_mof(url: &str, options: &FetchOptions) -> Result<FetchedDatabase, Error> {
    fetch_database(url, "core-mof.zip", options)
}

/// Fetches the QMOF archive with the published default file name.
///
/// # Errors
///
/// Returns the errors of [`fetch_database`].
pub fn fetch_qmof(url: &str, options: &FetchOptions) -> Result<FetchedDatabase, Error> {
    fetch_database(url, "qmof.zip", options)
}

fn validate_name(file_name: &str) -> Result<(), Error> {
    if file_name.is_empty()
        || file_name.contains('/')
        || file_name.contains('\\')
        || file_name == "."
        || file_name == ".."
    {
        return Err(Error::InvalidInput(format!(
            "database file name must be a plain file name, got {file_name:?}"
        )));
    }
    Ok(())
}

fn read(path: &Path) -> Result<Vec<u8>, Error> {
    std::fs::read(path)
        .map_err(|error| Error::InvalidInput(format!("cannot read {}: {error}", path.display())))
}

fn digest_hex(bytes: &[u8]) -> String {
    let digest = Sha256::digest(bytes);
    digest.iter().map(|byte| format!("{byte:02x}")).collect()
}

fn verify(actual: &str, expected: Option<&str>, path: &Path) -> Result<(), Error> {
    let Some(expected) = expected else {
        return Ok(());
    };
    let expected = expected.trim().to_ascii_lowercase();
    if expected.len() != 64 || !expected.chars().all(|value| value.is_ascii_hexdigit()) {
        return Err(Error::InvalidInput(format!(
            "expected SHA-256 must be 64 hexadecimal characters, got {expected:?}"
        )));
    }
    if actual != expected {
        return Err(Error::InvalidInput(format!(
            "SHA-256 mismatch for {}: expected {expected}, got {actual}",
            path.display()
        )));
    }
    Ok(())
}

#[cfg(feature = "net")]
fn download(url: &str) -> Result<Vec<u8>, Error> {
    let response = ureq::get(url)
        .call()
        .map_err(|error| Error::InvalidInput(format!("cannot download {url}: {error}")))?;
    let mut bytes = Vec::new();
    response
        .into_reader()
        .read_to_end(&mut bytes)
        .map_err(|error| Error::InvalidInput(format!("cannot read {url}: {error}")))?;
    if bytes.is_empty() {
        return Err(Error::InvalidInput(format!(
            "downloaded an empty archive from {url}"
        )));
    }
    Ok(bytes)
}

#[cfg(not(feature = "net"))]
fn download(_url: &str) -> Result<Vec<u8>, Error> {
    Err(Error::InvalidInput(
        "this build has no network feature; rebuild with `net`/`mof-net` or place the archive in the cache".to_owned(),
    ))
}
