//! Integration tests for the on-demand database cache (Task 17.1.4).

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::path::PathBuf;

use xtalforge_mof::platform::{FetchOptions, fetch_core_mof, fetch_database};

fn temp_dir(name: &str) -> PathBuf {
    let path = std::env::temp_dir().join(format!(
        "xtalforge-fetch-{name}-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&path).unwrap();
    path
}

fn digest(bytes: &[u8]) -> String {
    use sha2::{Digest, Sha256};
    Sha256::digest(bytes)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

#[test]
fn a_cached_archive_is_verified_and_reused_offline() {
    let dir = temp_dir("cache");
    let payload = b"not really a zip, but the cache only stores bytes";
    std::fs::write(dir.join("core-mof.zip"), payload).unwrap();
    let options = FetchOptions::new(&dir)
        .offline(true)
        .with_sha256(digest(payload));
    let fetched = fetch_core_mof("https://example.invalid/core.zip", &options).unwrap();
    assert!(fetched.from_cache);
    assert_eq!(fetched.bytes, payload.len());
    assert_eq!(fetched.sha256, digest(payload));
}

#[test]
fn offline_mode_refuses_a_missing_archive() {
    let dir = temp_dir("missing");
    let options = FetchOptions::new(&dir).offline(true);
    let error = fetch_core_mof("https://example.invalid/core.zip", &options).unwrap_err();
    assert!(error.to_string().contains("offline mode"));
}

#[test]
fn a_digest_mismatch_is_refused() {
    let dir = temp_dir("digest");
    std::fs::write(dir.join("core-mof.zip"), b"payload").unwrap();
    let options = FetchOptions::new(&dir)
        .offline(true)
        .with_sha256("0".repeat(64));
    let error = fetch_core_mof("https://example.invalid/core.zip", &options).unwrap_err();
    assert!(error.to_string().contains("SHA-256 mismatch"));
}

#[test]
fn a_path_in_the_file_name_is_refused() {
    let dir = temp_dir("name");
    let options = FetchOptions::new(&dir).offline(true);
    for name in ["../escape.zip", "a/b.zip", "..", ""] {
        let error = fetch_database("https://example.invalid", name, &options).unwrap_err();
        assert!(error.to_string().contains("plain file name"), "{name}");
    }
}

#[test]
fn a_non_http_url_is_refused() {
    let dir = temp_dir("url");
    let options = FetchOptions::new(&dir).offline(true);
    let error = fetch_database("file:///tmp/core.zip", "core.zip", &options).unwrap_err();
    assert!(error.to_string().contains("http or https"));
}

#[cfg(feature = "net")]
#[test]
fn a_downloaded_archive_is_written_and_verified() {
    use std::io::{Read, Write};
    use std::net::TcpListener;

    let payload = b"downloaded database bytes".to_vec();
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let address = listener.local_addr().unwrap();
    let body = payload.clone();
    let server = std::thread::spawn(move || {
        let (mut stream, _) = listener.accept().unwrap();
        let mut request = [0u8; 1024];
        let _ = stream.read(&mut request);
        write!(
            stream,
            "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
            body.len()
        )
        .unwrap();
        stream.write_all(&body).unwrap();
    });
    let dir = temp_dir("download");
    let options = FetchOptions::new(&dir).with_sha256(digest(&payload));
    let fetched = fetch_core_mof(&format!("http://{address}/core.zip"), &options).unwrap();
    server.join().unwrap();
    assert!(!fetched.from_cache);
    assert_eq!(std::fs::read(&fetched.path).unwrap(), payload);
    let cached = fetch_core_mof(&format!("http://{address}/core.zip"), &options).unwrap();
    assert!(cached.from_cache);
}
