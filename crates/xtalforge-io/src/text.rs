//! Text decoding shared by the file readers (Task 1.15.1).

/// Decodes UTF-8 with or without a byte-order mark, falling back to Latin-1.
///
/// A UTF-8 BOM is stripped. Invalid UTF-8 is decoded byte by byte as Latin-1,
/// which never fails; the reader can still report a warning at a higher level
/// when it cares about the fallback.
#[must_use]
pub fn decode_text(bytes: &[u8]) -> String {
    let bytes = bytes.strip_prefix(&[0xEF, 0xBB, 0xBF]).unwrap_or(bytes);
    match std::str::from_utf8(bytes) {
        Ok(text) => text.to_owned(),
        Err(_) => bytes.iter().map(|byte| char::from(*byte)).collect(),
    }
}

/// Reads a file and decodes it as text.
///
/// # Errors
///
/// Returns [`xtalforge_core::Error::InvalidInput`] when the file cannot be
/// read.
pub fn read_text(path: &std::path::Path) -> Result<String, xtalforge_core::Error> {
    let bytes = std::fs::read(path).map_err(|error| {
        xtalforge_core::Error::InvalidInput(format!("could not read {}: {error}", path.display()))
    })?;
    Ok(decode_text(&bytes))
}
