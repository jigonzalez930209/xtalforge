//! Provenance embedding and extraction (Task 16.2.2).
//!
//! Structures carry a [`Provenance`] record in the `provenance` metadata key
//! (see `xtalforge_core::provenance`). This module embeds and reads that record
//! in the formats that have room for it:
//!
//! | Surface | Carrier |
//! |---|---|
//! | CIF | `_xtalforge_provenance` tag with a semicolon text block |
//! | CIF comment | `# xtalforge-provenance: <json>` (read only) |
//! | Extended XYZ | `xtalforge_provenance=<base64 JSON>` comment pair |
//! | Project file | `provenance` key of the `session.json` part |
//! | Structure JSON | `metadata.provenance` of the structure object |
//!
//! [`provenance`] accepts a structure, a text or the bytes of a project
//! archive and returns the record, or `None` when the source carries none.
//! Malformed payloads are errors; the functions never panic and never guess a
//! different schema version.

use std::collections::BTreeMap;

use xtalforge_core::{Error, Provenance, Structure};

use crate::cif::CifDocument;
use crate::project::read_project;
use crate::trajectory::parse_info;

/// CIF tag carrying a provenance record as a semicolon text block.
pub const CIF_PROVENANCE_TAG: &str = "_xtalforge_provenance";

/// CIF comment prefix accepted when reading (the writer emits the tag).
pub const CIF_PROVENANCE_COMMENT: &str = "# xtalforge-provenance:";

/// Extended XYZ comment key carrying a base64-encoded provenance record.
pub const XYZ_PROVENANCE_KEY: &str = "xtalforge_provenance";

/// Key of the `session.json` project part that carries the record.
pub const PROJECT_PROVENANCE_KEY: &str = "provenance";

/// A value that may carry a provenance record.
///
/// Implemented for [`Structure`], `str`, `String`, `[u8]` and `Vec<u8>`. Byte
/// sources that start with the ZIP magic are read as project files; any other
/// bytes are decoded as text.
pub trait ProvenanceSource {
    /// Returns the record carried by this source, if any.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the source is malformed or its
    /// embedded provenance record is invalid.
    fn provenance(&self) -> Result<Option<Provenance>, Error>;
}

/// Extracts the provenance record from a structure, text or project bytes.
///
/// This is the surface used by the `xtalforge provenance` command.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the source cannot be recognised or its
/// embedded record is invalid.
pub fn provenance<S>(source: &S) -> Result<Option<Provenance>, Error>
where
    S: ProvenanceSource + ?Sized,
{
    source.provenance()
}

impl ProvenanceSource for Structure {
    fn provenance(&self) -> Result<Option<Provenance>, Error> {
        Provenance::from_structure(self)
    }
}

impl ProvenanceSource for str {
    fn provenance(&self) -> Result<Option<Provenance>, Error> {
        let text = self.strip_prefix('\u{feff}').unwrap_or(self);
        let trimmed = text.trim_start();
        if trimmed.starts_with('{') {
            return json_provenance(trimmed);
        }
        if text.contains(CIF_PROVENANCE_TAG)
            || text.contains(CIF_PROVENANCE_COMMENT)
            || looks_like_cif(text)
        {
            return cif_provenance(text);
        }
        if let Some(record) = xyz_provenance(text)? {
            return Ok(Some(record));
        }
        if looks_like_xyz(text) {
            return Ok(None);
        }
        Err(Error::InvalidInput(
            "provenance: could not detect the source format; pass a structure, CIF, extended XYZ, structure JSON or project archive"
                .to_owned(),
        ))
    }
}

impl ProvenanceSource for String {
    fn provenance(&self) -> Result<Option<Provenance>, Error> {
        self.as_str().provenance()
    }
}

impl ProvenanceSource for [u8] {
    fn provenance(&self) -> Result<Option<Provenance>, Error> {
        if self.starts_with(b"PK\x03\x04") {
            return read_project(self)?.provenance();
        }
        crate::text::decode_text(self).provenance()
    }
}

impl ProvenanceSource for Vec<u8> {
    fn provenance(&self) -> Result<Option<Provenance>, Error> {
        self.as_slice().provenance()
    }
}

fn json_provenance(text: &str) -> Result<Option<Provenance>, Error> {
    if text.contains("\"xtalforge_version\"") {
        return Provenance::from_json(text).map(Some);
    }
    let structure = if text.contains("\"schema_version\"") {
        xtalforge_core::StructureDocument::from_json(text)?.structure
    } else {
        serde_json::from_str::<Structure>(text).map_err(|error| {
            Error::InvalidInput(format!("provenance: invalid structure JSON: {error}"))
        })?
    };
    Provenance::from_structure(&structure)
}

fn looks_like_cif(text: &str) -> bool {
    text.contains("_cell_length_")
        || text
            .lines()
            .next()
            .is_some_and(|line| line.trim_start().starts_with("data_"))
}

fn cif_provenance(text: &str) -> Result<Option<Provenance>, Error> {
    if let Some(payload) = cif_comment_payload(text) {
        return Provenance::from_json(payload).map(Some);
    }
    let document = CifDocument::parse(text)?;
    for block in &document.blocks {
        if let Some(payload) = block.item(CIF_PROVENANCE_TAG) {
            return Provenance::from_json(payload).map(Some);
        }
    }
    Ok(None)
}

fn cif_comment_payload(text: &str) -> Option<&str> {
    let mut in_block = false;
    for line in text.lines() {
        if line.starts_with(';') {
            in_block = !in_block;
            continue;
        }
        if in_block {
            continue;
        }
        let line = line.trim();
        for prefix in [CIF_PROVENANCE_COMMENT, "# xtalforge_provenance:"] {
            if let Some(payload) = line.strip_prefix(prefix) {
                return Some(payload.trim());
            }
        }
    }
    None
}

fn looks_like_xyz(text: &str) -> bool {
    let mut lines = text.lines();
    lines
        .next()
        .is_some_and(|line| line.trim().parse::<usize>().is_ok())
        && lines.next().is_some()
}

fn xyz_provenance(text: &str) -> Result<Option<Provenance>, Error> {
    let mut lines = text.lines();
    let mut count = match lines.next() {
        Some(line) => match line.trim().parse::<usize>() {
            Ok(count) => count,
            Err(_) => return Ok(None),
        },
        None => return Ok(None),
    };
    loop {
        let Some(comment) = lines.next() else {
            return Ok(None);
        };
        if let Some(record) = provenance_from_xyz_comment(comment)? {
            return Ok(Some(record));
        }
        for _ in 0..count {
            if lines.next().is_none() {
                return Ok(None);
            }
        }
        loop {
            match lines.next() {
                Some(line) if line.trim().is_empty() => continue,
                Some(line) => match line.trim().parse::<usize>() {
                    Ok(next) => {
                        count = next;
                        break;
                    }
                    Err(_) => return Ok(None),
                },
                None => return Ok(None),
            }
        }
    }
}

pub(crate) fn provenance_from_xyz_comment(comment: &str) -> Result<Option<Provenance>, Error> {
    let info = parse_info(comment);
    provenance_from_xyz_info(&info)
}

pub(crate) fn provenance_from_xyz_info(
    info: &BTreeMap<String, String>,
) -> Result<Option<Provenance>, Error> {
    let Some(encoded) = info.get(XYZ_PROVENANCE_KEY) else {
        return Ok(None);
    };
    let bytes = base64_decode(encoded)?;
    let json = String::from_utf8(bytes).map_err(|error| {
        Error::InvalidInput(format!(
            "provenance: extended XYZ payload is not UTF-8: {error}"
        ))
    })?;
    Provenance::from_json(&json).map(Some)
}

pub(crate) fn xyz_provenance_pair(record: &Provenance) -> Result<String, Error> {
    let json = record.to_json()?;
    Ok(format!(
        "{XYZ_PROVENANCE_KEY}={}",
        base64_encode(json.as_bytes())
    ))
}

fn base64_encode(bytes: &[u8]) -> String {
    const ALPHABET: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for chunk in bytes.chunks(3) {
        let first = chunk.first().copied().unwrap_or(0) as u32;
        let second = chunk.get(1).copied().unwrap_or(0) as u32;
        let third = chunk.get(2).copied().unwrap_or(0) as u32;
        let triple = (first << 16) | (second << 8) | third;
        let used = chunk.len() + 1;
        for (position, shift) in [18, 12, 6, 0].into_iter().enumerate() {
            if position < used {
                let index = ((triple >> shift) & 0x3f) as usize;
                out.push(char::from(ALPHABET.get(index).copied().unwrap_or(b'=')));
            } else {
                out.push('=');
            }
        }
    }
    out
}

fn base64_decode(text: &str) -> Result<Vec<u8>, Error> {
    let text = text.trim();
    let payload = text.trim_end_matches('=').as_bytes();
    if payload.len() % 4 == 1 {
        return Err(base64_error(text));
    }
    let mut values = Vec::with_capacity(payload.len());
    for &byte in payload {
        values.push(base64_value(byte).ok_or_else(|| base64_error(text))?);
    }
    let mut out = Vec::with_capacity(values.len() / 4 * 3 + 3);
    for chunk in values.chunks(4) {
        let first = *chunk.first().unwrap_or(&0) as u32;
        let second = *chunk.get(1).unwrap_or(&0) as u32;
        let third = *chunk.get(2).unwrap_or(&0) as u32;
        let fourth = *chunk.get(3).unwrap_or(&0) as u32;
        let triple = (first << 18) | (second << 12) | (third << 6) | fourth;
        out.push((triple >> 16) as u8);
        if chunk.len() > 2 {
            out.push((triple >> 8) as u8);
        }
        if chunk.len() > 3 {
            out.push(triple as u8);
        }
    }
    Ok(out)
}

fn base64_value(byte: u8) -> Option<u8> {
    match byte {
        b'A'..=b'Z' => Some(byte - b'A'),
        b'a'..=b'z' => Some(byte - b'a' + 26),
        b'0'..=b'9' => Some(byte - b'0' + 52),
        b'+' => Some(62),
        b'/' => Some(63),
        _ => None,
    }
}

fn base64_error(text: &str) -> Error {
    Error::InvalidInput(format!(
        "provenance: extended XYZ payload is not valid base64: {text:?}"
    ))
}
