//! Minimal CIF reader and writer (Task 1.15).
//!
//! Supports data blocks, tags, loops, quoted and semicolon strings, comments,
//! standard-uncertainty numbers, atom sites, symmetry operations and
//! anisotropic displacement parameters. It is a pure Rust reader for the
//! corpus; the full crystal-builder priority rules are documented in the
//! module comments and the API page.

use std::collections::BTreeMap;

use xtalforge_core::{Error, Structure};

mod coincident;
mod stream;
mod stream_atoms;
mod structure;
mod values;
mod writer;

pub use stream::read_cif_stream;
use structure::structure_from_block as extract_structure;
pub use writer::write_cif;

/// Extracts a structure from one parsed block (shared with the mmCIF reader).
pub(crate) fn structure_from_block(block: &CifBlock) -> Result<Structure, Error> {
    extract_structure(block)
}

/// One parsed CIF data block.
#[derive(Debug, Clone, PartialEq)]
pub struct CifBlock {
    /// Block name without the `data_` prefix.
    pub name: String,
    /// Single-valued items.
    pub items: BTreeMap<String, String>,
    /// Loop tables.
    pub loops: Vec<CifLoop>,
}

/// One CIF loop.
#[derive(Debug, Clone, PartialEq)]
pub struct CifLoop {
    /// Column tags.
    pub tags: Vec<String>,
    /// Rows of values.
    pub rows: Vec<Vec<String>>,
}

impl CifBlock {
    /// Single item value.
    #[must_use]
    pub fn item(&self, tag: &str) -> Option<&str> {
        self.items
            .iter()
            .find(|(candidate, _)| tag_matches(candidate, tag))
            .map(|(_, value)| value.as_str())
    }

    /// First loop containing `tag`.
    #[must_use]
    pub fn loop_with(&self, tag: &str) -> Option<&CifLoop> {
        self.loops.iter().find(|loop_| {
            loop_
                .tags
                .iter()
                .any(|candidate| tag_matches(candidate, tag))
        })
    }
}

impl CifLoop {
    /// Column index of a tag.
    #[must_use]
    pub fn column(&self, tag: &str) -> Option<usize> {
        self.tags
            .iter()
            .position(|candidate| tag_matches(candidate, tag))
    }

    /// Value at `row` and `tag`.
    #[must_use]
    pub fn value(&self, row: usize, tag: &str) -> Option<&str> {
        self.column(tag)
            .and_then(|column| self.rows.get(row)?.get(column))
            .map(String::as_str)
    }
}

/// Matches CIF and mmCIF spellings: case-insensitive, and `.` equals `_`.
fn tag_matches(candidate: &str, wanted: &str) -> bool {
    let candidate = candidate.as_bytes();
    let wanted = wanted.as_bytes();
    if candidate.len() != wanted.len() {
        return false;
    }
    candidate.iter().zip(wanted).all(|(left, right)| {
        let left = left.to_ascii_lowercase();
        let right = right.to_ascii_lowercase();
        let left = if left == b'.' { b'_' } else { left };
        let right = if right == b'.' { b'_' } else { right };
        left == right
    })
}

/// A parsed CIF document.
#[derive(Debug, Clone, PartialEq)]
pub struct CifDocument {
    /// Data blocks in file order.
    pub blocks: Vec<CifBlock>,
}

impl CifDocument {
    /// Parses a CIF text.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an unbalanced loop or an empty
    /// document.
    pub fn parse(text: &str) -> Result<Self, Error> {
        let tokens = tokenize(text);
        let mut blocks = Vec::new();
        let mut index = 0;
        while index < tokens.len() {
            let token = tokens.get(index).cloned().unwrap_or_default();
            if let Some(name) = token.strip_prefix("data_") {
                blocks.push(CifBlock {
                    name: name.to_owned(),
                    items: BTreeMap::new(),
                    loops: Vec::new(),
                });
                index += 1;
                continue;
            }
            if token == "loop_" {
                let mut tags = Vec::new();
                index += 1;
                while let Some(next) = tokens.get(index) {
                    if next.starts_with('_') {
                        tags.push(next.clone());
                        index += 1;
                    } else {
                        break;
                    }
                }
                let mut values = Vec::new();
                while let Some(next) = tokens.get(index) {
                    if next.starts_with('_')
                        || next == "loop_"
                        || next.starts_with("data_")
                        || next.starts_with("save_")
                    {
                        break;
                    }
                    values.push(next.clone());
                    index += 1;
                }
                if tags.is_empty() {
                    return Err(Error::InvalidInput("CIF loop has no tags".to_owned()));
                }
                if values.len() % tags.len() != 0 {
                    return Err(Error::InvalidInput(format!(
                        "CIF loop has {} values for {} tags",
                        values.len(),
                        tags.len()
                    )));
                }
                let rows = values.chunks(tags.len()).map(<[String]>::to_vec).collect();
                if let Some(block) = blocks.last_mut() {
                    block.loops.push(CifLoop { tags, rows });
                }
                continue;
            }
            if token.starts_with('_') {
                let value = tokens.get(index + 1).cloned().unwrap_or_default();
                if let Some(block) = blocks.last_mut() {
                    block.items.insert(token, value);
                }
                index += 2;
                continue;
            }
            index += 1;
        }
        if blocks.is_empty() {
            return Err(Error::InvalidInput("CIF has no data block".to_owned()));
        }
        Ok(Self { blocks })
    }
}

/// Reads the first data block as a structure.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an invalid cell, missing sites or an
/// unsupported atom loop.
pub fn read_cif(text: &str) -> Result<Structure, Error> {
    let document = CifDocument::parse(text)?;
    let block = document
        .blocks
        .first()
        .ok_or_else(|| Error::InvalidInput("CIF has no data block".to_owned()))?;
    structure_from_block(block)
}

/// Reads a named block as a structure.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the block is missing.
pub fn read_cif_block(text: &str, name: &str) -> Result<Structure, Error> {
    let document = CifDocument::parse(text)?;
    let block = document
        .blocks
        .iter()
        .find(|block| block.name == name)
        .ok_or_else(|| Error::InvalidInput(format!("CIF block not found: {name}")))?;
    structure_from_block(block)
}

/// Reads every block as a structure.
///
/// # Errors
///
/// Returns the error of the first block that cannot become a structure.
pub fn read_cif_all(text: &str) -> Result<Vec<Structure>, Error> {
    let document = CifDocument::parse(text)?;
    document.blocks.iter().map(structure_from_block).collect()
}

fn tokenize(text: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    let mut text_block: Option<String> = None;
    for line in text.lines() {
        if let Some(mut block) = text_block.take() {
            if line.starts_with(';') {
                tokens.push(block);
            } else {
                block.push_str(line);
                block.push('\n');
                text_block = Some(block);
            }
            continue;
        }
        if line.starts_with(';') {
            text_block = Some(line.chars().skip(1).collect());
            continue;
        }
        let mut current = String::new();
        let mut quote: Option<char> = None;
        for character in line.chars() {
            if let Some(active) = quote {
                if character == active {
                    tokens.push(std::mem::take(&mut current));
                    quote = None;
                } else {
                    current.push(character);
                }
                continue;
            }
            match character {
                '\'' | '"' => {
                    if !current.is_empty() {
                        tokens.push(std::mem::take(&mut current));
                    }
                    quote = Some(character);
                }
                '#' => break,
                value if value.is_whitespace() => {
                    if !current.is_empty() {
                        tokens.push(std::mem::take(&mut current));
                    }
                }
                value => current.push(value),
            }
        }
        if !current.is_empty() {
            tokens.push(current);
        }
    }
    if let Some(block) = text_block {
        tokens.push(block);
    }
    tokens
}
