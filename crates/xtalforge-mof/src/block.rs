//! Curated building blocks (Tasks 6.11 and 6.12).
//!
//! A block is a small molecule in Cartesian Angstrom coordinates. Its origin
//! is the origin of the coordinate list, not the centroid. Connection points
//! are outward vectors from that origin; they are not atoms. The `cell` matrix
//! is a bounding box used for display and validation, not for placement.

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};
use xtalforge_core::elements::element_by_symbol;

use crate::Error;

/// Whether a block occupies vertices or edges.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum BlockKind {
    /// A node block, placed on every topology vertex.
    Node,
    /// A linker block, placed on every topology edge.
    Linker,
}

/// One building block.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Block {
    /// Block name, for example `Zn4O`.
    pub name: String,
    /// Node or linker.
    pub kind: BlockKind,
    /// 3x3 bounding-cell matrix in Angstrom, rows as vectors.
    pub cell: [[f64; 3]; 3],
    /// Cartesian atom positions in Angstrom.
    pub points: Vec<[f64; 3]>,
    /// Cartesian connection-point positions in Angstrom.
    pub connection_points: Vec<[f64; 3]>,
    /// Element symbol for each entry of `points`.
    pub elements: Vec<String>,
}

impl Block {
    /// Number of connection points.
    #[must_use]
    pub fn connection_count(&self) -> usize {
        self.connection_points.len()
    }

    /// Number of atoms.
    #[must_use]
    pub fn atom_count(&self) -> usize {
        self.points.len()
    }

    /// Checks the name, atom list, elements, connection points and cell.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty name, mismatched element
    /// and point lists, an unknown element, non-finite coordinates, a singular
    /// cell, a node with no connection point, or a linker with fewer than two.
    pub fn validate(&self) -> Result<(), Error> {
        if self.name.trim().is_empty() {
            return Err(Error::InvalidInput(
                "block name must not be empty".to_owned(),
            ));
        }
        if self.points.is_empty() {
            return Err(Error::InvalidInput(format!(
                "block `{}` has no atoms",
                self.name
            )));
        }
        if self.points.len() != self.elements.len() {
            return Err(Error::InvalidInput(format!(
                "block `{}` has {} points but {} elements",
                self.name,
                self.points.len(),
                self.elements.len()
            )));
        }
        for element in &self.elements {
            element_by_symbol(element)?;
        }
        if self
            .points
            .iter()
            .chain(self.connection_points.iter())
            .flatten()
            .any(|coordinate| !coordinate.is_finite())
        {
            return Err(Error::InvalidInput(format!(
                "block `{}` has a non-finite coordinate",
                self.name
            )));
        }
        let required = match self.kind {
            BlockKind::Node => 1,
            BlockKind::Linker => 2,
        };
        if self.connection_points.len() < required {
            return Err(Error::InvalidInput(format!(
                "block `{}` is a {:?} and needs at least {required} connection point(s)",
                self.name, self.kind
            )));
        }
        let determinant = determinant(self.cell);
        if !determinant.is_finite() || determinant.abs() < 1e-9 {
            return Err(Error::InvalidInput(format!(
                "block `{}` has a singular bounding cell",
                self.name
            )));
        }
        Ok(())
    }
}

fn determinant(matrix: [[f64; 3]; 3]) -> f64 {
    let [first, second, third] = matrix;
    dot(first, cross(second, third))
}

fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    let [x1, y1, z1] = first;
    let [x2, y2, z2] = second;
    x1 * x2 + y1 * y2 + z1 * z2
}

fn cross(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    let [x1, y1, z1] = first;
    let [x2, y2, z2] = second;
    [y1 * z2 - z1 * y2, z1 * x2 - x1 * z2, x1 * y2 - y1 * x2]
}

/// All embedded building blocks, validated.
///
/// # Errors
///
/// Returns [`Error::Data`] when the embedded JSON is invalid, has a duplicate
/// name, or fails [`Block::validate`].
pub fn blocks() -> Result<Vec<Block>, Error> {
    #[cfg(feature = "mof-data")]
    let entries = {
        let mut entries = full_blocks()?.to_vec();
        let curated: Vec<Block> = serde_json::from_str(include_str!("../data/blocks.json"))
            .map_err(|error| Error::Data(format!("block data: {error}")))?;
        for block in curated {
            if !entries
                .iter()
                .any(|entry| entry.name.eq_ignore_ascii_case(&block.name))
            {
                entries.push(block);
            }
        }
        entries
    };
    #[cfg(not(feature = "mof-data"))]
    let entries: Vec<Block> = serde_json::from_str(include_str!("../data/blocks.json"))
        .map_err(|error| Error::Data(format!("block data: {error}")))?;
    validate_blocks(&entries)?;
    Ok(entries)
}

/// The always-embedded curated eight-block subset, independent of `mof-data`.
///
/// # Errors
///
/// Returns [`Error::Data`] when the embedded JSON is invalid.
pub fn curated_blocks() -> Result<Vec<Block>, Error> {
    let entries: Vec<Block> = serde_json::from_str(include_str!("../data/blocks.json"))
        .map_err(|error| Error::Data(format!("block data: {error}")))?;
    validate_blocks(&entries)?;
    Ok(entries)
}

fn validate_blocks(entries: &[Block]) -> Result<(), Error> {
    let mut names = BTreeSet::new();
    for block in entries {
        block
            .validate()
            .map_err(|error| Error::Data(format!("block data: {error}")))?;
        if !names.insert(block.name.to_ascii_lowercase()) {
            return Err(Error::Data(format!(
                "block data repeats name `{}`",
                block.name
            )));
        }
    }
    Ok(())
}

#[cfg(feature = "mof-data")]
fn full_blocks() -> Result<&'static [Block], Error> {
    use std::io::Read;
    use std::sync::OnceLock;

    static FULL: OnceLock<Result<Vec<Block>, String>> = OnceLock::new();
    let loaded = FULL.get_or_init(|| {
        let bytes = include_bytes!("../data/pormake-blocks.json.gz");
        let mut decoder = flate2::read::GzDecoder::new(&bytes[..]);
        let mut text = String::new();
        decoder
            .read_to_string(&mut text)
            .map_err(|error| format!("block bundle: {error}"))?;
        let entries: Vec<Block> =
            serde_json::from_str(&text).map_err(|error| format!("block bundle: {error}"))?;
        Ok(entries)
    });
    match loaded {
        Ok(entries) => Ok(entries),
        Err(message) => Err(Error::Data(message.clone())),
    }
}

/// Finds a block by case-insensitive name.
///
/// # Errors
///
/// Returns the errors of [`blocks`].
pub fn find_block(name: &str) -> Result<Option<Block>, Error> {
    if let Some(block) = curated_blocks()?
        .into_iter()
        .find(|block| block.name.eq_ignore_ascii_case(name))
    {
        return Ok(Some(block));
    }
    #[cfg(feature = "mof-data")]
    if let Some(block) = full_blocks()?
        .iter()
        .find(|block| block.name.eq_ignore_ascii_case(name))
    {
        return Ok(Some(block.clone()));
    }
    Ok(None)
}
