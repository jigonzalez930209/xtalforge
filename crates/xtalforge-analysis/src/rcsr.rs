//! RCSR identification against the full 2019-06-01 catalogue (Tasks 5.6.1,
//! 5.6.4 and 5.6.5).
//!
//! The runtime data is `data/rcsr_index.bin`: one fingerprint hash per
//! catalogue net, generated from Crystal-Builder's `rcsr-2019-06-01.json`
//! (which stores the canonical key of every net) by
//! `examples/build_rcsr_index.rs`. The full snapshots live in
//! `data/rcsr/` at the repository root for provenance and regeneration.
//!
//! Identification uses the coordination sequence and Schläfli point symbol
//! of every distinct vertex. Those invariants name 2,665 of the 2,679
//! three-periodic nets uniquely; the nets that share both are returned as an
//! ambiguous set, exactly as the reference reports them. The canonical keys
//! are in the source index but are not recomputed from an input graph, so a
//! tie between nets with identical invariants is reported rather than
//! silently broken.

#![allow(clippy::indexing_slicing)] // bounded local edge and index arrays

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

use crate::QuotientGraph;
use crate::topology::{NetFingerprint, fingerprint};

/// Licence and provenance note for RCSR data.
///
/// The RCSR (Reticular Chemistry Structure Resource, <https://rcsr.ou.edu>)
/// is a curated database of periodic nets. The full 2019-06-01 snapshot is
/// bundled in `data/rcsr/` with the maintainer's authorisation (2026-11-09);
/// the resource itself is freely available for use by the general public.
/// Cite M. O'Keeffe, M. A. Peskov, S. J. Ramsden and O. M. Yaghi,
/// *Acc. Chem. Res.* **41** (2008) 1782.
pub const RCSR_LICENSE_NOTE: &str = "RCSR data (Reticular Chemistry Structure \
Resource, rcsr.ou.edu) are bundled as the full 2019-06-01 snapshot with the \
maintainer's authorisation (2026-11-09). The resource is freely available for \
use by the general public; cite O'Keeffe, Peskov, Ramsden and Yaghi, Acc. \
Chem. Res. 41 (2008) 1782. See data/rcsr/ and THIRD_PARTY.md.";

const MAGIC: &[u8; 8] = b"XFRCS01\n";
const INDEX: &[u8] = include_bytes!("../data/rcsr_index.bin");
const NEAREST_LIMIT: usize = 8;

/// One catalogue entry as the runtime index stores it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RcsrEntry {
    /// RCSR name (`pcu`, `dia`, ...).
    pub name: String,
    /// Periodicity dimension.
    pub periodicity: usize,
    /// Distinct coordination numbers, sorted.
    pub coordination: Vec<usize>,
}

/// What one fingerprint matched in the catalogue.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RcsrIdentification {
    /// Invariants computed from the input graph.
    pub fingerprint: NetFingerprint,
    /// Catalogue names with exactly these invariants (one is unique).
    pub names: Vec<String>,
    /// Names worth looking at when `names` is empty.
    pub nearest: Vec<String>,
    /// What the nearest names have in common (`coordination number` or empty).
    pub nearest_by: String,
}

impl RcsrIdentification {
    /// `named`, `ambiguous` or `unknown`.
    #[must_use]
    pub fn verdict(&self) -> &'static str {
        if self.names.len() == 1 {
            "named"
        } else if self.names.is_empty() {
            "unknown"
        } else {
            "ambiguous"
        }
    }

    /// The single name, when the identification is unique.
    #[must_use]
    pub fn name(&self) -> Option<&str> {
        if self.names.len() == 1 {
            self.names.first().map(String::as_str)
        } else {
            None
        }
    }
}

/// The full RCSR fingerprint index.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct RcsrIndex {
    entries: Vec<RcsrEntry>,
    /// Sorted by fingerprint hash; each row is `(hash, entry)`.
    hashes: Vec<(u32, usize)>,
    /// Salt the embedded index was built with.
    salt: u32,
}

impl RcsrIndex {
    /// Loads the index embedded in the crate.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the embedded index is corrupt.
    pub fn from_embedded() -> Result<Self, Error> {
        Self::parse(INDEX)
    }

    fn parse(bytes: &[u8]) -> Result<Self, Error> {
        if bytes.len() < 16 || bytes.get(..8) != Some(MAGIC) {
            return Err(Error::InvalidInput(
                "the RCSR index has no header".to_owned(),
            ));
        }
        let count = u32::from_le_bytes(
            bytes
                .get(8..12)
                .and_then(|slice| slice.try_into().ok())
                .ok_or_else(|| Error::InvalidInput("the RCSR index count is missing".to_owned()))?,
        ) as usize;
        let salt = u32::from_le_bytes(
            bytes
                .get(12..16)
                .and_then(|slice| slice.try_into().ok())
                .ok_or_else(|| Error::InvalidInput("the RCSR index salt is missing".to_owned()))?,
        );
        let mut cursor = 16usize;
        let mut entries = Vec::with_capacity(count);
        let mut hashes = Vec::with_capacity(count);
        for _ in 0..count {
            let name_len = *bytes
                .get(cursor)
                .ok_or_else(|| Error::InvalidInput("a short RCSR index".to_owned()))?
                as usize;
            cursor += 1;
            let name = std::str::from_utf8(
                bytes
                    .get(cursor..cursor + name_len)
                    .ok_or_else(|| Error::InvalidInput("a short RCSR name".to_owned()))?,
            )
            .map_err(|_| Error::InvalidInput("an RCSR name is not UTF-8".to_owned()))?
            .to_owned();
            cursor += name_len;
            let periodicity = *bytes
                .get(cursor)
                .ok_or_else(|| Error::InvalidInput("a short RCSR index".to_owned()))?;
            cursor += 1;
            let coordination_len = *bytes
                .get(cursor)
                .ok_or_else(|| Error::InvalidInput("a short RCSR index".to_owned()))?
                as usize;
            cursor += 1;
            let coordination: Vec<usize> = bytes
                .get(cursor..cursor + coordination_len)
                .ok_or_else(|| Error::InvalidInput("a short RCSR index".to_owned()))?
                .iter()
                .map(|value| *value as usize)
                .collect();
            cursor += coordination_len;
            let hash = u32::from_le_bytes(
                bytes
                    .get(cursor..cursor + 4)
                    .and_then(|slice| slice.try_into().ok())
                    .ok_or_else(|| Error::InvalidInput("a short RCSR index".to_owned()))?,
            );
            cursor += 4;
            let entry = RcsrEntry {
                name,
                periodicity: usize::from(periodicity),
                coordination,
            };
            let index = entries.len();
            entries.push(entry);
            hashes.push((hash, index));
        }
        hashes.sort_unstable();
        Ok(Self {
            entries,
            hashes,
            salt,
        })
    }

    /// Number of catalogue nets.
    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the index holds no nets.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Catalogue names, in index order.
    #[must_use]
    pub fn names(&self) -> Vec<&str> {
        self.entries
            .iter()
            .map(|entry| entry.name.as_str())
            .collect()
    }

    /// Catalogue entries, in index order.
    #[must_use]
    pub fn entries(&self) -> &[RcsrEntry] {
        &self.entries
    }

    /// Identifies a quotient graph against the catalogue.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`fingerprint`].
    pub fn identify(&self, graph: &QuotientGraph) -> Result<RcsrIdentification, Error> {
        let fingerprint = fingerprint(graph)?;
        let hash = fingerprint.hash_with_salt(self.salt);
        let names: Vec<String> = self
            .hashes
            .iter()
            .filter(|(candidate, _)| *candidate == hash)
            .filter_map(|(_, index)| self.entries.get(*index))
            .map(|entry| entry.name.clone())
            .collect();
        if !names.is_empty() {
            return Ok(RcsrIdentification {
                fingerprint,
                names,
                nearest: Vec::new(),
                nearest_by: String::new(),
            });
        }
        let coordination = fingerprint.coordination();
        let mut nearest: Vec<String> = self
            .entries
            .iter()
            .filter(|entry| {
                entry.periodicity == fingerprint.periodicity && entry.coordination == coordination
            })
            .map(|entry| entry.name.clone())
            .collect();
        let nearest_by = "coordination number";
        nearest.truncate(NEAREST_LIMIT);
        let nearest_by = if nearest.is_empty() {
            String::new()
        } else {
            nearest_by.to_owned()
        };
        Ok(RcsrIdentification {
            fingerprint,
            names: Vec::new(),
            nearest,
            nearest_by,
        })
    }

    /// Checks the embedded index.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty index, a repeated name, a bad
    /// periodicity or unsorted coordination numbers.
    pub fn validate(&self) -> Result<usize, Error> {
        if self.entries.is_empty() {
            return Err(Error::InvalidInput("the RCSR index is empty".to_owned()));
        }
        let mut names = std::collections::BTreeSet::new();
        for entry in &self.entries {
            if entry.name.trim().is_empty() {
                return Err(Error::InvalidInput("an RCSR entry has no name".to_owned()));
            }
            if !names.insert(&entry.name) {
                return Err(Error::InvalidInput(format!(
                    "the RCSR index repeats the name {}",
                    entry.name
                )));
            }
            if entry.periodicity > 3 {
                return Err(Error::InvalidInput(format!(
                    "{} has periodicity {}",
                    entry.name, entry.periodicity
                )));
            }
            if entry.coordination.is_empty()
                || entry.coordination.windows(2).any(|pair| pair[0] > pair[1])
            {
                return Err(Error::InvalidInput(format!(
                    "{} has invalid coordination numbers",
                    entry.name
                )));
            }
        }
        if self.hashes.len() != self.entries.len() {
            return Err(Error::InvalidInput(
                "the RCSR hash table and entries disagree".to_owned(),
            ));
        }
        Ok(self.entries.len())
    }
}

/// Number of names per fingerprint hash.
#[must_use]
pub fn shared_fingerprints(index: &RcsrIndex) -> BTreeMap<u32, usize> {
    let mut counts = BTreeMap::new();
    for (hash, _) in &index.hashes {
        *counts.entry(*hash).or_insert(0) += 1;
    }
    counts.retain(|_, count| *count > 1);
    counts
}
