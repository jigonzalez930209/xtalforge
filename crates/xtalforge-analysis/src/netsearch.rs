//! Searching the RCSR catalogue by space group and transitivity
//! (Task 14.4.2).
//!
//! The embedded fingerprint index ([`crate::RcsrIndex`]) names a net but does
//! not hold the facts a chemist searches by: the Hermann-Mauguin group, its
//! ITA number and the Blatov transitivity `p q r s` (kinds of vertex, edge,
//! face and tile). Those live in `data/rcsr_facts.tsv`, generated from the
//! bundled RCSR snapshots in `data/rcsr/` (see `data/README.md` in this crate).
//!
//! # Deliberate subset
//!
//! * A 2-D net answers to its **plane-group** number (1 to 17): **hcb** is 17,
//!   not the space-group number of the cell it is built in. The 3-D and 2-D
//!   numberings are kept apart, exactly as the RCSR prints them.
//! * The transitivity is the RCSR's own; `null` (unknown) matches nothing, so
//!   a query for `faces = 2` never answers with a net whose face count is not
//!   known. Layers have no tiles, so their `s` is always unknown.
//! * A group symbol the bundled symmetry catalogue cannot resolve still
//!   matches by symbol but never by number. `RcsrFactsIndex::unresolved`
//!   reports those rows rather than pretending they have a number.

#![allow(clippy::indexing_slicing)] // bounded local index and parse arrays

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

use std::collections::{BTreeMap, BTreeSet};
use std::sync::OnceLock;

use crate::netsearch_query::{NetQuery, normalize};

const FACTS: &str = include_str!("../data/rcsr_facts.tsv");

static EMBEDDED: OnceLock<Result<Vec<NetFacts>, String>> = OnceLock::new();
pub(crate) const TRANSITIVITY_FIELDS: usize = 4;

/// The 17 plane groups in ITA order, with the spellings the RCSR uses.
const PLANE_GROUPS: [(&str, u16); 17] = [
    ("p1", 1),
    ("p2", 2),
    ("pm", 3),
    ("pg", 4),
    ("cm", 5),
    ("p2mm", 6),
    ("p2mg", 7),
    ("p2gg", 8),
    ("c2mm", 9),
    ("p4", 10),
    ("p4mm", 11),
    ("p4gm", 12),
    ("p3", 13),
    ("p3m1", 14),
    ("p31m", 15),
    ("p6", 16),
    ("p6mm", 17),
];

/// Aliases for the plane-group spellings the RCSR file uses.
const PLANE_ALIASES: [(&str, &str); 7] = [
    ("pmm", "p2mm"),
    ("pmg", "p2mg"),
    ("pgg", "p2gg"),
    ("cmm", "c2mm"),
    ("p4m", "p4mm"),
    ("p4g", "p4gm"),
    ("p6m", "p6mm"),
];

/// One catalogue net and the facts it is searched by.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct NetFacts {
    /// RCSR name.
    pub name: String,
    /// Periodicity dimension (2 or 3).
    pub periodicity: u8,
    /// Distinct coordination numbers, sorted.
    pub coordination: Vec<u8>,
    /// Hermann-Mauguin group symbol as the RCSR writes it.
    pub group: String,
    /// ITA space-group number, or the plane-group number for a layer.
    pub number: Option<u16>,
    /// Blatov transitivity `p q r s`; `None` where the RCSR does not know.
    pub transitivity: [Option<u8>; TRANSITIVITY_FIELDS],
}

impl NetFacts {
    /// The transitivity as a slice of four optional counts.
    #[must_use]
    pub fn transitivity(&self) -> [Option<u8>; TRANSITIVITY_FIELDS] {
        self.transitivity
    }

    /// The one-line summary both net lists show: coordination, group and
    /// number, and the known transitivity.
    #[must_use]
    pub fn summary(&self) -> String {
        let counts = self
            .coordination
            .iter()
            .map(|value| format!("{value}-c"))
            .collect::<Vec<_>>()
            .join(", ");
        let group = if let Some(number) = self.number {
            format!("{} ({number})", self.group)
        } else {
            self.group.clone()
        };
        let mut values: Vec<String> = self
            .transitivity
            .iter()
            .map(|value| value.map_or_else(|| "?".to_owned(), |count| count.to_string()))
            .collect();
        while values.len() > 2 && values.last().is_some_and(|value| value == "?") {
            values.pop();
        }
        format!("{counts}  ·  {group}  ·  [{}]", values.join(" "))
    }
}

/// The space-group number of a symbol, or `None` when it is not recognised.
///
/// A setting suffix such as `I41/amd:2` is accepted: the origin choice is not
/// part of which group it is. A symbol that is not a three-dimensional space
/// group falls back to the 17 plane groups; a spelling shared by both, such as
/// `p6mm` (183 and 17), must be disambiguated with [`plane_group_number`].
#[must_use]
pub fn space_group_number(symbol: &str) -> Option<u16> {
    space_group_numbers()
        .get(&normalize(symbol_base(symbol)))
        .copied()
        .or_else(|| plane_group_number(symbol))
}

/// The plane-group number of a symbol (1 to 17), or `None`.
#[must_use]
pub fn plane_group_number(symbol: &str) -> Option<u16> {
    let normalized = normalize(symbol_base(symbol));
    let plane = PLANE_ALIASES
        .iter()
        .find_map(|(alias, canonical)| (*alias == normalized).then_some(*canonical))
        .unwrap_or(normalized.as_str());
    PLANE_GROUPS
        .iter()
        .find(|(name, _)| *name == plane)
        .map(|(_, number)| *number)
}

fn number_for(symbol: &str, periodicity: u8) -> Option<u16> {
    if periodicity == 2 {
        plane_group_number(symbol)
    } else {
        space_group_number(symbol)
    }
}

fn symbol_base(symbol: &str) -> &str {
    symbol.split(':').next().unwrap_or(symbol)
}

fn space_group_numbers() -> &'static BTreeMap<String, u16> {
    static NUMBERS: OnceLock<BTreeMap<String, u16>> = OnceLock::new();
    NUMBERS.get_or_init(|| {
        let mut map = BTreeMap::new();
        for group in xtalforge_symmetry::all_settings().unwrap_or_default() {
            let Ok(number) = u16::try_from(group.number) else {
                continue;
            };
            for key in [&group.hm_short, &group.hm_full, &group.hall_symbol] {
                map.entry(normalize(symbol_base(key))).or_insert(number);
            }
        }
        map
    })
}

/// The embedded RCSR facts index.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct RcsrFactsIndex {
    entries: Vec<NetFacts>,
}

impl RcsrFactsIndex {
    /// Loads the facts embedded in the crate.
    ///
    /// The embedded JSON is parsed once per process and cloned afterwards.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the embedded file is corrupt.
    pub fn from_embedded() -> Result<Self, Error> {
        let cached = EMBEDDED.get_or_init(|| parse_rows(FACTS).map_err(|error| error.to_string()));
        cached
            .clone()
            .map(|entries| Self { entries })
            .map_err(Error::InvalidInput)
    }

    /// Parses a facts file.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for malformed JSON, an unknown format,
    /// an empty or repeated name, a periodicity outside 2 to 3, or empty
    /// coordination numbers.
    pub fn parse(text: &str) -> Result<Self, Error> {
        parse_rows(text).map(|entries| Self { entries })
    }
}

fn parse_rows(text: &str) -> Result<Vec<NetFacts>, Error> {
    let mut entries = Vec::new();
    let mut names = BTreeSet::new();
    for (offset, line) in text.lines().enumerate() {
        let line = line.trim_end_matches('\r');
        if line.trim().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 8 {
            return Err(Error::InvalidInput(format!(
                "RCSR facts line {} does not have eight columns",
                offset + 1
            )));
        }
        let name = fields[0].trim();
        if name.is_empty() || !names.insert(name.to_owned()) {
            return Err(Error::InvalidInput(
                "the RCSR facts repeat or omit a name".to_owned(),
            ));
        }
        let periodicity = fields[1]
            .parse::<u8>()
            .map_err(|_| Error::InvalidInput(format!("{name} has an invalid periodicity")))?;
        if !(2..=3).contains(&periodicity) {
            return Err(Error::InvalidInput(format!(
                "{name} has periodicity {periodicity}"
            )));
        }
        let mut coordination = Vec::new();
        for value in fields[2].split(',') {
            if value.is_empty() {
                continue;
            }
            coordination.push(value.parse::<u8>().map_err(|_| {
                Error::InvalidInput(format!("{name} has an invalid coordination number"))
            })?);
        }
        if coordination.is_empty() {
            return Err(Error::InvalidInput(format!(
                "{name} has no coordination number"
            )));
        }
        coordination.sort_unstable();
        coordination.dedup();
        let group = fields[3].trim().to_owned();
        if group.is_empty() {
            return Err(Error::InvalidInput(format!("{name} has no group")));
        }
        let number = number_for(&group, periodicity);
        let mut transitivity = [None; TRANSITIVITY_FIELDS];
        for (slot, value) in transitivity.iter_mut().zip(&fields[4..8]) {
            if value.trim().is_empty() {
                continue;
            }
            let count = value
                .trim()
                .parse::<u8>()
                .map_err(|_| Error::InvalidInput(format!("{name} has an invalid transitivity")))?;
            if count == 0 {
                return Err(Error::InvalidInput(format!(
                    "{name} has a zero transitivity"
                )));
            }
            *slot = Some(count);
        }
        entries.push(NetFacts {
            name: name.to_owned(),
            periodicity,
            coordination,
            group,
            number,
            transitivity,
        });
    }
    if entries.is_empty() {
        return Err(Error::InvalidInput("the RCSR facts are empty".to_owned()));
    }
    Ok(entries)
}

impl RcsrFactsIndex {
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

    /// Every row, in catalogue order.
    #[must_use]
    pub fn entries(&self) -> &[NetFacts] {
        &self.entries
    }

    /// The row of one net.
    #[must_use]
    pub fn facts(&self, name: &str) -> Option<&NetFacts> {
        self.entries.iter().find(|entry| entry.name == name)
    }

    /// Rows whose group symbol has no resolvable number.
    #[must_use]
    pub fn unresolved(&self) -> Vec<&str> {
        self.entries
            .iter()
            .filter(|entry| entry.number.is_none())
            .map(|entry| entry.name.as_str())
            .collect()
    }

    /// Every row matching a query, in catalogue order.
    #[must_use]
    pub fn search(&self, query: &NetQuery) -> Vec<&NetFacts> {
        self.entries
            .iter()
            .filter(|entry| query.matches(entry))
            .collect()
    }

    /// The names of every row matching a query, in catalogue order.
    #[must_use]
    pub fn search_names(&self, query: &NetQuery) -> Vec<&str> {
        self.search(query)
            .into_iter()
            .map(|facts| facts.name.as_str())
            .collect()
    }
}
