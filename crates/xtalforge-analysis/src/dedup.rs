//! Dataset deduplication (Task 16.3.2).
//!
//! Duplicate structures are found with the structure matcher of Task 9.1.1.
//! A cheap fingerprint pre-filter (Task 9.3.4) runs first: records are bucketed
//! by their sorted element set, the fingerprint composition fractions must
//! agree, and a supercell relationship must have an integer atom-count and
//! volume ratio. Only candidates that survive the pre-filter reach
//! [`match_structures`].
//!
//! The pre-filter only removes pairs that the matcher cannot accept, so it
//! never changes the result. The optional `min_fingerprint_similarity` gate
//! is an extra accelerator; it is disabled by default because a large
//! displacement tolerance can lower the similarity of a true duplicate.
//!
//! The first record of each group is kept and later members are reported as
//! removed. Groups and kept records follow the input order, and the result is
//! independent of how many times the operation is run. Tolerance matching is
//! not transitive, so a later record joins the first kept representative that
//! matches it and the input order decides the representative.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, RdfOptions, StructureFingerprint, fingerprint};
use xtalforge_symmetry::{MatcherOptions, SpeciesComparator, match_structures};

use crate::dataset::Dataset;

/// Options for [`deduplicate`].
#[derive(Debug, Clone, Default, PartialEq)]
pub struct DedupOptions {
    /// Structure matcher tolerances and species comparator.
    pub matcher: MatcherOptions,
    /// Fingerprint options for the pre-filter.
    pub fingerprint: RdfOptions,
    /// Optional minimum fingerprint similarity for a matcher candidate.
    ///
    /// `None` (the default) keeps every candidate that passes the exact
    /// composition and supercell pre-filter. A value in `[0, 1]` enables the
    /// extra similarity gate; `0.0` is equivalent to `None`.
    pub min_fingerprint_similarity: Option<f64>,
}

/// One equality group of the deduplication report.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DedupGroup {
    /// Input index of the record kept as the representative.
    pub kept: usize,
    /// Input indices of the removed duplicates, in input order.
    pub removed: Vec<usize>,
}

/// Deterministic kept and removed report.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DedupReport {
    /// Number of input records.
    pub n_input: usize,
    /// Input indices of all kept records, in input order.
    pub kept: Vec<usize>,
    /// One group per kept record, in input order.
    pub groups: Vec<DedupGroup>,
}

impl DedupReport {
    /// Number of kept records.
    #[must_use]
    pub fn n_kept(&self) -> usize {
        self.kept.len()
    }

    /// Number of removed records.
    #[must_use]
    pub fn n_removed(&self) -> usize {
        self.n_input.saturating_sub(self.kept.len())
    }

    /// Sorted input indices of every removed record.
    #[must_use]
    pub fn removed(&self) -> Vec<usize> {
        let mut removed: Vec<usize> = self
            .groups
            .iter()
            .flat_map(|group| group.removed.iter().copied())
            .collect();
        removed.sort_unstable();
        removed
    }

    /// Whether no duplicate was found.
    #[must_use]
    pub fn is_duplicate_free(&self) -> bool {
        self.groups.iter().all(|group| group.removed.is_empty())
    }

    /// Applies the report to a dataset, preserving metadata and record order.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an out-of-range index or a duplicate
    /// kept index.
    pub fn apply(&self, dataset: &Dataset) -> Result<Dataset, Error> {
        let mut seen = BTreeMap::new();
        let mut output = Dataset::new(dataset.metadata().name.clone());
        *output.metadata_mut() = dataset.metadata().clone();
        for &index in &self.kept {
            if seen.insert(index, ()).is_some() {
                return Err(Error::InvalidInput(format!(
                    "deduplication report lists record {index} twice"
                )));
            }
            let Some(record) = dataset.get(index) else {
                return Err(Error::InvalidInput(format!(
                    "deduplication report index {index} is outside {} records",
                    dataset.len()
                )));
            };
            output.push(record.clone())?;
        }
        Ok(output)
    }

    /// Serializes the report as pretty JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the serializer fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string_pretty(self)
            .map_err(|error| Error::InvalidInput(format!("deduplication report encoding: {error}")))
    }
}

/// Result of [`deduplicate`].
#[derive(Debug, Clone, PartialEq)]
pub struct DedupResult {
    /// Dataset containing only the kept records, in input order.
    pub dataset: Dataset,
    /// Kept and removed groups.
    pub report: DedupReport,
}

/// Removes duplicate structures from a dataset.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, an empty structure in
/// the dataset, or a record that cannot be fingerprinted. Matcher failures are
/// reported with the record index and identifier.
pub fn deduplicate(dataset: &Dataset, options: &DedupOptions) -> Result<DedupResult, Error> {
    validate(options)?;
    let mut entries = Vec::with_capacity(dataset.len());
    for (index, record) in dataset.iter().enumerate() {
        let entry = Entry::new(&record.structure, options)
            .map_err(|error| with_context(index, &record.id, error))?;
        entries.push(entry);
    }

    let mut buckets: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    let mut report = DedupReport {
        n_input: dataset.len(),
        kept: Vec::new(),
        groups: Vec::new(),
    };
    let mut group_of: Vec<Option<usize>> = vec![None; entries.len()];

    for index in 0..entries.len() {
        let Some(entry) = entries.get(index) else {
            continue;
        };
        let key = bucket_key(entry, options.matcher.comparator);
        let candidates = buckets.get(&key).cloned().unwrap_or_default();
        let mut found = None;
        for representative in candidates {
            let (Some(candidate), Some(record)) =
                (entries.get(representative), dataset.get(representative))
            else {
                continue;
            };
            if !prefilter(candidate, entry, options) {
                continue;
            }
            let Some(current) = dataset.get(index) else {
                continue;
            };
            let outcome = match_structures(&record.structure, &current.structure, &options.matcher)
                .map_err(|error| with_context(index, &current.id, error))?;
            if outcome.matched {
                found = Some(representative);
                break;
            }
        }
        match found {
            Some(representative) => {
                if let Some(group) = group_of
                    .get(representative)
                    .copied()
                    .flatten()
                    .and_then(|group| report.groups.get_mut(group))
                {
                    group.removed.push(index);
                }
            }
            None => {
                buckets.entry(key).or_default().push(index);
                report.kept.push(index);
                if let Some(slot) = group_of.get_mut(index) {
                    *slot = Some(report.groups.len());
                }
                report.groups.push(DedupGroup {
                    kept: index,
                    removed: Vec::new(),
                });
            }
        }
    }

    let output = report.apply(dataset)?;
    Ok(DedupResult {
        dataset: output,
        report,
    })
}

struct Entry {
    fingerprint: StructureFingerprint,
    atom_count: usize,
    volume: f64,
    composition: BTreeMap<String, f64>,
}

impl Entry {
    fn new(structure: &xtalforge_core::Structure, options: &DedupOptions) -> Result<Self, Error> {
        let p1 = structure.p1()?;
        if p1.atoms.is_empty() {
            return Err(Error::InvalidInput(
                "deduplication needs at least one atom per record".to_owned(),
            ));
        }
        let fingerprint = fingerprint(structure, &options.fingerprint)?;
        let mut composition: BTreeMap<String, f64> = BTreeMap::new();
        for (symbol, fraction) in &fingerprint.composition {
            let symbol = xtalforge_core::elements::element_by_symbol(symbol)
                .map(|element| element.symbol.clone())
                .unwrap_or_else(|_| symbol.clone());
            *composition.entry(symbol).or_insert(0.0) += fraction;
        }
        Ok(Self {
            fingerprint,
            atom_count: p1.atoms.len(),
            volume: p1.lattice.volume(),
            composition,
        })
    }
}

fn validate(options: &DedupOptions) -> Result<(), Error> {
    if let Some(threshold) = options.min_fingerprint_similarity
        && (!threshold.is_finite() || !(0.0..=1.0).contains(&threshold))
    {
        return Err(Error::InvalidInput(format!(
            "minimum fingerprint similarity must be in [0, 1], got {threshold}"
        )));
    }
    if !options.fingerprint.r_max.is_finite() || options.fingerprint.r_max <= 0.0 {
        return Err(Error::InvalidInput(
            "fingerprint range must be positive and finite".to_owned(),
        ));
    }
    if !options.fingerprint.bin_width.is_finite() || options.fingerprint.bin_width <= 0.0 {
        return Err(Error::InvalidInput(
            "fingerprint bin width must be positive and finite".to_owned(),
        ));
    }
    Ok(())
}

fn bucket_key(entry: &Entry, comparator: SpeciesComparator) -> String {
    if comparator == SpeciesComparator::Ignore {
        return "*".to_owned();
    }
    entry
        .composition
        .keys()
        .cloned()
        .collect::<Vec<String>>()
        .join("\u{1f}")
}

fn prefilter(first: &Entry, second: &Entry, options: &DedupOptions) -> bool {
    let (small, large) = if first.atom_count <= second.atom_count {
        (first, second)
    } else {
        (second, first)
    };
    if small.atom_count == 0 {
        return false;
    }
    if large.atom_count != small.atom_count {
        if !options.matcher.allow_supercell {
            return false;
        }
        let ratio = large.atom_count as f64 / small.atom_count as f64;
        let factor = ratio.round();
        if (ratio - factor).abs() > 1e-9 || factor < 1.0 || factor > 64.0 {
            return false;
        }
        let volume_ratio = large.volume / small.volume;
        if (volume_ratio - factor).abs() > 0.5 * factor {
            return false;
        }
    }
    if options.matcher.comparator != SpeciesComparator::Ignore
        && !composition_close(&first.composition, &second.composition, 1e-9)
    {
        return false;
    }
    if let Some(threshold) = options.min_fingerprint_similarity
        && threshold > 0.0
        && first.fingerprint.similarity(&second.fingerprint) < threshold
    {
        return false;
    }
    true
}

fn composition_close(
    first: &BTreeMap<String, f64>,
    second: &BTreeMap<String, f64>,
    tolerance: f64,
) -> bool {
    if first.len() != second.len() {
        return false;
    }
    first.iter().all(|(symbol, fraction)| {
        second
            .get(symbol)
            .is_some_and(|other| (fraction - other).abs() <= tolerance)
    })
}

fn with_context(index: usize, id: &str, error: Error) -> Error {
    Error::InvalidInput(format!("dataset record {index} ({id:?}): {error}"))
}
