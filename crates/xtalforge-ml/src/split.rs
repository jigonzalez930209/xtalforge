#![allow(clippy::indexing_slicing)] // bounded indices from validated loops and fixed-size vectors

//! Leakage-safe dataset splits (Task 21.2.2).
//!
//! Four strategies are available:
//!
//! * `random` — a seeded shuffle of records.
//! * `composition` — whole composition groups go to one split, so no
//!   composition appears on both sides.
//! * `cluster` — structures are grouped with the Phase 9 matcher and whole
//!   clusters go to one split.
//! * `prototype` — the space group and the sorted Wyckoff signature form the
//!   group key.
//!
//! Fractions are turned into integer counts with largest-remainder rounding,
//! so they always add up to the dataset size. Groups are shuffled with the
//! explicit seed; the assignment is deterministic.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_analysis::dataset::Dataset;
use xtalforge_core::Error;
use xtalforge_symmetry::{
    MatcherOptions, SymmetryOptions, detect_symmetry, match_structures, wyckoff_assignments,
};

use crate::rng::SplitMix64;

/// Which grouping the split uses.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SplitStrategy {
    /// Independent records.
    Random,
    /// Whole composition groups.
    Composition,
    /// Whole structure-matcher clusters.
    Cluster {
        /// Matcher tolerances.
        matcher: MatcherOptions,
    },
    /// Whole space-group/Wyckoff prototypes.
    Prototype,
}

/// Options for [`split_dataset`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SplitOptions {
    /// Training fraction.
    pub train: f64,
    /// Validation fraction.
    pub validation: f64,
    /// Test fraction.
    pub test: f64,
    /// Explicit seed.
    pub seed: u64,
    /// Grouping strategy.
    pub strategy: SplitStrategy,
}

impl Default for SplitOptions {
    fn default() -> Self {
        Self {
            train: 0.6,
            validation: 0.2,
            test: 0.2,
            seed: 0,
            strategy: SplitStrategy::Random,
        }
    }
}

/// Split label for one record.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SplitLabel {
    /// Training.
    Train,
    /// Validation.
    Validation,
    /// Test.
    Test,
}

/// Assignment of every input index.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SplitAssignment {
    /// Label for each dataset index.
    pub labels: Vec<SplitLabel>,
    /// Group key for each index (empty for random splits).
    pub groups: Vec<String>,
    /// Number of groups.
    pub n_groups: usize,
}

/// A summary of the split sizes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct SplitSummary {
    /// Training records.
    pub train: usize,
    /// Validation records.
    pub validation: usize,
    /// Test records.
    pub test: usize,
}

impl SplitAssignment {
    /// Counts per split.
    #[must_use]
    pub fn summary(&self) -> SplitSummary {
        let mut summary = SplitSummary {
            train: 0,
            validation: 0,
            test: 0,
        };
        for label in &self.labels {
            match label {
                SplitLabel::Train => summary.train += 1,
                SplitLabel::Validation => summary.validation += 1,
                SplitLabel::Test => summary.test += 1,
            }
        }
        summary
    }

    /// Indices assigned to one split.
    #[must_use]
    pub fn indices(&self, label: SplitLabel) -> Vec<usize> {
        self.labels
            .iter()
            .enumerate()
            .filter(|(_, value)| **value == label)
            .map(|(index, _)| index)
            .collect()
    }
}

/// Splits a dataset.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid fractions, an empty dataset,
/// and the symmetry or matcher errors of the grouping strategies.
pub fn split_dataset(dataset: &Dataset, options: &SplitOptions) -> Result<SplitAssignment, Error> {
    validate(options)?;
    let n = dataset.len();
    if n == 0 {
        return Err(Error::InvalidInput(
            "cannot split an empty dataset".to_owned(),
        ));
    }
    let groups = group_keys(dataset, &options.strategy)?;
    let mut unique: Vec<String> = groups
        .iter()
        .cloned()
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    unique.sort();
    let mut rng = SplitMix64::new(options.seed);
    rng.shuffle(&mut unique);
    let targets = split_counts(n, options);
    let n_groups = unique.len();
    let preferred: Vec<SplitLabel> = (0..n_groups)
        .map(|index| {
            if index * 3 < n_groups {
                SplitLabel::Train
            } else if index * 3 < 2 * n_groups {
                SplitLabel::Validation
            } else {
                SplitLabel::Test
            }
        })
        .collect();
    let mut group_label: BTreeMap<String, SplitLabel> = BTreeMap::new();
    let mut counts = SplitSummary {
        train: 0,
        validation: 0,
        test: 0,
    };
    for (group, label) in unique.iter().zip(preferred) {
        let size = groups.iter().filter(|key| *key == group).count();
        let fits = match label {
            SplitLabel::Train => counts.train + size <= targets.train,
            SplitLabel::Validation => counts.validation + size <= targets.validation,
            SplitLabel::Test => counts.test + size <= targets.test,
        };
        let label = if fits {
            label
        } else {
            // Place the whole group in the split with the most room.
            let room = [
                (
                    SplitLabel::Train,
                    targets.train.saturating_sub(counts.train),
                ),
                (
                    SplitLabel::Validation,
                    targets.validation.saturating_sub(counts.validation),
                ),
                (SplitLabel::Test, targets.test.saturating_sub(counts.test)),
            ];
            room.iter()
                .max_by_key(|(_, room)| *room)
                .map(|(label, _)| *label)
                .unwrap_or(SplitLabel::Train)
        };
        match label {
            SplitLabel::Train => counts.train += size,
            SplitLabel::Validation => counts.validation += size,
            SplitLabel::Test => counts.test += size,
        }
        group_label.insert(group.clone(), label);
    }
    let labels: Vec<SplitLabel> = groups
        .iter()
        .map(|group| group_label.get(group).copied().unwrap_or(SplitLabel::Train))
        .collect();
    Ok(SplitAssignment {
        labels,
        groups,
        n_groups,
    })
}

/// Serializes a split summary.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when serialization fails.
pub fn split_summary(assignment: &SplitAssignment) -> Result<String, Error> {
    serde_json::to_string(&assignment.summary())
        .map_err(|error| Error::InvalidInput(format!("could not serialize split: {error}")))
}

fn split_counts(n: usize, options: &SplitOptions) -> SplitSummary {
    let total = options.train + options.validation + options.test;
    let raw = [
        options.train / total * n as f64,
        options.validation / total * n as f64,
        options.test / total * n as f64,
    ];
    let mut counts = raw.map(|value| value.floor() as usize);
    let mut remainder = n.saturating_sub(counts.iter().sum());
    let mut order: Vec<usize> = (0..3).collect();
    order.sort_by(|left, right| {
        let left_fraction = raw[*left] - raw[*left].floor();
        let right_fraction = raw[*right] - raw[*right].floor();
        right_fraction
            .total_cmp(&left_fraction)
            .then_with(|| left.cmp(right))
    });
    for index in order {
        if remainder == 0 {
            break;
        }
        counts[index] += 1;
        remainder -= 1;
    }
    SplitSummary {
        train: counts[0],
        validation: counts[1],
        test: counts[2],
    }
}

fn group_keys(dataset: &Dataset, strategy: &SplitStrategy) -> Result<Vec<String>, Error> {
    match strategy {
        SplitStrategy::Random => Ok((0..dataset.len()).map(|index| index.to_string()).collect()),
        SplitStrategy::Composition => {
            let mut keys = Vec::with_capacity(dataset.len());
            for record in dataset.iter() {
                let p1 = record.structure.p1()?;
                let mut counts: BTreeMap<String, usize> = BTreeMap::new();
                for atom in &p1.atoms {
                    *counts.entry(atom.element.clone()).or_default() += 1;
                }
                keys.push(
                    counts
                        .into_iter()
                        .map(|(element, count)| format!("{element}{count}"))
                        .collect::<Vec<_>>()
                        .join(" "),
                );
            }
            Ok(keys)
        }
        SplitStrategy::Cluster { matcher } => {
            let records: Vec<&xtalforge_core::Structure> =
                dataset.iter().map(|record| &record.structure).collect();
            let mut keys = Vec::with_capacity(records.len());
            let mut representatives: Vec<usize> = Vec::new();
            for (index, structure) in records.iter().enumerate() {
                let mut assigned = None;
                for representative in &representatives {
                    if let Some(other) = records.get(*representative)
                        && match_structures(structure, other, matcher)
                            .map(|result| result.matched)
                            .unwrap_or(false)
                    {
                        assigned = Some(*representative);
                        break;
                    }
                }
                let key = assigned.unwrap_or_else(|| {
                    representatives.push(index);
                    index
                });
                keys.push(format!("cluster-{key}"));
            }
            Ok(keys)
        }
        SplitStrategy::Prototype => {
            let mut keys = Vec::with_capacity(dataset.len());
            for record in dataset.iter() {
                let dataset = detect_symmetry(&record.structure, SymmetryOptions::default())?;
                let assignments =
                    wyckoff_assignments(&record.structure, SymmetryOptions::default())?;
                let mut signature: Vec<String> = assignments
                    .iter()
                    .map(|assignment| format!("{}{}", assignment.multiplicity, assignment.letter))
                    .collect();
                signature.sort();
                signature.dedup();
                keys.push(format!("{}-{}", dataset.number, signature.join(" ")));
            }
            Ok(keys)
        }
    }
}

fn validate(options: &SplitOptions) -> Result<(), Error> {
    if !options.train.is_finite()
        || !options.validation.is_finite()
        || !options.test.is_finite()
        || options.train < 0.0
        || options.validation < 0.0
        || options.test < 0.0
    {
        return Err(Error::InvalidInput(
            "split fractions must be finite and non-negative".to_owned(),
        ));
    }
    let sum = options.train + options.validation + options.test;
    if (sum - 1.0).abs() > 1e-9 || sum <= 0.0 {
        return Err(Error::InvalidInput(
            "split fractions must add up to one".to_owned(),
        ));
    }
    Ok(())
}
