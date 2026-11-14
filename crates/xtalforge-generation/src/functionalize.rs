//! Functionalisation sweeps (Task 20.4.3).

#![allow(clippy::indexing_slicing)] // fixed three-component vectors

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};
use xtalforge_mof::platform::descriptors::{FrameworkDescriptors, framework_descriptors};
use xtalforge_mof::platform::modify::{FunctionalGroup, FunctionalizeOptions, functionalize};
use xtalforge_mof::platform::porosity::{PoreDescriptors, PoreOptions, pore_descriptors};

/// One functionalisation variant.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FunctionalizationRecord {
    /// Group name.
    pub group: String,
    /// Fraction of eligible positions.
    pub fraction: f64,
    /// Functionalised framework.
    pub structure: Option<Structure>,
    /// Framework descriptors.
    pub descriptors: Option<FrameworkDescriptors>,
    /// Pore descriptors.
    pub pores: Option<PoreDescriptors>,
    /// Error message when the modification failed.
    pub error: Option<String>,
}

/// Result of a functionalisation sweep.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FunctionalizationSweep {
    /// Variants in group x fraction order.
    pub records: Vec<FunctionalizationRecord>,
    /// Variants built.
    pub built: usize,
    /// Variants that failed.
    pub failed: usize,
}

impl FunctionalizationSweep {
    /// Serializes the sweep to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self).map_err(|error| {
            Error::InvalidInput(format!(
                "could not serialize functionalization sweep: {error}"
            ))
        })
    }
}

/// Applies every group and fraction to a framework.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty group or fraction list, a
/// fraction outside `[0, 1]`, or invalid pore spacing.
pub fn functionalization_sweep(
    structure: &Structure,
    groups: &[FunctionalGroup],
    fractions: &[f64],
    seed: u64,
    descriptors: bool,
) -> Result<FunctionalizationSweep, Error> {
    if groups.is_empty() || fractions.is_empty() {
        return Err(Error::InvalidInput(
            "the sweep needs at least one group and one fraction".to_owned(),
        ));
    }
    for fraction in fractions {
        if !fraction.is_finite() || !(0.0..=1.0).contains(fraction) {
            return Err(Error::InvalidInput(format!(
                "functionalisation fractions must be in [0, 1], got {fraction}"
            )));
        }
    }
    let mut records = Vec::new();
    let mut built = 0usize;
    let mut failed = 0usize;
    for (group_index, group) in groups.iter().enumerate() {
        for (fraction_index, fraction) in fractions.iter().enumerate() {
            let options = FunctionalizeOptions {
                group: *group,
                positions: Vec::new(),
                seed: seed
                    .wrapping_add(group_index as u64)
                    .wrapping_add(fraction_index as u64 * 7919),
                fraction: *fraction,
            };
            match functionalize(structure, &options) {
                Ok(variant) => {
                    let framework = if descriptors {
                        framework_descriptors(&variant).ok()
                    } else {
                        None
                    };
                    let pore = if descriptors {
                        pore_descriptors(&variant, &PoreOptions::default())
                            .ok()
                            .map(|report| report.descriptors)
                    } else {
                        None
                    };
                    built += 1;
                    records.push(FunctionalizationRecord {
                        group: format!("{group:?}").to_lowercase(),
                        fraction: *fraction,
                        structure: Some(variant),
                        descriptors: framework,
                        pores: pore,
                        error: None,
                    });
                }
                Err(error) => {
                    failed += 1;
                    records.push(FunctionalizationRecord {
                        group: format!("{group:?}").to_lowercase(),
                        fraction: *fraction,
                        structure: None,
                        descriptors: None,
                        pores: None,
                        error: Some(error.to_string()),
                    });
                }
            }
        }
    }
    Ok(FunctionalizationSweep {
        records,
        built,
        failed,
    })
}
