#![allow(clippy::indexing_slicing)] // bounded 3-vectors and six-element permutations

//! Hypothetical MOFs at scale (Task 20.4).
//!
//! * [`global_orientation_fit`] is the public form of the Kabsch orientation
//!   used by the builder: every one-to-one assignment of connection points to
//!   slot directions is tried and the smallest RMSD wins.
//! * [`generate_mof_dataset`] enumerates topology x node x linker
//!   combinations from the PORMAKE bundle (the curated subset by default, the
//!   complete bundle with the `mof-data` feature), builds each framework,
//!   optionally relaxes it with UFF4MOF and records descriptors.
//! * [`crate::functionalization_sweep`] applies the Phase 10 post-synthetic
//!   modification to every group and fraction and records the descriptors of
//!   each variant.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};
use xtalforge_ff::{CalculatorOptions, Minimizer, OptimizeOptions, UffCalculator, minimize};
use xtalforge_mof::platform::descriptors::FrameworkDescriptors;
use xtalforge_mof::platform::descriptors::framework_descriptors;
use xtalforge_mof::platform::porosity::PoreDescriptors;
use xtalforge_mof::platform::porosity::{PoreOptions, pore_descriptors};
use xtalforge_mof::{BlockKind, MofRequest, blocks, build_mof, topologies};

use crate::rng::SplitMix64;

/// A global orientation fit.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OrientationFit {
    /// Proper rotation matrix.
    pub rotation: [[f64; 3]; 3],
    /// Connection index assigned to each target direction.
    pub assignment: Vec<usize>,
    /// Root-mean-square distance after the fit, in Angstrom.
    pub rmsd: f64,
}

/// Finds the rotation and assignment that minimise the RMSD between unit
/// connection directions and unit target directions.
///
/// The search is exhaustive for up to six connections (720 permutations) and
/// deterministic: ties keep the first lexicographic assignment. This is the
/// global fit the MOF builder uses (Task 20.4.1).
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the two lists differ in length, are
/// empty, have more than six entries, or contain a zero-length direction.
pub fn global_orientation_fit(
    connections: &[[f64; 3]],
    targets: &[[f64; 3]],
) -> Result<OrientationFit, Error> {
    let count = connections.len();
    if count == 0 || count != targets.len() || count > 6 {
        return Err(Error::InvalidInput(
            "the orientation fit needs one to six connection and target vectors".to_owned(),
        ));
    }
    let sources: Vec<[f64; 3]> = connections
        .iter()
        .map(|point| normalize(*point))
        .collect::<Option<Vec<_>>>()
        .ok_or_else(|| Error::InvalidInput("a connection point has zero length".to_owned()))?;
    let targets: Vec<[f64; 3]> = targets
        .iter()
        .map(|point| normalize(*point))
        .collect::<Option<Vec<_>>>()
        .ok_or_else(|| Error::InvalidInput("a target direction has zero length".to_owned()))?;
    let mut permutation: Vec<usize> = (0..count).collect();
    let mut best: Option<(f64, [[f64; 3]; 3], Vec<usize>)> = None;
    loop {
        let rotation = horn(&sources, &targets, &permutation);
        let mut sum = 0.0;
        for (target, source_index) in targets.iter().zip(&permutation) {
            let Some(source) = sources.get(*source_index) else {
                continue;
            };
            let moved = apply(rotation, *source);
            sum += (0..3)
                .map(|axis| (moved[axis] - target[axis]).powi(2))
                .sum::<f64>();
        }
        let rmsd = (sum / count as f64).sqrt();
        if best
            .as_ref()
            .is_none_or(|(value, _, _)| rmsd < *value - 1e-12)
        {
            best = Some((rmsd, rotation, permutation.clone()));
        }
        if !next_permutation(&mut permutation) {
            break;
        }
    }
    let (rmsd, rotation, assignment) = best.ok_or_else(|| {
        Error::InvalidInput("the orientation fit found no valid assignment".to_owned())
    })?;
    Ok(OrientationFit {
        rotation,
        assignment,
        rmsd,
    })
}

fn next_permutation(values: &mut [usize]) -> bool {
    if values.len() < 2 {
        return false;
    }
    let mut pivot = values.len() - 2;
    loop {
        if values[pivot] < values[pivot + 1] {
            break;
        }
        if pivot == 0 {
            return false;
        }
        pivot -= 1;
    }
    let mut successor = values.len() - 1;
    while values[successor] <= values[pivot] {
        successor -= 1;
    }
    values.swap(pivot, successor);
    values[pivot + 1..].reverse();
    true
}

fn normalize([x, y, z]: [f64; 3]) -> Option<[f64; 3]> {
    let norm = (x * x + y * y + z * z).sqrt();
    if norm < 1e-12 {
        None
    } else {
        Some([x / norm, y / norm, z / norm])
    }
}

fn apply(matrix: [[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    matrix.map(|row| row[0] * vector[0] + row[1] * vector[1] + row[2] * vector[2])
}

/// Horn's quaternion form of the Kabsch fit for one permutation.
fn horn(sources: &[[f64; 3]], targets: &[[f64; 3]], permutation: &[usize]) -> [[f64; 3]; 3] {
    let mut s = [[0.0; 3]; 3];
    for (target, source) in targets.iter().zip(permutation) {
        if let Some(source) = sources.get(*source) {
            for a in 0..3 {
                for b in 0..3 {
                    s[a][b] += source[a] * target[b];
                }
            }
        }
    }
    let [[sxx, sxy, sxz], [syx, syy, syz], [szx, szy, szz]] = s;
    let n = [
        [sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
        [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
        [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
        [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
    ];
    let shift = n.iter().flatten().map(|value| value.abs()).sum::<f64>() + 1.0;
    let mut q = [1.0, 0.31, 0.17, 0.07];
    for _ in 0..500 {
        let mut next = [0.0; 4];
        for (row, value) in next.iter_mut().enumerate() {
            *value = shift * q[row] + (0..4).map(|column| n[row][column] * q[column]).sum::<f64>();
        }
        let length = next.iter().map(|value| value * value).sum::<f64>().sqrt();
        if length < 1e-300 {
            break;
        }
        q = next.map(|value| value / length);
    }
    let [w, x, y, z] = q;
    [
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y - z * w),
            2.0 * (x * z + y * w),
        ],
        [
            2.0 * (x * y + z * w),
            1.0 - 2.0 * (x * x + z * z),
            2.0 * (y * z - x * w),
        ],
        [
            2.0 * (x * z - y * w),
            2.0 * (y * z + x * w),
            1.0 - 2.0 * (x * x + y * y),
        ],
    ]
}

/// Options for [`generate_mof_dataset`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MofDatasetOptions {
    /// Topology names; empty means every topology in the bundle.
    pub topologies: Vec<String>,
    /// Node block names; empty means every node.
    pub nodes: Vec<String>,
    /// Linker block names; empty means every linker.
    pub linkers: Vec<String>,
    /// Maximum number of frameworks to build.
    pub limit: usize,
    /// Whether to relax each framework with UFF4MOF.
    pub relax: bool,
    /// Maximum optimizer steps.
    pub max_steps: usize,
    /// Whether to compute framework and pore descriptors.
    pub descriptors: bool,
    /// Pore-grid spacing in Angstrom.
    pub pore_spacing: f64,
    /// Explicit seed for the combination order.
    pub seed: u64,
}

impl Default for MofDatasetOptions {
    fn default() -> Self {
        Self {
            topologies: Vec::new(),
            nodes: Vec::new(),
            linkers: Vec::new(),
            limit: 32,
            relax: false,
            max_steps: 100,
            descriptors: true,
            pore_spacing: 1.0,
            seed: 0,
        }
    }
}

/// One generated MOF record.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MofGenerationRecord {
    /// Topology name.
    pub topology: String,
    /// Node name.
    pub node: String,
    /// Linker name.
    pub linker: String,
    /// The built framework, when assembly succeeded.
    pub structure: Option<Structure>,
    /// UFF energy in eV, when evaluated.
    pub energy_ev: Option<f64>,
    /// Energy per atom in eV.
    pub energy_per_atom: Option<f64>,
    /// Framework descriptors.
    pub descriptors: Option<FrameworkDescriptors>,
    /// Pore descriptors.
    pub pores: Option<PoreDescriptors>,
    /// Error message when assembly failed.
    pub error: Option<String>,
    /// Whether this is a duplicate of an earlier record.
    pub duplicate: bool,
}

/// Result of a combinatorial MOF dataset run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MofDatasetReport {
    /// Records in enumeration order.
    pub records: Vec<MofGenerationRecord>,
    /// Combinations enumerated before the limit.
    pub combinations: usize,
    /// Frameworks built.
    pub built: usize,
    /// Combinations that failed.
    pub failed: usize,
    /// Duplicates detected by MOFkey.
    pub duplicates: usize,
    /// Whether the complete PORMAKE bundle was active.
    pub full_bundle: bool,
}

impl MofDatasetReport {
    /// Serializes the report to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self).map_err(|error| {
            Error::InvalidInput(format!("could not serialize MOF dataset report: {error}"))
        })
    }
}

/// Enumerates and builds MOF combinations.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options and the errors of the
/// catalogue. A framework that cannot be assembled is recorded with its error
/// instead of aborting the dataset.
pub fn generate_mof_dataset(options: &MofDatasetOptions) -> Result<MofDatasetReport, Error> {
    if options.limit == 0 {
        return Err(Error::InvalidInput(
            "the dataset limit must be at least one".to_owned(),
        ));
    }
    if !options.pore_spacing.is_finite() || options.pore_spacing <= 0.0 {
        return Err(Error::InvalidInput(
            "pore spacing must be positive and finite".to_owned(),
        ));
    }
    let topologies = topologies().map_err(|error| Error::InvalidInput(error.to_string()))?;
    let blocks = blocks().map_err(|error| Error::InvalidInput(error.to_string()))?;
    let mut topology_names: Vec<String> = if options.topologies.is_empty() {
        topologies
            .iter()
            .map(|topology| topology.name.clone())
            .collect()
    } else {
        options.topologies.clone()
    };
    let mut node_names: Vec<String> = if options.nodes.is_empty() {
        blocks
            .iter()
            .filter(|block| block.kind == BlockKind::Node)
            .map(|block| block.name.clone())
            .collect()
    } else {
        options.nodes.clone()
    };
    let mut linker_names: Vec<String> = if options.linkers.is_empty() {
        blocks
            .iter()
            .filter(|block| block.kind == BlockKind::Linker)
            .map(|block| block.name.clone())
            .collect()
    } else {
        options.linkers.clone()
    };
    let mut rng = SplitMix64::new(options.seed);
    rng.shuffle(&mut topology_names);
    rng.shuffle(&mut node_names);
    rng.shuffle(&mut linker_names);
    let mut records = Vec::new();
    let mut built = 0usize;
    let mut failed = 0usize;
    let mut combinations = 0usize;
    let calculator = if options.relax {
        Some(UffCalculator::new(CalculatorOptions::default()).map_err(crate::ff_error)?)
    } else {
        None
    };
    let optimizer = OptimizeOptions {
        max_steps: options.max_steps,
        ..OptimizeOptions::default()
    };
    let mut seen_keys: BTreeMap<String, usize> = BTreeMap::new();
    'outer: for topology in &topology_names {
        for node in &node_names {
            for linker in &linker_names {
                if records.len() >= options.limit {
                    break 'outer;
                }
                combinations += 1;
                let request = MofRequest::new(topology.clone(), node.clone(), linker.clone());
                let built_mof = match build_mof(&request) {
                    Ok(built_mof) => built_mof,
                    Err(error) => {
                        failed += 1;
                        records.push(MofGenerationRecord {
                            topology: topology.clone(),
                            node: node.clone(),
                            linker: linker.clone(),
                            structure: None,
                            energy_ev: None,
                            energy_per_atom: None,
                            descriptors: None,
                            pores: None,
                            error: Some(error.to_string()),
                            duplicate: false,
                        });
                        continue;
                    }
                };
                let mut structure = built_mof.structure;
                let mut energy_ev = None;
                let mut energy_per_atom = None;
                if let Some(calculator) = &calculator {
                    match minimize(calculator, &structure, &optimizer, Minimizer::Fire) {
                        Ok(result) => {
                            energy_per_atom =
                                Some(result.energy / structure.p1()?.atoms.len().max(1) as f64);
                            energy_ev = Some(result.energy);
                            structure = result.structure;
                        }
                        Err(error) => {
                            failed += 1;
                            records.push(MofGenerationRecord {
                                topology: topology.clone(),
                                node: node.clone(),
                                linker: linker.clone(),
                                structure: Some(structure),
                                energy_ev: None,
                                energy_per_atom: None,
                                descriptors: None,
                                pores: None,
                                error: Some(error.to_string()),
                                duplicate: false,
                            });
                            continue;
                        }
                    }
                }
                let descriptors = if options.descriptors {
                    framework_descriptors(&structure).ok()
                } else {
                    None
                };
                let pores = if options.descriptors {
                    pore_descriptors(
                        &structure,
                        &PoreOptions {
                            spacing: options.pore_spacing,
                            ..PoreOptions::default()
                        },
                    )
                    .ok()
                    .map(|report| report.descriptors)
                } else {
                    None
                };
                let key = xtalforge_mof::mofid(&structure)
                    .map(|id| id.mofkey)
                    .unwrap_or_else(|_| format!("{topology}/{node}/{linker}"));
                let duplicate = seen_keys.contains_key(&key);
                if !duplicate {
                    seen_keys.insert(key, records.len());
                }
                built += 1;
                records.push(MofGenerationRecord {
                    topology: topology.clone(),
                    node: node.clone(),
                    linker: linker.clone(),
                    structure: Some(structure),
                    energy_ev,
                    energy_per_atom,
                    descriptors,
                    pores,
                    error: None,
                    duplicate,
                });
            }
        }
    }
    let duplicates = records.iter().filter(|record| record.duplicate).count();
    Ok(MofDatasetReport {
        records,
        combinations,
        built,
        failed,
        duplicates,
        full_bundle: cfg!(feature = "mof-data"),
    })
}
