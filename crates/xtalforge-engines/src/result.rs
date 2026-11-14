//! Typed calculation results shared by every engine parser (Task 18.2.1).
//!
//! A [`CalculationResult`] keeps the final structure, the per-ionic-step
//! trajectory, energies, forces, stress, magnetic moments, charges,
//! convergence flags, electronic bands, densities of states, orbital
//! projections, timings and parser warnings. Every parser in [`crate::parse`]
//! returns this model, so consumers do not need to know which code produced
//! the output.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::Structure;

use crate::Error;

/// Convergence flags reported by an engine.
///
/// The `*_reported` flags distinguish "the run says it did not converge" from
/// "the output never states it". Parsers set a warning when a flag is missing.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Convergence {
    /// Self-consistent field converged.
    pub scf: bool,
    /// Ionic relaxation or optimisation converged.
    pub ionic: bool,
    /// The output explicitly reported the SCF state.
    #[serde(default)]
    pub scf_reported: bool,
    /// The output explicitly reported the ionic state.
    #[serde(default)]
    pub ionic_reported: bool,
}

/// One electronic or ionic step energy record.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ElectronicStep {
    /// Energy in eV.
    pub energy: f64,
    /// Fermi energy in eV, when reported.
    pub fermi: Option<f64>,
    /// Maximum force in eV/Angstrom, when reported.
    pub max_force: Option<f64>,
}

/// One timing entry.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Timing {
    /// Stage name.
    pub stage: String,
    /// Seconds.
    pub seconds: f64,
}

/// One parser warning.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Warning {
    /// Stable code.
    pub code: String,
    /// Message.
    pub message: String,
}

/// One ionic or molecular-dynamics snapshot.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Snapshot {
    /// Zero-based step index in file order.
    pub step: u64,
    /// Structure at this step.
    pub structure: Structure,
    /// Potential energy in eV, when reported.
    #[serde(default)]
    pub energy: Option<f64>,
    /// Forces in eV/Angstrom in site order, when reported.
    #[serde(default)]
    pub forces: Option<Vec<[f64; 3]>>,
    /// Stress tensor in GPa, when reported.
    #[serde(default)]
    pub stress: Option<[[f64; 3]; 3]>,
}

/// Band energies along a k-point path.
///
/// `energies[k][b]` is band `b` at k-point `k` in eV, `occupations[k][b]` its
/// occupation. `kpoints` and `weights` may be empty when the source file does
/// not print the k-point coordinates (VASP `vasprun.xml` eigenvalues, for
/// example).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct BandData {
    /// K-point coordinates in crystal coordinates when printed.
    #[serde(default)]
    pub kpoints: Vec<[f64; 3]>,
    /// K-point weights, parallel to `kpoints`.
    #[serde(default)]
    pub weights: Vec<f64>,
    /// Band energies in eV, one row per k-point.
    #[serde(default)]
    pub energies: Vec<Vec<f64>>,
    /// Occupations, one row per k-point.
    #[serde(default)]
    pub occupations: Vec<Vec<f64>>,
    /// Fermi energy in eV, when reported.
    #[serde(default)]
    pub fermi: Option<f64>,
}

/// One density-of-states channel.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DosChannel {
    /// Channel label (`total`, `spin_up`, `spin_down`, an element or orbital).
    pub label: String,
    /// Values on the energy grid.
    pub values: Vec<f64>,
}

/// Total and integrated density of states.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct DosData {
    /// Energy grid in eV.
    #[serde(default)]
    pub energies: Vec<f64>,
    /// Total DOS values.
    #[serde(default)]
    pub total: Vec<f64>,
    /// Integrated total DOS values.
    #[serde(default)]
    pub integrated: Vec<f64>,
    /// Spin or projected channels.
    #[serde(default)]
    pub channels: Vec<DosChannel>,
    /// Fermi energy in eV, when reported.
    #[serde(default)]
    pub fermi: Option<f64>,
}

/// One atom/orbital projection weight.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OrbitalWeight {
    /// Zero-based site index.
    pub atom: usize,
    /// Orbital label as printed by the engine (`s`, `p`, `dxy`, ...).
    pub orbital: String,
    /// Projection weight in electrons.
    pub weight: f64,
}

/// Projected weights for one band at one k-point.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BandProjection {
    /// K-point coordinates in crystal coordinates.
    pub kpoint: [f64; 3],
    /// Zero-based band index in file order.
    pub band: usize,
    /// Band energy in eV.
    pub energy: f64,
    /// Occupation.
    pub occupation: f64,
    /// Per-site, per-orbital weights.
    pub weights: Vec<OrbitalWeight>,
}

/// Orbital projections read from `PROCAR` or `projwfc`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ProjectionData {
    /// One record per k-point and band.
    #[serde(default)]
    pub bands: Vec<BandProjection>,
}

/// A typed result shared by every parser.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct CalculationResult {
    /// Final structure, when present.
    pub structure: Option<Structure>,
    /// Energies in eV, in file order.
    pub energies: Vec<f64>,
    /// Electronic steps.
    pub steps: Vec<ElectronicStep>,
    /// Convergence.
    pub convergence: Convergence,
    /// Timings.
    pub timings: Vec<Timing>,
    /// Warnings.
    pub warnings: Vec<Warning>,
    /// Extra scalar fields.
    pub scalars: BTreeMap<String, f64>,
    /// Per-step structures and their forces and stress.
    #[serde(default)]
    pub trajectory: Vec<Snapshot>,
    /// Final forces in eV/Angstrom in site order, when reported.
    #[serde(default)]
    pub forces: Option<Vec<[f64; 3]>>,
    /// Final stress tensor in GPa, when reported.
    #[serde(default)]
    pub stress: Option<[[f64; 3]; 3]>,
    /// Per-site magnetic moments in Bohr magnetons, when reported.
    #[serde(default)]
    pub magnetic_moments: Vec<f64>,
    /// Per-site charges in electrons, when reported (for example Bader).
    #[serde(default)]
    pub charges: Vec<f64>,
    /// Fermi energy in eV, when reported.
    #[serde(default)]
    pub fermi: Option<f64>,
    /// Electronic bands, when the file contains them.
    #[serde(default)]
    pub bands: Option<BandData>,
    /// Total and projected density of states, when the file contains it.
    #[serde(default)]
    pub dos: Option<DosData>,
    /// Orbital projections, when the file contains them.
    #[serde(default)]
    pub projections: Option<ProjectionData>,
    /// Producing code (`vasp`, `quantum_espresso`, `cp2k`, `lammps`, `gulp`).
    #[serde(default)]
    pub engine: Option<String>,
}

impl CalculationResult {
    /// Creates an empty result tagged with the producing engine.
    #[must_use]
    pub fn new(engine: &str) -> Self {
        Self {
            engine: Some(engine.to_owned()),
            ..Self::default()
        }
    }

    /// Appends a warning.
    pub fn push_warning(&mut self, code: &str, message: impl Into<String>) {
        self.warnings.push(Warning {
            code: code.to_owned(),
            message: message.into(),
        });
    }

    /// Appends a timing entry, ignoring negative or non-finite values.
    pub fn push_timing(&mut self, stage: &str, seconds: f64) {
        if seconds.is_finite() && seconds >= 0.0 {
            self.timings.push(Timing {
                stage: stage.to_owned(),
                seconds,
            });
        }
    }

    /// Last energy in eV, when any energy was parsed.
    #[must_use]
    pub fn final_energy(&self) -> Option<f64> {
        self.energies.last().copied()
    }

    /// Copies the final trajectory entry into the top-level structure, forces
    /// and stress fields when those are still empty.
    pub fn close_trajectory(&mut self) {
        if let Some(last) = self.trajectory.last() {
            if self.structure.is_none() {
                self.structure = Some(last.structure.clone());
            }
            if self.forces.is_none() {
                self.forces = last.forces.clone();
            }
            if self.stress.is_none() {
                self.stress = last.stress;
            }
        }
    }
}

/// One Bader charge entry from an `ACF.dat` file.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BaderEntry {
    /// One-based site index printed by the Bader program.
    pub site: usize,
    /// Partial charge in electrons.
    pub charge: f64,
    /// Distance to the nearest Bader maximum in Angstrom, when printed.
    pub min_distance: Option<f64>,
    /// Atomic volume in cubic Angstrom, when printed.
    pub volume: Option<f64>,
}

/// Reads a Bader `ACF.dat` file into its per-site table.
///
/// Rows are returned sorted by the one-based `site` index. The reader rejects
/// indices that are not contiguous from 1, so the returned vector can be
/// mapped directly onto a structure in site order:
/// `charges[i]` belongs to site `i`.
///
/// # Errors
///
/// Returns [`Error::Parse`] when the file has no charge rows or the indices
/// are not contiguous from 1.
pub fn bader_table(text: &str) -> Result<Vec<BaderEntry>, Error> {
    let mut entries = Vec::new();
    for line in text.lines() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        let Some(first) = fields.first() else {
            continue;
        };
        let Ok(site) = first.parse::<usize>() else {
            continue;
        };
        let Some(charge) = fields.get(4).and_then(|value| number(value)) else {
            continue;
        };
        entries.push(BaderEntry {
            site,
            charge,
            min_distance: fields.get(5).and_then(|value| number(value)),
            volume: fields.get(6).and_then(|value| number(value)),
        });
    }
    if entries.is_empty() {
        return Err(Error::Parse(
            "Bader ACF file has no per-site charge rows".to_owned(),
        ));
    }
    entries.sort_by_key(|entry| entry.site);
    for (index, entry) in entries.iter().enumerate() {
        if entry.site != index + 1 {
            return Err(Error::Parse(format!(
                "Bader ACF site indices are not contiguous from 1 (found {})",
                entry.site
            )));
        }
    }
    Ok(entries)
}

/// Reads Bader charges from an `ACF.dat` file.
///
/// The returned vector follows the one-based site index of the file, so
/// `charges[i]` belongs to site `i` of the structure.
///
/// # Errors
///
/// Returns the errors of [`bader_table`].
pub fn bader_charges(text: &str) -> Result<Vec<f64>, Error> {
    Ok(bader_table(text)?
        .into_iter()
        .map(|entry| entry.charge)
        .collect())
}

/// One DDEC charge entry.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DdecEntry {
    /// One-based row index in file order.
    pub site: usize,
    /// Element symbol printed on the row, when present.
    pub element: Option<String>,
    /// Net atomic charge in electrons.
    pub charge: f64,
    /// Spin moment in Bohr magnetons, when printed.
    pub spin: Option<f64>,
}

/// Reads a DDEC charge file into its per-site table.
///
/// Two layouts produced by Chargemol are accepted:
///
/// 1. the `DDEC6_even_tempered_net_atomic_charges.xyz` file (an atom count, a
///    comment line, then `element x y z charge [spin]` rows);
/// 2. a columnar table whose header names a net charge column, with
///    `index element charge [spin]` rows.
///
/// Rows are returned in file order with a one-based `site` index, so
/// `charges[i]` belongs to site `i` of the structure when the file follows
/// the structure atom order (the DDEC convention).
///
/// # Errors
///
/// Returns [`Error::Parse`] when neither layout can be read.
pub fn ddec_table(text: &str) -> Result<Vec<DdecEntry>, Error> {
    let lines: Vec<&str> = text
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect();
    if let Some(count) = lines.first().and_then(|line| line.parse::<usize>().ok())
        && let Some(entries) = ddec_xyz(&lines, count)
    {
        return Ok(entries);
    }
    if let Some(entries) = ddec_columns(&lines) {
        return Ok(entries);
    }
    Err(Error::Parse(
        "DDEC file has no recognisable charge rows".to_owned(),
    ))
}

/// Reads DDEC charges in file order (one-based site `i` is vector index `i`).
///
/// # Errors
///
/// Returns the errors of [`ddec_table`].
pub fn ddec_charges(text: &str) -> Result<Vec<f64>, Error> {
    Ok(ddec_table(text)?
        .into_iter()
        .map(|entry| entry.charge)
        .collect())
}

fn ddec_xyz(lines: &[&str], count: usize) -> Option<Vec<DdecEntry>> {
    let start = 2usize;
    let rows = lines.get(start..)?;
    if rows.len() < count {
        return None;
    }
    let mut entries = Vec::with_capacity(count);
    for (index, line) in rows.iter().take(count).enumerate() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        let element = fields.first().map(|value| (*value).to_owned());
        let charge = fields.get(4).and_then(|value| number(value))?;
        entries.push(DdecEntry {
            site: index + 1,
            element,
            charge,
            spin: fields.get(5).and_then(|value| number(value)),
        });
    }
    Some(entries)
}

fn ddec_columns(lines: &[&str]) -> Option<Vec<DdecEntry>> {
    let header = lines.iter().position(|line| {
        let lower = line.to_ascii_lowercase();
        lower.contains("charge") && (lower.contains("atom") || lower.contains("net"))
    })?;
    let mut entries = Vec::new();
    for line in lines.get(header + 1..)? {
        let fields: Vec<&str> = line.split_whitespace().collect();
        let Some(first) = fields.first() else {
            continue;
        };
        let Ok(index) = first.parse::<usize>() else {
            if entries.is_empty() {
                continue;
            }
            break;
        };
        let Some(charge) = fields.get(2).and_then(|value| number(value)) else {
            continue;
        };
        entries.push(DdecEntry {
            site: index,
            element: fields.get(1).map(|value| (*value).to_owned()),
            charge,
            spin: fields.get(3).and_then(|value| number(value)),
        });
    }
    (!entries.is_empty()).then_some(entries)
}

fn number(text: &str) -> Option<f64> {
    let cleaned: String = text
        .trim()
        .trim_matches(',')
        .chars()
        .map(|character| match character {
            'D' | 'd' => 'E',
            other => other,
        })
        .collect();
    cleaned.parse().ok()
}
