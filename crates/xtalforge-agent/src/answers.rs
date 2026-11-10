//! Typed answers (Task 8.1.5).
//!
//! Every verb returns one [`Answer`]. The enum is the printable and
//! serializable form: `Answer::to_json` produces the object the CLI prints in
//! `--json` mode and the MCP server returns to a client.

use serde::Serialize;
use serde_json::Value;

use crate::capabilities::Capabilities;
use crate::changes::ChangeBurst;
use crate::diagnostics::DiagnosticReport;
use crate::inspect::Inspection;

/// One document as the session describes it.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct DocumentSummary {
    /// Session name.
    pub name: String,
    /// Structure name.
    pub title: String,
    /// Formula.
    pub formula: String,
    /// Number of asymmetric sites.
    pub n_sites: usize,
    /// Cell volume in cubic Angstrom.
    pub volume: f64,
    /// Density in g/cm3.
    pub density: f64,
    /// Where it came from, when it was loaded from a file.
    pub source: Option<String>,
    /// Whether it has unsaved edits.
    pub dirty: bool,
}

/// One perceived bond.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct BondSummary {
    /// First atom in the P1 cell.
    pub i: usize,
    /// Second atom in the P1 cell.
    pub j: usize,
    /// Lattice image of the second atom.
    pub image: [i32; 3],
    /// Distance in Angstrom.
    pub distance: f64,
    /// Bond order.
    pub order: f64,
    /// Stored explicitly rather than perceived.
    pub explicit: bool,
}

/// One entry of the command log.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct HistoryEntry {
    /// Description of the edit.
    pub description: String,
    /// `user`, `script` or `agent`.
    pub origin: String,
}

/// Ring counts by size.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct RingCount {
    /// Ring size in atoms.
    pub size: usize,
    /// Number of primitive rings of this size.
    pub count: usize,
}

/// A sampled powder pattern point.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PatternSample {
    /// `2 theta` in degrees.
    pub two_theta: f64,
    /// Normalised intensity, maximum 100.
    pub intensity: f64,
}

/// One symmetry summary.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SymmetrySummary {
    /// ITA number.
    pub number: u32,
    /// Hermann-Mauguin symbol.
    pub hm_symbol: String,
    /// Pearson symbol.
    pub pearson_symbol: String,
    /// Number of symmetry operations in the input cell.
    pub n_operations: usize,
    /// Number of atoms in the input cell.
    pub n_atoms: usize,
    /// Wyckoff letter of each input atom.
    pub wyckoffs: Vec<char>,
    /// Site-symmetry symbol of each input atom.
    pub site_symmetries: Vec<String>,
}

/// Porosity numbers read off the distance grid.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PorositySummary {
    /// Grid spacing in Angstrom.
    pub spacing: f64,
    /// Probe radius in Angstrom.
    pub probe_radius: f64,
    /// Cell volume in cubic Angstrom.
    pub cell_volume: f64,
    /// Volume reachable by the probe.
    pub accessible_volume: f64,
    /// Accessible fraction of the cell.
    pub accessible_fraction: f64,
    /// Accessible surface area in square Angstrom.
    pub surface_area: f64,
    /// Largest included sphere diameter in Angstrom.
    pub largest_included_diameter: f64,
}

/// Result of a geometry optimization.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct OptimizationSummary {
    /// Final energy in eV.
    pub energy_ev: f64,
    /// Whether the run met its tolerances.
    pub converged: bool,
    /// Steps taken.
    pub steps: usize,
    /// Largest remaining force in eV/Angstrom.
    pub max_force_ev_per_angstrom: f64,
}

/// One verb answer.
#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Answer {
    /// `capabilities`.
    Capabilities {
        /// The probed report.
        capabilities: Capabilities,
    },
    /// `list_documents` and `close`.
    Documents {
        /// Open documents.
        documents: Vec<DocumentSummary>,
        /// Active document name.
        active: Option<String>,
    },
    /// `load`.
    Loaded {
        /// The new document.
        document: DocumentSummary,
    },
    /// `save`.
    Saved {
        /// Document name.
        document: String,
        /// Written path.
        path: String,
        /// Format name.
        format: String,
        /// Bytes written.
        bytes: usize,
    },
    /// `info`.
    Info {
        /// Document name.
        document: String,
        /// Summary, cell, composition and site table.
        inspection: Inspection,
        /// What is wrong with it.
        diagnostics: DiagnosticReport,
    },
    /// `diagnostics`.
    Diagnostics {
        /// Document name.
        document: String,
        /// Findings.
        report: DiagnosticReport,
    },
    /// A mutation that leaves a structure.
    Structure {
        /// Document name.
        document: String,
        /// Its new summary.
        summary: DocumentSummary,
        /// Description of the edit.
        description: String,
    },
    /// `symmetry`.
    Symmetry {
        /// Document name.
        document: String,
        /// Detected symmetry.
        symmetry: SymmetrySummary,
    },
    /// `bonds`.
    Bonds {
        /// Document name.
        document: String,
        /// Atoms in the P1 cell.
        n_atoms: usize,
        /// Perceived bonds.
        bonds: Vec<BondSummary>,
    },
    /// `rings`.
    Rings {
        /// Document name.
        document: String,
        /// Largest ring searched.
        max_size: usize,
        /// Counts by size, smallest first.
        counts: Vec<RingCount>,
        /// Total number of primitive rings.
        total: usize,
    },
    /// `pxrd`.
    Pxrd {
        /// Document name.
        document: String,
        /// Radiation name.
        radiation: String,
        /// Number of allowed reflections.
        n_reflections: usize,
        /// First reflections with relative intensity.
        reflections: Vec<[f64; 2]>,
        /// Normalised profile, decimated for transport.
        points: Vec<PatternSample>,
    },
    /// `porosity`.
    Porosity {
        /// Document name.
        document: String,
        /// The numbers.
        porosity: PorositySummary,
    },
    /// `uff_energy`.
    Energy {
        /// Document name.
        document: String,
        /// Energy in eV.
        energy_ev: f64,
        /// Atoms in the P1 cell.
        n_atoms: usize,
        /// Largest force in eV/Angstrom.
        max_force_ev_per_angstrom: f64,
        /// Stress in GPa, when periodic.
        stress_gpa: Option<[[f64; 3]; 3]>,
    },
    /// `optimize`.
    Optimization {
        /// Document name.
        document: String,
        /// The run.
        optimization: OptimizationSummary,
        /// New structure summary.
        summary: crate::inspect::StructureSummaryView,
    },
    /// `render`.
    Picture {
        /// Document name.
        document: String,
        /// `svg` or `png`.
        format: String,
        /// Width in pixels.
        width: u32,
        /// Height in pixels.
        height: u32,
        /// Written path, when one was asked for.
        path: Option<String>,
        /// Encoded bytes, when no path was asked for.
        bytes: Option<usize>,
    },
    /// `undo`, `redo` and `history`.
    History {
        /// Document name.
        document: String,
        /// Whether an undo or redo actually changed anything.
        applied: bool,
        /// The command log.
        entries: Vec<HistoryEntry>,
        /// Whether an undo is available.
        can_undo: bool,
        /// Whether a redo is available.
        can_redo: bool,
    },
    /// `changes`.
    Changes {
        /// Document name.
        document: String,
        /// Bursts of agent edits.
        bursts: Vec<ChangeBurst>,
    },
}

impl Answer {
    /// Serializes the answer to one JSON value.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::Json`] if a number is not finite.
    pub fn to_json(&self) -> Result<Value, crate::Error> {
        serde_json::to_value(self).map_err(|error| crate::Error::Json(error.to_string()))
    }

    /// The `kind` tag of the answer.
    #[must_use]
    pub fn kind(&self) -> &'static str {
        match self {
            Self::Capabilities { .. } => "capabilities",
            Self::Documents { .. } => "documents",
            Self::Loaded { .. } => "loaded",
            Self::Saved { .. } => "saved",
            Self::Info { .. } => "info",
            Self::Diagnostics { .. } => "diagnostics",
            Self::Structure { .. } => "structure",
            Self::Symmetry { .. } => "symmetry",
            Self::Bonds { .. } => "bonds",
            Self::Rings { .. } => "rings",
            Self::Pxrd { .. } => "pxrd",
            Self::Porosity { .. } => "porosity",
            Self::Energy { .. } => "energy",
            Self::Optimization { .. } => "optimization",
            Self::Picture { .. } => "picture",
            Self::History { .. } => "history",
            Self::Changes { .. } => "changes",
        }
    }
}
