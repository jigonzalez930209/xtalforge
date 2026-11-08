//! Command-line interface for `xtalforge`.
//!
//! Phase 1 provides structure input and conversion, supercells, slabs,
//! connectivity reports and geometry queries. Every command accepts the global
//! flags, emits one JSON value in `--json` mode and reports runtime failures
//! with exit code 3.

mod analyze;
mod basic;
mod build;
mod convert;
mod dispatch;
mod format;
mod input;
mod measure;
mod mesh;
mod output;
mod spacegroup;
mod symmetry;
mod symmetry_analysis;

use clap::{Parser, Subcommand, ValueEnum};
use std::process::ExitCode;

use format::StructureFormat;

/// Exit code used for every runtime failure.
const ERROR_EXIT_CODE: u8 = 3;

#[derive(Debug, Parser)]
#[command(
    name = "xtalforge",
    version,
    about = "Build, manipulate, analyse and export crystal structures",
    disable_help_subcommand = true
)]
struct Cli {
    /// Emit machine-readable JSON instead of human-readable text.
    #[arg(long, global = true)]
    json: bool,
    /// Suppress headers and non-essential output.
    #[arg(long, global = true)]
    quiet: bool,
    /// Turn loss warnings into errors.
    #[arg(long, global = true)]
    strict: bool,
    /// Worker threads (reserved for the parallel algorithms).
    #[arg(long, global = true, value_name = "N")]
    threads: Option<usize>,
    /// Limit policy used by the size checks.
    #[arg(long, global = true, value_enum, default_value_t = Limits::Standard)]
    limits: Limits,
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Print the xtalforge version.
    Version,
    /// Print reference data for an element.
    Element {
        /// Element symbol or atomic number.
        input: String,
    },
    /// Print a summary of a structure or structure JSON document.
    Info {
        /// Input path, or `-` for standard input.
        path: String,
    },
    /// Build a cell from its six parameters.
    Cell {
        /// a, b, c in Angstrom and alpha, beta, gamma in degrees.
        #[arg(long, num_args = 6, value_names = ["a", "b", "c", "alpha", "beta", "gamma"])]
        params: Vec<f64>,
    },
    /// Compute the volume of a cell from its six parameters.
    CellVolume {
        /// Cell length a in Angstrom.
        a: f64,
        /// Cell length b in Angstrom.
        b: f64,
        /// Cell length c in Angstrom.
        c: f64,
        /// Angle alpha in degrees.
        alpha: f64,
        /// Angle beta in degrees.
        beta: f64,
        /// Angle gamma in degrees.
        gamma: f64,
    },
    /// Convert a structure between CIF, POSCAR and XYZ.
    Convert {
        /// Input path, or `-` for standard input.
        path: String,
        /// Output format.
        #[arg(long, value_enum)]
        format: StructureFormat,
        /// Output path, or `-` for standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
    },
    /// Build a printable STL mesh.
    Mesh {
        /// Input path, or `-` for standard input.
        path: String,
        /// Draw style: ball-stick, space-fill, wireframe or sticks.
        #[arg(long, default_value = "ball-stick")]
        style: String,
        /// Atom radius multiplier.
        #[arg(long, default_value_t = 1.0)]
        scale: f64,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
    },
    /// Build an na x nb x nc supercell.
    Supercell {
        /// Input path, or `-` for standard input.
        path: String,
        /// Repeats along a.
        na: u32,
        /// Repeats along b.
        nb: u32,
        /// Repeats along c.
        nc: u32,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Cut a slab along the (h k l) plane.
    Slab {
        /// Input path, or `-` for standard input.
        path: String,
        /// First Miller index.
        h: i32,
        /// Second Miller index.
        k: i32,
        /// Third Miller index.
        l: i32,
        /// Number of layers.
        #[arg(long, value_name = "N")]
        layers: u32,
        /// Vacuum above the slab in Angstrom.
        #[arg(long, value_name = "A")]
        vacuum: f64,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// List the perceived bonds of the P1 cell.
    Bonds {
        /// Input path, or `-` for standard input.
        path: String,
    },
    /// List the chordless rings of the P1 cell.
    Rings {
        /// Input path, or `-` for standard input.
        path: String,
        /// Largest ring size in atoms.
        #[arg(long, value_name = "N", default_value_t = 6)]
        max_size: usize,
    },
    /// List the matched functional groups.
    Groups {
        /// Input path, or `-` for standard input.
        path: String,
    },
    /// List the periodic neighbors within a cutoff.
    Neighbors {
        /// Input path, or `-` for standard input.
        path: String,
        /// Cutoff radius in Angstrom.
        #[arg(long, value_name = "A")]
        cutoff: f64,
    },
    /// Minimum-image distance between two P1 atoms.
    Distance {
        /// Input path, or `-` for standard input.
        path: String,
        /// First atom index.
        i: usize,
        /// Second atom index.
        j: usize,
    },
    /// Angle i-j-k in degrees.
    Angle {
        /// Input path, or `-` for standard input.
        path: String,
        /// First atom index.
        i: usize,
        /// Pivot atom index.
        j: usize,
        /// Third atom index.
        k: usize,
    },
    /// Browse the space-group catalogue.
    Spacegroup {
        #[command(subcommand)]
        action: SpacegroupAction,
    },
    /// Detect the space group and Wyckoff positions.
    Symmetry {
        #[command(subcommand)]
        action: SymmetryAction,
    },
    /// Torsion i-j-k-l in degrees.
    Torsion {
        /// Input path, or `-` for standard input.
        path: String,
        /// First atom index.
        i: usize,
        /// Second atom index.
        j: usize,
        /// Third atom index.
        k: usize,
        /// Fourth atom index.
        l: usize,
    },
}

#[derive(Debug, Subcommand)]
enum SpacegroupAction {
    /// Show one space group by symbol, Hall symbol or number.
    Show {
        /// H-M symbol, Hall symbol or number.
        symbol: String,
    },
    /// List catalogue entries.
    List {
        /// Filter by crystal system.
        #[arg(long)]
        system: Option<String>,
        /// Filter by centering letter.
        #[arg(long)]
        centering: Option<char>,
        /// Filter by point-group symbol.
        #[arg(long)]
        point_group: Option<String>,
    },
    /// Search by symbol or number substring.
    Search {
        /// Text to search for.
        text: String,
    },
}

#[derive(Debug, Subcommand)]
enum SymmetryAction {
    /// Detect the symmetry of a structure.
    Detect {
        /// Input path, or `-` for standard input.
        path: String,
        /// Position tolerance in Angstrom.
        #[arg(long, value_name = "A", default_value_t = 1e-3)]
        symprec: f64,
        /// Angle tolerance in degrees.
        #[arg(long, value_name = "DEG")]
        angle_tolerance: Option<f64>,
        /// Preferred Hall number (1 to 530).
        #[arg(long, value_name = "N")]
        hall_number: Option<u32>,
    },
    /// Expand a structure to P1.
    P1 {
        /// Input path, or `-` for standard input.
        path: String,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Reduce a structure to its detected asymmetric unit.
    Asymmetrize {
        /// Input path, or `-` for standard input.
        path: String,
        /// Position tolerance in Angstrom.
        #[arg(long, value_name = "A", default_value_t = 1e-3)]
        symprec: f64,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// List the translationengleiche subgroups of a Hall setting.
    Subgroups {
        /// Hall number (1 to 530).
        hall_number: u32,
        /// Comparison tolerance.
        #[arg(long, default_value_t = 1e-5)]
        epsilon: f64,
    },
    /// Lower a structure to a subgroup number.
    Lower {
        /// Input path, or `-` for standard input.
        path: String,
        /// Target space-group number.
        subgroup: u32,
        /// Matching tolerance in Angstrom.
        #[arg(long, default_value_t = 0.5)]
        tolerance: f64,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Build the conventional standardized cell.
    Standardize {
        /// Input path, or `-` for standard input.
        path: String,
        /// Position tolerance in Angstrom.
        #[arg(long, default_value_t = 1e-3)]
        symprec: f64,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Build the standardized primitive cell.
    Primitive {
        /// Input path, or `-` for standard input.
        path: String,
        /// Position tolerance in Angstrom.
        #[arg(long, default_value_t = 1e-3)]
        symprec: f64,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Reduce the lattice with Niggli, Delaunay or Minkowski reduction.
    Reduce {
        /// Input path, or `-` for standard input.
        path: String,
        /// Reduction algorithm.
        #[arg(long, default_value = "niggli")]
        algorithm: String,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Print Wyckoff assignments.
    Wyckoff {
        /// Input path, or `-` for standard input.
        path: String,
        /// Position tolerance in Angstrom.
        #[arg(long, default_value_t = 1e-3)]
        symprec: f64,
    },
    /// Detect the molecular point group (finite molecules).
    Pointgroup {
        /// Input path, or `-` for standard input.
        path: String,
        /// Position tolerance in Angstrom.
        #[arg(long, default_value_t = 1e-3)]
        tolerance: f64,
    },
    /// Print the band path of a structure.
    Kpath {
        /// Input path, or `-` for standard input.
        path: String,
    },
    /// Change the setting of a structure to a Hall number.
    Setting {
        /// Input path, or `-` for standard input.
        path: String,
        /// Hall number (1 to 530).
        hall_number: u32,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Detect the collinear magnetic space group (experimental).
    Magnetic {
        /// Input path, or `-` for standard input.
        path: String,
        /// Position tolerance in Angstrom.
        #[arg(long, default_value_t = 1e-3)]
        symprec: f64,
    },
    /// Detect the layer group of a slab.
    Layer {
        /// Input path, or `-` for standard input.
        path: String,
        /// Position tolerance in Angstrom.
        #[arg(long, default_value_t = 1e-3)]
        symprec: f64,
    },
    /// Print the reciprocal-space Brillouin zone.
    Brillouin {
        /// Input path, or `-` for standard input.
        path: String,
    },
    /// Detect across a tolerance ladder.
    Scan {
        /// Input path, or `-` for standard input.
        path: String,
    },
    /// Report per-atom deviations from a space group.
    Diagnose {
        /// Input path, or `-` for standard input.
        path: String,
        /// H-M symbol, Hall symbol or number.
        symbol: String,
    },
    /// Apply a space group to a structure.
    Set {
        /// Input path, or `-` for standard input.
        path: String,
        /// H-M symbol, Hall symbol or number.
        symbol: String,
        /// `reinterpret` or `impose`.
        #[arg(long, default_value = "reinterpret")]
        mode: String,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Limits {
    Standard,
    Generous,
    WarnOnly,
}

impl Limits {
    fn profile(self) -> xtalforge::core::Profile {
        match self {
            Self::Standard => xtalforge::core::Profile::Standard,
            Self::Generous => xtalforge::core::Profile::Generous,
            Self::WarnOnly => xtalforge::core::Profile::WarnOnly,
        }
    }
}

#[derive(Debug)]
pub(crate) struct CliError {
    pub(crate) code: &'static str,
    pub(crate) message: String,
}

impl CliError {
    pub(crate) fn new(code: &'static str, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
        }
    }

    pub(crate) fn input(message: impl Into<String>) -> Self {
        Self::new("invalid_input", message)
    }

    pub(crate) fn core(error: xtalforge::core::Error) -> Self {
        Self::new(error.code(), error.message())
    }
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match dispatch::run(&cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            dispatch::report_error(&cli, &error);
            ExitCode::from(ERROR_EXIT_CODE)
        }
    }
}
