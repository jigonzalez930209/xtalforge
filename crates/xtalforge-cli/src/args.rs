//! The clap command tree for the CLI.

use clap::Subcommand;

use crate::analysis_commands::AnalysisAction;
use crate::ff_commands::FfAction;
use crate::format::StructureFormat;
use crate::spacegroup;
use crate::symmetry::SymmetryAction;

use crate::build_structures;

#[derive(Debug, Subcommand)]
pub(crate) enum ProjectAction {
    /// Pack a structure into a `.xtalproj` archive.
    Pack {
        /// Input structure path.
        path: String,
        /// Output archive path.
        out: String,
    },
    /// Print the parts of a project archive.
    Info {
        /// Archive path.
        path: String,
    },
    /// Extract a project into a directory.
    Extract {
        /// Archive path.
        path: String,
        /// Destination directory.
        dir: String,
    },
}

#[derive(Debug, Subcommand)]
pub(crate) enum Command {
    /// Print the xtalforge version.
    Version,
    /// Add the hydrogens an X-ray structure is missing.
    Hydrogens {
        /// Input path, or `-` for standard input.
        path: String,
        /// X-H distance in Angstrom.
        #[arg(long, default_value_t = 1.09)]
        distance: f64,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Porosity, nets and Zeo++ analysis.
    Analysis {
        #[command(subcommand)]
        action: AnalysisAction,
    },
    /// Force-field energies, optimizers, charges and vibrations.
    Ff {
        #[command(subcommand)]
        action: FfAction,
    },
    /// List the format registry.
    Formats,
    /// List parameter presets or print one schema.
    Params {
        /// Preset name; omit to list the presets.
        name: Option<String>,
    },
    /// Pack, inspect or extract a project archive.
    Project {
        #[command(subcommand)]
        action: ProjectAction,
    },
    /// Inspect or convert a multi-frame trajectory.
    Trajectory {
        /// Input path, or `-` for standard input.
        path: String,
        /// Frame index to write as a structure.
        #[arg(long)]
        frame: Option<usize>,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Apply a JSON command list to a structure.
    Script {
        /// Script path.
        path: String,
    },
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
    /// Calculate a PXRD pattern.
    Pxrd {
        /// Input path, or `-` for standard input.
        path: String,
        /// Radiation: cu-ka, cu-ka1, mo-ka, co-ka, fe-ka, cr-ka or ag-ka.
        #[arg(long, default_value = "cu-ka")]
        radiation: String,
        /// Highest 2-theta in degrees.
        #[arg(long, default_value_t = 80.0)]
        max_two_theta: f64,
        /// Step in degrees.
        #[arg(long, default_value_t = 0.02)]
        step: f64,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
    },
    /// Build graphene or a nanotube.
    Build {
        #[command(subcommand)]
        action: build_structures::BuildAction,
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
        action: spacegroup::SpacegroupAction,
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
