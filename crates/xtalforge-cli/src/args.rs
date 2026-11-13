//! The clap command tree for the CLI.

use clap::Subcommand;

use crate::analysis_commands::AnalysisAction;
use crate::feature_args::FeatureAction;
use crate::fetch_commands::FetchAction;
use crate::ff_commands::FfAction;
use crate::format::StructureFormat;
use crate::mof_commands::MofAction;
use crate::powder_commands::PowderAction;
use crate::spacegroup;
use crate::symmetry::SymmetryAction;

use crate::build_structures;

pub(crate) use crate::actions::{ModelAction, ProjectAction, SkillAction};

#[derive(Debug, Subcommand)]
pub(crate) enum Command {
    /// Print the xtalforge version.
    Version,
    /// Report the features, engines, data bundles and verbs of this installation.
    Capabilities,
    /// Diagnose the installation and optionally one structure.
    Doctor {
        /// Structure to check, when one is supplied.
        path: Option<String>,
    },
    /// Draw a structure to SVG or PNG.
    Render {
        /// Input path, or `-` for standard input.
        path: String,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Picture format.
        #[arg(long, value_enum, default_value_t = crate::format::PictureFormat::Svg)]
        format: crate::format::PictureFormat,
        /// Width in pixels.
        #[arg(long, default_value_t = 800)]
        width: u32,
        /// Height in pixels.
        #[arg(long, default_value_t = 600)]
        height: u32,
    },
    /// Run a JSON list of agent verbs against one session.
    Agent {
        /// Script path, or `-` for standard input.
        path: String,
    },
    /// Run a declarative workflow over many inputs.
    Run(crate::workflow_commands::RunArgs),
    /// Serve the Model Context Protocol over standard input and output.
    Mcp,
    /// Skill documents for an assistant: path, list, show or install.
    Skill {
        #[command(subcommand)]
        action: SkillAction,
    },
    /// Manage the optional machine-learning model cache.
    Models {
        #[command(subcommand)]
        action: ModelAction,
    },
    /// Print a shell completion script.
    Completions {
        /// Shell to generate for.
        #[arg(value_enum)]
        shell: clap_complete::Shell,
    },
    /// Print the roff manual page.
    Man,
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
    /// Powder refinement: indexing, Le Bail/Pawley, Rietveld and cells.
    Powder {
        #[command(subcommand)]
        action: PowderAction,
    },
    /// MOF topologies, blocks and framework assembly.
    Mof {
        #[command(subcommand)]
        action: MofAction,
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
    /// Fetch a structure from OPTIMADE, COD or the Materials Project.
    Fetch {
        #[command(subcommand)]
        action: FetchAction,
    },
    /// Print the xtalforge provenance embedded in a file.
    Provenance {
        /// Input path, or `-` for standard input.
        path: String,
    },
    /// List the module registry.
    Modules,
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
    /// Decide whether two structures are the same.
    Compare {
        /// First structure.
        first: String,
        /// Second structure.
        second: String,
        /// Fractional length tolerance.
        #[arg(long, default_value_t = 0.2)]
        ltol: f64,
        /// Site tolerance in Angstrom.
        #[arg(long, default_value_t = 0.3)]
        stol: f64,
        /// Angle tolerance in degrees.
        #[arg(long, default_value_t = 5.0)]
        angle_tol: f64,
        /// Do not accept a supercell relationship.
        #[arg(long)]
        no_supercell: bool,
    },
    /// Group structures into equality classes.
    Group {
        /// Structure paths.
        paths: Vec<String>,
        /// Fractional length tolerance.
        #[arg(long, default_value_t = 0.2)]
        ltol: f64,
        /// Site tolerance in Angstrom.
        #[arg(long, default_value_t = 0.3)]
        stol: f64,
    },
    /// Apply a strain, deformation gradient, volume or density.
    Deform {
        /// Input structure.
        path: String,
        /// Six engineering Voigt components: xx,yy,zz,yz,xz,xy.
        #[arg(long)]
        voigt: Option<String>,
        /// Nine deformation-gradient components.
        #[arg(long)]
        gradient: Option<String>,
        /// Target cell volume in cubic Angstrom.
        #[arg(long)]
        volume: Option<f64>,
        /// Target density in g/cm3.
        #[arg(long)]
        density: Option<f64>,
        /// Output path.
        #[arg(long)]
        out: Option<String>,
        /// Output format.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Core crystal features: matcher edits, local environments, surfaces,
    /// defects, orderings, SQS, interfaces, elasticity and thermodynamics.
    Features {
        #[command(subcommand)]
        action: FeatureAction,
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
    /// Convert a structure between the registered formats: cif, poscar, xyz, pmg-json, cssr, gen, lammps-data, shelx and pdb.
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
        /// Caglioti U V W width coefficients.
        #[arg(long, num_args = 3, value_names = ["U", "V", "W"])]
        caglioti: Option<Vec<f64>>,
        /// Pseudo-Voigt mixing in [0, 1].
        #[arg(long)]
        eta: Option<f64>,
        /// Profile shape: gaussian, lorentzian or pseudo-voigt.
        #[arg(long)]
        shape: Option<String>,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
    },
    /// Parse a SMILES string and embed it in 3D.
    Molecule {
        /// SMILES string.
        #[arg(long)]
        smiles: String,
        /// Deterministic embedder seed.
        #[arg(long, default_value_t = 7)]
        seed: u64,
        /// Embedder gradient steps.
        #[arg(long, default_value_t = 1500)]
        iterations: usize,
        /// Output path; defaults to standard output.
        #[arg(long)]
        out: Option<String>,
        /// Output format.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// List the fragment library or embed one entry.
    Fragment {
        /// Fragment name; omit to list the library.
        name: Option<String>,
        /// Output path for an embedded fragment.
        #[arg(long)]
        out: Option<String>,
        /// Output format.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Round-trip a SMILES through the sketch model.
    Sketch {
        /// SMILES string.
        smiles: String,
    },
    /// Replace one site after a clearance check.
    Substitute {
        /// Input structure.
        path: String,
        /// Asymmetric-unit site index.
        index: usize,
        /// New element.
        element: String,
        /// Output path.
        #[arg(long)]
        out: Option<String>,
        /// Output format.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Place a guest molecule after a clearance check.
    Guest {
        /// Host structure.
        host: String,
        /// Guest SMILES string.
        #[arg(long)]
        smiles: String,
        /// Cartesian offset x,y,z.
        #[arg(long, default_value = "0,0,0")]
        offset: String,
        /// Minimum host-guest clearance in Angstrom.
        #[arg(long, default_value_t = 1.0)]
        clearance: f64,
        /// Output path.
        #[arg(long)]
        out: Option<String>,
        /// Output format.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Build graphene, nanotube, C20, C60, twisted bilayer, nanoribbon, ZTC, schwarzite, disordered carbon or polymer.
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
        /// Emit the full scene geometry as JSON instead of a mesh.
        #[arg(long)]
        scene: bool,
        /// Mesh format to write.
        #[arg(long, value_enum, default_value_t = crate::format::MeshFormat::Stl)]
        mesh_format: crate::format::MeshFormat,
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
