//! Arguments for the Phase 20 `generate` command.

use clap::{Args, Subcommand};

use crate::format::StructureFormat;

/// Structure generation and search.
#[derive(Debug, Subcommand)]
pub(crate) enum GenerateAction {
    /// Generate a Wyckoff-constrained random crystal.
    Crystal(CrystalArgs),
    /// Generate one crystal per space group in a range.
    CrystalRange(CrystalRangeArgs),
    /// Place a rigid molecule on Wyckoff positions.
    Molecule(MoleculeArgs),
    /// List the layer groups.
    LayerGroups,
    /// Generate a layer crystal.
    Layer(LayerArgs),
    /// List the rod groups.
    RodGroups,
    /// Generate a rod crystal.
    Rod(RodArgs),
    /// List the prototype library.
    Prototypes,
    /// Build a prototype by label.
    Prototype(PrototypeArgs),
    /// Identify the prototype label of a structure.
    Identify(IdentifyArgs),
    /// Propose chemically similar substitutions.
    Substitute(SubstituteArgs),
    /// Run an AIRSS-style random search.
    Search(SearchArgs),
    /// Run an evolutionary search.
    Evolve(EvolveArgs),
    /// Emit a checkpointed search workflow document.
    Workflow(WorkflowArgs),
    /// Enumerate hypothetical MOF combinations.
    MofDataset(MofDatasetArgs),
    /// Functionalise a framework over groups and fractions.
    Functionalize(FunctionalizeArgs),
}

/// Shared crystal-generation options.
#[derive(Debug, Args)]
pub(crate) struct CrystalArgs {
    /// Space-group number (1-230).
    #[arg(long, default_value_t = 225)]
    pub space_group: u32,
    /// Composition as `Element:count,Element:count`.
    #[arg(long)]
    pub composition: String,
    /// Seed.
    #[arg(long, default_value_t = 0)]
    pub seed: u64,
    /// Volume per atom in cubic Angstrom.
    #[arg(long, default_value_t = 20.0)]
    pub volume_per_atom: f64,
    /// Minimum distance in Angstrom.
    #[arg(long, default_value_t = 1.0)]
    pub min_distance: f64,
    /// Attempts.
    #[arg(long, default_value_t = 64)]
    pub attempts: usize,
    /// Output path; prints a summary without it.
    #[arg(long)]
    pub out: Option<String>,
    /// Output format for `--out`.
    #[arg(long)]
    pub format: Option<StructureFormat>,
}

/// Generation over a space-group range.
#[derive(Debug, Args)]
pub(crate) struct CrystalRangeArgs {
    /// Groups as `1-10` or `225,227`.
    #[arg(long)]
    pub groups: String,
    /// Composition as `Element:count,Element:count`.
    #[arg(long)]
    pub composition: String,
    /// Seed.
    #[arg(long, default_value_t = 0)]
    pub seed: u64,
    /// Volume per atom in cubic Angstrom.
    #[arg(long, default_value_t = 20.0)]
    pub volume_per_atom: f64,
    /// Minimum distance in Angstrom.
    #[arg(long, default_value_t = 1.0)]
    pub min_distance: f64,
    /// Attempts per group.
    #[arg(long, default_value_t = 64)]
    pub attempts: usize,
}

/// Molecular crystal options.
#[derive(Debug, Args)]
pub(crate) struct MoleculeArgs {
    /// Space-group number.
    #[arg(long, default_value_t = 1)]
    pub space_group: u32,
    /// Molecule file (any registered structure format).
    #[arg(long)]
    pub molecule: String,
    /// Independent molecules (Z').
    #[arg(long, default_value_t = 1)]
    pub copies: u32,
    /// Seed.
    #[arg(long, default_value_t = 0)]
    pub seed: u64,
    /// Volume per atom in cubic Angstrom.
    #[arg(long, default_value_t = 20.0)]
    pub volume_per_atom: f64,
    /// Minimum intermolecular distance in Angstrom.
    #[arg(long, default_value_t = 1.0)]
    pub min_distance: f64,
    /// Attempts.
    #[arg(long, default_value_t = 64)]
    pub attempts: usize,
    /// Output path.
    #[arg(long)]
    pub out: Option<String>,
    /// Output format.
    #[arg(long)]
    pub format: Option<StructureFormat>,
}

/// Layer options.
#[derive(Debug, Args)]
pub(crate) struct LayerArgs {
    /// Layer-group number (1-80).
    #[arg(long, default_value_t = 1)]
    pub layer_group: u32,
    /// Composition as `Element:count,Element:count`.
    #[arg(long)]
    pub composition: String,
    /// Area per atom in square Angstrom.
    #[arg(long, default_value_t = 10.0)]
    pub area_per_atom: f64,
    /// Minimum distance in Angstrom.
    #[arg(long, default_value_t = 1.0)]
    pub min_distance: f64,
    /// Seed.
    #[arg(long, default_value_t = 0)]
    pub seed: u64,
    /// Attempts.
    #[arg(long, default_value_t = 64)]
    pub attempts: usize,
    /// Output path.
    #[arg(long)]
    pub out: Option<String>,
    /// Output format.
    #[arg(long)]
    pub format: Option<StructureFormat>,
}

/// Rod options.
#[derive(Debug, Args)]
pub(crate) struct RodArgs {
    /// Rod-group catalogue number (1-14).
    #[arg(long, default_value_t = 1)]
    pub rod_group: u32,
    /// Composition as `Element:count,Element:count`.
    #[arg(long)]
    pub composition: String,
    /// Repeat length in Angstrom.
    #[arg(long, default_value_t = 5.0)]
    pub repeat: f64,
    /// Minimum distance in Angstrom.
    #[arg(long, default_value_t = 1.0)]
    pub min_distance: f64,
    /// Seed.
    #[arg(long, default_value_t = 0)]
    pub seed: u64,
    /// Attempts.
    #[arg(long, default_value_t = 64)]
    pub attempts: usize,
    /// Output path.
    #[arg(long)]
    pub out: Option<String>,
    /// Output format.
    #[arg(long)]
    pub format: Option<StructureFormat>,
}

/// Prototype build options.
#[derive(Debug, Args)]
pub(crate) struct PrototypeArgs {
    /// Prototype label, for example `B1`.
    pub label: String,
    /// Uniform lattice scale.
    #[arg(long, default_value_t = 1.0)]
    pub scale: f64,
    /// Output path.
    #[arg(long)]
    pub out: Option<String>,
    /// Output format.
    #[arg(long)]
    pub format: Option<StructureFormat>,
}

/// Prototype identification options.
#[derive(Debug, Args)]
pub(crate) struct IdentifyArgs {
    /// Input structure.
    pub path: String,
    /// Matcher site tolerance in Angstrom.
    #[arg(long, default_value_t = 0.3)]
    pub stol: f64,
}

/// Substitution options.
#[derive(Debug, Args)]
pub(crate) struct SubstituteArgs {
    /// Input structure.
    pub path: String,
    /// Maximum relative ionic-radius difference.
    #[arg(long, default_value_t = 0.15)]
    pub radius_tolerance: f64,
    /// Maximum candidates.
    #[arg(long, default_value_t = 32)]
    pub limit: usize,
    /// Allow a different oxidation state.
    #[arg(long)]
    pub allow_charge_change: bool,
}

/// Random search options.
#[derive(Debug, Args)]
pub(crate) struct SearchArgs {
    /// Composition as `Element:count,Element:count`.
    #[arg(long)]
    pub composition: String,
    /// Space groups as `1-10` or `225,227`.
    #[arg(long, default_value = "225")]
    pub space_groups: String,
    /// Number of structures.
    #[arg(long, default_value_t = 8)]
    pub structures: usize,
    /// Seed.
    #[arg(long, default_value_t = 0)]
    pub seed: u64,
    /// Skip relaxation.
    #[arg(long)]
    pub no_relax: bool,
    /// Relax the cell too.
    #[arg(long)]
    pub cell: bool,
    /// Maximum optimizer steps.
    #[arg(long, default_value_t = 200)]
    pub steps: usize,
    /// Write the best structure here.
    #[arg(long)]
    pub out: Option<String>,
    /// Format for `--out`.
    #[arg(long)]
    pub format: Option<StructureFormat>,
    /// Write the JSON report here instead of standard output.
    #[arg(long)]
    pub report: Option<String>,
}

/// Evolutionary search options.
#[derive(Debug, Args)]
pub(crate) struct EvolveArgs {
    /// Composition as `Element:count,Element:count`.
    #[arg(long)]
    pub composition: String,
    /// Space groups as `1-10` or `225,227`.
    #[arg(long, default_value = "225")]
    pub space_groups: String,
    /// Population size.
    #[arg(long, default_value_t = 8)]
    pub population: usize,
    /// Generations.
    #[arg(long, default_value_t = 5)]
    pub generations: usize,
    /// Seed.
    #[arg(long, default_value_t = 0)]
    pub seed: u64,
    /// Relax the cell too.
    #[arg(long)]
    pub cell: bool,
    /// Maximum optimizer steps.
    #[arg(long, default_value_t = 150)]
    pub steps: usize,
    /// Write the best structure here.
    #[arg(long)]
    pub out: Option<String>,
    /// Format for `--out`.
    #[arg(long)]
    pub format: Option<StructureFormat>,
}

/// Search workflow options.
#[derive(Debug, Args)]
pub(crate) struct WorkflowArgs {
    /// `random` or `evolutionary`.
    #[arg(long, default_value = "random")]
    pub kind: String,
    /// Seed structure for the workflow's read step.
    #[arg(long)]
    pub seed_path: String,
    /// Search options passed as JSON.
    #[arg(long, default_value = "{}")]
    pub options: String,
    /// Checkpoint path.
    #[arg(long)]
    pub checkpoint: Option<String>,
    /// Output report path.
    #[arg(long, default_value = "search-results.json")]
    pub output: String,
    /// Workflow name.
    #[arg(long, default_value = "search")]
    pub name: String,
}

/// MOF dataset options.
#[derive(Debug, Args)]
pub(crate) struct MofDatasetArgs {
    /// Topologies as a comma list; empty means every topology.
    #[arg(long, default_value = "")]
    pub topologies: String,
    /// Node blocks as a comma list; empty means every node.
    #[arg(long, default_value = "")]
    pub nodes: String,
    /// Linkers as a comma list; empty means every linker.
    #[arg(long, default_value = "")]
    pub linkers: String,
    /// Maximum frameworks.
    #[arg(long, default_value_t = 32)]
    pub limit: usize,
    /// Skip UFF4MOF relaxation.
    #[arg(long)]
    pub no_relax: bool,
    /// Skip descriptors.
    #[arg(long)]
    pub no_descriptors: bool,
    /// Seed.
    #[arg(long, default_value_t = 0)]
    pub seed: u64,
    /// Write the JSON report here.
    #[arg(long)]
    pub out: Option<String>,
}

/// Functionalisation sweep options.
#[derive(Debug, Args)]
pub(crate) struct FunctionalizeArgs {
    /// Input framework.
    pub path: String,
    /// Groups as a comma list: amino, hydroxyl, methyl, fluorine.
    #[arg(long, default_value = "amino")]
    pub groups: String,
    /// Fractions as a comma list.
    #[arg(long, default_value = "0.5,1.0")]
    pub fractions: String,
    /// Seed.
    #[arg(long, default_value_t = 0)]
    pub seed: u64,
    /// Skip descriptors.
    #[arg(long)]
    pub no_descriptors: bool,
    /// Write the JSON report here.
    #[arg(long)]
    pub out: Option<String>,
}
