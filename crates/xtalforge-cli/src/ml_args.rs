//! Arguments for the Phase 21 `ml` command.

use clap::{Args, Subcommand};

use crate::format::StructureFormat;

/// Machine-learning representations, datasets and inference.
#[derive(Debug, Subcommand)]
pub(crate) enum MlAction {
    /// Build a periodic crystal graph.
    Graph(GraphArgs),
    /// Calculate the exact SOAP power spectrum.
    Soap(SoapArgs),
    /// Calculate ACSF features.
    Acsf(AcsfArgs),
    /// Calculate an MBTR descriptor.
    Mbtr(MbtrArgs),
    /// Calculate MOF-specific features.
    MofFeatures(MofFeaturesArgs),
    /// Split a dataset without leakage.
    Split(SplitArgs),
    /// Import an extxyz, JSON or OPTIMADE file into a dataset.
    Import(ImportArgs),
    /// Predict a property for every dataset record.
    Predict(PredictArgs),
    /// Select a batch from a feature table.
    Select(SelectArgs),
    /// Emit an active-learning workflow document.
    ActiveWorkflow(ActiveWorkflowArgs),
    /// Write a JSON Lines dataset as Parquet.
    ParquetWrite(ParquetWriteArgs),
    /// Read a Parquet dataset as JSON Lines.
    ParquetRead(ParquetReadArgs),
    /// Run a single-input ONNX model.
    OnnxPredict(OnnxPredictArgs),
    /// Compute the ONNX Cartesian energy of a structure.
    OnnxEnergy(OnnxEnergyArgs),
}

/// Graph options.
#[derive(Debug, Args)]
pub(crate) struct GraphArgs {
    /// Input structure.
    pub path: String,
    /// Cutoff in Angstrom.
    #[arg(long, default_value_t = 5.0)]
    pub cutoff: f64,
    /// Keep at most this many outgoing edges per node.
    #[arg(long)]
    pub k: Option<usize>,
    /// Write the JSON graph here.
    #[arg(long)]
    pub out: Option<String>,
}

/// SOAP options.
#[derive(Debug, Args)]
pub(crate) struct SoapArgs {
    /// Input structure.
    pub path: String,
    /// Cutoff in Angstrom.
    #[arg(long, default_value_t = 5.0)]
    pub cutoff: f64,
    /// Radial basis functions.
    #[arg(long, default_value_t = 4)]
    pub n_max: usize,
    /// Maximum spherical harmonic degree.
    #[arg(long, default_value_t = 4)]
    pub l_max: usize,
    /// Gaussian width.
    #[arg(long, default_value_t = 0.5)]
    pub sigma: f64,
    /// Average over centres.
    #[arg(long)]
    pub average: bool,
    /// Write the JSON descriptor here.
    #[arg(long)]
    pub out: Option<String>,
}

/// ACSF options.
#[derive(Debug, Args)]
pub(crate) struct AcsfArgs {
    /// Input structure.
    pub path: String,
    /// Cutoff in Angstrom.
    #[arg(long, default_value_t = 5.0)]
    pub cutoff: f64,
    /// Average over centres.
    #[arg(long)]
    pub average: bool,
    /// Write the JSON descriptor here.
    #[arg(long)]
    pub out: Option<String>,
}

/// MBTR options.
#[derive(Debug, Args)]
pub(crate) struct MbtrArgs {
    /// Input structure.
    pub path: String,
    /// Cutoff in Angstrom.
    #[arg(long, default_value_t = 5.0)]
    pub cutoff: f64,
    /// Grid points.
    #[arg(long, default_value_t = 64)]
    pub grid: usize,
    /// Gaussian width.
    #[arg(long, default_value_t = 0.1)]
    pub sigma: f64,
    /// Skip the angle block.
    #[arg(long)]
    pub no_angles: bool,
    /// Skip the distance block.
    #[arg(long)]
    pub no_distances: bool,
    /// Write the JSON descriptor here.
    #[arg(long)]
    pub out: Option<String>,
}

/// MOF feature options.
#[derive(Debug, Args)]
pub(crate) struct MofFeaturesArgs {
    /// Input framework.
    pub path: String,
    /// Longest bond path for RACs.
    #[arg(long, default_value_t = 3)]
    pub max_path: usize,
    /// Skip the pore block.
    #[arg(long)]
    pub no_pores: bool,
    /// Write the JSON feature vector here.
    #[arg(long)]
    pub out: Option<String>,
}

/// Split options.
#[derive(Debug, Args)]
pub(crate) struct SplitArgs {
    /// Dataset JSON Lines file.
    pub path: String,
    /// `random`, `composition`, `cluster` or `prototype`.
    #[arg(long, default_value = "random")]
    pub strategy: String,
    /// Training fraction.
    #[arg(long, default_value_t = 0.6)]
    pub train: f64,
    /// Validation fraction.
    #[arg(long, default_value_t = 0.2)]
    pub validation: f64,
    /// Test fraction.
    #[arg(long, default_value_t = 0.2)]
    pub test: f64,
    /// Seed.
    #[arg(long, default_value_t = 0)]
    pub seed: u64,
    /// Write the assignment JSON here.
    #[arg(long)]
    pub out: Option<String>,
}

/// Import options.
#[derive(Debug, Args)]
pub(crate) struct ImportArgs {
    /// Input file.
    pub path: String,
    /// `extxyz`, `json` or `optimade`.
    #[arg(long, default_value = "extxyz")]
    pub format: String,
    /// Dataset name.
    #[arg(long, default_value = "imported")]
    pub name: String,
    /// Write the JSON Lines dataset here.
    #[arg(long)]
    pub out: Option<String>,
}

/// Prediction options.
#[derive(Debug, Args)]
pub(crate) struct PredictArgs {
    /// Dataset JSON Lines file.
    pub path: String,
    /// Property model JSON file.
    #[arg(long)]
    pub model: String,
    /// Feature kind: `soap`, `acsf` or `mbtr`.
    #[arg(long, default_value = "soap")]
    pub features: String,
    /// Property name written to the dataset.
    #[arg(long, default_value = "predicted")]
    pub target: String,
    /// Write the updated dataset here.
    #[arg(long)]
    pub out: Option<String>,
}

/// Selection options.
#[derive(Debug, Args)]
pub(crate) struct SelectArgs {
    /// Feature table JSON file (`FeatureTable` or `{features, uncertainties}`).
    pub path: String,
    /// Batch size.
    #[arg(long, default_value_t = 8)]
    pub batch: usize,
    /// `uncertainty`, `farthest_point` or `diversity`.
    #[arg(long, default_value = "uncertainty")]
    pub strategy: String,
    /// Minimum distance for the diversity strategy.
    #[arg(long, default_value_t = 0.0)]
    pub diversity_radius: f64,
    /// Write the selected indices here.
    #[arg(long)]
    pub out: Option<String>,
}

/// Active-learning workflow options.
#[derive(Debug, Args)]
pub(crate) struct ActiveWorkflowArgs {
    /// Workflow name.
    #[arg(long, default_value = "active-learning")]
    pub name: String,
    /// Dataset path read by the workflow.
    #[arg(long, default_value = "dataset.jsonl")]
    pub dataset: String,
    /// Model path.
    #[arg(long, default_value = "model.json")]
    pub model: String,
    /// Batch size.
    #[arg(long, default_value_t = 8)]
    pub batch: usize,
    /// Rounds.
    #[arg(long, default_value_t = 4)]
    pub rounds: usize,
    /// Checkpoint path.
    #[arg(long)]
    pub checkpoint: Option<String>,
    /// Write the workflow JSON here.
    #[arg(long)]
    pub out: Option<String>,
}

/// Keeps the structure-format import used by future options.
pub(crate) type _StructureFormat = StructureFormat;

/// Parquet write options.
#[derive(Debug, Args)]
pub(crate) struct ParquetWriteArgs {
    /// Dataset JSON Lines file.
    pub path: String,
    /// Output Parquet file.
    #[arg(long)]
    pub out: String,
}

/// Parquet read options.
#[derive(Debug, Args)]
pub(crate) struct ParquetReadArgs {
    /// Input Parquet file.
    pub path: String,
    /// Output JSON Lines file.
    #[arg(long)]
    pub out: Option<String>,
}

/// ONNX prediction options.
#[derive(Debug, Args)]
pub(crate) struct OnnxPredictArgs {
    /// ONNX model file.
    pub model: String,
    /// Feature vector as JSON.
    #[arg(long)]
    pub features: String,
    /// Input width.
    #[arg(long)]
    pub n_features: usize,
}

/// ONNX energy options.
#[derive(Debug, Args)]
pub(crate) struct OnnxEnergyArgs {
    /// Input structure.
    pub path: String,
    /// ONNX model file.
    #[arg(long)]
    pub model: String,
    /// Input width.
    #[arg(long)]
    pub n_features: usize,
    /// Finite-difference step.
    #[arg(long, default_value_t = 0.001)]
    pub displacement: f64,
}
