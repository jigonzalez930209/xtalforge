//! The `build` command family (graphene and nanotubes).

use clap::Subcommand;

use crate::format::StructureFormat;
use crate::{Cli, CliError, output};

#[derive(Debug, Subcommand)]
pub(crate) enum BuildAction {
    /// Build an AA-stacked graphene sheet.
    Graphene {
        /// Repeats along a.
        nx: u32,
        /// Repeats along b.
        ny: u32,
        /// Number of layers.
        #[arg(long, default_value_t = 1)]
        layers: u32,
        /// Interlayer spacing in Angstrom.
        #[arg(long, default_value_t = 3.35)]
        interlayer: f64,
        /// Vacuum above the stack in Angstrom.
        #[arg(long, default_value_t = 15.0)]
        vacuum: f64,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Build a chiral nanotube.
    Nanotube {
        /// Chiral index n.
        n: i32,
        /// Chiral index m.
        m: i32,
        /// Number of translation cells.
        #[arg(long, default_value_t = 1)]
        length: u32,
        /// Vacuum around the tube in Angstrom.
        #[arg(long, default_value_t = 10.0)]
        vacuum: f64,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
}

/// Options for the graphene builder.
pub(crate) struct GrapheneOptions<'a> {
    /// Repeats along a.
    pub nx: u32,
    /// Repeats along b.
    pub ny: u32,
    /// Number of layers.
    pub layers: u32,
    /// Interlayer spacing.
    pub interlayer: f64,
    /// Vacuum above the stack.
    pub vacuum: f64,
    /// Output path.
    pub out: Option<&'a str>,
    /// Output format.
    pub format: Option<StructureFormat>,
}

/// Builds a graphene sheet or stack.
///
/// # Errors
///
/// Returns [`CliError`] for a builder or output failure.
pub(crate) fn graphene(cli: &Cli, options: GrapheneOptions<'_>) -> Result<(), CliError> {
    let structure = xtalforge::build::graphene(
        options.nx,
        options.ny,
        options.layers,
        options.interlayer,
        options.vacuum,
    )
    .map_err(CliError::core)?;
    output::emit_structure(cli, &structure, options.out, options.format)
}

/// Builds a nanotube.
///
/// # Errors
///
/// Returns [`CliError`] for a builder or output failure.
pub(crate) fn nanotube(
    cli: &Cli,
    n: i32,
    m: i32,
    length: u32,
    vacuum: f64,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure = xtalforge::build::nanotube(n, m, length, vacuum).map_err(CliError::core)?;
    output::emit_structure(cli, &structure, out, format)
}

/// Adds the missing hydrogens.
///
/// # Errors
///
/// Returns [`CliError`] for a read or builder failure.
pub(crate) fn hydrogens(
    cli: &Cli,
    path: &str,
    distance: f64,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    use xtalforge::build::hydrogens::{HydrogenOptions, add_hydrogens};
    let structure = crate::input::read_structure_path(path).map_err(CliError::core)?;
    let filled = add_hydrogens(
        &structure,
        &HydrogenOptions {
            distance,
            ..Default::default()
        },
    )
    .map_err(CliError::core)?;
    output::emit_structure(cli, &filled, out, format)
}
