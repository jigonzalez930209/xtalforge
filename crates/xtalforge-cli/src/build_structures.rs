//! The `build` command family (graphene, nanotubes and fullerenes).

use clap::Subcommand;
use serde_json::json;

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
    /// Build a C20 dodecahedron.
    C20 {
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Build a C60 fullerene.
    C60 {
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Build two graphene layers with a relative twist.
    Twist {
        /// Repeats along a.
        n: u32,
        /// Repeats along b.
        m: u32,
        /// Requested twist angle in degrees.
        theta: f64,
        /// Number of layers.
        #[arg(long, default_value_t = 2)]
        layers: u32,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Build a honeycomb nanoribbon with vacuum.
    Ribbon {
        /// Width in repeat units.
        width: u32,
        /// Length in repeat units.
        length: u32,
        /// Vacuum on both sides in Angstrom.
        #[arg(long, default_value_t = 12.0)]
        vacuum: f64,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Build a gyroid ZTC carbon model.
    Ztc {
        /// Samples per cube axis.
        #[arg(long, default_value_t = 16)]
        samples: usize,
        /// Cubic cell side in Angstrom.
        #[arg(long, default_value_t = 12.0)]
        cell: f64,
        /// Level-set value.
        #[arg(long, default_value_t = 0.0)]
        isovalue: f64,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Build a Schwarz P carbon model.
    Schwarzite {
        /// Samples per cube axis.
        #[arg(long, default_value_t = 16)]
        samples: usize,
        /// Cubic cell side in Angstrom.
        #[arg(long, default_value_t = 12.0)]
        cell: f64,
        /// Level-set value.
        #[arg(long, default_value_t = 0.0)]
        isovalue: f64,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Build a disordered carbon from a recipe (dia or pcu net).
    Carbon {
        /// Net name (`dia` or `pcu`).
        #[arg(long, default_value = "dia")]
        net: String,
        /// Repeats along a, b and c.
        #[arg(long, default_value_t = 1)]
        repeat_x: usize,
        #[arg(long, default_value_t = 1)]
        repeat_y: usize,
        #[arg(long, default_value_t = 1)]
        repeat_z: usize,
        /// Carbon density in g/cm3.
        #[arg(long, default_value_t = 0.42)]
        density: f64,
        /// Strut radius over the mean net edge.
        #[arg(long, default_value_t = 0.25)]
        radius_ratio: f64,
        /// Fraction of the closed sheet kept.
        #[arg(long, default_value_t = 0.38)]
        coverage: f64,
        /// Deterministic seed.
        #[arg(long, default_value_t = 0)]
        seed: u64,
        /// Hydrogen per carbon.
        #[arg(long, default_value_t = 0.07)]
        hydrogen: f64,
        /// Fluorine per carbon.
        #[arg(long, default_value_t = 0.29)]
        fluorine: f64,
        /// Oxygen per carbon.
        #[arg(long, default_value_t = 0.044)]
        oxygen: f64,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Build a polymer chain from a repeat-unit SMILES.
    Polymer {
        /// Repeat-unit SMILES (`C=C` or `C=Cc1ccccc1`).
        smiles: String,
        /// Number of monomers.
        monomers: u32,
        /// Target density in gram per cubic centimetre.
        #[arg(long, default_value_t = 1.0)]
        density: f64,
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

/// Builds a C20 dodecahedron.
///
/// # Errors
///
/// Returns [`CliError`] for a builder or output failure.
pub(crate) fn c20(
    cli: &Cli,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure = xtalforge::build::face_spiral(20).map_err(CliError::core)?;
    output::emit_structure(cli, &structure, out, format)
}

/// Builds two twisted graphene layers.
///
/// # Errors
///
/// Returns [`CliError`] for a builder or output failure.
pub(crate) fn twist(
    cli: &Cli,
    n: u32,
    m: u32,
    theta: f64,
    layers: u32,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure =
        xtalforge::build::commensurate_twist(n, m, theta, layers).map_err(CliError::core)?;
    output::emit_structure(cli, &structure, out, format)
}

/// Builds a honeycomb nanoribbon.
///
/// # Errors
///
/// Returns [`CliError`] for a builder or output failure.
pub(crate) fn ribbon(
    cli: &Cli,
    width: u32,
    length: u32,
    vacuum: f64,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure = xtalforge::build::nanoribbon(width, length, vacuum).map_err(CliError::core)?;
    output::emit_structure(cli, &structure, out, format)
}

fn ztc_options(samples: usize, cell: f64, isovalue: f64) -> xtalforge::build::ZtcOptions {
    xtalforge::build::ZtcOptions {
        samples,
        cell,
        isovalue,
    }
}

/// Builds a gyroid ZTC carbon model.
///
/// # Errors
///
/// Returns [`CliError`] for a builder or output failure.
pub(crate) fn ztc(
    cli: &Cli,
    samples: usize,
    cell: f64,
    isovalue: f64,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure = xtalforge::build::build_ztc(&ztc_options(samples, cell, isovalue))
        .map_err(CliError::core)?;
    output::emit_structure(cli, &structure, out, format)
}

/// Builds a Schwarz P carbon model.
///
/// # Errors
///
/// Returns [`CliError`] for a builder or output failure.
pub(crate) fn schwarzite(
    cli: &Cli,
    samples: usize,
    cell: f64,
    isovalue: f64,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure = xtalforge::build::build_schwarzite(&ztc_options(samples, cell, isovalue))
        .map_err(CliError::core)?;
    output::emit_structure(cli, &structure, out, format)
}

/// Builds a polymer chain and records its steps.
///
/// # Errors
///
/// Returns [`CliError`] for a builder or output failure.
pub(crate) fn polymer(
    cli: &Cli,
    smiles: &str,
    monomers: u32,
    density: f64,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let build = xtalforge::build::build_polymer(
        smiles,
        monomers,
        &xtalforge::build::PolymerOptions { density },
    )
    .map_err(CliError::core)?;
    if out.is_none() && cli.json {
        println!(
            "{}",
            serde_json::json!({
                "structure": build.structure,
                "steps": build.steps,
            })
        );
        return Ok(());
    }
    output::emit_structure(cli, &build.structure, out, format)?;
    if !cli.quiet && !cli.json {
        for step in &build.steps {
            println!("step: {step}");
        }
    }
    Ok(())
}

/// Builds a C60 fullerene.
///
/// # Errors
///
/// Returns [`CliError`] for a builder or output failure.
pub(crate) fn c60(
    cli: &Cli,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure = xtalforge::build::fullerene(60).map_err(CliError::core)?;
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

/// Builds a disordered carbon from a recipe.
///
/// # Errors
///
/// Returns [`CliError`] for a builder or output failure.
#[allow(clippy::too_many_arguments)] // one flag per recipe field
pub(crate) fn carbon(
    cli: &Cli,
    net: &str,
    repeat: [usize; 3],
    density: f64,
    radius_ratio: f64,
    coverage: f64,
    seed: u64,
    hydrogen: f64,
    fluorine: f64,
    oxygen: f64,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let recipe = xtalforge::build::CarbonRecipe {
        net: net.to_owned(),
        repeat,
        density,
        radius_ratio,
        coverage,
        seed,
        hydrogen,
        fluorine,
        oxygen,
        ..xtalforge::build::CarbonRecipe::default()
    };
    let build = xtalforge::build::build_carbon(&recipe).map_err(CliError::core)?;
    if cli.json {
        println!("{}", json!(build));
        return Ok(());
    }
    println!(
        "{} atoms, {:.3} g/cm3, mean C-C {:.3} A, {} piece(s), periodicity {}",
        build.structure.sites().len(),
        build.density,
        build.mean_bond,
        build.pieces,
        build.periodicity
    );
    for note in &build.notes {
        println!("note: {note}");
    }
    output::emit_structure(cli, &build.structure, out, format)
}
