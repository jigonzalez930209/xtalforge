//! The `symmetry` command family.

use xtalforge::symmetry::{
    SetMode, SettingChoice, SpaceGroup, SymmetryOptions, detect_symmetry, set_space_group, to_p1,
};

use clap::Subcommand;

use crate::format::StructureFormat;
use crate::{Cli, CliError, input, output};

/// Detects the space group of a structure and prints the dataset.
///
/// # Errors
///
/// Returns [`CliError`] for a read or detection failure.
pub(crate) fn detect(
    cli: &Cli,
    path: &str,
    symprec: f64,
    angle_tolerance: Option<f64>,
    hall_number: Option<u32>,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let options = SymmetryOptions {
        symprec,
        angle_tolerance,
        setting: hall_number.map_or(SettingChoice::Spglib, SettingChoice::HallNumber),
    };
    let dataset = detect_symmetry(&structure, options).map_err(CliError::core)?;
    if cli.json {
        let payload = serde_json::to_string(&dataset).map_err(|error| {
            CliError::new(
                "json_error",
                format!("could not serialize symmetry: {error}"),
            )
        })?;
        println!("{payload}");
    } else {
        println!(
            "{} (No. {}, Hall {}, Pearson {})",
            dataset.hm_symbol, dataset.number, dataset.hall_number, dataset.pearson_symbol
        );
        println!(
            "operations={} atoms={} symprec={} A",
            dataset.operations.len(),
            dataset.n_atoms(),
            symprec
        );
        for (index, (wyckoff, site)) in dataset
            .wyckoffs
            .iter()
            .zip(dataset.site_symmetry_symbols.iter())
            .enumerate()
        {
            println!("  atom {index}: Wyckoff {wyckoff}, site symmetry {site}");
        }
    }
    Ok(())
}

/// Expands a structure to P1.
///
/// # Errors
///
/// Returns [`CliError`] for a read or expansion failure.
pub(crate) fn p1(
    cli: &Cli,
    path: &str,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let expanded = to_p1(&structure).map_err(CliError::core)?;
    output::emit_structure(cli, &expanded, out, format)
}

/// Reduces a structure to its detected asymmetric unit.
///
/// # Errors
///
/// Returns [`CliError`] for a read or detection failure.
pub(crate) fn asymmetrize(
    cli: &Cli,
    path: &str,
    symprec: f64,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let reduced = xtalforge::symmetry::asymmetrize(
        &structure,
        SymmetryOptions {
            symprec,
            ..SymmetryOptions::default()
        },
    )
    .map_err(CliError::core)?;
    output::emit_structure(cli, &reduced, out, format)
}

/// Applies a space group in reinterpret or impose mode.
///
/// # Errors
///
/// Returns [`CliError`] for a read, parse or averaging failure.
pub(crate) fn set(
    cli: &Cli,
    path: &str,
    symbol: &str,
    mode: &str,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let group = SpaceGroup::parse(symbol).map_err(CliError::core)?;
    let mode = match mode {
        "reinterpret" => SetMode::Reinterpret,
        "impose" => SetMode::Impose,
        other => {
            return Err(CliError::input(format!(
                "unknown set mode {other:?}; use reinterpret or impose"
            )));
        }
    };
    let report = set_space_group(&structure, &group, mode).map_err(CliError::core)?;
    if cli.json && out.is_none() {
        println!(
            "{}",
            serde_json::to_string(&report.structure).map_err(|error| {
                CliError::new(
                    "json_error",
                    format!("could not serialize structure: {error}"),
                )
            })?
        );
        return Ok(());
    }
    output::emit_structure(cli, &report.structure, out, format)
}

/// Lists the translationengleiche subgroups of a Hall setting.
///
/// # Errors
///
/// Returns [`CliError`] when enumeration fails.
pub(crate) fn subgroups(cli: &Cli, hall_number: u32, epsilon: f64) -> Result<(), CliError> {
    let groups = xtalforge::symmetry::translationengleiche_subgroups(hall_number, epsilon)
        .map_err(CliError::core)?;
    if cli.json {
        println!(
            "{}",
            serde_json::to_string(&groups).map_err(|error| {
                CliError::new(
                    "json_error",
                    format!("could not serialize subgroups: {error}"),
                )
            })?
        );
    } else {
        for group in groups {
            println!(
                "index {:>3} depth {} {:?} operations {}{}",
                group.index,
                group.depth,
                group.kind,
                group.operations.len(),
                group
                    .subgroup_number
                    .map(|number| format!(" -> No. {number}"))
                    .unwrap_or_default()
            );
        }
    }
    Ok(())
}

/// Lowers a structure to a subgroup number.
///
/// # Errors
///
/// Returns [`CliError`] for a read or lowering failure.
pub(crate) fn lower(
    cli: &Cli,
    path: &str,
    subgroup: u32,
    tolerance: f64,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let lowered = xtalforge::symmetry::lower_to_subgroup(&structure, subgroup, tolerance)
        .map_err(CliError::core)?;
    output::emit_structure(cli, &lowered, out, format)
}

#[derive(Debug, Subcommand)]
pub(crate) enum SymmetryAction {
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
