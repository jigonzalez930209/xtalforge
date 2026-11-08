//! The `symmetry` command family.

use xtalforge::symmetry::{
    SetMode, SettingChoice, SpaceGroup, SymmetryOptions, asymmetrize as reduce, detect_symmetry,
    set_space_group, to_p1,
};

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
    let reduced = reduce(
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
