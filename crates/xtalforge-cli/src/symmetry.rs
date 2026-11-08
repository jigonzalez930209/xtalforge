//! The `symmetry` command family.

use xtalforge::symmetry::{SettingChoice, SymmetryOptions, detect_symmetry};

use crate::{Cli, CliError, input};

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
