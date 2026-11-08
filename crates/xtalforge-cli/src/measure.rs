//! Geometry commands: `distance`, `angle` and `torsion`.

use serde_json::json;

use crate::{Cli, CliError, input};

/// Minimum-image distance between two P1 atoms.
///
/// # Errors
///
/// Returns [`CliError`] for a read error or an out-of-range atom index.
pub(crate) fn distance(cli: &Cli, path: &str, i: usize, j: usize) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let value = xtalforge::core::distance(&structure, i, j).map_err(CliError::core)?;
    if cli.json {
        println!("{}", json!({ "i": i, "j": j, "distance": value }));
    } else {
        println!("{value:.6}");
    }
    Ok(())
}

/// Angle `i-j-k` in degrees.
///
/// # Errors
///
/// Returns [`CliError`] for a read error or a degenerate angle.
pub(crate) fn angle(cli: &Cli, path: &str, i: usize, j: usize, k: usize) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let value = xtalforge::core::angle(&structure, i, j, k).map_err(CliError::core)?;
    if cli.json {
        println!("{}", json!({ "i": i, "j": j, "k": k, "angle": value }));
    } else {
        println!("{value:.6}");
    }
    Ok(())
}

/// Torsion `i-j-k-l` in degrees.
///
/// # Errors
///
/// Returns [`CliError`] for a read error or a collinear segment.
pub(crate) fn torsion(
    cli: &Cli,
    path: &str,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let value = xtalforge::core::torsion(&structure, i, j, k, l).map_err(CliError::core)?;
    if cli.json {
        println!(
            "{}",
            json!({ "i": i, "j": j, "k": k, "l": l, "torsion": value })
        );
    } else {
        println!("{value:.6}");
    }
    Ok(())
}
