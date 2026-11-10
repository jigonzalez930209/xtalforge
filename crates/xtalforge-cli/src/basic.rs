//! Foundation commands that predate the Phase 1 domain surface.

use crate::{Cli, CliError, input};

/// Prints a structure summary from JSON, CIF, POSCAR or XYZ.
///
/// # Errors
///
/// Returns [`CliError`] for an unreadable input, an unknown format or a
/// structure without a summary.
pub(crate) fn info(cli: &Cli, path: &str) -> Result<(), CliError> {
    let text = input::read_source(path).map_err(CliError::core)?;
    let is_json =
        path.to_ascii_lowercase().ends_with(".json") || text.trim_start().starts_with('{');
    let summary = if is_json {
        xtalforge::core::structure_info(&text)
    } else {
        input::structure_from_text(&text, path).and_then(|structure| {
            if !cli.quiet {
                for key in ["symmetry_warning", "coincident_atoms_warning"] {
                    if let Some(message) = structure
                        .metadata()
                        .get(key)
                        .and_then(serde_json::Value::as_str)
                    {
                        eprintln!("warning: {message}");
                    }
                }
            }
            structure.summary()
        })
    }
    .map_err(CliError::core)?;
    if cli.json {
        let payload = serde_json::to_string(&summary)
            .map_err(|error| CliError::new("json_error", error.to_string()))?;
        println!("{payload}");
    } else {
        println!(
            "{}: {} [{}] sites={} volume={:.3} A^3 density={:.4} g/cm^3 charge={:.3} p1={}",
            summary.name,
            summary.formula,
            summary.hall,
            summary.n_sites,
            summary.volume,
            summary.density,
            summary.total_charge,
            summary.is_p1
        );
    }
    Ok(())
}

/// Builds a lattice from six cell parameters and reports it.
///
/// # Errors
///
/// Returns [`CliError`] when the argument count is not six or the parameters
/// do not describe a valid cell.
pub(crate) fn cell(cli: &Cli, params: &[f64]) -> Result<(), CliError> {
    let [a, b, c, alpha, beta, gamma] = match params {
        [a, b, c, alpha, beta, gamma] => [*a, *b, *c, *alpha, *beta, *gamma],
        _ => {
            return Err(CliError::input(
                "cell needs six parameters: a b c alpha beta gamma",
            ));
        }
    };
    let lattice =
        xtalforge::core::lattice_info(a, b, c, alpha, beta, gamma).map_err(CliError::core)?;
    if cli.json {
        println!(
            "{}",
            serde_json::json!({
                "parameters": lattice.parameters(),
                "volume": lattice.volume(),
                "matrix": lattice.matrix(),
                "right_handed": lattice.is_right_handed(),
            })
        );
    } else {
        let [a, b, c, alpha, beta, gamma] = lattice.parameters();
        println!(
            "a={a:.6} b={b:.6} c={c:.6} alpha={alpha:.6} beta={beta:.6} gamma={gamma:.6} volume={:.6}",
            lattice.volume()
        );
        for row in lattice.matrix() {
            let [x, y, z] = row;
            println!("{x:12.6} {y:12.6} {z:12.6}");
        }
    }
    Ok(())
}

/// Computes the volume of a cell from six parameters.
///
/// # Errors
///
/// Returns [`CliError`] when the parameters do not describe a valid cell.
pub(crate) fn cell_volume(
    cli: &Cli,
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
) -> Result<(), CliError> {
    let volume =
        xtalforge::core::cell_volume(a, b, c, alpha, beta, gamma).map_err(CliError::core)?;
    if cli.json {
        println!("{}", serde_json::json!({ "volume": volume }));
    } else {
        println!("{volume}");
    }
    Ok(())
}
