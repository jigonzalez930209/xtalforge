//! The `pxrd` command: calculated powder patterns to XY.

use xtalforge::analysis::{PxrdOptions, Radiation, calculate_pxrd};
use xtalforge::io::write_xy;

use crate::{Cli, CliError, input, output};

/// Calculates a PXRD pattern and writes an XY file.
///
/// # Errors
///
/// Returns [`CliError`] for a read, calculation or output failure.
pub(crate) fn run(
    cli: &Cli,
    path: &str,
    radiation: &str,
    two_theta_max: f64,
    step: f64,
    out: Option<&str>,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let radiation = match radiation {
        "cu-ka" | "cuka" => Radiation::CuKa,
        "cu-ka1" | "cuka1" => Radiation::CuKa1,
        "mo-ka" | "moka" => Radiation::MoKa,
        "co-ka" | "coka" => Radiation::CoKa,
        "fe-ka" | "feka" => Radiation::FeKa,
        "cr-ka" | "crka" => Radiation::CrKa,
        "ag-ka" | "agka" => Radiation::AgKa,
        other => {
            return Err(CliError::input(format!(
                "unknown radiation {other:?}; use cu-ka, cu-ka1, mo-ka, co-ka, fe-ka, cr-ka or ag-ka"
            )));
        }
    };
    let pattern = calculate_pxrd(
        &structure,
        &PxrdOptions {
            radiation,
            two_theta_range: [5.0, two_theta_max],
            step,
            ..PxrdOptions::default()
        },
    )
    .map_err(CliError::core)?;
    let points: Vec<xtalforge::io::PatternPoint> = pattern
        .points
        .iter()
        .map(|point| xtalforge::io::PatternPoint {
            two_theta: point.two_theta,
            intensity: point.intensity,
            uncertainty: None,
        })
        .collect();
    let text = write_xy(&points);
    if let Some(destination) = out.filter(|value| *value != "-") {
        output::write_destination(destination, &text)?;
        if cli.json {
            println!(
                "{}",
                serde_json::json!({"out": destination, "reflections": pattern.reflections.len(), "points": pattern.points.len()})
            );
        } else if !cli.quiet {
            println!(
                "wrote {} points ({} reflections) to {destination}",
                pattern.points.len(),
                pattern.reflections.len()
            );
        }
        return Ok(());
    }
    if cli.json {
        println!(
            "{}",
            serde_json::json!({
                "reflections": pattern.reflections,
                "points": pattern.points.iter().map(|point| [point.two_theta, point.intensity]).collect::<Vec<_>>(),
            })
        );
    } else {
        print!("{text}");
    }
    Ok(())
}
