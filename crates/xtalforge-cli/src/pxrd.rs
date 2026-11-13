//! The `pxrd` command: calculated powder patterns to XY.

use xtalforge::analysis::{
    Caglioti, LineShape, PxrdOptions, PxrdProfileOptions, Radiation, calculate_pxrd,
    calculate_pxrd_profile,
};
use xtalforge::io::write_xy;

use crate::{Cli, CliError, input, output};

/// One `pxrd` invocation.
#[derive(Clone, Copy)]
pub(crate) struct PxrdRequest<'a> {
    /// Input path, or `-`.
    pub(crate) path: &'a str,
    /// Radiation name.
    pub(crate) radiation: &'a str,
    /// Highest 2-theta.
    pub(crate) two_theta_max: f64,
    /// Sampling step.
    pub(crate) step: f64,
    /// Optional Caglioti U, V, W.
    pub(crate) caglioti: Option<&'a [f64]>,
    /// Optional pseudo-Voigt mixing.
    pub(crate) eta: Option<f64>,
    /// Optional line shape name.
    pub(crate) shape: Option<&'a str>,
    /// Output path.
    pub(crate) out: Option<&'a str>,
}

/// Calculates a PXRD pattern and writes an XY file.
///
/// # Errors
///
/// Returns [`CliError`] for a read, calculation or output failure.
pub(crate) fn run(cli: &Cli, request: &PxrdRequest<'_>) -> Result<(), CliError> {
    let PxrdRequest {
        path,
        radiation,
        two_theta_max,
        step,
        caglioti,
        eta,
        shape,
        out,
    } = *request;
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
    let base = PxrdOptions {
        radiation,
        two_theta_range: [5.0, two_theta_max],
        step,
        ..PxrdOptions::default()
    };
    let custom = caglioti.is_some() || eta.is_some() || shape.is_some();
    let pattern = if custom {
        let mut options = PxrdProfileOptions::from_pxrd(&base).map_err(CliError::core)?;
        if let Some(values) = caglioti {
            let [u, v, w] = values else {
                return Err(CliError::input("--caglioti needs exactly three values"));
            };
            options.profile.caglioti = Some(Caglioti {
                u: *u,
                v: *v,
                w: *w,
            });
        }
        if let Some(eta) = eta {
            options.profile.eta = eta;
        }
        if let Some(shape) = shape {
            options.profile.shape = match shape {
                "gaussian" => LineShape::Gaussian,
                "lorentzian" => LineShape::Lorentzian,
                "pseudo-voigt" | "pseudovoigt" => LineShape::PseudoVoigt,
                other => {
                    return Err(CliError::input(format!(
                        "unknown profile shape {other:?}; use gaussian, lorentzian or pseudo-voigt"
                    )));
                }
            };
        }
        calculate_pxrd_profile(&structure, &options).map_err(CliError::core)?
    } else {
        calculate_pxrd(&structure, &base).map_err(CliError::core)?
    };
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
