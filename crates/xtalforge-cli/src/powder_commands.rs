//! The `powder` command family (Phase 7).

use clap::Subcommand;
use serde_json::json;
use xtalforge::analysis::{
    AutomaticOptions, IndexingOptions, LebailOptions, PawleyOptions, PeakOptions, RefineOptions,
    RietveldOptions, extinction_symbol, find_peaks, free_cell_parameters, index_cubic,
    index_pattern, le_bail, pawley_refine, refine_phase, rietx_inventory, run_automatic,
};

use crate::{Cli, CliError, input};

/// Subcommands of `xtalforge powder`.
#[derive(Debug, Subcommand)]
pub(crate) enum PowderAction {
    /// Find peaks and index the pattern.
    Index {
        /// Measured XY pattern.
        pattern: String,
        /// Wavelength in Angstrom.
        #[arg(long, default_value_t = 1.54184)]
        wavelength: f64,
        /// Largest hkl index.
        #[arg(long, default_value_t = 6)]
        max_index: i32,
        /// Position tolerance in degrees.
        #[arg(long, default_value_t = 0.1)]
        tolerance: f64,
        /// Run the general metric-family search and rank its candidates.
        #[arg(long)]
        general: bool,
        /// Subtract a SNIP background before finding peaks.
        #[arg(long)]
        background: bool,
        /// Fit every detected line in its own window.
        #[arg(long)]
        refine_peaks: bool,
    },
    /// Le Bail intensity extraction.
    Lebail {
        /// Structure path.
        structure: String,
        /// Measured XY pattern.
        pattern: String,
    },
    /// Pawley intensity extraction.
    Pawley {
        /// Structure path.
        structure: String,
        /// Measured XY pattern.
        pattern: String,
        /// Profile shape: gaussian, lorentzian or pseudo-voigt.
        #[arg(long, default_value = "pseudo-voigt")]
        profile: String,
    },
    /// Rietveld scale/background or structural refinement.
    Rietveld {
        /// Structure path.
        structure: String,
        /// Measured XY pattern.
        pattern: String,
        /// Refine positions, ADPs and occupancies.
        #[arg(long)]
        refine: bool,
    },
    /// Automatic peaks, indexing, Le Bail and Rietveld run.
    Automatic {
        /// Structure path.
        structure: String,
        /// Measured XY pattern.
        pattern: String,
    },
    /// Free cell parameters and extinction symbol.
    Cells {
        /// Space-group number (1 to 230).
        number: u32,
    },
    /// Print the rietx bridge inventory.
    Rietx,
}

/// Runs one `powder` action.
///
/// # Errors
///
/// Returns [`CliError`] for read, analysis or output failures.
pub(crate) fn run(cli: &Cli, action: &PowderAction) -> Result<(), CliError> {
    match action {
        PowderAction::Index {
            pattern,
            wavelength,
            max_index,
            tolerance,
            general,
            background,
            refine_peaks,
        } => index(
            cli,
            pattern,
            *wavelength,
            *max_index,
            *tolerance,
            *general,
            *background,
            *refine_peaks,
        ),
        PowderAction::Lebail { structure, pattern } => lebail(cli, structure, pattern, false, None),
        PowderAction::Pawley {
            structure,
            pattern,
            profile,
        } => lebail(cli, structure, pattern, true, Some(profile)),
        PowderAction::Rietveld {
            structure,
            pattern,
            refine,
        } => rietveld(cli, structure, pattern, *refine),
        PowderAction::Automatic { structure, pattern } => automatic(cli, structure, pattern),
        PowderAction::Cells { number } => cells(cli, *number),
        PowderAction::Rietx => {
            let bridge = rietx_inventory();
            if cli.json {
                println!(
                    "{}",
                    json!({
                        "name": bridge.name,
                        "capabilities": bridge.capabilities,
                        "license": bridge.license,
                        "source": bridge.source,
                    })
                );
            } else {
                println!("{} ({})", bridge.name, bridge.license);
                println!("source: {}", bridge.source);
                for capability in bridge.capabilities {
                    println!("  {capability}");
                }
            }
            Ok(())
        }
    }
}

fn pattern(path: &str) -> Result<Vec<xtalforge::io::PatternPoint>, CliError> {
    let text = input::read_source(path).map_err(CliError::core)?;
    xtalforge::io::read_xy(&text).map_err(CliError::core)
}

#[allow(clippy::too_many_arguments)] // one flag per documented search switch
fn index(
    cli: &Cli,
    path: &str,
    wavelength: f64,
    max_index: i32,
    tolerance: f64,
    general: bool,
    background: bool,
    refine_peaks: bool,
) -> Result<(), CliError> {
    let points = pattern(path)?;
    let peak_options = PeakOptions {
        background,
        refine: refine_peaks,
        ..PeakOptions::default()
    };
    let peaks = find_peaks(&points, &peak_options).map_err(CliError::analysis)?;
    let cubic =
        index_cubic(&peaks, wavelength, max_index, tolerance).map_err(CliError::analysis)?;
    let candidates = if general {
        Some(
            index_pattern(
                &peaks,
                wavelength,
                &IndexingOptions {
                    max_index,
                    tolerance,
                    ..IndexingOptions::default()
                },
            )
            .map_err(CliError::analysis)?,
        )
    } else {
        None
    };
    if cli.json {
        println!(
            "{}",
            json!({"peaks": peaks, "index": cubic, "candidates": candidates})
        );
        return Ok(());
    }
    println!("{} peaks", peaks.len());
    match cubic {
        Some(index) => {
            println!(
                "a = {:.6} A, completeness {:.2}, mean error {:.4} deg",
                index.a, index.completeness, index.mean_error
            );
        }
        None => println!("no cubic candidate"),
    }
    if let Some(candidates) = candidates {
        for (rank, candidate) in candidates.candidates.iter().take(5).enumerate() {
            println!(
                "{}. {} a={:.4} b={:.4} c={:.4} beta={:.2} M20={:.1} F_N={:.1} coverage={:.2}",
                rank + 1,
                candidate.system.name(),
                candidate.cell.a,
                candidate.cell.b,
                candidate.cell.c,
                candidate.cell.beta,
                candidate.m20,
                candidate.f_n,
                candidate.completeness
            );
        }
    }
    Ok(())
}

fn lebail(
    cli: &Cli,
    structure_path: &str,
    pattern_path: &str,
    use_pawley: bool,
    profile: Option<&str>,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(structure_path).map_err(CliError::core)?;
    let points = pattern(pattern_path)?;
    if use_pawley {
        let kind = match profile
            .unwrap_or("pseudo-voigt")
            .to_ascii_lowercase()
            .as_str()
        {
            "gaussian" => xtalforge::analysis::ProfileKind::Gaussian,
            "lorentzian" => xtalforge::analysis::ProfileKind::Lorentzian,
            _ => xtalforge::analysis::ProfileKind::PseudoVoigt,
        };
        let result = pawley_refine(
            &structure,
            &points,
            &PawleyOptions {
                kind,
                ..PawleyOptions::default()
            },
        )
        .map_err(CliError::analysis)?;
        if cli.json {
            println!("{}", json!(result));
        } else {
            println!(
                "{} reflections, Rp = {:.4}, Rwp = {:.4}, Rexp = {:.4}, chi2 = {:.4}",
                result.reflections.len(),
                result.rp,
                result.rwp,
                result.rexp,
                result.chi_squared
            );
            println!(
                "a = {:.6} (esd {:.6}), zero = {:.5} (esd {:.5})",
                result.cell[0], result.cell_esds[0], result.zero_shift, result.zero_shift_esd
            );
        }
        return Ok(());
    }
    let result =
        le_bail(&structure, &points, &LebailOptions::default()).map_err(CliError::analysis)?;
    if cli.json {
        println!(
            "{}",
            json!({
                "intensities": result.intensities,
                "positions": result.positions,
                "r_factor": result.r_factor,
            })
        );
    } else {
        println!(
            "{} reflections, R = {:.4}",
            result.intensities.len(),
            result.r_factor
        );
    }
    Ok(())
}

fn rietveld(
    cli: &Cli,
    structure_path: &str,
    pattern_path: &str,
    refine: bool,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(structure_path).map_err(CliError::core)?;
    let points = pattern(pattern_path)?;
    let options = RietveldOptions::default();
    if refine {
        let result = refine_phase(
            &structure,
            &points,
            &options,
            &RefineOptions {
                max_iterations: 20,
                refine_positions: true,
                ..RefineOptions::default()
            },
        )
        .map_err(CliError::analysis)?;
        if cli.json {
            println!(
                "{}",
                json!({
                    "rwp": result.rwp,
                    "converged": result.converged,
                    "iterations": result.iterations,
                    "parameters": result.parameters,
                })
            );
        } else {
            println!(
                "Rwp {:.4}, converged {}, {} iterations",
                result.rwp, result.converged, result.iterations
            );
            for (name, value) in &result.parameters {
                println!("  {name} = {value:.6}");
            }
        }
    } else {
        let phase = xtalforge::analysis::fit_phase(&structure, &points, &options)
            .map_err(CliError::analysis)?;
        if cli.json {
            println!(
                "{}",
                json!({"scale": phase.scale, "background": phase.background, "rwp": phase.rwp})
            );
        } else {
            println!(
                "scale {:.6}, background {:.6}, Rwp {:.4}",
                phase.scale, phase.background, phase.rwp
            );
        }
    }
    Ok(())
}

fn automatic(cli: &Cli, structure_path: &str, pattern_path: &str) -> Result<(), CliError> {
    let structure = input::read_structure_path(structure_path).map_err(CliError::core)?;
    let points = pattern(pattern_path)?;
    let result = run_automatic(&structure, &points, &AutomaticOptions::default())
        .map_err(CliError::analysis)?;
    if cli.json {
        println!(
            "{}",
            json!({
                "peaks": result.peaks.len(),
                "index": result.index,
                "pawley_rwp": result.pawley.as_ref().map(|pawley| pawley.rwp),
                "lebail": result.lebail.as_ref().map(|lebail| lebail.r_factor),
                "rwp": result.phase.as_ref().map(|phase| phase.rwp),
                "steps": result.steps,
            })
        );
    } else {
        for step in &result.steps {
            println!("{step}");
        }
    }
    Ok(())
}

fn cells(cli: &Cli, number: u32) -> Result<(), CliError> {
    let parameters = free_cell_parameters(number).map_err(CliError::analysis)?;
    let extinction = extinction_symbol(number).map_err(CliError::analysis)?;
    if cli.json {
        println!(
            "{}",
            json!({"number": number, "free_parameters": parameters, "extinction": extinction})
        );
    } else {
        println!(
            "space group {number}: free parameters {}, extinction {extinction}",
            parameters.join(", ")
        );
    }
    Ok(())
}
