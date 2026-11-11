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
        /// Refine positions, ADPs, zero shift, cell and width.
        #[arg(long)]
        refine: bool,
        /// Also refine one occupancy per asymmetric site.
        #[arg(long, requires = "refine")]
        occupancies: bool,
        /// Also refine the March-Dollase texture about the c axis.
        #[arg(long, requires = "refine")]
        texture: bool,
        /// Model one wavelength instead of the Cu K alpha 1/K alpha 2 doublet.
        #[arg(long, requires = "refine")]
        no_doublet: bool,
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
            occupancies,
            texture,
            no_doublet,
        } => rietveld(
            cli,
            structure,
            pattern,
            *refine,
            *occupancies,
            *texture,
            !*no_doublet,
        ),
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
    // `--background` also applies the noise threshold, so a noisy pattern on a
    // raised baseline is not indexed on its noise maxima.
    let peak_options = PeakOptions {
        refine: refine_peaks,
        ..if background {
            PeakOptions::noise_limited(&points)
        } else {
            PeakOptions::default()
        }
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
                "intensity_esds": result.intensity_esds,
                "positions": result.positions,
                "r_factor": result.r_factor,
                "rwp": result.rwp,
                "chi_squared": result.chi_squared,
            })
        );
    } else {
        println!(
            "{} reflections, R = {:.4}, Rwp = {:.4}, chi2 = {:.4}",
            result.intensities.len(),
            result.r_factor,
            result.rwp,
            result.chi_squared
        );
    }
    Ok(())
}

fn rietveld(
    cli: &Cli,
    structure_path: &str,
    pattern_path: &str,
    refine: bool,
    occupancies: bool,
    texture: bool,
    doublet: bool,
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
                refine_positions: true,
                refine_adps: true,
                refine_occupancies: occupancies,
                refine_zero: true,
                refine_cell: true,
                refine_width: true,
                refine_texture: texture,
                doublet,
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
                    "stages": result.stages,
                    "parameters": result.parameters,
                    "parameter_esds": result.parameter_esds,
                })
            );
        } else {
            println!(
                "Rwp {:.4}, converged {}, {} iterations, stages {}",
                result.rwp,
                result.converged,
                result.iterations,
                result.stages.join(" -> ")
            );
            for (name, value) in &result.parameters {
                let esd = result
                    .parameter_esds
                    .iter()
                    .find(|(key, _)| key == name)
                    .map_or(f64::NAN, |(_, value)| *value);
                println!("  {name} = {value:.6} (esd {esd:.6})");
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
    // The Pawley cell of the leading candidate, or its indexed cell when the
    // Pawley fit was skipped (the same fallback as Python and WASM).
    let best_cell = result
        .pawley
        .as_ref()
        .map(|pawley| pawley.cell)
        .or_else(|| {
            result.leading_cells.first().map(|row| {
                let cell = &row.index.cell;
                [cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma]
            })
        });
    if cli.json {
        println!(
            "{}",
            json!({
                "peaks": result.peaks.len(),
                "index": result.index,
                "pawley_rwp": result.pawley.as_ref().map(|pawley| pawley.rwp),
                "best_cell": best_cell,
                "leading_cells": result.leading_cells,
                "lebail": result.lebail.as_ref().map(|lebail| lebail.r_factor),
                "rwp": result.phase.as_ref().map(|phase| phase.rwp),
                "steps": result.steps,
            })
        );
    } else {
        for step in &result.steps {
            println!("{step}");
        }
        if let Some(pawley) = &result.pawley {
            let [a, b, c, alpha, beta, gamma] = pawley.cell;
            println!(
                "best cell (Pawley): a={a:.4} b={b:.4} c={c:.4} alpha={alpha:.2} beta={beta:.2} gamma={gamma:.2} Rwp={:.4}",
                pawley.rwp
            );
        }
        for (rank, row) in result.leading_cells.iter().enumerate() {
            println!(
                "{}. cell {} class {} {} dBIC={:.2} {} a={:.4} b={:.4} c={:.4} beta={:.2} M20={:.1} F_N={:.1} Rwp={}",
                rank + 1,
                row.cell_rank,
                row.class_rank,
                row.class_symbol,
                row.delta_bic,
                row.index.system.name(),
                row.index.cell.a,
                row.index.cell.b,
                row.index.cell.c,
                row.index.cell.beta,
                row.index.m20,
                row.index.f_n,
                row.pawley_rwp
                    .map_or_else(|| "not fitted".to_owned(), |value| format!("{value:.4}"))
            );
            if !row.space_groups.is_empty() {
                let groups = row
                    .space_groups
                    .iter()
                    .map(|group| {
                        format!(
                            "{} {} ({}, order {})",
                            group.number, group.symbol, group.centering, group.order
                        )
                    })
                    .collect::<Vec<_>>()
                    .join("; ");
                println!("   space groups: {groups}");
            }
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
