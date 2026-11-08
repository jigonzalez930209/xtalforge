//! The `ff` command family: energies, optimizers, charges and vibrations.

use clap::Subcommand;
use serde_json::json;
use xtalforge::ff::{
    Calculator, CalculatorOptions, ChargeOptions, HessianOptions, IrOptions, Minimizer,
    OptimizeOptions, ThermoOptions, UffCalculator, UffModel, eqeq_charges, hessian, ir_spectrum,
    minimize, qeq_charges, thermochemistry, vibrational_modes,
};

use crate::engine_commands::{self, ExternalAction};
use crate::format::StructureFormat;
use crate::{Cli, CliError, input, output};

/// Subcommands of `xtalforge ff`.
#[derive(Debug, Subcommand)]
pub(crate) enum FfAction {
    /// Print the UFF energy of a structure.
    Energy {
        /// Input path, or `-` for standard input.
        path: String,
    },
    /// Relax a structure with FIRE or L-BFGS.
    Optimize {
        /// Input path, or `-` for standard input.
        path: String,
        /// Minimizer: fire or lbfgs.
        #[arg(long, default_value = "fire")]
        minimizer: String,
        /// Maximum number of outer steps.
        #[arg(long, default_value_t = 500)]
        max_steps: usize,
        /// Force tolerance in eV/Angstrom.
        #[arg(long, default_value_t = 0.05)]
        force_tolerance: f64,
        /// Whether the cell shape and volume are variables.
        #[arg(long)]
        cell: bool,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Compute QEq or EQeq charges.
    Charges {
        /// Input path, or `-` for standard input.
        path: String,
        /// Method: qeq, eqeq or eeq.
        #[arg(long, default_value = "qeq")]
        method: String,
        /// Total charge of the cell.
        #[arg(long, default_value_t = 0.0)]
        total_charge: f64,
    },
    /// Gamma-point vibrational modes.
    Vibrations {
        /// Input path, or `-` for standard input.
        path: String,
        /// Central-difference step in Angstrom.
        #[arg(long, default_value_t = 0.005)]
        step: f64,
    },
    /// Gamma-point IR spectrum.
    Ir {
        /// Input path, or `-` for standard input.
        path: String,
        /// Line width in cm^-1.
        #[arg(long, default_value_t = 10.0)]
        fwhm: f64,
    },
    /// Harmonic thermochemistry at a temperature.
    Thermo {
        /// Input path, or `-` for standard input.
        path: String,
        /// Temperature in kelvin.
        #[arg(long, default_value_t = 298.15)]
        temperature: f64,
        /// Pressure in bar.
        #[arg(long, default_value_t = 1.0)]
        pressure: f64,
    },
    /// Print the UFF atom types with their confidence.
    Types {
        /// Input path, or `-` for standard input.
        path: String,
    },
    /// Run the preparation pipeline: relax, charges, Ewald, vibrations.
    Prepare {
        /// Input path, or `-` for standard input.
        path: String,
        /// Skip the geometry relaxation.
        #[arg(long)]
        no_minimize: bool,
        /// Skip the charge equilibration.
        #[arg(long)]
        no_charges: bool,
        /// Evaluate the periodic Coulomb energy.
        #[arg(long)]
        ewald: bool,
        /// Compute the Gamma-point frequencies.
        #[arg(long)]
        vibrations: bool,
    },
    /// Probe or write inputs for an external engine.
    External {
        #[command(subcommand)]
        action: ExternalAction,
    },
    /// List the registered calculators and their availability.
    Engines,
}

/// Runs one `ff` action.
///
/// # Errors
///
/// Returns [`CliError`] for read, calculator or output failures.
pub(crate) fn run(cli: &Cli, action: &FfAction) -> Result<(), CliError> {
    match action {
        FfAction::Energy { path } => energy(cli, path),
        FfAction::Optimize {
            path,
            minimizer,
            max_steps,
            force_tolerance,
            cell,
            out,
            format,
        } => optimize(
            cli,
            path,
            minimizer,
            *max_steps,
            *force_tolerance,
            *cell,
            out.as_deref(),
            *format,
        ),
        FfAction::Charges {
            path,
            method,
            total_charge,
        } => charges(cli, path, method, *total_charge),
        FfAction::Vibrations { path, step } => vibrations(cli, path, *step),
        FfAction::Ir { path, fwhm } => ir(cli, path, *fwhm),
        FfAction::Thermo {
            path,
            temperature,
            pressure,
        } => thermo(cli, path, *temperature, *pressure),
        FfAction::Types { path } => types(cli, path),
        FfAction::Prepare {
            path,
            no_minimize,
            no_charges,
            ewald,
            vibrations,
        } => engine_commands::prepare(cli, path, !*no_minimize, !*no_charges, *ewald, *vibrations),
        FfAction::External { action } => engine_commands::external(cli, action),
        FfAction::Engines => engines(cli),
    }
}

pub(crate) fn calculator() -> Result<UffCalculator, CliError> {
    UffCalculator::with_defaults().map_err(CliError::ff)
}

fn energy(cli: &Cli, path: &str) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let energy = calculator()?.energy(&structure).map_err(CliError::ff)?;
    if cli.json {
        println!("{}", json!({"calculator": "uff", "energy": energy}));
    } else {
        println!("UFF energy {energy:.9} eV");
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn optimize(
    cli: &Cli,
    path: &str,
    minimizer: &str,
    max_steps: usize,
    force_tolerance: f64,
    cell: bool,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let minimizer = match minimizer {
        "fire" => Minimizer::Fire,
        "lbfgs" => Minimizer::Lbfgs,
        other => {
            return Err(CliError::input(format!(
                "unknown minimizer {other:?}; use fire or lbfgs"
            )));
        }
    };
    let options = OptimizeOptions {
        max_steps,
        force_tolerance,
        cell,
        ..OptimizeOptions::default()
    };
    let result = minimize(&calculator()?, &structure, &options, minimizer).map_err(CliError::ff)?;
    if out.is_some() {
        return output::emit_structure(cli, &result.structure, out, format);
    }
    if cli.json {
        println!(
            "{}",
            json!({
                "energy": result.energy,
                "max_force": result.max_force,
                "steps": result.steps,
                "converged": result.converged,
            })
        );
    } else {
        println!(
            "{} after {} steps: E = {:.9} eV, max force {:.6} eV/A, converged {}",
            if cell { "variable-cell" } else { "fixed-cell" },
            result.steps,
            result.energy,
            result.max_force,
            result.converged
        );
        print!(
            "{}",
            xtalforge::io::write_structure(&result.structure, "cif").map_err(CliError::core)?
        );
    }
    Ok(())
}

fn charges(cli: &Cli, path: &str, method: &str, total_charge: f64) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let options = ChargeOptions {
        total_charge,
        ..ChargeOptions::default()
    };
    let charges = match method {
        "qeq" => qeq_charges(&structure, &options),
        "eqeq" => eqeq_charges(&structure, &options),
        "eeq" => xtalforge::ff::eeq_charges(&structure, &options),
        other => {
            return Err(CliError::input(format!(
                "unknown charge method {other:?}; use qeq, eqeq or eeq"
            )));
        }
    }
    .map_err(CliError::ff)?;
    if cli.json {
        println!(
            "{}",
            json!({
                "method": method,
                "total": charges.total,
                "per_atom": charges.per_atom,
            })
        );
    } else {
        println!("{method}: total charge {:.6}", charges.total);
        for (index, charge) in charges.per_atom.iter().enumerate() {
            println!("  atom {index}: {charge:.6}");
        }
    }
    Ok(())
}

fn modes_for(
    path: &str,
    step: f64,
) -> Result<(xtalforge::core::Structure, xtalforge::ff::VibrationalModes), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let p1 = xtalforge::symmetry::to_p1(&structure).map_err(CliError::core)?;
    let options = HessianOptions {
        step,
        ..HessianOptions::default()
    };
    let matrix = hessian(&calculator()?, &p1, &options).map_err(CliError::ff)?;
    let modes = vibrational_modes(&matrix, &p1, &options).map_err(CliError::ff)?;
    Ok((p1, modes))
}

fn vibrations(cli: &Cli, path: &str, step: f64) -> Result<(), CliError> {
    let (_, modes) = modes_for(path, step)?;
    if cli.json {
        println!(
            "{}",
            json!({
                "frequencies": modes.frequencies,
                "zero_modes": modes.zero_modes,
                "imaginary": modes.imaginary,
            })
        );
    } else {
        for (index, frequency) in modes.frequencies.iter().enumerate() {
            println!("mode {index}: {frequency:.3} cm^-1");
        }
    }
    Ok(())
}

fn ir(cli: &Cli, path: &str, fwhm: f64) -> Result<(), CliError> {
    let (p1, modes) = modes_for(path, 0.005)?;
    let spectrum = ir_spectrum(
        &p1,
        &modes,
        &IrOptions {
            fwhm,
            ..IrOptions::default()
        },
    )
    .map_err(CliError::ff)?;
    if cli.json {
        println!(
            "{}",
            json!({
                "peaks": spectrum.peaks.iter().map(|peak| json!({
                    "mode": peak.mode,
                    "frequency": peak.frequency,
                    "intensity": peak.intensity,
                })).collect::<Vec<_>>(),
                "points": spectrum.points,
            })
        );
    } else {
        for peak in &spectrum.peaks {
            println!(
                "mode {}: {:.3} cm^-1, {:.3} km/mol",
                peak.mode, peak.frequency, peak.intensity
            );
        }
    }
    Ok(())
}

fn thermo(cli: &Cli, path: &str, temperature: f64, pressure: f64) -> Result<(), CliError> {
    let (_, modes) = modes_for(path, 0.005)?;
    let result = thermochemistry(
        &modes.frequencies,
        &ThermoOptions {
            temperature,
            pressure,
        },
    )
    .map_err(CliError::ff)?;
    if cli.json {
        println!(
            "{}",
            json!({
                "temperature": temperature,
                "pressure": pressure,
                "zpe": result.zpe,
                "thermal_u": result.thermal_u,
                "thermal_h": result.thermal_h,
                "entropy": result.entropy,
                "gibbs": result.gibbs,
            })
        );
    } else {
        println!("T = {temperature} K, p = {pressure} bar");
        println!("ZPE {:.6} eV", result.zpe);
        println!("U {:.6} eV, H {:.6} eV", result.thermal_u, result.thermal_h);
        println!("S {:.6} eV/K, G {:.6} eV", result.entropy, result.gibbs);
    }
    Ok(())
}

fn types(cli: &Cli, path: &str) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let model = UffModel::build(&structure, CalculatorOptions::default()).map_err(CliError::ff)?;
    if cli.json {
        let types: Vec<serde_json::Value> = model
            .typing()
            .types
            .iter()
            .enumerate()
            .map(|(atom, entry)| {
                json!({
                    "atom": atom,
                    "type": entry.name,
                    "confidence": format!("{:?}", entry.confidence),
                    "reason": entry.reason,
                })
            })
            .collect();
        println!("{}", json!(types));
    } else {
        for (atom, typing) in model.typing().types.iter().enumerate() {
            println!("atom {atom}: {} ({})", typing.name, typing.reason);
        }
    }
    Ok(())
}

fn engines(cli: &Cli) -> Result<(), CliError> {
    let registry = xtalforge::ff::builtin_registry();
    let infos: Vec<serde_json::Value> = registry
        .infos()
        .map(|info| {
            json!({
                "name": info.name,
                "label": info.label,
                "description": info.description,
                "availability": format!("{:?}", info.availability),
            })
        })
        .collect();
    if cli.json {
        println!("{}", json!(infos));
    } else {
        for info in registry.infos() {
            println!(
                "{:<10} {} — {} ({:?})",
                info.name, info.label, info.description, info.availability
            );
        }
    }
    Ok(())
}
