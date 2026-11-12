//! MOF adsorption commands: TraPPE guests, Widom, GCMC, grids and RASPA.

use crate::feature_commands::{load, parse_numbers};
use crate::mof_commands::AdsorptionAction;
use crate::{Cli, CliError, output};

fn emit_json<T: serde::Serialize>(value: &T) -> Result<(), CliError> {
    let text = serde_json::to_string(value)
        .map_err(|error| CliError::new("json_error", error.to_string()))?;
    println!("{text}");
    Ok(())
}

pub(crate) fn run(cli: &Cli, action: &AdsorptionAction) -> Result<(), CliError> {
    match action {
        AdsorptionAction::Guests => {
            let guests = xtalforge::mof::guests().map_err(CliError::core)?;
            if cli.json {
                emit_json(&guests)
            } else {
                for guest in guests {
                    println!("{}: {} atoms", guest.name, guest.atoms.len());
                }
                Ok(())
            }
        }
        AdsorptionAction::ForceField { path } => {
            let force_field = xtalforge::mof::force_field(&load(path)?).map_err(CliError::core)?;
            if cli.json {
                emit_json(&force_field)
            } else {
                println!(
                    "{} host atoms, guests: {}",
                    force_field.hosts.len(),
                    force_field
                        .guests
                        .keys()
                        .cloned()
                        .collect::<Vec<_>>()
                        .join(", ")
                );
                Ok(())
            }
        }
        AdsorptionAction::Widom {
            path,
            guest,
            temperature,
            samples,
            seed,
        } => {
            let report = xtalforge::mof::widom_insertion(
                &load(path)?,
                guest,
                &xtalforge::mof::WidomOptions {
                    temperature: *temperature,
                    samples: *samples,
                    seed: *seed,
                },
            )
            .map_err(CliError::core)?;
            if cli.json {
                emit_json(&report)
            } else {
                println!(
                    "{} at {} K: Henry={:.6} E={:.6} eV Qst={:.6} eV Emin={:.6} eV",
                    report.guest,
                    report.temperature,
                    report.henry,
                    report.mean_energy,
                    report.heat_of_adsorption,
                    report.minimum_energy
                );
                Ok(())
            }
        }
        AdsorptionAction::Gcmc {
            path,
            guest,
            temperature,
            pressures,
            steps,
            seed,
        } => {
            let pressures = parse_numbers(pressures)?;
            let report = xtalforge::mof::gcmc_isotherm(
                &load(path)?,
                guest,
                &xtalforge::mof::GcmcOptions {
                    temperature: *temperature,
                    pressures,
                    steps: *steps,
                    seed: *seed,
                },
            )
            .map_err(CliError::core)?;
            if cli.json {
                emit_json(&report)
            } else {
                for point in &report.points {
                    println!(
                        "P={:.3} bar loading={:.6} acceptance={:.3}",
                        point.pressure, point.loading, point.acceptance
                    );
                }
                Ok(())
            }
        }
        AdsorptionAction::Grid {
            path,
            guest,
            spacing,
        } => {
            let grid = xtalforge::mof::energy_grid(&load(path)?, guest, *spacing)
                .map_err(CliError::core)?;
            if cli.json {
                emit_json(&serde_json::json!({
                    "guest": grid.guest,
                    "shape": grid.grid.shape,
                    "minimum": grid.minimum,
                    "maximum": grid.maximum,
                }))
            } else {
                println!(
                    "{} grid {:?} min={:.6} eV max={:.6} eV",
                    grid.guest, grid.grid.shape, grid.minimum, grid.maximum
                );
                Ok(())
            }
        }
        AdsorptionAction::Raspa {
            path,
            guest,
            temperature,
            pressure,
            dir,
        } => {
            let files = xtalforge::mof::raspa_input(&load(path)?, guest, *temperature, *pressure)
                .map_err(CliError::core)?;
            if let Some(dir) = dir {
                let directory = std::path::Path::new(dir);
                std::fs::create_dir_all(directory)
                    .map_err(|error| CliError::new("io_error", error.to_string()))?;
                output::write_destination(
                    &directory.join("simulation.input").to_string_lossy(),
                    &files.simulation_input,
                )?;
                output::write_destination(
                    &directory.join("framework.cif").to_string_lossy(),
                    &files.framework_cif,
                )?;
                output::write_destination(
                    &directory.join("force_field.def").to_string_lossy(),
                    &files.force_field,
                )?;
                output::write_destination(
                    &directory.join("mixing_rules.def").to_string_lossy(),
                    &files.mixing_rules,
                )?;
                if !cli.quiet {
                    println!("wrote RASPA input to {dir}");
                }
                Ok(())
            } else if cli.json {
                emit_json(&files)
            } else {
                println!("{}", files.simulation_input);
                Ok(())
            }
        }
    }
}
