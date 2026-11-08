//! The `ff prepare` and `ff external` commands.

use clap::Subcommand;
use serde_json::json;

use crate::{Cli, CliError, input, output};

use crate::ff_commands::calculator;

/// Subcommands of `xtalforge ff external`.
#[derive(Debug, Subcommand)]
pub(crate) enum ExternalAction {
    /// Check whether an engine binary runs.
    Probe {
        /// Engine: dftb+, xtb or orca.
        engine: String,
    },
    /// Write a DFTB+ `dftb_in.hsd`.
    Dftb {
        /// Input path.
        path: String,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
    },
    /// Write an xTB input.
    Xtb {
        /// Input path.
        path: String,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
    },
    /// Write an ORCA input.
    Orca {
        /// Input path.
        path: String,
        /// ORCA method.
        #[arg(long, default_value = "B3LYP")]
        method: String,
        /// ORCA basis.
        #[arg(long, default_value = "def2-SVP")]
        basis: String,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
    },
}

pub(crate) fn prepare(
    cli: &Cli,
    path: &str,
    minimize: bool,
    charges: bool,
    ewald: bool,
    vibrations: bool,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let options = xtalforge::ff::PrepareOptions {
        minimize,
        charges,
        ewald,
        vibrations,
        ..xtalforge::ff::PrepareOptions::default()
    };
    let report =
        xtalforge::ff::prepare(&calculator()?, &structure, &options).map_err(CliError::ff)?;
    if cli.json {
        println!(
            "{}",
            json!({
                "energy": report.energy,
                "total_charge": report.total_charge,
                "ewald_energy": report.ewald_energy,
                "frequencies": report.frequencies,
                "steps": report.steps.iter().map(|step| json!({
                    "name": step.name,
                    "done": step.done,
                    "detail": step.detail,
                })).collect::<Vec<_>>(),
            })
        );
    } else {
        for step in &report.steps {
            println!(
                "{}: {} ({})",
                step.name,
                if step.done { "done" } else { "skipped" },
                step.detail
            );
        }
        println!(
            "energy {:.9} eV, total charge {:.6}",
            report.energy, report.total_charge
        );
        if let Some(value) = report.ewald_energy {
            println!("Ewald {value:.9} eV");
        }
    }
    Ok(())
}

pub(crate) fn external(cli: &Cli, action: &ExternalAction) -> Result<(), CliError> {
    match action {
        ExternalAction::Probe { engine } => {
            let kind = match engine.as_str() {
                "dftb+" | "dftb" => xtalforge::ff::EngineKind::DftbPlus,
                "xtb" => xtalforge::ff::EngineKind::Xtb,
                "orca" => xtalforge::ff::EngineKind::Orca,
                other => {
                    return Err(CliError::input(format!(
                        "unknown engine {other:?}; use dftb+, xtb or orca"
                    )));
                }
            };
            let status = xtalforge::ff::probe(kind);
            if cli.json {
                println!(
                    "{}",
                    json!({"engine": engine, "status": format!("{status:?}")})
                );
            } else {
                println!("{engine}: {status:?}");
            }
            Ok(())
        }
        ExternalAction::Dftb { path, out } => {
            let structure = input::read_structure_path(path).map_err(CliError::core)?;
            let text = xtalforge::ff::write_dftb_hsd(&structure, &Default::default())
                .map_err(CliError::ff)?;
            write_text(cli, out.as_deref(), &text)
        }
        ExternalAction::Xtb { path, out } => {
            let structure = input::read_structure_path(path).map_err(CliError::core)?;
            let text = xtalforge::ff::write_xtb_input(&structure, &Default::default())
                .map_err(CliError::ff)?;
            write_text(cli, out.as_deref(), &text)
        }
        ExternalAction::Orca {
            path,
            method,
            basis,
            out,
        } => {
            let structure = input::read_structure_path(path).map_err(CliError::core)?;
            let options = xtalforge::ff::OrcaOptions {
                method: method.clone(),
                basis: basis.clone(),
                ..Default::default()
            };
            let text =
                xtalforge::ff::write_orca_input(&structure, &options).map_err(CliError::ff)?;
            write_text(cli, out.as_deref(), &text)
        }
    }
}

fn write_text(cli: &Cli, out: Option<&str>, text: &str) -> Result<(), CliError> {
    match out {
        Some(destination) => output::write_destination(destination, text),
        None => {
            if cli.json {
                println!("{}", json!({"text": text}));
            } else {
                print!("{text}");
            }
            Ok(())
        }
    }
}
