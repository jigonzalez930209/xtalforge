//! The `engines` command family: input sets, parsers, electronic structure and
//! scheduler scripts (Phase 18).

use std::path::Path;

use clap::Subcommand;
use serde_json::json;

use crate::{Cli, CliError, input};

/// Subcommands of `xtalforge engines`.
#[derive(Debug, Subcommand)]
pub(crate) enum EngineAction {
    /// List the input-set presets and their provenance.
    Presets,
    /// Write one input set into a directory.
    Inputs {
        /// Preset name, for example `vasp/relax`.
        preset: String,
        /// Input structure.
        path: String,
        /// Output directory.
        #[arg(long, default_value = "engine-input")]
        out: String,
        /// JSON object merged over the preset settings.
        #[arg(long)]
        overrides: Option<String>,
    },
    /// Parse one engine output into the typed result model.
    Parse {
        /// Output code: vasp-vasprun, vasp-outcar, vasp-oszicar,
        /// vasp-eigenval, vasp-doscar, vasp-procar, qe, qe-xml, qe-bands,
        /// qe-projwfc, cp2k, lammps, lammps-dump or gulp.
        code: String,
        /// Output file.
        path: String,
    },
    /// Parse a band file and report the gap and effective masses.
    Bands {
        /// Output code (vasp-eigenval, qe-bands or qe-xml).
        code: String,
        /// Output file.
        path: String,
    },
    /// Parse a DOS file and report the total and projected curves.
    Dos {
        /// Output code (vasp-doscar).
        code: String,
        /// Output file.
        path: String,
        /// Gaussian smearing in eV.
        #[arg(long, default_value_t = 0.1)]
        sigma: f64,
    },
    /// Write a SLURM or PBS submission script.
    Script {
        /// Scheduler: slurm or pbs.
        scheduler: String,
        /// Command executed by the scheduler.
        command: String,
        /// Output directory.
        #[arg(long, default_value = "job")]
        out: String,
        /// Queue or partition.
        #[arg(long, default_value = "main")]
        queue: String,
        /// Number of nodes.
        #[arg(long, default_value_t = 1)]
        nodes: usize,
        /// Walltime.
        #[arg(long, default_value = "01:00:00")]
        walltime: String,
        /// Write the script only.
        #[arg(long)]
        dry_run: bool,
    },
    /// List the known failure patterns and their corrective actions.
    Failures,
}

pub(crate) fn run(cli: &Cli, action: &EngineAction) -> Result<(), CliError> {
    match action {
        EngineAction::Presets => {
            let presets = xtalforge::engines::list_presets()
                .into_iter()
                .map(|name| {
                    (
                        name,
                        xtalforge::engines::preset_description(name).unwrap_or(""),
                    )
                })
                .collect::<Vec<_>>();
            if cli.json {
                let values = presets
                    .iter()
                    .map(|(name, description)| json!({"name": name, "description": description}))
                    .collect::<Vec<_>>();
                println!("{}", json!(values));
            } else {
                for (name, description) in presets {
                    println!("{name}: {description}");
                }
            }
            Ok(())
        }
        EngineAction::Inputs {
            preset,
            path,
            out,
            overrides,
        } => {
            let structure = input::read_structure_path(path).map_err(CliError::core)?;
            let overrides = match overrides {
                Some(text) => serde_json::from_str(text).map_err(|error| {
                    CliError::input(format!("overrides are not a JSON object: {error}"))
                })?,
                None => serde_json::Map::new(),
            };
            let summary = xtalforge::engines::write_input_set(
                &structure,
                &xtalforge::engines::InputSetRequest {
                    preset: preset.clone(),
                    overrides,
                },
            )
            .map_err(engine_error)?;
            let directory = Path::new(out);
            std::fs::create_dir_all(directory).map_err(|error| {
                CliError::input(format!("cannot create {}: {error}", directory.display()))
            })?;
            for file in &summary.files {
                let destination = directory.join(&file.name);
                if let Some(parent) = destination.parent() {
                    std::fs::create_dir_all(parent).map_err(|error| {
                        CliError::input(format!("cannot create {}: {error}", parent.display()))
                    })?;
                }
                std::fs::write(&destination, &file.text).map_err(|error| {
                    CliError::input(format!("cannot write {}: {error}", destination.display()))
                })?;
            }
            if cli.json {
                println!("{}", serde_json::to_string(&summary).unwrap_or_default());
            } else {
                println!("{} files in {}", summary.files.len(), directory.display());
                for file in &summary.files {
                    println!("  {} ({} bytes)", file.name, file.text.len());
                }
            }
            Ok(())
        }
        EngineAction::Parse { code, path } => {
            let text = input::read_source(path).map_err(CliError::core)?;
            let result = parse(code, &text)?;
            if cli.json {
                println!("{}", serde_json::to_string(&result).unwrap_or_default());
            } else {
                print_result(&result);
            }
            Ok(())
        }
        EngineAction::Bands { code, path } => {
            let text = input::read_source(path).map_err(CliError::core)?;
            let result = parse(code, &text)?;
            let Some(bands) = &result.bands else {
                return Err(CliError::input(format!("{code} output has no band data")));
            };
            let electronic = electronic_from_bands(bands);
            let structure =
                xtalforge::engines::band_structure(&electronic).map_err(CliError::core)?;
            let gap = xtalforge::engines::band_gap(&electronic).map_err(CliError::core)?;
            if cli.json {
                println!("{}", json!({"bands": structure, "gap": gap}));
            } else {
                println!(
                    "{} bands, {} k-points, gap {:.4} eV ({})",
                    structure.bands.len(),
                    structure.kpoints.len(),
                    gap.energy,
                    if gap.direct { "direct" } else { "indirect" }
                );
            }
            Ok(())
        }
        EngineAction::Dos { code, path, sigma } => {
            let text = input::read_source(path).map_err(CliError::core)?;
            let result = parse(code, &text)?;
            let Some(dos) = &result.dos else {
                return Err(CliError::input(format!("{code} output has no DOS data")));
            };
            let smeared = smear(&dos.energies, &dos.total, *sigma)?;
            let centre = xtalforge::engines::d_band_centre(
                &dos.energies,
                dos.channels
                    .iter()
                    .find(|channel| channel.label == "d")
                    .map_or(&dos.total, |channel| &channel.values),
            )
            .map_err(CliError::core)?;
            if cli.json {
                println!(
                    "{}",
                    json!({
                        "energies": dos.energies,
                        "total": dos.total,
                        "smeared": smeared,
                        "channels": dos.channels,
                        "fermi": dos.fermi,
                        "d_band_centre": centre,
                    })
                );
            } else {
                println!(
                    "{} DOS samples, {} channels, d-band centre {:.4} eV",
                    dos.energies.len(),
                    dos.channels.len(),
                    centre
                );
            }
            Ok(())
        }
        EngineAction::Script {
            scheduler,
            command,
            out,
            queue,
            nodes,
            walltime,
            dry_run,
        } => {
            let scheduler = match scheduler.to_ascii_lowercase().as_str() {
                "slurm" => xtalforge::engines::Scheduler::Slurm,
                "pbs" => xtalforge::engines::Scheduler::Pbs,
                other => {
                    return Err(CliError::input(format!(
                        "unknown scheduler {other:?}; use slurm or pbs"
                    )));
                }
            };
            let options = xtalforge::engines::JobOptions {
                scheduler,
                queue: queue.clone(),
                nodes: *nodes,
                walltime: walltime.clone(),
                command: command.clone(),
                dry_run: *dry_run,
                ..xtalforge::engines::JobOptions::new(scheduler, command.clone())
            };
            let path =
                xtalforge::engines::write_scheduler_script(&options, out).map_err(engine_error)?;
            if cli.json {
                println!(
                    "{}",
                    json!({"path": path.display().to_string(), "dry_run": dry_run})
                );
            } else {
                println!("wrote {}", path.display());
            }
            Ok(())
        }
        EngineAction::Failures => {
            let patterns = xtalforge::engines::failure_catalogue();
            if cli.json {
                println!("{}", serde_json::to_string(&patterns).unwrap_or_default());
            } else {
                for pattern in patterns {
                    println!(
                        "{}: {} -> {}",
                        pattern.code, pattern.pattern, pattern.action
                    );
                }
            }
            Ok(())
        }
    }
}

fn parse(code: &str, text: &str) -> Result<xtalforge::engines::CalculationResult, CliError> {
    let result = match code {
        "vasp-vasprun" | "vasprun" => xtalforge::engines::parse_vasp_vasprun(text),
        "vasp-outcar" | "outcar" => xtalforge::engines::parse_vasp_outcar(text),
        "vasp-oszicar" | "oszicar" => xtalforge::engines::parse_vasp_oszicar(text),
        "vasp-eigenval" | "eigenval" => xtalforge::engines::parse_vasp_eigenval(text),
        "vasp-doscar" | "doscar" => xtalforge::engines::parse_vasp_doscar(text),
        "vasp-procar" | "procar" => xtalforge::engines::parse_vasp_procar(text),
        "qe" | "pw" => xtalforge::engines::parse_qe(text),
        "qe-xml" => xtalforge::engines::parse_qe_xml(text),
        "qe-bands" => xtalforge::engines::parse_qe_bands(text),
        "qe-projwfc" => xtalforge::engines::parse_qe_projwfc(text),
        "cp2k" => xtalforge::engines::parse_cp2k(text),
        "lammps" => xtalforge::engines::parse_lammps(text),
        "lammps-dump" => xtalforge::engines::parse_lammps_dump(text),
        "gulp" => xtalforge::engines::parse_gulp(text),
        other => {
            return Err(CliError::input(format!(
                "unknown engine output code {other:?}"
            )));
        }
    };
    result.map_err(engine_error)
}

fn electronic_from_bands(
    bands: &xtalforge::engines::BandData,
) -> xtalforge::engines::ElectronicStructure {
    let n_bands = bands.energies.first().map_or(0, Vec::len);
    let mut energies = Vec::with_capacity(n_bands);
    for band in 0..n_bands {
        energies.push(
            bands
                .energies
                .iter()
                .map(|row| row.get(band).copied().unwrap_or(f64::NAN))
                .collect(),
        );
    }
    let mut structure =
        xtalforge::engines::ElectronicStructure::new(bands.kpoints.clone(), energies);
    if !bands.occupations.is_empty() {
        let mut occupations = Vec::with_capacity(n_bands);
        for band in 0..n_bands {
            occupations.push(
                bands
                    .occupations
                    .iter()
                    .map(|row| row.get(band).copied().unwrap_or(0.0))
                    .collect(),
            );
        }
        structure = structure.with_occupations(occupations);
    }
    if let Some(fermi) = bands.fermi {
        structure = structure.with_fermi(fermi);
    }
    structure
}

fn smear(energies: &[f64], values: &[f64], sigma: f64) -> Result<Vec<f64>, CliError> {
    xtalforge::engines::smear_dos(energies, values, sigma).map_err(CliError::core)
}

fn print_result(result: &xtalforge::engines::CalculationResult) {
    println!(
        "engine {}: {} energies, {} steps, SCF={} ionic={}",
        result.engine.as_deref().unwrap_or("unknown"),
        result.energies.len(),
        result.steps.len(),
        result.convergence.scf,
        result.convergence.ionic
    );
    if let Some(last) = result.energies.last() {
        println!("final energy {last:.9} eV");
    }
    for warning in &result.warnings {
        println!("warning {}: {}", warning.code, warning.message);
    }
}

fn engine_error(error: xtalforge::engines::Error) -> CliError {
    CliError::new(error.code(), error.to_string())
}
