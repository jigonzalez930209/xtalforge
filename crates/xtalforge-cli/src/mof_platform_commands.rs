//! MOF platform commands: checking, deconstruction, descriptors, porosity,
//! modification and database screening.

use serde_json::json;

use crate::feature_commands::{load, parse_numbers};
use crate::mof_commands::{DatabaseAction, ModifyAction};
use crate::{Cli, CliError, output};

fn emit_json<T: serde::Serialize>(value: &T) -> Result<(), CliError> {
    let text = serde_json::to_string(value)
        .map_err(|error| CliError::new("json_error", error.to_string()))?;
    println!("{text}");
    Ok(())
}

pub(crate) fn check(cli: &Cli, path: &str) -> Result<(), CliError> {
    let report = xtalforge::mof::check_mof(&load(path)?).map_err(CliError::core)?;
    if cli.json {
        emit_json(&report)
    } else {
        println!(
            "clean={} errors={} warnings={}",
            report.clean,
            report.counts.get("error").copied().unwrap_or(0),
            report.counts.get("warning").copied().unwrap_or(0)
        );
        for finding in &report.findings {
            println!(
                "  [{}] {}: {} atoms={:?} fix={:?}",
                finding.severity.as_str(),
                finding.code,
                finding.message,
                finding.atoms,
                finding.fix
            );
        }
        Ok(())
    }
}

pub(crate) fn solvent_remove(
    cli: &Cli,
    path: &str,
    remove_ions: bool,
    out: Option<&str>,
    format: Option<crate::format::StructureFormat>,
) -> Result<(), CliError> {
    let report = xtalforge::mof::remove_solvent(
        &load(path)?,
        &xtalforge::mof::SolventOptions {
            remove_ions,
            ..xtalforge::mof::SolventOptions::default()
        },
    )
    .map_err(CliError::core)?;
    if cli.json {
        emit_json(&report)?;
    } else {
        println!(
            "removed {} fragments, kept {}",
            report.removed.len(),
            report.kept
        );
        for fragment in &report.removed {
            println!(
                "  {} ({:?}) atoms={:?} ion={}",
                fragment.formula, fragment.solvent, fragment.atoms, fragment.ion
            );
        }
    }
    output::emit_structure(cli, &report.structure, out, format)
}

pub(crate) fn disorder(
    cli: &Cli,
    path: &str,
    out: Option<&str>,
    format: Option<crate::format::StructureFormat>,
) -> Result<(), CliError> {
    let structure =
        xtalforge::mof::resolve_disorder(&load(path)?, &xtalforge::mof::DisorderOptions::default())
            .map_err(CliError::core)?;
    if cli.json {
        emit_json(&structure.summary().map_err(CliError::core)?)?;
    }
    output::emit_structure(cli, &structure, out, format)
}

pub(crate) fn charge_balance(cli: &Cli, path: &str) -> Result<(), CliError> {
    let report = xtalforge::mof::charge_balance(&load(path)?).map_err(CliError::core)?;
    if cli.json {
        return emit_json(&report);
    }
    for metal in &report.metals {
        println!(
            "{}: state {:+} atoms={:?} BVS={:.3}",
            metal.element, metal.state, metal.atoms, metal.bvs
        );
    }
    println!(
        "linker charge {:.3}, framework {:.3}, total {:.3}, neutral={}",
        report.linker_charge, report.framework_charge, report.total_charge, report.neutral
    );
    Ok(())
}

pub(crate) fn deconstruct(cli: &Cli, path: &str) -> Result<(), CliError> {
    let report = xtalforge::mof::deconstruct(&load(path)?).map_err(CliError::core)?;
    if cli.json {
        return emit_json(&report);
    }
    println!(
        "{} units, topology={:?} nearest={:?} interpenetration={}",
        report.units.len(),
        report.topology,
        report.nearest_topologies,
        report.interpenetration
    );
    for (index, unit) in report.units.iter().enumerate() {
        println!(
            "  {index} {:?} {} atoms={} connections={} span={:.3}",
            unit.kind,
            unit.formula,
            unit.atoms.len(),
            unit.connections.len(),
            unit.span
        );
    }
    println!(
        "net: {} nodes, {} edges",
        report.net.nodes.len(),
        report.net.edges.len()
    );
    Ok(())
}

pub(crate) fn mofid(cli: &Cli, path: &str) -> Result<(), CliError> {
    let report = xtalforge::mof::mofid(&load(path)?).map_err(CliError::core)?;
    if cli.json {
        emit_json(&report)
    } else {
        println!("MOFid  {}", report.mofid);
        println!("MOFkey {}", report.mofkey);
        Ok(())
    }
}

pub(crate) fn open_metal_sites(cli: &Cli, path: &str) -> Result<(), CliError> {
    let report = xtalforge::mof::open_metal_sites(
        &load(path)?,
        &xtalforge::mof::OpenMetalOptions::default(),
    )
    .map_err(CliError::core)?;
    if cli.json {
        emit_json(&report)
    } else {
        for site in &report {
            println!(
                "atom {} {} CN={}/{} missing={} accessible={} clearance={:.3}",
                site.atom,
                site.element,
                site.coordination,
                site.expected_coordination,
                site.missing,
                site.accessible,
                site.clearance
            );
        }
        Ok(())
    }
}

pub(crate) fn descriptors(cli: &Cli, path: &str) -> Result<(), CliError> {
    let report = xtalforge::mof::framework_descriptors(&load(path)?).map_err(CliError::core)?;
    if cli.json {
        emit_json(&report)
    } else {
        println!(
            "{} atoms={} volume={:.3} density={:.4} topology={:?} interpenetration={}",
            report.formula,
            report.n_atoms,
            report.cell_volume,
            report.density,
            report.topology,
            report.interpenetration
        );
        for node in &report.nodes {
            println!(
                "  node {} metals={:?} nuclearity={} connectivity={}",
                node.formula, node.metals, node.nuclearity, node.connectivity
            );
        }
        for linker in &report.linkers {
            println!(
                "  linker {} connections={} span={:.3} groups={:?}",
                linker.formula, linker.connections, linker.span, linker.groups
            );
        }
        Ok(())
    }
}

pub(crate) fn pores(cli: &Cli, path: &str, spacing: f64, probe: f64) -> Result<(), CliError> {
    let report = xtalforge::mof::pore_descriptors(
        &load(path)?,
        &xtalforge::mof::PoreOptions {
            spacing,
            probe_radius: probe,
            channel_probe: probe,
        },
    )
    .map_err(CliError::core)?;
    if cli.json {
        emit_json(&report)
    } else {
        let pores = &report.descriptors;
        println!(
            "PLD={:.3} LCD={:.3} ASA={:.3} m2/g ({:.3} m2/cm3) AV={:.3} VF={:.4} channels={}",
            pores.pld,
            pores.lcd,
            pores.asa_m2_g,
            pores.asa_m2_cm3,
            pores.accessible_volume,
            pores.void_fraction,
            pores.channel_dimension
        );
        Ok(())
    }
}

pub(crate) fn modify(cli: &Cli, action: &ModifyAction) -> Result<(), CliError> {
    match action {
        ModifyAction::Functionalize {
            path,
            group,
            positions,
            fraction,
            seed,
            out,
            format,
        } => {
            let group = match group.as_str() {
                "amino" | "nh2" => xtalforge::mof::FunctionalGroup::Amino,
                "hydroxyl" | "oh" => xtalforge::mof::FunctionalGroup::Hydroxyl,
                "methyl" | "ch3" => xtalforge::mof::FunctionalGroup::Methyl,
                "fluoro" | "f" => xtalforge::mof::FunctionalGroup::Fluorine,
                other => {
                    return Err(CliError::input(format!(
                        "unknown functional group {other:?}"
                    )));
                }
            };
            let structure = xtalforge::mof::functionalize(
                &load(path)?,
                &xtalforge::mof::FunctionalizeOptions {
                    group,
                    positions: positions.clone(),
                    seed: *seed,
                    fraction: *fraction,
                },
            )
            .map_err(CliError::core)?;
            output::emit_structure(cli, &structure, out.as_deref(), *format)
        }
        ModifyAction::Defects {
            path,
            linker_fraction,
            cluster_fraction,
            cap,
            seed,
            out,
            format,
        } => {
            let cap = match cap.as_str() {
                "formate" => xtalforge::mof::CapSpecies::Formate,
                "water" => xtalforge::mof::CapSpecies::Water,
                "hydroxide" | "oh" => xtalforge::mof::CapSpecies::Hydroxide,
                other => return Err(CliError::input(format!("unknown cap {other:?}"))),
            };
            let report = xtalforge::mof::missing_defects(
                &load(path)?,
                &xtalforge::mof::DefectOptions {
                    linker_fraction: *linker_fraction,
                    cluster_fraction: *cluster_fraction,
                    cap,
                    seed: *seed,
                },
            )
            .map_err(CliError::core)?;
            if cli.json {
                emit_json(&report)?;
            } else {
                println!(
                    "removed {} linkers, {} clusters, added {} caps",
                    report.removed_linkers, report.removed_clusters, report.caps
                );
            }
            output::emit_structure(cli, &report.structure, out.as_deref(), *format)
        }
        ModifyAction::Interpenetrate {
            path,
            folds,
            out,
            format,
        } => {
            let structure =
                xtalforge::mof::interpenetrate(&load(path)?, *folds).map_err(CliError::core)?;
            output::emit_structure(cli, &structure, out.as_deref(), *format)
        }
        ModifyAction::StackLayers {
            path,
            layers,
            offset,
            out,
            format,
        } => {
            let values = parse_numbers(offset)?;
            let [u, v] = match values.as_slice() {
                [u, v] => [*u, *v],
                _ => return Err(CliError::input("--offset needs two numbers")),
            };
            let structure = xtalforge::mof::stack_layers(&load(path)?, *layers, [u, v])
                .map_err(CliError::core)?;
            output::emit_structure(cli, &structure, out.as_deref(), *format)
        }
        ModifyAction::Mixed {
            path,
            map,
            ordered,
            fraction,
            seed,
            out,
            format,
        } => {
            let mut replacements = std::collections::BTreeMap::new();
            for entry in map {
                let (from, to) = entry
                    .split_once(':')
                    .ok_or_else(|| CliError::input("--map needs old:new"))?;
                replacements.insert(from.trim().to_owned(), to.trim().to_owned());
            }
            let structure = xtalforge::mof::mixed_framework(
                &load(path)?,
                &xtalforge::mof::MixtureOptions {
                    map: replacements,
                    ordered: *ordered,
                    fraction: *fraction,
                    seed: *seed,
                },
            )
            .map_err(CliError::core)?;
            output::emit_structure(cli, &structure, out.as_deref(), *format)
        }
    }
}

pub(crate) fn database(cli: &Cli, action: &DatabaseAction) -> Result<(), CliError> {
    match action {
        DatabaseAction::Read { path, kind } => {
            let path = std::path::Path::new(path);
            let (entries, info) = if kind.eq_ignore_ascii_case("qmof") {
                xtalforge::mof::read_qmof(path)
            } else {
                xtalforge::mof::read_core_mof(path)
            }
            .map_err(CliError::core)?;
            if cli.json {
                let summaries: Vec<serde_json::Value> = entries
                    .iter()
                    .map(|entry| {
                        json!({
                            "id": entry.id,
                            "formula": entry.structure.formula(),
                            "properties": entry.properties,
                        })
                    })
                    .collect();
                emit_json(&json!({ "info": info, "entries": summaries }))
            } else {
                println!(
                    "{:?}: {} entries; {}",
                    info.source, info.entries, info.license
                );
                for entry in &entries {
                    println!("  {} {}", entry.id, entry.structure.formula());
                }
                Ok(())
            }
        }
        DatabaseAction::Screen {
            path,
            guest,
            spacing,
            samples,
            limit,
        } => {
            let (entries, _) = xtalforge::mof::read_core_mof(std::path::Path::new(path))
                .map_err(CliError::core)?;
            let records = xtalforge::mof::screen_database(
                &entries,
                &xtalforge::mof::ScreeningOptions {
                    spacing: *spacing,
                    guest: guest.clone(),
                    samples: *samples,
                    limit: *limit,
                    ..xtalforge::mof::ScreeningOptions::default()
                },
            )
            .map_err(CliError::core)?;
            if cli.json {
                emit_json(&records)
            } else {
                println!("{} records", records.len());
                for record in &records {
                    println!(
                        "  {} MOFkey={} PLD={:.3} LCD={:.3} Henry={:?}",
                        record.id, record.mofkey, record.pores.pld, record.pores.lcd, record.henry
                    );
                }
                Ok(())
            }
        }
        DatabaseAction::Fetch {
            url,
            cache,
            sha256,
            offline,
            kind,
        } => {
            let mut options = xtalforge::mof::FetchOptions::new(cache).offline(*offline);
            if let Some(sha256) = sha256 {
                options = options.with_sha256(sha256);
            }
            let fetched = if kind.eq_ignore_ascii_case("qmof") {
                xtalforge::mof::fetch_qmof(url, &options)
            } else {
                xtalforge::mof::fetch_core_mof(url, &options)
            }
            .map_err(CliError::core)?;
            if cli.json {
                emit_json(&fetched)
            } else {
                println!(
                    "{} {} ({} bytes, sha256 {})",
                    if fetched.from_cache {
                        "cached"
                    } else {
                        "downloaded"
                    },
                    fetched.path.display(),
                    fetched.bytes,
                    fetched.sha256
                );
                Ok(())
            }
        }
        DatabaseAction::Deduplicate { paths } => {
            let mut structures = Vec::with_capacity(paths.len());
            for path in paths {
                structures.push(load(path)?);
            }
            let groups = xtalforge::mof::deduplicate(&structures).map_err(CliError::core)?;
            if cli.json {
                emit_json(&groups)
            } else {
                println!(
                    "{} unique MOFs in {} groups",
                    structures.len(),
                    groups.len()
                );
                Ok(())
            }
        }
    }
}
