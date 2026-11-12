//! Phase 9 defect, ordering and SQS commands.

use xtalforge::analysis::{
    OrderingOptions, SqsOptions, interstitials, orderings, special_quasirandom, substitutions,
    vacancies,
};

use crate::feature_args::FeatureAction;
use crate::feature_commands::{analysis_error, emit_structure, load};
use crate::{Cli, CliError};

pub(crate) fn run(cli: &Cli, action: &FeatureAction) -> Result<(), CliError> {
    match action {
        FeatureAction::Vacancies {
            path,
            min_image,
            index,
            out,
            format,
        } => {
            let defects = vacancies(&load(path)?, *min_image).map_err(analysis_error)?;
            emit_defects(cli, &defects, *index, out.as_deref(), *format)
        }
        FeatureAction::Substitutions {
            path,
            species,
            min_image,
            index,
            out,
            format,
        } => {
            let defects =
                substitutions(&load(path)?, species, *min_image).map_err(analysis_error)?;
            emit_defects(cli, &defects, *index, out.as_deref(), *format)
        }
        FeatureAction::Interstitials {
            path,
            species,
            clearance,
            cutoff,
            min_image,
        } => {
            let defects = interstitials(&load(path)?, species, *min_image, *cutoff, *clearance)
                .map_err(analysis_error)?;
            if cli.json {
                println!("{}", serde_json::to_string(&defects).map_err(json_error)?);
            } else {
                for defect in &defects {
                    println!(
                        "interstitial {} at {:?} supercell {:?} image distance {:.3}",
                        defect.species, defect.frac, defect.supercell, defect.image_distance
                    );
                }
            }
            Ok(())
        }
        FeatureAction::Orderings {
            path,
            supercell,
            max_orderings,
            index,
            out,
            format,
        } => {
            let supercell = parse_supercell(supercell)?;
            let results = orderings(
                &load(path)?,
                &OrderingOptions {
                    supercell,
                    max_orderings: *max_orderings,
                    ..OrderingOptions::default()
                },
            )
            .map_err(analysis_error)?;
            if cli.json {
                let summaries: Vec<serde_json::Value> = results
                    .iter()
                    .enumerate()
                    .map(|(position, ordering)| {
                        serde_json::json!({
                            "index": position,
                            "energy": ordering.energy,
                            "assignment": ordering.assignment,
                            "n_atoms": ordering.structure.sites().len(),
                        })
                    })
                    .collect();
                println!(
                    "{}",
                    serde_json::json!({
                        "n_orderings": results.len(),
                        "orderings": summaries,
                    })
                );
            } else {
                println!("{} orderings", results.len());
                for (position, ordering) in results.iter().enumerate() {
                    println!(
                        "  {} energy={:?} assignment={:?}",
                        position, ordering.energy, ordering.assignment
                    );
                }
            }
            let selected = index.unwrap_or(0);
            if let Some(ordering) = results.get(selected) {
                emit_structure(cli, &ordering.structure, out.as_deref(), *format)?;
            } else if out.is_some() {
                return Err(CliError::input(format!(
                    "ordering index {selected} is out of range"
                )));
            }
            Ok(())
        }
        FeatureAction::Sqs {
            path,
            supercell,
            seed,
            cutoff,
            iterations,
            out,
            format,
        } => {
            let result = special_quasirandom(
                &load(path)?,
                &SqsOptions {
                    supercell: parse_supercell(supercell)?,
                    seed: *seed,
                    cutoff: *cutoff,
                    iterations: *iterations,
                    ..SqsOptions::default()
                },
            )
            .map_err(analysis_error)?;
            if cli.json {
                println!(
                    "{}",
                    serde_json::json!({
                        "assignment": result.assignment,
                        "shells": result.shells,
                        "triplet_correlation": result.triplet_correlation,
                        "objective": result.objective,
                        "accepted": result.accepted,
                        "attempted": result.attempted,
                        "n_atoms": result.structure.sites().len(),
                    })
                );
            } else {
                println!(
                    "SQS: objective={:.6} triplet={:.4} accepted={}/{}",
                    result.objective, result.triplet_correlation, result.accepted, result.attempted
                );
                for shell in &result.shells {
                    println!(
                        "  r={:.3} pairs={} correlation={:.4}",
                        shell.distance, shell.pairs, shell.correlation
                    );
                }
            }
            emit_structure(cli, &result.structure, out.as_deref(), *format)
        }
        _ => Err(CliError::input("unsupported defect action")),
    }
}

fn emit_defects(
    cli: &Cli,
    defects: &[xtalforge::analysis::PointDefect],
    index: Option<usize>,
    out: Option<&str>,
    format: Option<crate::format::StructureFormat>,
) -> Result<(), CliError> {
    if cli.json {
        let summaries: Vec<serde_json::Value> = defects
            .iter()
            .enumerate()
            .map(|(position, defect)| {
                serde_json::json!({
                    "index": position,
                    "kind": defect.kind,
                    "site_index": defect.site_index,
                    "original": defect.original,
                    "species": defect.species,
                    "frac": defect.frac,
                    "supercell": defect.supercell,
                    "image_distance": defect.image_distance,
                    "n_atoms": defect.structure.sites().len(),
                })
            })
            .collect();
        println!(
            "{}",
            serde_json::json!({
                "n_defects": defects.len(),
                "defects": summaries,
            })
        );
    } else {
        println!("{} defects", defects.len());
        for (position, defect) in defects.iter().enumerate() {
            println!(
                "  {} {:?} {:?} -> {} supercell {:?} image distance {:.3}",
                position,
                defect.kind,
                defect.original,
                defect.species,
                defect.supercell,
                defect.image_distance
            );
        }
    }
    let selected = index.unwrap_or(0);
    if let Some(defect) = defects.get(selected) {
        emit_structure(cli, &defect.structure, out, format)?;
    } else if out.is_some() {
        return Err(CliError::input(format!(
            "defect index {selected} is out of range"
        )));
    }
    Ok(())
}

fn parse_supercell(values: &str) -> Result<[u32; 3], CliError> {
    let values: Vec<u32> = values
        .split(',')
        .map(|value| {
            value
                .trim()
                .parse()
                .map_err(|_| CliError::input("--supercell needs three integers: a,b,c"))
        })
        .collect::<Result<_, _>>()?;
    match values.as_slice() {
        [a, b, c] => Ok([*a, *b, *c]),
        _ => Err(CliError::input("--supercell needs three integers: a,b,c")),
    }
}

pub(crate) fn json_error(error: serde_json::Error) -> CliError {
    CliError::new("json_error", error.to_string())
}
