//! Phase 9 local-environment commands.

use xtalforge::analysis::{BvsOptions, bond_valence_sums, oxidation_states};
use xtalforge::core::{
    CoordinationStrategy, RdfOptions, bond_angle_distribution, coordination, fingerprint,
    radial_distribution, voronoi,
};

use crate::feature_args::FeatureAction;
use crate::feature_commands::{analysis_error, core_error, load};
use crate::{Cli, CliError};

pub(crate) fn run(cli: &Cli, action: &FeatureAction) -> Result<(), CliError> {
    match action {
        FeatureAction::Voronoi { path, cutoff } => {
            let structure = load(path)?;
            let report = voronoi(&structure, *cutoff).map_err(core_error)?;
            if cli.json {
                println!(
                    "{}",
                    serde_json::json!({
                        "total_volume": report.total_volume,
                        "cells": report.cells,
                    })
                );
            } else {
                println!("{} Voronoi cells", report.cells.len());
                for cell in &report.cells {
                    println!(
                        "  atom {}: volume={:.4} A^3 neighbors={}",
                        cell.index,
                        cell.volume,
                        cell.neighbors.len()
                    );
                }
            }
            Ok(())
        }
        FeatureAction::Coordination {
            path,
            index,
            strategy,
            cutoff,
            solid_angle,
        } => {
            let structure = load(path)?;
            let strategy = match strategy.as_str() {
                "cutoff" => CoordinationStrategy::Cutoff(*cutoff),
                "voronoi" => CoordinationStrategy::Voronoi {
                    cutoff: *cutoff,
                    solid_angle_min: *solid_angle,
                },
                "crystalnn" | "crystal_nn" => CoordinationStrategy::CrystalNn { cutoff: *cutoff },
                other => {
                    return Err(CliError::input(format!(
                        "unknown coordination strategy '{other}'"
                    )));
                }
            };
            let shell = coordination(&structure, *index, strategy).map_err(core_error)?;
            if cli.json {
                println!("{}", serde_json::to_string(&shell).map_err(json_error)?);
            } else {
                println!(
                    "atom {}: CN={} effective={:.2}",
                    shell.index, shell.number, shell.effective_number
                );
                for neighbor in &shell.neighbors {
                    println!(
                        "  atom {} image {:?} distance={:.4} weight={:.4}",
                        neighbor.index, neighbor.image, neighbor.distance, neighbor.weight
                    );
                }
                for shape in &shell.shapes {
                    println!("  shape {} cshm={:.4}", shape.shape, shape.cshm);
                }
            }
            Ok(())
        }
        FeatureAction::Bvs { path, cutoff } => {
            let structure = load(path)?;
            let report = bond_valence_sums(
                &structure,
                &BvsOptions {
                    cutoff: *cutoff,
                    ..BvsOptions::default()
                },
            )
            .map_err(analysis_error)?;
            if cli.json {
                println!("{}", serde_json::to_string(&report).map_err(json_error)?);
            } else {
                for site in &report.sites {
                    println!(
                        "atom {} {}: BVS={:.4} bonds={}",
                        site.index, site.element, site.bvs, site.bonds
                    );
                }
            }
            Ok(())
        }
        FeatureAction::Oxidation { path, cutoff } => {
            let structure = load(path)?;
            let report = oxidation_states(
                &structure,
                &BvsOptions {
                    cutoff: *cutoff,
                    ..BvsOptions::default()
                },
            )
            .map_err(analysis_error)?;
            if cli.json {
                println!("{}", serde_json::to_string(&report).map_err(json_error)?);
            } else {
                for (element, state) in &report.by_element {
                    println!("{element}: {state:+}");
                }
                println!(
                    "total charge={:.6} neutral={} score={:.4}",
                    report.total_charge, report.neutral, report.score
                );
            }
            Ok(())
        }
        FeatureAction::Rdf {
            path,
            rmax,
            bin,
            smoothing,
        } => {
            let structure = load(path)?;
            let options = RdfOptions {
                r_max: *rmax,
                bin_width: *bin,
                smoothing: *smoothing,
            };
            let distribution = radial_distribution(&structure, &options).map_err(core_error)?;
            let angles = bond_angle_distribution(&structure, rmax.min(4.0), 2.0).ok();
            if cli.json {
                println!(
                    "{}",
                    serde_json::json!({
                        "rdf": distribution,
                        "angles": angles,
                    })
                );
            } else {
                println!(
                    "RDF: {} bins to {:.2} A, density={:.6} atoms/A^3",
                    distribution.r.len(),
                    rmax,
                    distribution.density
                );
                println!("partial pairs: {}", distribution.partial.len());
                if let Some(angles) = angles {
                    println!("angle bins: {}", angles.counts.len());
                }
            }
            Ok(())
        }
        FeatureAction::Fingerprint {
            first,
            second,
            rmax,
            bin,
        } => {
            let options = RdfOptions {
                r_max: *rmax,
                bin_width: *bin,
                smoothing: 1,
            };
            let first = fingerprint(&load(first)?, &options).map_err(core_error)?;
            let second = fingerprint(&load(second)?, &options).map_err(core_error)?;
            let similarity = first.similarity(&second);
            if cli.json {
                println!(
                    "{}",
                    serde_json::json!({
                        "similarity": similarity,
                        "first": first,
                        "second": second,
                    })
                );
            } else {
                println!("similarity={similarity:.6}");
            }
            Ok(())
        }
        _ => Err(CliError::input("unsupported local-environment action")),
    }
}

pub(crate) fn json_error(error: serde_json::Error) -> CliError {
    CliError::new("json_error", error.to_string())
}
