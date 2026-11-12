//! Phase 9 edit commands: rattle, site edits, interpolation and optimal cells.

use xtalforge::core::{RattleOptions, StructureSummary};

use crate::feature_args::FeatureAction;
use crate::feature_commands::{core_error, edit, emit_structure, load, write_frames};
use crate::{Cli, CliError};

pub(crate) fn run(cli: &Cli, action: &FeatureAction) -> Result<(), CliError> {
    match action {
        FeatureAction::Rattle {
            path,
            amplitude,
            seed,
            bounded,
            out,
            format,
        } => {
            let structure = load(path)?;
            let report = xtalforge::core::rattle(
                &structure,
                &RattleOptions {
                    amplitude: *amplitude,
                    seed: *seed,
                    bounded: *bounded,
                    ..RattleOptions::default()
                },
            )
            .map_err(core_error)?;
            if cli.json {
                println!(
                    "{}",
                    serde_json::json!({
                        "n_atoms": report.n_atoms,
                        "max_displacement": report.max_displacement,
                        "rms_displacement": report.rms_displacement,
                    })
                );
            } else if !cli.quiet {
                println!(
                    "rattled {} atoms: max={:.4} A rms={:.4} A",
                    report.n_atoms, report.max_displacement, report.rms_displacement
                );
            }
            emit_structure(cli, &report.structure, out.as_deref(), *format)
        }
        FeatureAction::Edit {
            path,
            add,
            remove,
            replace,
            move_site,
            species,
            fraction,
            seed,
            sort,
            wrap,
            out,
            format,
        } => edit(
            cli,
            path,
            add,
            remove,
            replace,
            move_site,
            species,
            *fraction,
            *seed,
            sort.as_deref(),
            *wrap,
            out.as_deref(),
            *format,
        ),
        FeatureAction::Interpolate {
            first,
            second,
            images,
            out,
            format: _,
        } => {
            let start = load(first)?;
            let end = load(second)?;
            let options = xtalforge::symmetry::MatcherOptions {
                allow_supercell: true,
                ..xtalforge::symmetry::MatcherOptions::default()
            };
            let frames = xtalforge::symmetry::interpolate_matched(&start, &end, *images, &options)
                .map_err(core_error)?;
            if cli.json {
                let summaries: Vec<StructureSummary> = frames
                    .iter()
                    .map(|frame| frame.summary())
                    .collect::<Result<_, _>>()
                    .map_err(core_error)?;
                println!(
                    "{}",
                    serde_json::json!({
                        "n_frames": frames.len(),
                        "frames": summaries,
                    })
                );
            } else if !cli.quiet {
                println!("interpolated {} frames", frames.len());
            }
            write_frames(cli, &frames, out.as_deref())
        }
        FeatureAction::OptimalSupercell {
            path,
            cells,
            out,
            format,
        } => {
            let structure = load(path)?;
            let choice =
                xtalforge::core::optimal_supercell(&structure, *cells).map_err(core_error)?;
            if cli.json {
                println!(
                    "{}",
                    serde_json::json!({
                        "matrix": choice.matrix,
                        "cells": choice.cells,
                        "score": choice.score,
                        "n_atoms": choice.structure.sites().len(),
                    })
                );
            } else if !cli.quiet {
                println!(
                    "supercell det={} score={:.6} sites={}",
                    choice.cells,
                    choice.score,
                    choice.structure.sites().len()
                );
            }
            emit_structure(cli, &choice.structure, out.as_deref(), *format)
        }
        _ => Err(CliError::input("unsupported edit action")),
    }
}
