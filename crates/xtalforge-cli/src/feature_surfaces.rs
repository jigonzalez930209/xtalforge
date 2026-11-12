//! Phase 9 surface commands: Miller indices, terminations, Wulff and
//! adsorption sites.

use std::collections::BTreeSet;

use xtalforge::analysis::{slab_terminations, unique_miller_indices};
use xtalforge::core::{SurfaceEnergy, adsorption_sites, place_molecule, wulff};
use xtalforge::symmetry::{SymmetryOptions, detect_symmetry};

use crate::feature_args::FeatureAction;
use crate::feature_commands::{analysis_error, emit_structure, load};
use crate::{Cli, CliError};

pub(crate) fn run(cli: &Cli, action: &FeatureAction) -> Result<(), CliError> {
    match action {
        FeatureAction::Miller {
            path,
            max_index,
            symprec,
        } => {
            let structure = load(path)?;
            let indices =
                unique_miller_indices(&structure, *max_index, *symprec).map_err(analysis_error)?;
            if cli.json {
                println!("{}", serde_json::to_string(&indices).map_err(json_error)?);
            } else {
                for index in &indices {
                    println!(
                        "  ({:2} {:2} {:2}) d={:.4} multiplicity={}",
                        index.hkl[0],
                        index.hkl[1],
                        index.hkl[2],
                        index.d_spacing,
                        index.multiplicity
                    );
                }
            }
            Ok(())
        }
        FeatureAction::Terminations {
            path,
            hkl,
            layers,
            vacuum,
            symmetric,
            index,
            out,
            format,
        } => {
            let structure = load(path)?;
            let hkl = parse_hkl(hkl)?;
            let report = slab_terminations(&structure, hkl, *layers, *vacuum, *symmetric)
                .map_err(analysis_error)?;
            if cli.json {
                let summaries: Vec<serde_json::Value> = report
                    .terminations
                    .iter()
                    .map(|termination| {
                        serde_json::json!({
                            "index": termination.index,
                            "shift": termination.shift,
                            "n_atoms": termination.n_atoms,
                            "thickness": termination.thickness,
                            "dipole": termination.dipole,
                            "polar": termination.polar,
                        })
                    })
                    .collect();
                println!(
                    "{}",
                    serde_json::json!({
                        "hkl": report.hkl,
                        "n_terminations": report.n_terminations,
                        "terminations": summaries,
                    })
                );
            } else {
                println!(
                    "({} {} {}): {} terminations",
                    report.hkl[0], report.hkl[1], report.hkl[2], report.n_terminations
                );
                for termination in &report.terminations {
                    println!(
                        "  {} shift={:.4} atoms={} thickness={:.3} dipole={:.4} polar={}",
                        termination.index,
                        termination.shift,
                        termination.n_atoms,
                        termination.thickness,
                        termination.dipole,
                        termination.polar
                    );
                }
            }
            let selected = index.unwrap_or(0);
            if let Some(termination) = report.terminations.get(selected) {
                emit_structure(cli, &termination.structure, out.as_deref(), *format)?;
            } else if out.is_some() {
                return Err(CliError::input(format!(
                    "termination index {selected} is out of range"
                )));
            }
            Ok(())
        }
        FeatureAction::Wulff { path, surfaces } => {
            let structure = load(path)?;
            if surfaces.is_empty() {
                return Err(CliError::input("--surface needs four values: h,k,l,energy"));
            }
            let operations = detect_symmetry(
                &structure,
                SymmetryOptions {
                    symprec: 1e-4,
                    ..SymmetryOptions::default()
                },
            )
            .map(|dataset| dataset.operations)
            .unwrap_or_default();
            let mut directions: BTreeSet<[i32; 3]> = BTreeSet::new();
            let mut energies = Vec::new();
            for entry in surfaces {
                let parts: Vec<&str> = entry.split(',').collect();
                if parts.len() != 4 {
                    return Err(CliError::input("--surface needs four values: h,k,l,energy"));
                }
                let h = parse_index(parts.first(), "h")?;
                let k = parse_index(parts.get(1), "k")?;
                let l = parse_index(parts.get(2), "l")?;
                let energy: f64 = parts
                    .get(3)
                    .ok_or_else(|| CliError::input("missing surface energy"))?
                    .trim()
                    .parse()
                    .map_err(|_| CliError::input("surface energy must be a number"))?;
                let mut family = BTreeSet::new();
                family.insert([h, k, l]);
                family.insert([-h, -k, -l]);
                for operation in &operations {
                    let rotation = operation.rotation();
                    let transformed = [
                        rotation[0][0] * h + rotation[0][1] * k + rotation[0][2] * l,
                        rotation[1][0] * h + rotation[1][1] * k + rotation[1][2] * l,
                        rotation[2][0] * h + rotation[2][1] * k + rotation[2][2] * l,
                    ];
                    family.insert(transformed);
                    family.insert([-transformed[0], -transformed[1], -transformed[2]]);
                }
                for direction in family {
                    if directions.insert(direction) {
                        energies.push(SurfaceEnergy {
                            hkl: direction,
                            energy,
                        });
                    }
                }
            }
            let shape = wulff(structure.lattice(), &energies).map_err(analysis_error)?;
            if cli.json {
                println!("{}", serde_json::to_string(&shape).map_err(json_error)?);
            } else {
                println!(
                    "Wulff: volume={:.4} surface={:.4}",
                    shape.volume, shape.surface_area
                );
                for facet in &shape.facets {
                    println!(
                        "  ({:2} {:2} {:2}) area={:.4} fraction={:.4} energy={:.4}",
                        facet.hkl[0],
                        facet.hkl[1],
                        facet.hkl[2],
                        facet.area,
                        facet.fraction,
                        facet.energy
                    );
                }
            }
            Ok(())
        }
        FeatureAction::Adsorption {
            path,
            tolerance,
            molecule,
            offset,
            out,
            format,
        } => {
            let slab = load(path)?;
            let sites = adsorption_sites(&slab, *tolerance).map_err(analysis_error)?;
            if cli.json {
                println!("{}", serde_json::to_string(&sites).map_err(json_error)?);
            } else {
                println!("{} adsorption sites", sites.len());
                for (index, site) in sites.iter().enumerate() {
                    println!(
                        "  {} {:?} atoms={:?} position={:?}",
                        index, site.kind, site.atoms, site.position
                    );
                }
            }
            if let Some(molecule) = molecule {
                let molecule = load(molecule)?;
                let site = sites
                    .first()
                    .ok_or_else(|| CliError::input("the slab has no adsorption site"))?;
                let placed =
                    place_molecule(&slab, &molecule, site, *offset).map_err(analysis_error)?;
                emit_structure(cli, &placed, out.as_deref(), *format)?;
            } else if out.is_some() {
                emit_structure(cli, &slab, out.as_deref(), *format)?;
            }
            Ok(())
        }
        _ => Err(CliError::input("unsupported surface action")),
    }
}

fn parse_hkl(values: &str) -> Result<[i32; 3], CliError> {
    let parsed: Vec<i32> = values
        .split(',')
        .map(|value| {
            value
                .trim()
                .parse()
                .map_err(|_| CliError::input("--hkl needs three integers: h,k,l"))
        })
        .collect::<Result<_, _>>()?;
    match parsed.as_slice() {
        [h, k, l] => Ok([*h, *k, *l]),
        _ => Err(CliError::input("--hkl needs three integers: h,k,l")),
    }
}

fn parse_index(value: Option<&&str>, name: &str) -> Result<i32, CliError> {
    value
        .ok_or_else(|| CliError::input(format!("missing {name}")))?
        .trim()
        .parse()
        .map_err(|_| CliError::input(format!("{name} must be an integer")))
}

pub(crate) fn json_error(error: serde_json::Error) -> CliError {
    CliError::new("json_error", error.to_string())
}
