//! Phase 9 interface, elastic and thermodynamic commands.

use xtalforge::analysis::{InterfaceOptions, lattice_matches, stack_slabs};
use xtalforge::core::{EosKind, EosPoint, PhaseEntry, convex_hull, fit_eos};
use xtalforge::ff::{CalculatorOptions, ElasticOptions, UffCalculator};

use crate::feature_args::FeatureAction;
use crate::feature_commands::{
    analysis_error, core_error, emit_structure, ff_error, load, parse_numbers,
};
use crate::{Cli, CliError, input};

pub(crate) fn run(cli: &Cli, action: &FeatureAction) -> Result<(), CliError> {
    match action {
        FeatureAction::LatticeMatch {
            film,
            substrate,
            max_index,
            max_area,
            tolerance,
        } => {
            let matches = lattice_matches(
                &load(film)?,
                &load(substrate)?,
                &InterfaceOptions {
                    max_index: *max_index,
                    max_area: *max_area,
                    length_tolerance: *tolerance,
                    ..InterfaceOptions::default()
                },
            )
            .map_err(analysis_error)?;
            if cli.json {
                println!("{}", serde_json::to_string(&matches).map_err(json_error)?);
            } else {
                println!("{} lattice matches", matches.len());
                for entry in matches.iter().take(10) {
                    println!(
                        "  strain={:.4} angle={:.3} area={:.3} film={:?} substrate={:?}",
                        entry.strain,
                        entry.angle_mismatch,
                        entry.area,
                        entry.film_matrix,
                        entry.substrate_matrix
                    );
                }
            }
            Ok(())
        }
        FeatureAction::Stack {
            film,
            substrate,
            gap,
            vacuum,
            shift,
            out,
            format,
        } => {
            let values = parse_numbers(shift)?;
            let shift = match values.as_slice() {
                [u, v] => [*u, *v, 0.0],
                [u, v, w] => [*u, *v, *w],
                _ => return Err(CliError::input("--shift needs two or three numbers")),
            };
            let result = stack_slabs(&load(film)?, &load(substrate)?, *gap, shift, *vacuum)
                .map_err(analysis_error)?;
            if cli.json {
                println!(
                    "{}",
                    serde_json::json!({
                        "strain": result.strain,
                        "gap": result.gap,
                        "vacuum": result.vacuum,
                        "n_atoms": result.structure.sites().len(),
                    })
                );
            } else {
                println!(
                    "stacked interface: strain={:.4} gap={} vacuum={}",
                    result.strain, result.gap, result.vacuum
                );
            }
            emit_structure(cli, &result.structure, out.as_deref(), *format)
        }
        FeatureAction::Elastic { path, strain } => {
            let structure = load(path)?;
            let calculator = UffCalculator::new(CalculatorOptions::default()).map_err(ff_error)?;
            let tensor = xtalforge::ff::elastic_tensor(
                &structure,
                &calculator,
                &ElasticOptions { strain: *strain },
            )
            .map_err(core_error)?;
            if cli.json {
                println!("{}", serde_json::to_string(&tensor).map_err(json_error)?);
            } else {
                println!(
                    "bulk={:.3} GPa shear={:.3} GPa young={:.3} GPa poisson={:.4} stable={}",
                    tensor.bulk_modulus_hill,
                    tensor.shear_modulus_hill,
                    tensor.young_modulus,
                    tensor.poisson_ratio,
                    tensor.stable
                );
                for row in &tensor.stiffness {
                    println!(
                        "{:10.4} {:10.4} {:10.4} {:10.4} {:10.4} {:10.4}",
                        row[0], row[1], row[2], row[3], row[4], row[5]
                    );
                }
            }
            Ok(())
        }
        FeatureAction::Eos { points, model } => {
            let mut data = Vec::new();
            for point in points {
                let values = parse_numbers(point)?;
                let [volume, energy] = match values.as_slice() {
                    [volume, energy] => [*volume, *energy],
                    _ => return Err(CliError::input("--point needs volume,energy")),
                };
                data.push(EosPoint { volume, energy });
            }
            let kind = match model.as_str() {
                "birch" | "birch_murnaghan" => EosKind::BirchMurnaghan,
                "vinet" => EosKind::Vinet,
                "murnaghan" => EosKind::Murnaghan,
                other => {
                    return Err(CliError::input(format!("unknown EOS model '{other}'")));
                }
            };
            let fit = fit_eos(&data, kind).map_err(core_error)?;
            if cli.json {
                println!("{}", serde_json::to_string(&fit).map_err(json_error)?);
            } else {
                println!(
                    "V0={:.5} A^3 E0={:.6} eV B0={:.3} GPa B0'={:.4} rms={:.6}",
                    fit.v0, fit.e0, fit.b0, fit.b0_prime, fit.rms
                );
                println!(
                    "errors: V0={:.5} E0={:.6} B0={:.3} B0'={:.4}",
                    fit.errors[0], fit.errors[1], fit.errors[2], fit.errors[3]
                );
            }
            Ok(())
        }
        FeatureAction::Hull { path } => {
            let text = input::read_source(path).map_err(CliError::core)?;
            let entries: Vec<PhaseEntry> = serde_json::from_str(&text).map_err(json_error)?;
            let hull = convex_hull(&entries, None).map_err(core_error)?;
            if cli.json {
                println!("{}", serde_json::to_string(&hull).map_err(json_error)?);
            } else {
                println!(
                    "{} phases, {} elements",
                    hull.points.len(),
                    hull.elements.len()
                );
                for point in &hull.points {
                    println!(
                        "  {} Ef={:.6} hull={:.6} distance={:.6} stable={}",
                        point.name,
                        point.formation_energy,
                        point.hull_energy,
                        point.distance,
                        point.stable
                    );
                }
            }
            Ok(())
        }
        _ => Err(CliError::input("unsupported mechanics action")),
    }
}

pub(crate) fn json_error(error: serde_json::Error) -> CliError {
    CliError::new("json_error", error.to_string())
}
