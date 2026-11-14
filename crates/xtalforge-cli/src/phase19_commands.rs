//! Phase 19 commands: magnetism, disorder, diffraction and spectroscopy.

use serde_json::json;

use crate::feature_args::FeatureAction;
use crate::{Cli, CliError, input};

pub(crate) fn feature(cli: &Cli, action: &FeatureAction) -> Result<(), CliError> {
    match action {
        FeatureAction::MagneticOrderings {
            path,
            supercell,
            max,
        } => {
            let structure = input::read_structure_path(path).map_err(CliError::core)?;
            let options = xtalforge::analysis::MagneticOrderingOptions {
                supercell: [*supercell, *supercell, *supercell],
                max_orderings: *max,
                ..xtalforge::analysis::MagneticOrderingOptions::default()
            };
            let orderings = xtalforge::analysis::magnetic_orderings(&structure, &options)
                .map_err(CliError::analysis)?;
            if cli.json {
                let values = orderings
                    .iter()
                    .map(|ordering| {
                        json!({
                            "kind": ordering.kind,
                            "moments": ordering.moments,
                            "net_moment": ordering.net_moment,
                            "uni_number": ordering.uni_number,
                            "operation_count": ordering.operation_count,
                            "magmom": ordering.magmom,
                        })
                    })
                    .collect::<Vec<_>>();
                println!("{}", json!(values));
            } else {
                for (index, ordering) in orderings.iter().enumerate() {
                    println!(
                        "{index}: {:?} net={:.3} operations={} MAGMOM={}",
                        ordering.kind,
                        ordering.net_moment,
                        ordering.operation_count,
                        ordering.magmom
                    );
                }
            }
            Ok(())
        }
        FeatureAction::Disorder { path } => {
            let structure = input::read_structure_path(path).map_err(CliError::core)?;
            let groups =
                xtalforge::analysis::disorder_groups(&structure).map_err(CliError::analysis)?;
            let composition = xtalforge::analysis::disorder_composition(&structure)
                .map_err(CliError::analysis)?;
            if cli.json {
                println!(
                    "{}",
                    json!({
                        "formula": composition.formula(),
                        "density": composition.density(structure.lattice().volume())
                            .map_err(CliError::analysis)?,
                        "groups": groups.iter().map(|group| json!({
                            "disorder_group": group.disorder_group,
                            "total_occupancy": group.total_occupancy(),
                            "positions": group.positions.iter().map(|position| json!({
                                "frac": position.frac,
                                "occupants": position.occupants.iter().map(|occupant| json!({
                                    "element": occupant.element,
                                    "occupancy": occupant.occupancy,
                                })).collect::<Vec<_>>(),
                            })).collect::<Vec<_>>(),
                        })).collect::<Vec<_>>(),
                    })
                );
            } else {
                println!("{} groups, formula {}", groups.len(), composition.formula());
                for group in &groups {
                    println!(
                        "  group {:?} total occupancy {:.3}",
                        group.disorder_group,
                        group.total_occupancy()
                    );
                }
            }
            Ok(())
        }
        FeatureAction::Derivatives {
            path,
            max_index,
            max,
        } => {
            let structure = input::read_structure_path(path).map_err(CliError::core)?;
            let options = xtalforge::analysis::DerivativeOptions {
                max_index: *max_index,
                max_structures: *max,
                ..xtalforge::analysis::DerivativeOptions::default()
            };
            let derivatives = xtalforge::analysis::derivative_structures(&structure, &options)
                .map_err(CliError::analysis)?;
            if cli.json {
                let values = derivatives
                    .iter()
                    .map(|derivative| {
                        json!({
                            "matrix": derivative.matrix,
                            "decoration": derivative.decoration,
                            "symmetry_count": derivative.symmetry_count,
                        })
                    })
                    .collect::<Vec<_>>();
                println!("{}", json!(values));
            } else {
                println!("{} derivative structures", derivatives.len());
                for derivative in &derivatives {
                    println!(
                        "  {:?} symmetry={} decoration={}",
                        derivative.matrix,
                        derivative.symmetry_count,
                        derivative.decoration.join("")
                    );
                }
            }
            Ok(())
        }
        _ => Err(CliError::input(
            "this action is not a Phase 19 feature".to_owned(),
        )),
    }
}

pub(crate) fn neutron(cli: &Cli, path: &str, wavelength: f64, step: f64) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let options = xtalforge::analysis::NeutronOptions {
        geometry: xtalforge::analysis::NeutronGeometry::ConstantWavelength { wavelength },
        step,
        ..xtalforge::analysis::NeutronOptions::default()
    };
    let pattern =
        xtalforge::analysis::neutron_powder(&structure, &options).map_err(CliError::analysis)?;
    if cli.json {
        println!("{}", serde_json::to_string(&pattern).unwrap_or_default());
    } else {
        println!(
            "{} reflections, {} points, max {:.3}",
            pattern.reflections.len(),
            pattern.points.len(),
            pattern.max_intensity()
        );
    }
    Ok(())
}

pub(crate) fn electron(
    cli: &Cli,
    path: &str,
    zone: &str,
    voltage: f64,
    image: Option<&str>,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let zone_axis = parse_zone(zone)?;
    let options = xtalforge::analysis::SaedOptions {
        acceleration_kv: voltage,
        zone_axis,
        image_size: 256,
        ..xtalforge::analysis::SaedOptions::default()
    };
    let pattern =
        xtalforge::analysis::zone_axis_pattern(&structure, &options).map_err(CliError::analysis)?;
    if let Some(destination) = image {
        let buffer = pattern.image().map_err(CliError::analysis)?;
        let mut text = format!("P5\n{} {}\n255\n", buffer.width, buffer.height).into_bytes();
        text.extend_from_slice(&buffer.to_u8());
        std::fs::write(destination, text)
            .map_err(|error| CliError::input(format!("cannot write {destination}: {error}")))?;
    }
    if cli.json {
        println!(
            "{}",
            json!({
                "spots": pattern.spots,
                "max_intensity": pattern.max_intensity(),
            })
        );
    } else {
        println!(
            "{} spots, max intensity {:.3}",
            pattern.spots.len(),
            pattern.max_intensity()
        );
    }
    Ok(())
}

pub(crate) fn raman(cli: &Cli, path: &str, tensors: &str, fwhm: f64) -> Result<(), CliError> {
    let (structure, modes) = modes_for(path)?;
    let text = std::fs::read_to_string(tensors)
        .map_err(|error| CliError::input(format!("cannot read {tensors}: {error}")))?;
    let tensors = xtalforge::ff::RamanTensorSet::from_text(&text).map_err(CliError::ff)?;
    let spectrum = xtalforge::ff::raman_spectrum(
        &modes,
        &tensors,
        &xtalforge::ff::RamanOptions {
            fwhm,
            ..xtalforge::ff::RamanOptions::default()
        },
    )
    .map_err(CliError::ff)?;
    let _ = structure;
    if cli.json {
        println!(
            "{}",
            json!({
                "peaks": spectrum.peaks.iter().map(|peak| json!({
                    "mode": peak.mode,
                    "frequency": peak.frequency,
                    "activity": peak.activity,
                })).collect::<Vec<_>>(),
                "points": spectrum.points,
            })
        );
    } else {
        for peak in &spectrum.peaks {
            println!(
                "mode {}: {:.3} cm^-1, activity {:.6}",
                peak.mode, peak.frequency, peak.activity
            );
        }
    }
    Ok(())
}

pub(crate) fn ir_born(cli: &Cli, path: &str, born_path: &str, fwhm: f64) -> Result<(), CliError> {
    let (structure, modes) = modes_for(path)?;
    let text = std::fs::read_to_string(born_path)
        .map_err(|error| CliError::input(format!("cannot read {born_path}: {error}")))?;
    let born = xtalforge::ff::BornCharges::from_text(&text).map_err(CliError::ff)?;
    let spectrum = xtalforge::ff::ir_spectrum_with_born(
        &structure,
        &modes,
        &born,
        &xtalforge::ff::IrOptions {
            fwhm,
            ..xtalforge::ff::IrOptions::default()
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

fn modes_for(
    path: &str,
) -> Result<(xtalforge::core::Structure, xtalforge::ff::VibrationalModes), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let p1 = xtalforge::symmetry::to_p1(&structure).map_err(CliError::core)?;
    let options = xtalforge::ff::HessianOptions::default();
    let matrix = xtalforge::ff::hessian(&crate::ff_commands::calculator()?, &p1, &options)
        .map_err(CliError::ff)?;
    let modes = xtalforge::ff::vibrational_modes(&matrix, &p1, &options).map_err(CliError::ff)?;
    Ok((p1, modes))
}

fn parse_zone(text: &str) -> Result<[i32; 3], CliError> {
    let values = text
        .split(',')
        .map(|value| {
            value
                .trim()
                .parse::<i32>()
                .map_err(|_| CliError::input(format!("bad zone axis {text:?}")))
        })
        .collect::<Result<Vec<_>, _>>()?;
    match values.as_slice() {
        [u, v, w] => Ok([*u, *v, *w]),
        _ => Err(CliError::input(format!(
            "zone axis needs three integers, got {text:?}"
        ))),
    }
}
