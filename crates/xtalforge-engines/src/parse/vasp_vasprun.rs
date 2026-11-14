//! VASP `vasprun.xml` result parser (Task 18.2.2).
//!
//! The structure, per-ionic-step energies, forces and stress come from the
//! Phase 13 reader in `xtalforge_io::vasp`, so the trajectory semantics stay
//! identical. This module adds the result-model fields that the tree reader
//! does not expose: magnetic moments, the Fermi level, bands, DOS and timings.
//! Those extra sections are streamed tag by tag with [`super::xml`]; only one
//! `<calculation>`, `<dos>` or `<eigenvalues>` section is materialised at a
//! time.

use xtalforge_io::trajectory::Trajectory;
use xtalforge_io::vasp::{read_vasprun, read_vasprun_structure};

use crate::Error;
use crate::result::{BandData, CalculationResult, DosChannel, DosData, ElectronicStep, Snapshot};

use super::common::{max_force, warn};
use super::xml::{self, Event, Node, Reader};

/// Parses a `vasprun.xml` into the common result model.
///
/// # Errors
///
/// Returns [`Error::Parse`] for malformed or truncated XML, a calculation
/// with an inconsistent atom count, or a missing cell.
pub fn parse(text: &str) -> Result<CalculationResult, Error> {
    let trajectory =
        read_vasprun(text).map_err(|error| Error::Parse(format!("vasprun.xml: {error}")))?;
    let structure = read_vasprun_structure(text)
        .map_err(|error| Error::Parse(format!("vasprun.xml: {error}")))?;
    let mut result = CalculationResult::new("vasp");
    result.structure = Some(structure);
    snapshots(&trajectory, &mut result)?;
    let finalpos = text.contains("name=\"finalpos\"");
    result.convergence.scf = !result.trajectory.is_empty();
    result.convergence.ionic = finalpos;
    result.convergence.scf_reported = false;
    result.convergence.ionic_reported = finalpos;
    warn(
        &mut result,
        "vasp_scf_convergence_inferred",
        "vasprun.xml does not print an explicit SCF flag; a completed calculation was treated as converged",
    );
    if !finalpos {
        warn(
            &mut result,
            "vasp_ionic_convergence_inferred",
            "vasprun.xml has no finalpos structure; the ionic state was inferred from the last calculation",
        );
    }
    extras(text, &mut result)?;
    result.close_trajectory();
    Ok(result)
}

fn snapshots(trajectory: &Trajectory, result: &mut CalculationResult) -> Result<(), Error> {
    for (index, frame) in trajectory.frames.iter().enumerate() {
        let structure = trajectory
            .structure(index)
            .map_err(|error| Error::Parse(format!("vasprun.xml step {index}: {error}")))?;
        result.trajectory.push(Snapshot {
            step: frame.step.unwrap_or(index as u64),
            structure,
            energy: frame.energy,
            forces: frame.forces.clone(),
            stress: frame.stress,
        });
        if let Some(energy) = frame.energy {
            result.energies.push(energy);
            result.steps.push(ElectronicStep {
                energy,
                fermi: None,
                max_force: frame.forces.as_deref().map(max_force),
            });
        }
    }
    Ok(())
}

fn extras(text: &str, result: &mut CalculationResult) -> Result<(), Error> {
    let mut reader = Reader::new(text);
    let mut magnetic = None;
    let mut efermi = None;
    let mut bands = None;
    let mut dos = None;
    let mut kpoints: Vec<[f64; 3]> = Vec::new();
    let mut weights: Vec<f64> = Vec::new();
    let mut timings: Vec<(String, f64)> = Vec::new();
    loop {
        match reader.next()? {
            Event::Eof => break,
            Event::Start { name, attrs, empty } => match name {
                "calculation" => {
                    let node = xml::read_start(&mut reader, name, attrs, empty)?;
                    if let Some(rows) = node.varray("magnetization") {
                        magnetic = Some(
                            rows.iter()
                                .filter_map(|row| row.first().copied())
                                .collect::<Vec<f64>>(),
                        );
                    }
                    if let Some(values) = node.find("eigenvalues") {
                        bands = parse_eigenvalues(values).or(bands);
                    }
                }
                "eigenvalues" => {
                    let node = xml::read_start(&mut reader, name, attrs, empty)?;
                    bands = parse_eigenvalues(&node);
                }
                "dos" => {
                    let node = xml::read_start(&mut reader, name, attrs, empty)?;
                    efermi = efermi_from(&node).or(efermi);
                    dos = parse_dos(&node, efermi);
                }
                "varray" => {
                    let label = attr(&attrs, "name").unwrap_or("").to_owned();
                    let node = xml::read_start(&mut reader, name, attrs, empty)?;
                    let rows = node
                        .children_named("v")
                        .map(|row| row.numbers())
                        .collect::<Vec<Vec<f64>>>();
                    match label.as_str() {
                        "kpointlist" => {
                            kpoints = rows
                                .iter()
                                .filter_map(|row| three(row))
                                .collect::<Vec<[f64; 3]>>();
                        }
                        "weights" => {
                            weights = rows.iter().filter_map(|row| row.first().copied()).collect();
                        }
                        _ => {}
                    }
                }
                "timing" => {
                    let node = xml::read_start(&mut reader, name, attrs, empty)?;
                    for label in ["total_cpu_time", "total_wall_time"] {
                        if let Some(value) = named_value(&node, label) {
                            timings.push((label.to_owned(), value));
                        }
                    }
                }
                _ => {}
            },
            Event::End { .. } | Event::Text(_) => {}
        }
    }
    if let Some(moments) = magnetic {
        let sites = result
            .structure
            .as_ref()
            .map_or(0, |value| value.sites().len());
        if sites == 0 || moments.len() == sites {
            result.magnetic_moments = moments;
        } else {
            warn(
                result,
                "vasp_magnetization_mismatch",
                "vasprun.xml magnetization rows do not match the number of sites",
            );
        }
    }
    if let Some(mut data) = bands {
        if data.energies.len() == kpoints.len() {
            data.kpoints = kpoints;
            if weights.len() == data.kpoints.len() {
                data.weights = weights;
            }
        }
        data.fermi = efermi;
        result.bands = Some(data);
    }
    if let Some(mut data) = dos {
        data.fermi = efermi;
        result.dos = Some(data);
    }
    result.fermi = efermi;
    for (stage, seconds) in timings {
        result.push_timing(&stage, seconds);
    }
    Ok(())
}

fn attr<'a>(attrs: &'a [(&str, String)], name: &str) -> Option<&'a str> {
    attrs
        .iter()
        .find(|(key, _)| *key == name)
        .map(|(_, value)| value.as_str())
}

fn three(row: &[f64]) -> Option<[f64; 3]> {
    Some([*row.first()?, *row.get(1)?, *row.get(2)?])
}

fn named_value(node: &Node, name: &str) -> Option<f64> {
    let mut nodes = Vec::new();
    node.find_all("i", &mut nodes);
    nodes
        .into_iter()
        .find(|child| child.attr("name") == Some(name))
        .and_then(Node::value)
}

fn efermi_from(node: &Node) -> Option<f64> {
    let mut nodes = Vec::new();
    node.find_all("i", &mut nodes);
    nodes
        .into_iter()
        .find(|child| child.attr("name") == Some("efermi"))
        .and_then(Node::value)
}

fn parse_eigenvalues(node: &Node) -> Option<BandData> {
    let array = node.child("array").unwrap_or(node);
    let sets: Vec<&Node> = array.children_named("set").collect();
    let spin_resolved = sets
        .first()
        .and_then(|set| set.attr("comment"))
        .is_some_and(|comment| comment.starts_with("spin"));
    let kpoint_sets: Vec<&Node> = if spin_resolved {
        sets.first()
            .map(|set| set.children_named("set").collect())
            .unwrap_or_default()
    } else {
        sets
    };
    let mut energies = Vec::new();
    let mut occupations = Vec::new();
    for set in kpoint_sets {
        let mut row_energies = Vec::new();
        let mut row_occupations = Vec::new();
        for row in set.children_named("r") {
            let values = row.numbers();
            if let (Some(energy), Some(occupation)) = (values.first(), values.get(1)) {
                row_energies.push(*energy);
                row_occupations.push(*occupation);
            }
        }
        if !row_energies.is_empty() {
            energies.push(row_energies);
            occupations.push(row_occupations);
        }
    }
    if energies.is_empty() {
        return None;
    }
    Some(BandData {
        energies,
        occupations,
        ..BandData::default()
    })
}

fn parse_dos(node: &Node, fermi: Option<f64>) -> Option<DosData> {
    let total = node.child("total")?;
    let array = total.child("array")?;
    let sets: Vec<&Node> = array.children_named("set").collect();
    let mut energies = Vec::new();
    let mut first = Vec::new();
    let mut first_integrated = Vec::new();
    let mut second = Vec::new();
    let mut second_integrated = Vec::new();
    for (set_index, set) in sets.iter().enumerate() {
        for row in set.children_named("r") {
            let values = row.numbers();
            let (Some(energy), Some(value), Some(integrated)) =
                (values.first(), values.get(1), values.get(2))
            else {
                continue;
            };
            if set_index == 0 {
                energies.push(*energy);
                first.push(*value);
                first_integrated.push(*integrated);
            } else if set_index == 1 {
                second.push(*value);
                second_integrated.push(*integrated);
            }
        }
    }
    if energies.is_empty() {
        return None;
    }
    let spin = second.len() == first.len();
    let total_values = if spin {
        first
            .iter()
            .zip(second.iter())
            .map(|(up, down)| up + down)
            .collect::<Vec<f64>>()
    } else {
        first.clone()
    };
    let integrated = if spin {
        first_integrated
            .iter()
            .zip(second_integrated.iter())
            .map(|(up, down)| up + down)
            .collect::<Vec<f64>>()
    } else {
        first_integrated
    };
    let mut channels = Vec::new();
    if spin {
        channels.push(DosChannel {
            label: "spin_up".to_owned(),
            values: first,
        });
        channels.push(DosChannel {
            label: "spin_down".to_owned(),
            values: second,
        });
    }
    Some(DosData {
        energies,
        total: total_values,
        integrated,
        channels,
        fermi,
    })
}
