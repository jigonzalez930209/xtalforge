//! VASP `vasprun.xml` reader (Task 13.2.4).
//!
//! Every `<calculation>` becomes one [`Frame`] with the cell and fractional
//! positions of its `<structure>`, the `e_fr_energy` value in eV, the forces
//! in eV/Angstrom and the stress in kBar converted to GPa. The structure
//! reader prefers the `finalpos` structure over the last calculation.

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked

use std::collections::BTreeMap;

use xtalforge_core::elements::parse_symbol;
use xtalforge_core::{Error, Lattice, Structure};

use super::super::trajectory::{Frame, Trajectory};
use super::xml::{self, Node};

type StructureData = ([[f64; 3]; 3], Vec<[f64; 3]>);

/// Reads a `vasprun.xml` trajectory.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for malformed XML, a missing atom list, a
/// missing cell or a calculation with incomparable atom counts.
pub fn read_vasprun(text: &str) -> Result<Trajectory, Error> {
    let document = xml::parse(text)?;
    let root = document
        .children
        .first()
        .ok_or_else(|| Error::InvalidInput("vasprun.xml has no root element".to_owned()))?;
    let elements = atom_symbols(&document)?;
    let generator = generator_info(root);
    let mut calculations = Vec::new();
    document.find_all("calculation", &mut calculations);
    let mut frames = Vec::new();
    for calculation in calculations {
        let Some((matrix, positions)) = structure_of(calculation) else {
            continue;
        };
        Lattice::new(matrix)?;
        if positions.len() != elements.len() {
            return Err(Error::InvalidInput(format!(
                "vasprun.xml calculation has {} positions for {} atoms",
                positions.len(),
                elements.len()
            )));
        }
        let energy = calculation
            .child("energy")
            .and_then(energy_value)
            .or_else(|| calculation.find("energy").and_then(energy_value));
        let forces = find_varray(calculation, "forces")
            .or_else(|| calculation.varray("forces"))
            .and_then(|rows| rows_to_forces(&rows));
        let stress = find_varray(calculation, "stress")
            .or_else(|| calculation.varray("stress"))
            .and_then(|rows| rows_to_stress(&rows));
        let mut info = BTreeMap::new();
        info.insert("source".to_owned(), "vasprun".to_owned());
        info.extend(generator.clone());
        frames.push(Frame {
            frac: positions,
            lattice: Some(matrix),
            energy,
            step: Some(frames.len() as u64),
            info,
            forces,
            stress,
        });
    }
    if frames.is_empty() {
        if let Some((matrix, positions)) = final_structure(&document)? {
            Lattice::new(matrix)?;
            let mut info = BTreeMap::new();
            info.insert("source".to_owned(), "vasprun".to_owned());
            info.extend(generator);
            return Ok(Trajectory {
                elements,
                lattice: Some(matrix),
                frames: vec![Frame {
                    frac: positions.clone(),
                    lattice: None,
                    energy: None,
                    step: Some(0),
                    info,
                    forces: None,
                    stress: None,
                }],
            });
        }
        return Err(Error::InvalidInput(
            "vasprun.xml has no calculation or final structure".to_owned(),
        ));
    }
    let first = frames
        .first()
        .and_then(|frame| frame.lattice)
        .ok_or_else(|| Error::InvalidInput("vasprun.xml has no cell".to_owned()))?;
    let fixed = frames.iter().all(|frame| frame.lattice == Some(first));
    if fixed {
        for frame in &mut frames {
            frame.lattice = None;
        }
    }
    Ok(Trajectory {
        elements,
        lattice: fixed.then_some(first),
        frames,
    })
}

/// Reads the final structure of a `vasprun.xml`.
///
/// # Errors
///
/// Returns the errors of [`read_vasprun`] or of the frame-to-structure
/// conversion.
pub fn read_vasprun_structure(text: &str) -> Result<Structure, Error> {
    let document = xml::parse(text)?;
    if let Some((matrix, positions)) = final_structure(&document)? {
        let lattice = Lattice::new(matrix)?;
        let elements = atom_symbols(&document)?;
        if positions.len() != elements.len() {
            return Err(Error::InvalidInput(format!(
                "vasprun.xml final structure has {} positions for {} atoms",
                positions.len(),
                elements.len()
            )));
        }
        let mut builder = xtalforge_core::StructureBuilder::new("vasprun", lattice);
        for (index, (element, point)) in elements.iter().zip(positions.iter()).enumerate() {
            builder = builder.site(xtalforge_core::Site::new(
                format!("{element}{}", index + 1),
                element,
                *point,
            )?);
        }
        return builder.build();
    }
    let trajectory = read_vasprun(text)?;
    let last = trajectory
        .frames
        .len()
        .checked_sub(1)
        .ok_or_else(|| Error::InvalidInput("vasprun.xml has no frames".to_owned()))?;
    trajectory.structure(last)
}

fn atom_symbols(document: &Node) -> Result<Vec<String>, Error> {
    let mut arrays = Vec::new();
    document.find_all("array", &mut arrays);
    let array = arrays
        .into_iter()
        .find(|node| node.attr("name") == Some("atoms"))
        .ok_or_else(|| Error::InvalidInput("vasprun.xml has no atominfo atoms array".to_owned()))?;
    let mut symbols = Vec::new();
    for rc in array
        .children_named("set")
        .flat_map(|set| set.children_named("rc"))
    {
        let symbol = rc
            .children_named("c")
            .next()
            .map(Node::trimmed)
            .ok_or_else(|| Error::InvalidInput("vasprun.xml atom row has no symbol".to_owned()))?;
        symbols.push(parse_symbol(symbol)?);
    }
    if symbols.is_empty() {
        return Err(Error::InvalidInput(
            "vasprun.xml atom array is empty".to_owned(),
        ));
    }
    Ok(symbols)
}

fn generator_info(root: &Node) -> BTreeMap<String, String> {
    let mut info = BTreeMap::new();
    if let Some(generator) = root.find("generator") {
        for (key, name) in [("program", "program"), ("version", "version")] {
            if let Some(value) = generator
                .children_named("i")
                .find(|node| node.attr("name") == Some(name))
            {
                info.insert(key.to_owned(), value.trimmed().to_owned());
            }
        }
    }
    info
}

fn structure_of(node: &Node) -> Option<StructureData> {
    let structure = node
        .children_named("structure")
        .next()
        .or_else(|| node.find("structure"))?;
    structure_data(structure)
}

fn structure_data(structure: &Node) -> Option<StructureData> {
    let crystal = structure.child("crystal")?;
    let basis = crystal.varray("basis")?;
    let positions = structure.varray("positions")?;
    if basis.len() < 3 || positions.is_empty() {
        return None;
    }
    let matrix = [
        [
            *basis.first()?.first()?,
            *basis.first()?.get(1)?,
            *basis.first()?.get(2)?,
        ],
        [
            *basis.get(1)?.first()?,
            *basis.get(1)?.get(1)?,
            *basis.get(1)?.get(2)?,
        ],
        [
            *basis.get(2)?.first()?,
            *basis.get(2)?.get(1)?,
            *basis.get(2)?.get(2)?,
        ],
    ];
    let mut points = Vec::with_capacity(positions.len());
    for row in &positions {
        if row.len() < 3 {
            return None;
        }
        points.push([
            *row.first().unwrap_or(&0.0),
            *row.get(1).unwrap_or(&0.0),
            *row.get(2).unwrap_or(&0.0),
        ]);
    }
    Some((matrix, points))
}

fn final_structure(document: &Node) -> Result<Option<StructureData>, Error> {
    let mut structures = Vec::new();
    document.find_all("structure", &mut structures);
    let mut found = structures
        .iter()
        .filter(|node| node.attr("name") == Some("finalpos"))
        .find_map(|node| structure_data(node));
    if found.is_none() {
        found = structures
            .iter()
            .rev()
            .find_map(|node| structure_data(node));
    }
    Ok(found)
}

fn energy_value(energy: &Node) -> Option<f64> {
    ["e_fr_energy", "e_0_energy", "e_wo_entrp"]
        .iter()
        .find_map(|name| {
            energy
                .children_named("i")
                .find(|node| node.attr("name") == Some(name))
                .and_then(Node::number)
        })
}

fn rows_to_forces(rows: &[Vec<f64>]) -> Option<Vec<[f64; 3]>> {
    let mut forces = Vec::with_capacity(rows.len());
    for row in rows {
        if row.len() < 3 {
            return None;
        }
        forces.push([row[0], row[1], row[2]]);
    }
    Some(forces)
}

fn find_varray(node: &Node, name: &str) -> Option<Vec<Vec<f64>>> {
    if node.name == "varray" && node.attr("name") == Some(name) {
        return Some(
            node.children_named("v")
                .map(|row| {
                    row.trimmed()
                        .split_whitespace()
                        .filter_map(|value| value.parse().ok())
                        .collect()
                })
                .collect(),
        );
    }
    node.children
        .iter()
        .find_map(|child| find_varray(child, name))
}

fn rows_to_stress(rows: &[Vec<f64>]) -> Option<[[f64; 3]; 3]> {
    if rows.len() == 1 {
        let row = rows.first()?;
        if row.len() >= 6 {
            let value = |index: usize| row.get(index).copied().unwrap_or(0.0) * 0.1;
            return Some([
                [value(0), value(3), value(5)],
                [value(3), value(1), value(4)],
                [value(5), value(4), value(2)],
            ]);
        }
    }
    if rows.len() < 3 {
        return None;
    }
    let kbar = 0.1;
    let value = |row: usize, column: usize| -> f64 {
        rows.get(row)
            .and_then(|values| values.get(column))
            .copied()
            .unwrap_or(0.0)
            * kbar
    };
    Some([
        [value(0, 0), value(0, 1), value(0, 2)],
        [value(1, 0), value(1, 1), value(1, 2)],
        [value(2, 0), value(2, 1), value(2, 2)],
    ])
}
