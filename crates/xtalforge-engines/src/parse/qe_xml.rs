//! Quantum ESPRESSO `data-file-schema.xml` parser (Task 18.2.3).
//!
//! QE stores the total energy and eigenvalues in Hartree, the Fermi level in
//! eV, forces in Hartree/bohr and stress in the unit named by the `unit`
//! attribute (normally kbar). All values are converted to the common model
//! units and documented.

use xtalforge_core::{Lattice, Site, StructureBuilder};

use crate::Error;
use crate::result::{BandData, CalculationResult, ElectronicStep, Snapshot};

use super::common::{
    BOHR_ANGSTROM, HARTREE_BOHR_EV_ANGSTROM, HARTREE_EV, KBAR_GPA, flag, stress_from_rows, warn,
};
use super::xml::{self, Node, Reader};

/// Parses a QE `data-file-schema.xml` into the common result model.
///
/// # Errors
///
/// Returns [`Error::Parse`] for malformed XML, a missing `<output>`, a missing
/// total energy or an inconsistent atomic structure.
pub fn parse(text: &str) -> Result<CalculationResult, Error> {
    let mut reader = Reader::new(text);
    let root = xml::read_root(&mut reader)?;
    let output = root
        .find("output")
        .ok_or_else(|| Error::Parse("QE XML has no <output> element".to_owned()))?;
    let mut result = CalculationResult::new("quantum_espresso");
    if let Some(value) = output
        .find("convergence_achieved")
        .map(Node::trimmed)
        .and_then(flag)
    {
        result.convergence.scf = value;
        result.convergence.scf_reported = true;
    }
    if let Some(steps) = output.find("n_scf_steps").and_then(Node::value) {
        result.scalars.insert("n_scf_steps".to_owned(), steps);
    }
    let energy = output
        .find("total_energy")
        .and_then(|node| node.child("etot").or_else(|| node.child("etot_ry")))
        .and_then(Node::value)
        .ok_or_else(|| Error::Parse("QE XML has no total energy".to_owned()))?;
    let energy_ev = energy * HARTREE_EV;
    result.energies.push(energy_ev);
    result.fermi = output
        .find("fermi_energy")
        .or_else(|| output.find("highestOccupiedLevel"))
        .and_then(Node::value);
    result.bands = bands(output, result.fermi);
    if let Some(structure) = structure(output)? {
        result.trajectory.push(Snapshot {
            step: 0,
            structure: structure.clone(),
            energy: Some(energy_ev),
            forces: None,
            stress: None,
        });
        result.structure = Some(structure);
    }
    result.forces = forces(output);
    result.stress = stress(output)?;
    if let Some(last) = result.trajectory.last_mut() {
        last.forces = result.forces.clone();
        last.stress = result.stress;
    }
    result.steps.push(ElectronicStep {
        energy: energy_ev,
        fermi: result.fermi,
        max_force: result.forces.as_deref().map(super::common::max_force),
    });
    if let Some(timing) = output
        .find("timing_info")
        .and_then(|node| {
            node.child("total_seconds")
                .or_else(|| node.find("total_seconds"))
        })
        .and_then(Node::value)
    {
        result.push_timing("total", timing);
    }
    if !result.convergence.scf_reported {
        warn(
            &mut result,
            "qe_convergence_missing",
            "QE XML has no convergence_achieved flag",
        );
    }
    Ok(result)
}

fn bands(output: &Node, fermi: Option<f64>) -> Option<BandData> {
    let band_structure = output.find("band_structure")?;
    let mut kpoints = Vec::new();
    let mut weights = Vec::new();
    let mut energies = Vec::new();
    let mut occupations = Vec::new();
    let mut entries = Vec::new();
    band_structure.find_all("ks_energies", &mut entries);
    for entry in entries {
        let Some(point) = entry.child("k_point") else {
            continue;
        };
        let coordinates = point.numbers();
        let Some(point_values) = coordinates.get(..3) else {
            continue;
        };
        kpoints.push([
            point_values.first().copied().unwrap_or(0.0),
            point_values.get(1).copied().unwrap_or(0.0),
            point_values.get(2).copied().unwrap_or(0.0),
        ]);
        weights.push(
            point
                .attr("weight")
                .and_then(crate::parse::common::number)
                .unwrap_or(1.0),
        );
        let row = entry
            .child("eigenvalues")
            .map(|node| {
                node.numbers()
                    .iter()
                    .map(|value| value * HARTREE_EV)
                    .collect()
            })
            .unwrap_or_default();
        let occupation = entry
            .child("occupations")
            .map(Node::numbers)
            .unwrap_or_default();
        energies.push(row);
        occupations.push(occupation);
    }
    if energies.is_empty() {
        return None;
    }
    Some(BandData {
        kpoints,
        weights,
        energies,
        occupations,
        fermi,
    })
}

fn structure(output: &Node) -> Result<Option<xtalforge_core::Structure>, Error> {
    let Some(atomic) = output.find("atomic_structure") else {
        return Ok(None);
    };
    let cell = atomic
        .child("cell")
        .ok_or_else(|| Error::Parse("QE XML atomic structure has no cell".to_owned()))?;
    let scale = atomic
        .attr("alat")
        .and_then(crate::parse::common::number)
        .unwrap_or(1.0)
        * BOHR_ANGSTROM;
    let mut matrix = [[0.0; 3]; 3];
    for (index, name) in ["a1", "a2", "a3"].iter().enumerate() {
        let row = cell
            .child(name)
            .ok_or_else(|| Error::Parse(format!("QE XML cell has no {name} vector")))?
            .numbers();
        let target = matrix
            .get_mut(index)
            .ok_or_else(|| Error::Parse("QE XML cell has more than three vectors".to_owned()))?;
        for (column, slot) in target.iter_mut().enumerate() {
            *slot = row.get(column).copied().unwrap_or(0.0) * scale;
        }
    }
    let lattice =
        Lattice::new(matrix).map_err(|error| Error::Parse(format!("QE XML cell: {error}")))?;
    let positions = atomic.child("atomic_positions");
    let mut symbols = Vec::new();
    let mut cart = Vec::new();
    if let Some(positions) = positions {
        for atom in positions.children_named("atom") {
            let symbol = atom
                .attr("name")
                .ok_or_else(|| Error::Parse("QE XML atom has no name".to_owned()))?;
            let values = atom.numbers();
            symbols.push(symbol.to_owned());
            cart.push([
                values.first().copied().unwrap_or(0.0) * scale,
                values.get(1).copied().unwrap_or(0.0) * scale,
                values.get(2).copied().unwrap_or(0.0) * scale,
            ]);
        }
    }
    if symbols.is_empty() {
        return Err(Error::Parse(
            "QE XML atomic structure has no positions".to_owned(),
        ));
    }
    let nat = atomic
        .attr("nat")
        .and_then(|value| value.parse::<usize>().ok());
    if let Some(count) = nat
        && count != symbols.len()
    {
        return Err(Error::Parse(format!(
            "QE XML declares {count} atoms but prints {}",
            symbols.len()
        )));
    }
    let frac = cart
        .iter()
        .map(|point| lattice.to_frac(*point))
        .collect::<Vec<_>>();
    let mut builder = StructureBuilder::new("qe", lattice);
    for (index, (symbol, point)) in symbols.iter().zip(frac).enumerate() {
        builder = builder.site(Site::new(format!("{symbol}{}", index + 1), symbol, point)?);
    }
    builder
        .build()
        .map(Some)
        .map_err(|error| Error::Parse(format!("QE XML structure: {error}")))
}

fn forces(output: &Node) -> Option<Vec<[f64; 3]>> {
    let node = output.find("forces")?;
    let mut rows = Vec::new();
    for (index, force) in node.children_named("force").enumerate() {
        let values = force.numbers();
        if let (Some(x), Some(y), Some(z)) = (values.first(), values.get(1), values.get(2)) {
            let order = force
                .attr("index")
                .and_then(|value| value.parse::<usize>().ok())
                .unwrap_or(index + 1);
            rows.push((
                order,
                [
                    x * HARTREE_BOHR_EV_ANGSTROM,
                    y * HARTREE_BOHR_EV_ANGSTROM,
                    z * HARTREE_BOHR_EV_ANGSTROM,
                ],
            ));
        }
    }
    if rows.is_empty() {
        return None;
    }
    rows.sort_by_key(|(order, _)| *order);
    Some(rows.into_iter().map(|(_, force)| force).collect())
}

fn stress(output: &Node) -> Result<Option<[[f64; 3]; 3]>, Error> {
    let Some(node) = output.find("stress") else {
        return Ok(None);
    };
    let unit = node.attr("unit").unwrap_or("kbar").to_ascii_lowercase();
    let factor = match unit.as_str() {
        "kbar" => KBAR_GPA,
        "gpa" => 1.0,
        "ry/bohr^3" | "ry/bohr3" => super::common::RYDBERG_BOHR3_GPA,
        "hartree/bohr^3" | "hartree/bohr3" => 2.0 * super::common::RYDBERG_BOHR3_GPA,
        other => {
            return Err(Error::Parse(format!(
                "QE XML stress unit {other:?} is not supported"
            )));
        }
    };
    let rows = node
        .children
        .iter()
        .map(Node::numbers)
        .filter(|row| row.len() >= 3)
        .take(3)
        .collect::<Vec<Vec<f64>>>();
    Ok(stress_from_rows(&rows, factor))
}
