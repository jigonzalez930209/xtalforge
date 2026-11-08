//! ORCA simple input and whole-molecule reconstruction (Task 4.9.4).

use std::collections::VecDeque;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Structure, element_info};

use super::text::{push_line, validate_token};
use crate::Error;

/// ORCA input options.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OrcaOptions {
    /// Functional or composite method, for example `B3LYP`.
    pub method: String,
    /// Basis set, for example `def2-SVP`.
    pub basis: String,
    /// Total charge.
    pub charge: i32,
    /// Spin multiplicity; one for a closed-shell system.
    pub multiplicity: u32,
    /// Optional implicit-solvent name written as `CPCM(<name>)`.
    pub solvent: Option<String>,
}

impl Default for OrcaOptions {
    fn default() -> Self {
        Self {
            method: "B3LYP".to_owned(),
            basis: "def2-SVP".to_owned(),
            charge: 0,
            multiplicity: 1,
            solvent: None,
        }
    }
}

impl OrcaOptions {
    /// Checks the ORCA options.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an unsafe method, basis or solvent
    /// token, or a multiplicity below one.
    pub fn validate(&self) -> Result<(), Error> {
        validate_token("ORCA method", &self.method)?;
        validate_token("ORCA basis", &self.basis)?;
        if let Some(solvent) = &self.solvent {
            validate_token("ORCA solvent", solvent)?;
        }
        if self.multiplicity == 0 {
            return Err(Error::InvalidInput(
                "ORCA multiplicity must be at least one".to_owned(),
            ));
        }
        Ok(())
    }
}

/// Functionals and bases offered by the ORCA input UI, as `(name, kind)`.
///
/// The list is a curated catalogue, not a claim that every pair is valid: the
/// caller combines one functional with one basis.
#[must_use]
pub fn orca_catalogue() -> Vec<(&'static str, &'static str)> {
    vec![
        ("B3LYP", "functional"),
        ("PBE0", "functional"),
        ("M06-2X", "functional"),
        ("wB97X-D3", "functional"),
        ("TPSS", "functional"),
        ("def2-SVP", "basis"),
        ("def2-TZVP", "basis"),
        ("def2-TZVPP", "basis"),
        ("cc-pVTZ", "basis"),
        ("6-31G*", "basis"),
    ]
}

/// Writes an ORCA simple-input file for `structure`.
///
/// Each connected component of the perceived bond graph is made whole by
/// translating bonded images, then moved as a rigid unit so its centroid lies
/// in the cell.  Every P1 atom is written once in the `* xyz` block.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, an empty structure, a
/// dummy atom or a failed bond perception, and a core error when the P1
/// expansion fails.
pub fn write_orca_input(structure: &Structure, options: &OrcaOptions) -> Result<String, Error> {
    options.validate()?;
    let cell = structure.p1()?;
    if cell.n_atoms() == 0 {
        return Err(Error::InvalidInput(
            "cannot write an ORCA input for an empty structure".to_owned(),
        ));
    }
    for atom in &cell.atoms {
        if element_info(&atom.element)?.is_dummy() {
            return Err(Error::InvalidInput(
                "external engines do not accept the dummy X species".to_owned(),
            ));
        }
    }
    let positions = whole_molecule_positions(structure)?;
    let mut simple = format!("! {} {}", options.method, options.basis);
    if let Some(solvent) = &options.solvent {
        simple.push_str(&format!(" CPCM({solvent})"));
    }
    let mut out = String::new();
    push_line(&mut out, simple);
    push_line(
        &mut out,
        format!("* xyz {} {}", options.charge, options.multiplicity),
    );
    for (atom, position) in cell.atoms.iter().zip(&positions) {
        push_line(
            &mut out,
            format!(
                "{:<2} {:.10} {:.10} {:.10}",
                atom.element, position[0], position[1], position[2]
            ),
        );
    }
    push_line(&mut out, "*");
    Ok(out)
}

/// Cartesian positions with every bonded component translated into one piece.
fn whole_molecule_positions(structure: &Structure) -> Result<Vec<[f64; 3]>, Error> {
    let cell = structure.p1()?;
    let graph = structure.bond_graph()?;
    let n = cell.n_atoms();
    if graph.n_atoms() != n {
        return Err(Error::Topology(
            "bond graph does not match the P1 atom count".to_owned(),
        ));
    }
    let mut offsets = vec![[0_i32; 3]; n];
    let mut visited = vec![false; n];
    let mut components = Vec::new();
    for start in 0..n {
        if visited.get(start).copied().unwrap_or(true) {
            continue;
        }
        let mut component = Vec::new();
        let mut queue = VecDeque::from([start]);
        if let Some(slot) = visited.get_mut(start) {
            *slot = true;
        }
        while let Some(index) = queue.pop_front() {
            component.push(index);
            for bond_index in graph.bonds_of(index) {
                let Some(bond) = graph.bonds.get(*bond_index) else {
                    continue;
                };
                let (other, delta) = if bond.i == index {
                    (bond.j, bond.image)
                } else {
                    (bond.i, [-bond.image[0], -bond.image[1], -bond.image[2]])
                };
                if visited.get(other).copied().unwrap_or(true) {
                    continue;
                }
                if let Some(slot) = visited.get_mut(other) {
                    *slot = true;
                }
                let base = offsets.get(index).copied().unwrap_or([0; 3]);
                if let Some(slot) = offsets.get_mut(other) {
                    *slot = [base[0] + delta[0], base[1] + delta[1], base[2] + delta[2]];
                }
                queue.push_back(other);
            }
        }
        components.push(component);
    }
    let mut frac: Vec<[f64; 3]> = cell
        .atoms
        .iter()
        .enumerate()
        .map(|(index, atom)| {
            let offset = offsets.get(index).copied().unwrap_or([0; 3]);
            [
                atom.frac[0] + f64::from(offset[0]),
                atom.frac[1] + f64::from(offset[1]),
                atom.frac[2] + f64::from(offset[2]),
            ]
        })
        .collect();
    for component in components {
        if component.is_empty() {
            continue;
        }
        let mut mean = [0.0; 3];
        for index in &component {
            if let Some(point) = frac.get(*index) {
                mean[0] += point[0];
                mean[1] += point[1];
                mean[2] += point[2];
            }
        }
        let count = component.len() as f64;
        let shift = [
            (mean[0] / count).floor(),
            (mean[1] / count).floor(),
            (mean[2] / count).floor(),
        ];
        for index in &component {
            if let Some(point) = frac.get_mut(*index) {
                point[0] -= shift[0];
                point[1] -= shift[1];
                point[2] -= shift[2];
            }
        }
    }
    Ok(cell.lattice.to_cart_many(&frac))
}
