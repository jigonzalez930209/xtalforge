//! Symmetry workflow: detect, asymmetrize, reduce to P1 and set a group
//! (Task 2.3).

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Lattice, P1Cell, Site, Structure, StructureBuilder};

use crate::dataset::{SymmetryDataset, SymmetryOptions, detect_symmetry};
use crate::spacegroup::SpaceGroup;

/// How a new space group is applied to the existing sites.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SetMode {
    /// Keep the sites as the asymmetric unit of the new group.
    Reinterpret,
    /// Average the full cell onto the orbits of the new group.
    Impose,
}

/// Result of setting a space group.
#[derive(Debug, Clone)]
pub struct SetReport {
    /// The new structure.
    pub structure: Structure,
    /// Largest Cartesian displacement applied by [`SetMode::Impose`].
    pub max_displacement: f64,
    /// Number of sites before the operation.
    pub n_input: usize,
    /// Number of sites after the operation.
    pub n_output: usize,
}

/// Detection plus the symmetric structure in the detected setting.
#[derive(Debug, Clone)]
pub struct SymmetryResult {
    /// Asymmetric-unit structure carrying the detected operations.
    pub structure: Structure,
    /// The full dataset.
    pub dataset: SymmetryDataset,
}

/// Detects symmetry and returns the asymmetric-unit structure.
///
/// # Errors
///
/// Returns the errors of [`detect_symmetry`] or a builder error.
pub fn detect_and_symmetrize(
    structure: &Structure,
    options: SymmetryOptions,
) -> Result<SymmetryResult, Error> {
    let dataset = detect_symmetry(structure, options)?;
    let space_group = SpaceGroup::from_hall_number(dataset.hall_number)?;
    let asymmetric = asymmetric_structure(structure, &dataset, &space_group)?;
    Ok(SymmetryResult {
        structure: asymmetric,
        dataset,
    })
}

/// Reduces a structure to a representative asymmetric unit.
///
/// # Errors
///
/// Returns the errors of [`detect_and_symmetrize`].
pub fn asymmetrize(structure: &Structure, options: SymmetryOptions) -> Result<Structure, Error> {
    Ok(detect_and_symmetrize(structure, options)?.structure)
}

/// Expands a structure to an explicit P1 structure.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty structure.
pub fn to_p1(structure: &Structure) -> Result<Structure, Error> {
    let p1 = structure.p1()?;
    p1_structure(structure, &p1)
}

/// Applies a space group to a structure.
///
/// # Errors
///
/// Returns [`Error::Symmetry`] for an invalid mode result, or a builder error.
pub fn set_space_group(
    structure: &Structure,
    space_group: &SpaceGroup,
    mode: SetMode,
) -> Result<SetReport, Error> {
    let report = match mode {
        SetMode::Reinterpret => {
            let n_input = structure.sites().len();
            let mut builder = StructureBuilder::new(structure.name(), *structure.lattice())
                .hall(space_group.hall_symbol.clone());
            for operation in &space_group.operations {
                builder = builder.operation(*operation);
            }
            for site in structure.sites() {
                builder = builder.site(site.clone());
            }
            for (key, value) in structure.metadata() {
                builder = builder.metadata(key.clone(), value.clone());
            }
            let built = builder.build()?;
            Ok(SetReport {
                structure: built,
                max_displacement: 0.0,
                n_input,
                n_output: structure.sites().len(),
            })
        }
        SetMode::Impose => impose(structure, space_group),
    }?;
    validate(&report.structure, space_group)?;
    Ok(report)
}

fn validate(structure: &Structure, space_group: &SpaceGroup) -> Result<(), Error> {
    let dataset = detect_symmetry(
        structure,
        SymmetryOptions {
            symprec: 1e-3,
            ..SymmetryOptions::default()
        },
    )?;
    if dataset.number != space_group.number {
        return Err(Error::Symmetry(format!(
            "the result is space group {} but {} was requested",
            dataset.number, space_group.number
        )));
    }
    Ok(())
}

fn impose(structure: &Structure, space_group: &SpaceGroup) -> Result<SetReport, Error> {
    let p1 = structure.p1()?;
    let n_input = p1.atoms.len();
    if space_group.operations.is_empty() {
        return Err(Error::Symmetry(
            "an imposed group needs at least one operation".to_owned(),
        ));
    }
    let tolerance = 0.5;
    let mut assigned = vec![false; n_input];
    let mut averaged: Vec<([f64; 3], String, String, f64)> = Vec::new();
    let mut max_displacement: f64 = 0.0;
    for (index, atom) in p1.atoms.iter().enumerate() {
        if assigned.get(index).copied().unwrap_or(true) {
            continue;
        }
        if let Some(slot) = assigned.get_mut(index) {
            *slot = true;
        }
        let mut shift = [0.0; 3];
        let mut matched = 0u32;
        for operation in &space_group.operations {
            let image = operation.apply(atom.frac);
            let mut best: Option<(f64, usize, [f64; 3])> = None;
            for (other_index, other) in p1.atoms.iter().enumerate() {
                if other.element != atom.element {
                    continue;
                }
                let delta = p1.lattice.minimum_image([
                    other.frac[0] - image[0],
                    other.frac[1] - image[1],
                    other.frac[2] - image[2],
                ]);
                if delta.1 <= tolerance && best.is_none_or(|(distance, _, _)| delta.1 < distance) {
                    best = Some((delta.1, other_index, delta.0));
                }
            }
            if let Some((_, other_index, delta)) = best {
                shift[0] += delta[0];
                shift[1] += delta[1];
                shift[2] += delta[2];
                matched += 1;
                if let Some(slot) = assigned.get_mut(other_index) {
                    *slot = true;
                }
            }
        }
        let count = f64::from(matched.max(1));
        let position = Lattice::wrap([
            atom.frac[0] + shift[0] / count,
            atom.frac[1] + shift[1] / count,
            atom.frac[2] + shift[2] / count,
        ]);
        let displacement = p1.lattice.distance(atom.frac, position);
        max_displacement = max_displacement.max(displacement);
        let duplicate = averaged.iter().any(|(existing, element, _, _)| {
            element == &atom.element && p1.lattice.distance(*existing, position) <= 1e-5
        });
        if !duplicate {
            averaged.push((
                position,
                atom.element.clone(),
                atom.label.clone(),
                atom.occupancy,
            ));
        }
    }
    let mut builder =
        StructureBuilder::new(structure.name(), p1.lattice).hall(space_group.hall_symbol.clone());
    for operation in &space_group.operations {
        builder = builder.operation(*operation);
    }
    let n_output = averaged.len();
    for (index, (frac, element, label, occupancy)) in averaged.into_iter().enumerate() {
        let site = Site::new(
            if label.is_empty() {
                format!("{element}{}", index + 1)
            } else {
                label
            },
            &element,
            frac,
        )?
        .with_occupancy(occupancy)?;
        builder = builder.site(site);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    Ok(SetReport {
        structure: builder.build()?,
        max_displacement,
        n_input,
        n_output,
    })
}

/// Serialized report for the agent and CLI surfaces.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SetSummary {
    /// Sites before the operation.
    pub n_input: usize,
    /// Sites after the operation.
    pub n_output: usize,
    /// Largest Cartesian displacement in Angstrom.
    pub max_displacement: f64,
    /// Hall symbol of the applied group.
    pub hall: String,
    /// Space-group number.
    pub number: u32,
}

impl SetReport {
    /// Summary for the surfaces.
    #[must_use]
    pub fn summary(&self, space_group: &SpaceGroup) -> SetSummary {
        SetSummary {
            n_input: self.n_input,
            n_output: self.n_output,
            max_displacement: self.max_displacement,
            hall: space_group.hall_symbol.clone(),
            number: space_group.number,
        }
    }
}

fn asymmetric_structure(
    structure: &Structure,
    dataset: &SymmetryDataset,
    space_group: &SpaceGroup,
) -> Result<Structure, Error> {
    let p1 = structure.p1()?;
    let mut builder = StructureBuilder::new(structure.name(), *structure.lattice())
        .hall(space_group.hall_symbol.clone());
    for operation in &dataset.operations {
        builder = builder.operation(*operation);
    }
    for (index, atom) in p1.atoms.iter().enumerate() {
        if dataset.orbits.get(index).copied() != Some(index) {
            continue;
        }
        let parent = structure.sites().get(atom.site_idx);
        let mut site = Site::new(atom.label.clone(), &atom.element, atom.frac)?
            .with_occupancy(atom.occupancy)?;
        if let Some(parent) = parent {
            if let Some(u_iso) = parent.u_iso() {
                site = site.with_u_iso(u_iso)?;
            }
            if let Some(u_aniso) = parent.u_aniso() {
                site = site.with_u_aniso(u_aniso)?;
            }
            if let Some(charge) = parent.charge() {
                site = site.with_charge(charge)?;
            }
            if let Some(group) = parent.disorder_group() {
                site = site.with_disorder_group(group);
            }
            if let Some(moment) = parent.magnetic_moment() {
                site = site.with_magnetic_moment(moment);
            }
            for (key, value) in parent.props() {
                site = site.with_prop(key.clone(), value.clone());
            }
        }
        builder = builder.site(site);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

fn p1_structure(original: &Structure, p1: &P1Cell) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new(original.name(), p1.lattice);
    for (index, atom) in p1.atoms.iter().enumerate() {
        let mut site = Site::new(
            format!("{}{}", atom.element, index + 1),
            &atom.element,
            atom.frac,
        )?
        .with_occupancy(atom.occupancy)?;
        if let Some(parent) = original.sites().get(atom.site_idx)
            && let Some(moment) = parent.magnetic_moment()
        {
            site = site.with_magnetic_moment(moment);
        }
        builder = builder.site(site);
    }
    for (key, value) in original.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

/// Wyckoff assignment of one atom of the input cell.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WyckoffAssignment {
    /// Atom index in the input cell.
    pub atom: usize,
    /// Element symbol.
    pub element: String,
    /// Site label.
    pub label: String,
    /// Wyckoff letter.
    pub letter: char,
    /// Site-symmetry symbol.
    pub site_symmetry: String,
    /// Multiplicity in the conventional cell.
    pub multiplicity: usize,
    /// Coordinate template.
    pub coordinates: String,
    /// Free parameters present in the template with their current values.
    pub free_parameters: Vec<(char, f64)>,
}

/// Assigns a Wyckoff letter, multiplicity and free parameters to every atom.
///
/// # Errors
///
/// Returns the errors of [`detect_symmetry`] or the Hall lookup.
pub fn wyckoff_assignments(
    structure: &Structure,
    options: SymmetryOptions,
) -> Result<Vec<WyckoffAssignment>, Error> {
    let dataset = detect_symmetry(structure, options)?;
    let positions = crate::spacegroup::WyckoffPosition::list(dataset.hall_number)?;
    let p1 = structure.p1()?;
    let mut assignments = Vec::with_capacity(p1.atoms.len());
    for (index, atom) in p1.atoms.iter().enumerate() {
        let letter = dataset.wyckoffs.get(index).copied().unwrap_or('?');
        let Some(position) = positions.iter().find(|position| position.letter == letter) else {
            continue;
        };
        assignments.push(WyckoffAssignment {
            atom: index,
            element: atom.element.clone(),
            label: atom.label.clone(),
            letter,
            site_symmetry: dataset
                .site_symmetry_symbols
                .get(index)
                .cloned()
                .unwrap_or_default(),
            multiplicity: position.multiplicity,
            coordinates: position.coordinates.clone(),
            free_parameters: free_parameters(&position.coordinates, atom.frac),
        });
    }
    Ok(assignments)
}

/// Extracts `x`, `y` and `z` from a coordinate template with their values.
///
/// The template is split on commas; the n-th component supplies the value of
/// every variable it mentions from the n-th fractional coordinate. Templates
/// with a mixed term such as `x-y` still report `x` and `y` from their own
/// component.
#[must_use]
pub fn free_parameters(coordinates: &str, frac: [f64; 3]) -> Vec<(char, f64)> {
    let mut parameters: Vec<(char, f64)> = Vec::new();
    for (axis, component) in coordinates.split(',').enumerate() {
        for variable in ['x', 'y', 'z'] {
            if component.contains(variable)
                && !parameters.iter().any(|(existing, _)| *existing == variable)
            {
                parameters.push((variable, frac.get(axis).copied().unwrap_or(0.0)));
            }
        }
    }
    parameters
}

/// Lowers a structure to a subgroup identified by its space-group number.
///
/// The subgroup operations are taken in the default setting of
/// `subgroup_number`; atoms are grouped by matching the operation images to the
/// current P1 atoms and the first atom of every orbit is kept.
///
/// # Errors
///
/// Returns [`Error::Symmetry`] when the subgroup has no identified number or
/// when the tolerance is invalid.
pub fn lower_to_subgroup(
    structure: &Structure,
    subgroup_number: u32,
    tolerance: f64,
) -> Result<Structure, Error> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(Error::Symmetry(format!(
            "subgroup tolerance must be positive, got {tolerance}"
        )));
    }
    let group = SpaceGroup::from_number(subgroup_number)?;
    lower_with_group(structure, &group, tolerance)
}

fn lower_with_group(
    structure: &Structure,
    group: &SpaceGroup,
    tolerance: f64,
) -> Result<Structure, Error> {
    let p1 = structure.p1()?;
    let mut assigned = vec![false; p1.atoms.len()];
    let mut builder =
        StructureBuilder::new(structure.name(), p1.lattice).hall(group.hall_symbol.clone());
    for operation in &group.operations {
        builder = builder.operation(*operation);
    }
    for (index, atom) in p1.atoms.iter().enumerate() {
        if assigned.get(index).copied().unwrap_or(true) {
            continue;
        }
        if let Some(slot) = assigned.get_mut(index) {
            *slot = true;
        }
        for operation in &group.operations {
            let image = operation.apply(atom.frac);
            let mut best: Option<(f64, usize)> = None;
            for (other_index, other) in p1.atoms.iter().enumerate() {
                if other.element != atom.element {
                    continue;
                }
                let distance = p1.lattice.distance(image, other.frac);
                if distance <= tolerance && best.is_none_or(|(value, _)| distance < value) {
                    best = Some((distance, other_index));
                }
            }
            if let Some((_, other_index)) = best
                && let Some(slot) = assigned.get_mut(other_index)
            {
                *slot = true;
            }
        }
        let parent = structure.sites().get(atom.site_idx);
        let mut site = Site::new(atom.label.clone(), &atom.element, atom.frac)?
            .with_occupancy(atom.occupancy)?;
        if let Some(parent) = parent {
            if let Some(u_iso) = parent.u_iso() {
                site = site.with_u_iso(u_iso)?;
            }
            if let Some(charge) = parent.charge() {
                site = site.with_charge(charge)?;
            }
            if let Some(moment) = parent.magnetic_moment() {
                site = site.with_magnetic_moment(moment);
            }
            for (key, value) in parent.props() {
                site = site.with_prop(key.clone(), value.clone());
            }
        }
        builder = builder.site(site);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}
