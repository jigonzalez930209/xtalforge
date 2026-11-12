#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Orderings of partial occupancies (Task 9.5.2).
//!
//! A disordered site is a group of sites at the same fractional position with
//! partial occupancies (or an explicit `disorder_group`). The module expands
//! the group into a supercell, enumerates the symmetry-distinct decorations
//! and ranks them by the Ewald energy of the ordered cell.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use xtalforge_core::{Error, Site, Structure, StructureBuilder};
use xtalforge_symmetry::{MatcherOptions, match_structures};

use crate::bvs::{BvsOptions, oxidation_states};
use crate::ewald::{EwaldOptions, ewald_energy};

/// One ordered decoration.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Ordering {
    /// Ordered P1 structure.
    pub structure: Structure,
    /// Ewald energy in eV; `None` when the decoration is not charge neutral.
    pub energy: Option<f64>,
    /// Species chosen for each disordered position, in enumeration order.
    pub assignment: Vec<String>,
}

/// Options for [`orderings`].
#[derive(Debug, Clone, PartialEq)]
pub struct OrderingOptions {
    /// Supercell repeats.
    pub supercell: [u32; 3],
    /// Maximum number of decorations to evaluate.
    pub max_orderings: usize,
    /// Ewald options.
    pub ewald: EwaldOptions,
}

impl Default for OrderingOptions {
    fn default() -> Self {
        Self {
            supercell: [1, 1, 1],
            max_orderings: 256,
            ewald: EwaldOptions::default(),
        }
    }
}

/// Enumerates the decorations of the disordered sites of a structure.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a structure without partial occupancies,
/// a zero supercell, or a bad option.
pub fn orderings(structure: &Structure, options: &OrderingOptions) -> Result<Vec<Ordering>, Error> {
    if options.supercell.contains(&0) {
        return Err(Error::InvalidInput(
            "ordering supercell repeats must be positive".to_owned(),
        ));
    }
    if options.max_orderings == 0 {
        return Err(Error::InvalidInput(
            "max_orderings must be positive".to_owned(),
        ));
    }
    let p1 = structure.p1()?;
    let groups = disorder_groups(&p1)?;
    if groups.is_empty() {
        return Err(Error::InvalidInput(
            "the structure has no partial occupancies to order".to_owned(),
        ));
    }
    let [na, nb, nc] = options.supercell;
    let host = host_structure(&p1, &groups)?;
    let host = host.supercell(na, nb, nc)?;
    let charges = element_charges(&p1, &groups)?;
    let mut positions: Vec<([f64; 3], Vec<String>)> = Vec::new();
    for group in &groups {
        for cell_a in 0..na {
            for cell_b in 0..nb {
                for cell_c in 0..nc {
                    for position in &group.positions {
                        let frac = [
                            (position[0] + f64::from(cell_a)) / f64::from(na),
                            (position[1] + f64::from(cell_b)) / f64::from(nb),
                            (position[2] + f64::from(cell_c)) / f64::from(nc),
                        ];
                        positions.push((frac, group.species.clone()));
                    }
                }
            }
        }
    }
    let mut assignments: Vec<Vec<String>> = Vec::new();
    let mut current = Vec::new();
    enumerate(
        &positions,
        0,
        &mut current,
        &mut assignments,
        options.max_orderings,
    );
    let mut results: Vec<Ordering> = Vec::new();
    let matcher = MatcherOptions {
        ltol: 0.05,
        stol: 0.1,
        angle_tol: 1.0,
        ..MatcherOptions::default()
    };
    for assignment in assignments {
        let ordered = build_ordered(&host, &positions, &assignment)?;
        let mut duplicate = false;
        for existing in &results {
            if match_structures(&existing.structure, &ordered, &matcher)?.matched {
                duplicate = true;
                break;
            }
        }
        if duplicate {
            continue;
        }
        let charge_vector: Vec<f64> = ordered
            .sites()
            .iter()
            .map(|site| charges.get(site.element()).copied().unwrap_or(0.0))
            .collect();
        let energy = ewald_energy(&ordered, &charge_vector, &options.ewald).ok();
        results.push(Ordering {
            structure: ordered,
            energy,
            assignment,
        });
    }
    results.sort_by(|first, second| match (first.energy, second.energy) {
        (Some(a), Some(b)) => a.partial_cmp(&b).unwrap_or(std::cmp::Ordering::Equal),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => std::cmp::Ordering::Equal,
    });
    Ok(results)
}

#[derive(Debug, Clone)]
struct DisorderedGroup {
    positions: Vec<[f64; 3]>,
    species: Vec<String>,
}

fn disorder_groups(p1: &xtalforge_core::P1Cell) -> Result<Vec<DisorderedGroup>, Error> {
    let mut by_key: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (index, atom) in p1.atoms.iter().enumerate() {
        let key = atom.disorder_group.map_or_else(
            || {
                format!(
                    "{:.5}:{:.5}:{:.5}",
                    atom.frac[0], atom.frac[1], atom.frac[2]
                )
            },
            |group| format!("group:{group}"),
        );
        by_key.entry(key).or_default().push(index);
    }
    let mut groups = Vec::new();
    for members in by_key.into_values() {
        if members.len() < 2 {
            continue;
        }
        let mut positions: Vec<[f64; 3]> = Vec::new();
        let mut species: Vec<String> = Vec::new();
        for index in members {
            let Some(atom) = p1.atoms.get(index) else {
                continue;
            };
            if !positions.iter().any(|position| {
                (position[0] - atom.frac[0]).abs() < 1e-5
                    && (position[1] - atom.frac[1]).abs() < 1e-5
                    && (position[2] - atom.frac[2]).abs() < 1e-5
            }) {
                positions.push(atom.frac);
            }
            if !species.contains(&atom.element) {
                species.push(atom.element.clone());
            }
        }
        if positions.len() == 1 && species.len() >= 2 {
            species.sort();
            groups.push(DisorderedGroup { positions, species });
        }
    }
    Ok(groups)
}

fn host_structure(
    p1: &xtalforge_core::P1Cell,
    groups: &[DisorderedGroup],
) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new("host", p1.lattice);
    for atom in &p1.atoms {
        let is_group = groups.iter().any(|group| {
            group.positions.iter().any(|position| {
                (position[0] - atom.frac[0]).abs() < 1e-5
                    && (position[1] - atom.frac[1]).abs() < 1e-5
                    && (position[2] - atom.frac[2]).abs() < 1e-5
            })
        });
        if !is_group {
            builder = builder.site(
                Site::new(atom.label.clone(), &atom.element, atom.frac)?
                    .with_occupancy(atom.occupancy)?,
            );
        }
    }
    builder.build()
}

fn build_ordered(
    host: &Structure,
    positions: &[([f64; 3], Vec<String>)],
    assignment: &[String],
) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new(host.name(), *host.lattice());
    for site in host.sites() {
        builder = builder.site(site.clone());
    }
    for (index, ((frac, _), species)) in positions.iter().zip(assignment).enumerate() {
        builder = builder.site(Site::new(
            format!("{species}{}", index + 1),
            species,
            *frac,
        )?);
    }
    builder.build()
}

fn enumerate(
    positions: &[([f64; 3], Vec<String>)],
    index: usize,
    current: &mut Vec<String>,
    output: &mut Vec<Vec<String>>,
    limit: usize,
) {
    if output.len() >= limit {
        return;
    }
    if index >= positions.len() {
        output.push(current.clone());
        return;
    }
    let Some((_, species)) = positions.get(index) else {
        return;
    };
    for choice in species {
        current.push(choice.clone());
        enumerate(positions, index + 1, current, output, limit);
        current.pop();
        if output.len() >= limit {
            return;
        }
    }
}

fn element_charges(
    p1: &xtalforge_core::P1Cell,
    groups: &[DisorderedGroup],
) -> Result<BTreeMap<String, f64>, Error> {
    let mut builder = StructureBuilder::new("sample", p1.lattice);
    for atom in &p1.atoms {
        builder = builder.site(Site::new(atom.label.clone(), &atom.element, atom.frac)?);
    }
    for group in groups {
        if let (Some(position), Some(species)) = (group.positions.first(), group.species.first()) {
            builder = builder.site(Site::new(format!("{species}0"), species, *position)?);
        }
    }
    let structure = builder.build()?;
    let mut charges = BTreeMap::new();
    if let Ok(report) = oxidation_states(&structure, &BvsOptions::default()) {
        for (element, state) in report.by_element {
            charges.insert(element, f64::from(state));
        }
    }
    Ok(charges)
}
