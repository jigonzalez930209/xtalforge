//! Structural cleanup steps for Prepare for simulation (Task 4.7).
//!
//! These are the steps that change what atoms are present before the force
//! field sees them: merging sites written twice, choosing one disorder
//! component, removing known solvent molecules and checking the cell charge.

mod disorder;

use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, P1Cell, Site, Structure, StructureBuilder, element_info, neighbors};

pub(crate) use disorder::order_disorder;

/// One pair of atoms closer than the reference clash threshold.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CloseContact {
    /// First P1 atom.
    pub i: usize,
    /// Second P1 atom.
    pub j: usize,
    /// Minimum-image distance in Angstrom.
    pub distance: f64,
    /// Element symbols of the two atoms.
    pub elements: [String; 2],
}

/// Net charge of the prepared cell.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ChargeBalance {
    /// Sum of the P1 site charges in elementary charge units.
    pub total_charge: f64,
    /// Whether the total is zero within `1e-6 e`.
    pub neutral: bool,
}

/// Merges sites with the same element at the same periodic position.
///
/// The keeper's occupancy is the sum of the merged sites, clamped to one.
/// Returns the new structure and the number of sites removed.
///
/// # Errors
///
/// Returns the errors of the P1 expansion or the rebuilt structure.
pub(crate) fn merge_duplicates(structure: &Structure) -> Result<(Structure, usize), Error> {
    let p1 = structure.p1()?;
    let mut drop = BTreeSet::new();
    let mut occupancy = BTreeMap::new();
    for i in 0..structure.sites().len() {
        if drop.contains(&i) {
            continue;
        }
        let Some(first) = p1.images_of(i).first().copied() else {
            continue;
        };
        for j in (i + 1)..structure.sites().len() {
            if drop.contains(&j) {
                continue;
            }
            let Some(second) = p1.images_of(j).first().copied() else {
                continue;
            };
            if first.element != second.element
                || p1.lattice.distance(first.frac, second.frac) > 1e-3
            {
                continue;
            }
            let (Some(site_i), Some(site_j)) = (structure.sites().get(i), structure.sites().get(j))
            else {
                continue;
            };
            drop.insert(j);
            *occupancy.entry(i).or_insert(site_i.occupancy()) += site_j.occupancy();
        }
    }
    if drop.is_empty() {
        return Ok((structure.clone(), 0));
    }
    let mut sites = Vec::with_capacity(structure.sites().len() - drop.len());
    for (index, site) in structure.sites().iter().enumerate() {
        if drop.contains(&index) {
            continue;
        }
        let value = occupancy
            .get(&index)
            .copied()
            .unwrap_or_else(|| site.occupancy())
            .min(1.0);
        sites.push(site.clone().with_occupancy(value)?);
    }
    Ok((rebuild(structure, sites)?, drop.len()))
}

pub(crate) fn remove_solvent(structure: &Structure) -> Result<(Structure, usize), Error> {
    let p1 = structure.p1()?;
    let graph = structure.bond_graph()?;
    let mut drop = BTreeSet::new();
    let mut removed = 0;
    for component in graph.components() {
        let elements = component
            .iter()
            .filter_map(|index| p1.atoms.get(*index).map(|atom| atom.element.as_str()))
            .collect::<Vec<_>>();
        if elements.is_empty() || elements.iter().any(|element| is_metal(element)) {
            continue;
        }
        let heavy = elements
            .iter()
            .filter(|element| !matches!(**element, "H" | "D"))
            .copied()
            .collect::<Vec<_>>();
        if solvent_name(&heavy).is_some() {
            drop.extend(component);
            removed += 1;
        }
    }
    if drop.is_empty() {
        return Ok((structure.clone(), 0));
    }
    Ok((without_p1_atoms(structure, &p1, &drop)?, removed))
}

/// Sums the stored charges over the P1 cell.
///
/// # Errors
///
/// Returns the errors of the P1 expansion.
pub(crate) fn charge_balance(structure: &Structure) -> Result<ChargeBalance, Error> {
    let p1 = structure.p1()?;
    let mut total = 0.0;
    for atom in &p1.atoms {
        if let Some(site) = structure.sites().get(atom.site_idx) {
            total += site.charge().unwrap_or(0.0);
        }
    }
    Ok(ChargeBalance {
        total_charge: total,
        neutral: total.abs() <= 1e-6,
    })
}

/// Every P1 pair closer than the reference clash threshold.
///
/// H-H uses 1.4 A, O-O 2.0 A and halogen-halogen 1.9 A; other pairs use
/// `0.8 * (r_i + r_j)` with the covalent radii. The list is sorted by atom
/// indices and is empty for a clean structure.
///
/// # Errors
///
/// Returns the errors of the P1 expansion or the periodic neighbor search.
pub(crate) fn close_contacts(structure: &Structure) -> Result<Vec<CloseContact>, Error> {
    let p1 = structure.p1()?;
    let found = neighbors(structure, 2.1)?;
    let mut out = Vec::new();
    for pair in found {
        if pair.i >= pair.j {
            continue;
        }
        let (Some(first), Some(second)) = (p1.atoms.get(pair.i), p1.atoms.get(pair.j)) else {
            continue;
        };
        if !is_clash(&first.element, &second.element, pair.distance) {
            continue;
        }
        out.push(CloseContact {
            i: pair.i,
            j: pair.j,
            distance: pair.distance,
            elements: [first.element.clone(), second.element.clone()],
        });
    }
    out.sort_by_key(|contact| (contact.i, contact.j, contact.distance.to_bits()));
    Ok(out)
}

fn is_clash(first: &str, second: &str, distance: f64) -> bool {
    let both = |set: &[&str]| set.contains(&first) && set.contains(&second);
    if both(&["H", "D"]) {
        return distance < 1.4;
    }
    if first == "O" && second == "O" {
        return distance < 2.0;
    }
    let halogens = ["F", "Cl", "Br", "I"];
    let (Ok(a), Ok(b)) = (element_info(first), element_info(second)) else {
        return false;
    };
    let covalent = 0.8 * (a.covalent_radius + b.covalent_radius);
    if halogens.contains(&first) && halogens.contains(&second) {
        return distance < covalent.max(1.9);
    }
    distance < covalent
}

fn without_p1_atoms(
    structure: &Structure,
    p1: &P1Cell,
    drop: &BTreeSet<usize>,
) -> Result<Structure, Error> {
    let mut by_site: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (index, atom) in p1.atoms.iter().enumerate() {
        by_site.entry(atom.site_idx).or_default().push(index);
    }
    let whole = by_site.values().all(|atoms| {
        atoms.iter().all(|index| drop.contains(index))
            || atoms.iter().all(|index| !drop.contains(index))
    });
    if whole {
        let mut sites = Vec::new();
        for (index, site) in structure.sites().iter().enumerate() {
            let gone = by_site
                .get(&index)
                .is_some_and(|atoms| atoms.iter().all(|atom| drop.contains(atom)));
            if !gone {
                sites.push(site.clone().with_occupancy(1.0)?.without_disorder_group());
            }
        }
        rebuild(structure, sites)
    } else {
        let mut sites = Vec::new();
        for (index, atom) in p1.atoms.iter().enumerate() {
            if drop.contains(&index) {
                continue;
            }
            sites.push(
                Site::new(
                    format!("{}{}", atom.element, index + 1),
                    &atom.element,
                    atom.frac,
                )?
                .with_occupancy(1.0)?,
            );
        }
        rebuild_p1(structure, sites)
    }
}

pub(super) fn rebuild(structure: &Structure, sites: Vec<Site>) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new(structure.name().to_owned(), *structure.lattice())
        .hall(structure.hall().to_owned());
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    for site in sites {
        builder = builder.site(site);
    }
    builder.build()
}

fn rebuild_p1(structure: &Structure, sites: Vec<Site>) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new(structure.name().to_owned(), *structure.lattice());
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    for site in sites {
        builder = builder.site(site);
    }
    builder.build()
}

pub(super) fn is_metal(element: &str) -> bool {
    element_info(element).is_ok_and(|info| info.is_metal)
}

fn solvent_name(elements: &[&str]) -> Option<&'static str> {
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for element in elements {
        *counts.entry(element).or_default() += 1;
    }
    let formula = if counts.len() == 1 && counts.get("O") == Some(&1) {
        "O".to_owned()
    } else {
        counts
            .iter()
            .map(|(element, count)| format!("{element}{count}"))
            .collect::<String>()
    };
    match formula.as_str() {
        "O" => Some("water"),
        "C1O1" => Some("methanol"),
        "C2O1" => Some("ethanol"),
        "C3N1O1" => Some("DMF"),
        "C5N1O1" => Some("DEF"),
        "C2O1S1" => Some("DMSO"),
        "C3O1" => Some("acetone"),
        "C4O1" => Some("THF"),
        "C2N1" => Some("acetonitrile"),
        "C4O2" => Some("dioxane"),
        "C1Cl2" => Some("dichloromethane"),
        "C1Cl3" => Some("chloroform"),
        "C5N1" => Some("pyridine"),
        "C6" => Some("benzene"),
        "C7" => Some("toluene"),
        _ => None,
    }
}
