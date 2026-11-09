//! Structural cleanup steps for Prepare for simulation (Task 4.7).
//!
//! These are the steps that change what atoms are present before the force
//! field sees them: merging sites written twice, choosing one disorder
//! component, removing known solvent molecules and checking the cell charge.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, P1Cell, Site, Structure, StructureBuilder, element_info};

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

/// Chooses one component for every disordered group.
///
/// Sites with `_atom_site_disorder_group` are kept together as one component.
/// Ungrouped partial sites form their own component unless their periodic
/// positions are within 2 A, where they are alternatives. Overlapping
/// components are reduced to the one with the highest total occupancy (lowest
/// group number on a tie), and every survivor is set to full occupancy.
///
/// # Errors
///
/// Returns the errors of the P1 expansion or the rebuilt structure.
pub(crate) fn order_disorder(structure: &Structure) -> Result<(Structure, usize), Error> {
    let p1 = structure.p1()?;
    let partial = structure
        .sites()
        .iter()
        .enumerate()
        .filter(|(_, site)| site.occupancy() < 1.0 - 1e-6)
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    if partial.is_empty() {
        return Ok((structure.clone(), 0));
    }
    let mut groups: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for (offset, index) in partial.iter().enumerate() {
        let key = structure
            .sites()
            .get(*index)
            .and_then(Site::disorder_group)
            .unwrap_or(1_000_000 + offset as u32);
        groups.entry(key).or_default().push(*index);
    }
    let keys = groups.keys().copied().collect::<Vec<_>>();
    let mut adjacent: Vec<Vec<usize>> = vec![Vec::new(); keys.len()];
    for (a, first) in keys.iter().enumerate() {
        for (b, second) in keys.iter().enumerate().skip(a + 1) {
            if overlaps(&groups, *first, *second, &p1) {
                if let Some(row) = adjacent.get_mut(a) {
                    row.push(b);
                }
                if let Some(row) = adjacent.get_mut(b) {
                    row.push(a);
                }
            }
        }
    }
    let mut visited = vec![false; keys.len()];
    let mut drop = BTreeSet::new();
    for start in 0..keys.len() {
        if visited.get(start).copied().unwrap_or(true) {
            continue;
        }
        let mut component = Vec::new();
        let mut queue = VecDeque::from([start]);
        while let Some(index) = queue.pop_front() {
            if visited.get(index).copied().unwrap_or(true) {
                continue;
            }
            if let Some(slot) = visited.get_mut(index) {
                *slot = true;
            }
            component.push(index);
            for next in adjacent.get(index).into_iter().flatten() {
                queue.push_back(*next);
            }
        }
        let best = component
            .iter()
            .max_by(|first, second| {
                let a = group_occupancy(&groups, &keys, **first, structure);
                let b = group_occupancy(&groups, &keys, **second, structure);
                a.partial_cmp(&b)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| keys.get(**second).cmp(&keys.get(**first)))
            })
            .copied()
            .unwrap_or(start);
        for index in component {
            if index != best
                && let Some(group) = keys.get(index).and_then(|key| groups.get(key))
            {
                drop.extend(group.iter().copied());
            }
        }
    }
    let mut sites = Vec::with_capacity(structure.sites().len() - drop.len());
    for (index, site) in structure.sites().iter().enumerate() {
        if drop.contains(&index) {
            continue;
        }
        let occupancy = if site.occupancy() < 1.0 - 1e-6 {
            1.0
        } else {
            site.occupancy()
        };
        sites.push(
            site.clone()
                .with_occupancy(occupancy)?
                .without_disorder_group(),
        );
    }
    Ok((rebuild(structure, sites)?, drop.len()))
}

/// Removes isolated molecules whose heavy-atom formula is a known solvent.
///
/// A component with a metal is framework and stays. Anything else that is not
/// recognised is left in place and counted in the returned message. Returns
/// the new structure and the number of solvent molecules removed.
///
/// # Errors
///
/// Returns the errors of the P1 expansion, bond perception or rebuild.
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

fn overlaps(groups: &BTreeMap<u32, Vec<usize>>, first: u32, second: u32, p1: &P1Cell) -> bool {
    let (Some(a), Some(b)) = (groups.get(&first), groups.get(&second)) else {
        return false;
    };
    for site_a in a {
        for site_b in b {
            for atom_a in p1.images_of(*site_a) {
                for atom_b in p1.images_of(*site_b) {
                    if p1.lattice.distance(atom_a.frac, atom_b.frac) < 2.0 {
                        return true;
                    }
                }
            }
        }
    }
    false
}

fn group_occupancy(
    groups: &BTreeMap<u32, Vec<usize>>,
    keys: &[u32],
    index: usize,
    structure: &Structure,
) -> f64 {
    keys.get(index)
        .and_then(|key| groups.get(key))
        .map(|sites| {
            sites
                .iter()
                .filter_map(|site| structure.sites().get(*site))
                .map(Site::occupancy)
                .sum()
        })
        .unwrap_or(0.0)
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

fn rebuild(structure: &Structure, sites: Vec<Site>) -> Result<Structure, Error> {
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

fn is_metal(element: &str) -> bool {
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
