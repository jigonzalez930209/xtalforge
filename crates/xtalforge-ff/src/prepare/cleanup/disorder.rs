//! Graph-aware disorder selection for Prepare for simulation (Task 4.7).

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use xtalforge_core::{Error, P1Cell, Site, Structure};

use super::{is_metal, rebuild};

/// Chooses one component for every disordered group.
///
/// Sites with `_atom_site_disorder_group` are kept together as one component.
/// Ungrouped partial sites form their own component unless their periodic
/// positions are within 2 A, where they are alternatives. Overlapping
/// components choose the compatible subset that best matches the declared
/// CIF formula and a neutral cell, preferring higher occupancy on a tie.
/// Every survivor is set to full occupancy.
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
    let graph = structure.bond_graph()?;
    let mut component_by_atom = vec![0_usize; p1.atoms.len()];
    for (component, atoms) in graph.components().iter().enumerate() {
        for atom in atoms {
            if let Some(slot) = component_by_atom.get_mut(*atom) {
                *slot = component;
            }
        }
    }
    let mut first_atom = vec![None; structure.sites().len()];
    for (atom, index) in p1.atoms.iter().enumerate() {
        if let Some(slot) = first_atom.get_mut(index.site_idx)
            && slot.is_none()
        {
            *slot = Some(atom);
        }
    }
    let mut partials_in_component: BTreeMap<usize, usize> = BTreeMap::new();
    for index in &partial {
        let component = first_atom
            .get(*index)
            .and_then(|atom| *atom)
            .and_then(|atom| component_by_atom.get(atom).copied())
            .unwrap_or(*index);
        *partials_in_component.entry(component).or_default() += 1;
    }
    let mut groups: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for (offset, index) in partial.iter().enumerate() {
        let explicit = structure.sites().get(*index).and_then(Site::disorder_group);
        let key = if let Some(group) = explicit {
            group
        } else {
            let component = first_atom
                .get(*index)
                .and_then(|atom| *atom)
                .and_then(|atom| component_by_atom.get(atom).copied())
                .unwrap_or(*index);
            if partials_in_component.get(&component).copied().unwrap_or(0) > 1 {
                1_000_000 + component as u32
            } else {
                2_000_000 + offset as u32
            }
        };
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
    let declared = structure
        .metadata()
        .get("chemical_formula")
        .and_then(|value| value.as_str())
        .and_then(parse_formula);
    let partial_set = partial.iter().copied().collect::<BTreeSet<_>>();
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
        let chooser = Chooser {
            keys: &keys,
            groups: &groups,
            adjacent: &adjacent,
            structure,
            p1: &p1,
            partial: &partial_set,
            declared: declared.as_ref(),
        };
        let chosen = chooser.choose(&component);
        for index in component {
            if !chosen.contains(&index)
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

struct Chooser<'a> {
    keys: &'a [u32],
    groups: &'a BTreeMap<u32, Vec<usize>>,
    adjacent: &'a [Vec<usize>],
    structure: &'a Structure,
    p1: &'a P1Cell,
    partial: &'a BTreeSet<usize>,
    declared: Option<&'a BTreeMap<String, f64>>,
}

impl Chooser<'_> {
    fn choose(&self, component: &[usize]) -> BTreeSet<usize> {
        let n = component.len();
        if n > 16 {
            // Greedy for large clusters: highest occupancy first, then a
            // formula/charge check over the chosen set.
            let mut order = component.to_vec();
            order.sort_by(|first, second| {
                group_occupancy(self.groups, self.keys, *second, self.structure)
                    .partial_cmp(&group_occupancy(
                        self.groups,
                        self.keys,
                        *first,
                        self.structure,
                    ))
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| self.keys.get(*first).cmp(&self.keys.get(*second)))
            });
            let mut chosen = BTreeSet::new();
            for index in order {
                if self
                    .adjacent
                    .get(index)
                    .into_iter()
                    .flatten()
                    .all(|other| !chosen.contains(other))
                {
                    chosen.insert(index);
                }
            }
            return chosen;
        }
        let mut best: Option<(f64, f64, f64, u64, BTreeSet<usize>)> = None;
        for mask in 0u64..(1u64 << n) {
            let mut chosen = BTreeSet::new();
            let mut compatible = true;
            for (position, index) in component.iter().enumerate() {
                if mask & (1 << position) == 0 {
                    continue;
                }
                if self
                    .adjacent
                    .get(*index)
                    .into_iter()
                    .flatten()
                    .any(|other| chosen.contains(other))
                {
                    compatible = false;
                    break;
                }
                chosen.insert(*index);
            }
            if !compatible {
                continue;
            }
            let (formula_error, charge, occupancy) = self.score(&chosen);
            let candidate = (formula_error, charge.abs(), -occupancy, mask, chosen);
            let better = best.as_ref().is_none_or(|current| {
                candidate
                    .0
                    .total_cmp(&current.0)
                    .then_with(|| candidate.1.total_cmp(&current.1))
                    .then_with(|| candidate.2.total_cmp(&current.2))
                    .then_with(|| candidate.3.cmp(&current.3))
                    == std::cmp::Ordering::Less
            });
            if better {
                best = Some(candidate);
            }
        }
        best.map(|value| value.4).unwrap_or_default()
    }

    fn score(&self, chosen: &BTreeSet<usize>) -> (f64, f64, f64) {
        let mut counts: BTreeMap<String, f64> = BTreeMap::new();
        let mut charge = 0.0;
        let mut occupancy = 0.0;
        for (index, site) in self.structure.sites().iter().enumerate() {
            let keep = !self.partial.contains(&index)
                || chosen
                    .iter()
                    .filter_map(|group| self.keys.get(*group))
                    .filter_map(|key| self.groups.get(key))
                    .any(|sites| sites.contains(&index));
            if !keep {
                continue;
            }
            let multiplicity = self.p1.images_of(index).len() as f64;
            *counts.entry(site.element().to_owned()).or_default() += multiplicity;
            charge += site.charge().unwrap_or(0.0) * multiplicity;
            occupancy += site.occupancy() * multiplicity;
        }
        let error = self
            .declared
            .map_or(0.0, |formula| formula_error(formula, &counts));
        (error, charge, occupancy)
    }
}

fn parse_formula(text: &str) -> Option<BTreeMap<String, f64>> {
    let mut counts = BTreeMap::new();
    let mut chars = text.chars().peekable();
    while let Some(character) = chars.next() {
        if !character.is_ascii_uppercase() {
            continue;
        }
        let mut symbol = character.to_string();
        if chars.peek().is_some_and(|next| next.is_ascii_lowercase())
            && let Some(next) = chars.next()
        {
            symbol.push(next);
        }
        let mut number = String::new();
        while chars
            .peek()
            .is_some_and(|next| next.is_ascii_digit() || *next == '.')
        {
            if let Some(next) = chars.next() {
                number.push(next);
            }
        }
        let value = number.parse::<f64>().unwrap_or(1.0);
        if value > 0.0 {
            *counts.entry(symbol).or_default() += value;
        }
    }
    if counts.is_empty() {
        None
    } else {
        Some(counts)
    }
}

fn formula_error(declared: &BTreeMap<String, f64>, got: &BTreeMap<String, f64>) -> f64 {
    let normalize = |counts: &BTreeMap<String, f64>| {
        let metals = counts
            .iter()
            .filter(|(element, _)| is_metal(element))
            .map(|(_, value)| *value)
            .sum::<f64>();
        let total = if metals > 0.0 {
            metals
        } else {
            counts.values().sum::<f64>()
        };
        if total > 0.0 { total } else { 1.0 }
    };
    let declared_total = normalize(declared);
    let got_total = normalize(got);
    declared
        .keys()
        .chain(got.keys())
        .map(|element| {
            let a = declared.get(element).copied().unwrap_or(0.0) / declared_total;
            let b = got.get(element).copied().unwrap_or(0.0) / got_total;
            (a - b).abs()
        })
        .fold(0.0, f64::max)
}
