//! Derivative-structure enumeration (Task 19.2.2).
//!
//! The enumeration combines two deterministic levels:
//!
//! 1. every Hermite-normal-form (HNF) supercell in an index range, ordered by
//!    `(a, d, f, b, c, e)`;
//! 2. every species decoration of the disordered groups that matches the
//!    occupancy-derived concentration, ordered lexicographically.
//!
//! Decorations that are equivalent under the lattice symmetry are removed with
//! [`match_structures`]. No random seed is used, so repeating a call returns
//! the same list in the same order.

use serde::{Deserialize, Serialize};

use xtalforge_core::{Error, Site, Structure, StructureBuilder};
use xtalforge_symmetry::{MatcherOptions, match_structures};

use crate::disorder::{DisorderGroup, DisorderOccupant, disorder_groups};

/// Options for [`derivative_structures`].
#[derive(Debug, Clone, PartialEq)]
pub struct DerivativeOptions {
    /// Smallest supercell index, inclusive.
    pub min_index: u32,
    /// Largest supercell index, inclusive.
    pub max_index: u32,
    /// Maximum number of symmetry-distinct derivatives returned.
    pub max_structures: usize,
    /// Maximum number of candidate decorations generated overall.
    pub max_candidates: usize,
    /// Occupancy below which a species is treated as absent.
    pub concentration_tolerance: f64,
    /// Explicit transformation matrices; when non-empty, the index range is
    /// ignored and these matrices are used in the given order.
    pub matrices: Vec<[[i32; 3]; 3]>,
    /// Matcher used to remove symmetry-equivalent decorations.
    pub matcher: MatcherOptions,
}

impl Default for DerivativeOptions {
    fn default() -> Self {
        Self {
            min_index: 1,
            max_index: 4,
            max_structures: 64,
            max_candidates: 100_000,
            concentration_tolerance: 1e-6,
            matrices: Vec::new(),
            matcher: MatcherOptions {
                ltol: 0.05,
                stol: 0.1,
                angle_tol: 1.0,
                allow_supercell: false,
                ..MatcherOptions::default()
            },
        }
    }
}

/// One symmetry-distinct derivative structure.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DerivativeStructure {
    /// Ordered supercell.
    pub supercell: Structure,
    /// HNF matrix whose rows are the new lattice vectors in the parent basis.
    pub matrix: [[i32; 3]; 3],
    /// Species assigned to each disordered position, in the cell order.
    pub decoration: Vec<String>,
    /// Number of symmetry operations of the ordered cell, including identity.
    pub symmetry_count: usize,
}

mod canonical;
mod hnf;

use canonical::{CanonicalSignature, signature_and_symmetry};
use hnf::supercell_positions;
pub use hnf::{MAX_HNF_INDEX, hnf_supercells};

/// Enumerates symmetry-distinct derivative structures of a disordered cell.
///
/// The concentration of each disorder group comes from the summed occupancies
/// of its species. Counts in a supercell follow the largest-remainder rule;
/// a supercell too small for every species is skipped. A group must have one
/// fractional position (split disorder groups are rejected) and at least two
/// species with positive occupancy. When [`DerivativeOptions::matrices`] is
/// non-empty it replaces the HNF index range.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a structure without partial occupancy,
/// a split or single-species group, or invalid options, or the errors of the
/// supercell construction, the symmetry search or the matcher.
pub fn derivative_structures(
    structure: &Structure,
    options: &DerivativeOptions,
) -> Result<Vec<DerivativeStructure>, Error> {
    validate_options(options)?;
    let groups = disorder_groups(structure)?;
    if groups.is_empty() {
        return Err(Error::InvalidInput(
            "the structure has no partial occupancy to decorate".to_owned(),
        ));
    }
    for group in &groups {
        if group.positions.len() != 1 {
            return Err(Error::InvalidInput(
                "derivative enumeration needs one position per disorder group".to_owned(),
            ));
        }
    }
    let p1 = structure.as_p1()?;
    let host = host_without_groups(&p1, &groups)?;
    let mut results: Vec<DerivativeStructure> = Vec::new();
    let mut signatures: Vec<Option<CanonicalSignature>> = Vec::new();
    let mut generated = 0usize;
    let matrices = if options.matrices.is_empty() {
        let mut all = Vec::new();
        for index in options.min_index..=options.max_index {
            all.extend(hnf_supercells(index)?);
        }
        all
    } else {
        options.matrices.clone()
    };
    for matrix in matrices {
        if results.len() >= options.max_structures || generated >= options.max_candidates {
            break;
        }
        let Some(plans) = supercell_plans(&groups, matrix, options)? else {
            continue;
        };
        let host_supercell = host.supercell_matrix(matrix.map(|row| row.map(f64::from)))?;
        let mut positions = Vec::new();
        for plan in &plans {
            positions.extend(plan.positions.iter().copied());
        }
        let remaining = options.max_candidates - generated;
        let mut candidates = Vec::new();
        enumerate_groups(&plans, 0, &mut Vec::new(), &mut candidates, remaining);
        generated += candidates.len();
        for decoration in candidates {
            if results.len() >= options.max_structures {
                break;
            }
            let ordered = build_ordered(&host_supercell, &positions, &decoration)?;
            let mut duplicate = false;
            for existing in &results {
                if match_structures(&existing.supercell, &ordered, &options.matcher)?.matched {
                    duplicate = true;
                    break;
                }
            }
            if duplicate {
                continue;
            }
            let (signature, symmetry_count) = signature_and_symmetry(&ordered);
            if let Some(candidate) = &signature
                && signatures
                    .iter()
                    .flatten()
                    .any(|existing| existing == candidate)
            {
                continue;
            }
            results.push(DerivativeStructure {
                supercell: ordered,
                matrix,
                decoration,
                symmetry_count,
            });
            signatures.push(signature);
        }
    }
    Ok(results)
}

fn validate_options(options: &DerivativeOptions) -> Result<(), Error> {
    if options.matrices.is_empty() {
        if options.min_index == 0 || options.max_index < options.min_index {
            return Err(Error::InvalidInput(format!(
                "HNF index range must satisfy 1 <= min <= max, got {}..{}",
                options.min_index, options.max_index
            )));
        }
        if options.max_index > MAX_HNF_INDEX {
            return Err(Error::InvalidInput(format!(
                "HNF index must not exceed {MAX_HNF_INDEX}, got {}",
                options.max_index
            )));
        }
    }
    if options.max_structures == 0 || options.max_candidates == 0 {
        return Err(Error::InvalidInput(
            "derivative max_structures and max_candidates must be positive".to_owned(),
        ));
    }
    if !options.concentration_tolerance.is_finite() || options.concentration_tolerance < 0.0 {
        return Err(Error::InvalidInput(
            "the concentration tolerance must be non-negative and finite".to_owned(),
        ));
    }
    Ok(())
}

fn host_without_groups(p1: &Structure, groups: &[DisorderGroup]) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new(p1.name(), *p1.lattice());
    for site in p1.sites() {
        let disordered = groups.iter().any(|group| {
            group
                .positions
                .iter()
                .any(|position| same_position(site.frac(), position.frac))
        });
        if !disordered {
            builder = builder.site(site.clone());
        }
    }
    builder.build()
}

struct GroupPlan {
    positions: Vec<[f64; 3]>,
    assignments: Vec<Vec<String>>,
}

fn supercell_plans(
    groups: &[DisorderGroup],
    matrix: [[i32; 3]; 3],
    options: &DerivativeOptions,
) -> Result<Option<Vec<GroupPlan>>, Error> {
    let mut plans = Vec::with_capacity(groups.len());
    for group in groups {
        let occupants = positive_occupants(group, options.concentration_tolerance);
        if occupants.len() < 2 {
            return Err(Error::InvalidInput(
                "each disorder group needs at least two species with occupancy".to_owned(),
            ));
        }
        let total: f64 = occupants.iter().map(|occupant| occupant.occupancy).sum();
        if total <= 0.0 {
            return Err(Error::InvalidInput(
                "the disorder group has no occupancy".to_owned(),
            ));
        }
        let species: Vec<String> = occupants
            .iter()
            .map(|occupant| occupant.element.clone())
            .collect();
        let fractions: Vec<f64> = occupants
            .iter()
            .map(|occupant| occupant.occupancy / total)
            .collect();
        let Some(position) = group.positions.first() else {
            continue;
        };
        let positions = supercell_positions(position.frac, matrix)?;
        let Some(counts) = allocate_counts(&fractions, positions.len(), options) else {
            return Ok(None);
        };
        let assignments = enumerate_assignments(&species, &counts, options.max_candidates);
        plans.push(GroupPlan {
            positions,
            assignments,
        });
    }
    if plans.is_empty() {
        return Err(Error::InvalidInput(
            "the disorder groups have no position to decorate".to_owned(),
        ));
    }
    Ok(Some(plans))
}

fn positive_occupants(group: &DisorderGroup, tolerance: f64) -> Vec<DisorderOccupant> {
    group
        .aggregate_occupants()
        .into_iter()
        .filter(|occupant| occupant.occupancy > tolerance)
        .collect()
}

fn allocate_counts(
    fractions: &[f64],
    sites: usize,
    options: &DerivativeOptions,
) -> Option<Vec<usize>> {
    let mut counts = vec![0usize; fractions.len()];
    let mut remainders: Vec<(f64, usize)> = Vec::new();
    let mut assigned = 0usize;
    for (index, fraction) in fractions.iter().enumerate() {
        if *fraction <= options.concentration_tolerance {
            continue;
        }
        let exact = fraction * sites as f64;
        let floor = exact.floor() as usize;
        if let Some(slot) = counts.get_mut(index) {
            *slot = floor;
        }
        assigned += floor;
        remainders.push((exact - floor as f64, index));
    }
    remainders.sort_by(|first, second| {
        second
            .0
            .partial_cmp(&first.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| first.1.cmp(&second.1))
    });
    let mut remaining = sites.saturating_sub(assigned);
    for (_, index) in &remainders {
        if remaining == 0 {
            break;
        }
        if let Some(slot) = counts.get_mut(*index) {
            *slot += 1;
            remaining -= 1;
        }
    }
    for (fraction, count) in fractions.iter().zip(&counts) {
        if *fraction > options.concentration_tolerance && *count == 0 {
            return None;
        }
    }
    Some(counts)
}

fn enumerate_assignments(species: &[String], counts: &[usize], limit: usize) -> Vec<Vec<String>> {
    let mut output = Vec::new();
    let mut remaining = counts.to_vec();
    let mut current = Vec::new();
    enumerate_species(species, &mut remaining, &mut current, &mut output, limit);
    output
}

fn enumerate_species(
    species: &[String],
    remaining: &mut [usize],
    current: &mut Vec<String>,
    output: &mut Vec<Vec<String>>,
    limit: usize,
) {
    if output.len() >= limit {
        return;
    }
    let mut placed = false;
    for index in 0..remaining.len() {
        let available = remaining.get(index).copied().unwrap_or(0);
        if available == 0 {
            continue;
        }
        let Some(label) = species.get(index) else {
            continue;
        };
        let Some(slot) = remaining.get_mut(index) else {
            continue;
        };
        *slot -= 1;
        current.push(label.clone());
        placed = true;
        enumerate_species(species, remaining, current, output, limit);
        current.pop();
        if let Some(slot) = remaining.get_mut(index) {
            *slot += 1;
        }
        if output.len() >= limit {
            return;
        }
    }
    if !placed {
        output.push(current.clone());
    }
}

fn enumerate_groups(
    plans: &[GroupPlan],
    plan_index: usize,
    current: &mut Vec<String>,
    output: &mut Vec<Vec<String>>,
    limit: usize,
) {
    if output.len() >= limit {
        return;
    }
    let Some(plan) = plans.get(plan_index) else {
        output.push(current.clone());
        return;
    };
    for assignment in &plan.assignments {
        if output.len() >= limit {
            return;
        }
        let length = current.len();
        current.extend(assignment.iter().cloned());
        enumerate_groups(plans, plan_index + 1, current, output, limit);
        current.truncate(length);
    }
}

fn build_ordered(
    host: &Structure,
    positions: &[[f64; 3]],
    decoration: &[String],
) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new(host.name(), *host.lattice());
    for site in host.sites() {
        builder = builder.site(site.clone());
    }
    for (index, (frac, species)) in positions.iter().zip(decoration).enumerate() {
        builder = builder.site(Site::new(
            format!("{species}{}", index + 1),
            species,
            *frac,
        )?);
    }
    builder.build()
}

fn same_position(first: [f64; 3], second: [f64; 3]) -> bool {
    let delta = [
        first[0] - second[0],
        first[1] - second[1],
        first[2] - second[2],
    ];
    delta.iter().all(|value| {
        let wrapped = value - value.round();
        wrapped.abs() < 1e-5
    })
}
