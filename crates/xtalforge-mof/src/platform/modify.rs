#![allow(clippy::indexing_slicing)] // fixed 3x3 and bounded graph arrays

//! Post-synthetic modification and MOF variants (Tasks 10.5.3-10.5.6).

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Site, Structure, StructureBuilder};

use crate::platform::deconstruct::{UnitKind, deconstruct};
use crate::platform::rng::SplitMix64;
use crate::platform::util::{bond_adjacency, p1_and_bonds};

/// Functional group for post-synthetic modification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FunctionalGroup {
    /// Amino group.
    Amino,
    /// Hydroxyl group.
    Hydroxyl,
    /// Methyl group.
    Methyl,
    /// Fluorine atom.
    Fluorine,
}

/// Options for [`functionalize`].
#[derive(Debug, Clone, PartialEq)]
pub struct FunctionalizeOptions {
    /// Functional group.
    pub group: FunctionalGroup,
    /// P1 carbon indices to functionalise; empty means every eligible carbon.
    pub positions: Vec<usize>,
    /// Seed used when `positions` is empty and `fraction < 1`.
    pub seed: u64,
    /// Fraction of eligible positions to functionalise.
    pub fraction: f64,
}

impl Default for FunctionalizeOptions {
    fn default() -> Self {
        Self {
            group: FunctionalGroup::Amino,
            positions: Vec::new(),
            seed: 0,
            fraction: 1.0,
        }
    }
}

/// Functionalises linker carbons by replacing a hydrogen.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a bad fraction or non-finite geometry,
/// or a builder error.
pub fn functionalize(
    structure: &Structure,
    options: &FunctionalizeOptions,
) -> Result<Structure, Error> {
    if !options.fraction.is_finite() || !(0.0..=1.0).contains(&options.fraction) {
        return Err(Error::InvalidInput(
            "the functionalisation fraction must be in [0, 1]".to_owned(),
        ));
    }
    let (p1, graph) = p1_and_bonds(structure)?;
    let adjacency = bond_adjacency(&graph, p1.sites().len());
    let mut eligible: Vec<(usize, usize)> = Vec::new();
    for (carbon, site) in p1.sites().iter().enumerate() {
        if site.element() != "C" {
            continue;
        }
        let hydrogens: Vec<usize> = adjacency
            .get(carbon)
            .map(Vec::as_slice)
            .unwrap_or(&[])
            .iter()
            .copied()
            .filter(|index| {
                p1.sites()
                    .get(*index)
                    .is_some_and(|site| site.element() == "H")
            })
            .collect();
        for hydrogen in hydrogens {
            eligible.push((carbon, hydrogen));
        }
    }
    let mut selected: Vec<(usize, usize)> = if options.positions.is_empty() {
        eligible
    } else {
        eligible
            .into_iter()
            .filter(|(carbon, _)| options.positions.contains(carbon))
            .collect()
    };
    let mut rng = SplitMix64::new(options.seed);
    rng.shuffle(&mut selected);
    let count = ((selected.len() as f64) * options.fraction).round() as usize;
    selected.truncate(count);
    let mut removed: Vec<usize> = selected.iter().map(|(_, hydrogen)| *hydrogen).collect();
    removed.sort_unstable();
    let mut builder = StructureBuilder::new(structure.name(), *p1.lattice());
    for (index, site) in p1.sites().iter().enumerate() {
        if removed.binary_search(&index).is_ok() {
            continue;
        }
        builder = builder.site(site.clone());
    }
    for (serial, (carbon, hydrogen)) in selected.iter().enumerate() {
        let Some(carbon_site) = p1.sites().get(*carbon) else {
            continue;
        };
        let Some(hydrogen_site) = p1.sites().get(*hydrogen) else {
            continue;
        };
        let carbon_cart = p1.lattice().to_cart(carbon_site.frac());
        let hydrogen_cart = p1.lattice().to_cart(hydrogen_site.frac());
        let direction = normalize([
            hydrogen_cart[0] - carbon_cart[0],
            hydrogen_cart[1] - carbon_cart[1],
            hydrogen_cart[2] - carbon_cart[2],
        ]);
        builder = add_group(builder, &p1, carbon_cart, direction, options.group, serial)?;
    }
    builder.build()
}

fn add_group(
    builder: StructureBuilder,
    p1: &Structure,
    origin: [f64; 3],
    direction: [f64; 3],
    group: FunctionalGroup,
    serial: usize,
) -> Result<StructureBuilder, Error> {
    let place = |distance: f64| {
        p1.lattice().to_frac([
            origin[0] + direction[0] * distance,
            origin[1] + direction[1] * distance,
            origin[2] + direction[2] * distance,
        ])
    };
    let (element, distance) = match group {
        FunctionalGroup::Fluorine => ("F", 1.35),
        FunctionalGroup::Hydroxyl => ("O", 1.36),
        FunctionalGroup::Amino => ("N", 1.45),
        FunctionalGroup::Methyl => ("C", 1.50),
    };
    Ok(builder.site(Site::new(
        format!("{element}{}", serial + 1),
        element,
        place(distance),
    )?))
}

fn normalize(vector: [f64; 3]) -> [f64; 3] {
    let norm = (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]).sqrt();
    if norm < 1e-12 {
        [1.0, 0.0, 0.0]
    } else {
        [vector[0] / norm, vector[1] / norm, vector[2] / norm]
    }
}

/// Capping species for a missing linker or cluster defect.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CapSpecies {
    /// Formate.
    Formate,
    /// Water.
    Water,
    /// Hydroxide.
    Hydroxide,
}

/// Options for [`missing_defects`].
#[derive(Debug, Clone, PartialEq)]
pub struct DefectOptions {
    /// Fraction of symmetry-distinct linkers to remove.
    pub linker_fraction: f64,
    /// Fraction of symmetry-distinct nodes to remove.
    pub cluster_fraction: f64,
    /// Capping species.
    pub cap: CapSpecies,
    /// Deterministic seed.
    pub seed: u64,
}

impl Default for DefectOptions {
    fn default() -> Self {
        Self {
            linker_fraction: 0.0,
            cluster_fraction: 0.0,
            cap: CapSpecies::Formate,
            seed: 0,
        }
    }
}

/// Result of a defect generation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DefectReport {
    /// Defective P1 structure.
    pub structure: Structure,
    /// Linker units removed.
    pub removed_linkers: usize,
    /// Node units removed.
    pub removed_clusters: usize,
    /// Capping atoms added.
    pub caps: usize,
}

/// Removes linkers and clusters and caps the exposed connections.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for bad fractions, or the errors of
/// deconstruction and the builder.
pub fn missing_defects(
    structure: &Structure,
    options: &DefectOptions,
) -> Result<DefectReport, Error> {
    for value in [options.linker_fraction, options.cluster_fraction] {
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err(Error::InvalidInput(
                "defect fractions must be in [0, 1]".to_owned(),
            ));
        }
    }
    let deconstruction = deconstruct(structure)?;
    let p1 = &deconstruction.structure;
    let mut rng = SplitMix64::new(options.seed);
    let mut remove_atoms: Vec<usize> = Vec::new();
    let mut removed_linkers = 0;
    let mut removed_clusters = 0;
    let mut caps = 0;
    let mut linkers: Vec<_> = deconstruction
        .units
        .iter()
        .filter(|unit| unit.kind == UnitKind::Linker)
        .collect();
    rng.shuffle(&mut linkers);
    let linker_count = ((linkers.len() as f64) * options.linker_fraction).round() as usize;
    let mut builder = StructureBuilder::new(structure.name(), *p1.lattice());
    let mut nodes: Vec<_> = deconstruction
        .units
        .iter()
        .filter(|unit| unit.kind == UnitKind::Node)
        .collect();
    rng.shuffle(&mut nodes);
    let cluster_count = ((nodes.len() as f64) * options.cluster_fraction).round() as usize;
    for unit in nodes.iter().take(cluster_count) {
        removed_clusters += 1;
        remove_atoms.extend(unit.atoms.iter().copied());
    }
    for unit in linkers.iter().take(linker_count) {
        removed_linkers += 1;
        remove_atoms.extend(unit.atoms.iter().copied());
    }
    remove_atoms.sort_unstable();
    remove_atoms.dedup();
    let mut kept: Vec<usize> = Vec::new();
    for (index, site) in p1.sites().iter().enumerate() {
        if remove_atoms.binary_search(&index).is_err() {
            kept.push(index);
            builder = builder.site(site.clone());
        }
    }
    // Cap exposed node connections with one atom each.
    for unit in &deconstruction.units {
        if unit.kind != UnitKind::Node || remove_atoms.iter().any(|atom| unit.atoms.contains(atom))
        {
            continue;
        }
        for connection in &unit.connections {
            if !remove_atoms.contains(&connection.other) {
                continue;
            }
            let Some(site) = p1.sites().get(connection.atom) else {
                continue;
            };
            let here = site.frac();
            let here_cart = p1.lattice().to_cart(here);
            let other_cart = p1
                .sites()
                .get(connection.other)
                .map(|other| p1.lattice().to_cart(other.frac()));
            let direction = other_cart.map_or([0.0, 0.0, 1.0], |other| {
                normalize([
                    here_cart[0] - other[0],
                    here_cart[1] - other[1],
                    here_cart[2] - other[2],
                ])
            });
            let capped = p1.lattice().to_frac([
                here_cart[0] + direction[0] * 1.4,
                here_cart[1] + direction[1] * 1.4,
                here_cart[2] + direction[2] * 1.4,
            ]);
            let (element, label) = match options.cap {
                CapSpecies::Formate => ("C", "Ccap"),
                CapSpecies::Water => ("O", "Ocap"),
                CapSpecies::Hydroxide => ("O", "Ocap"),
            };
            builder = builder.site(Site::new(format!("{label}{caps}"), element, capped)?);
            caps += 1;
        }
    }
    Ok(DefectReport {
        structure: builder.build()?,
        removed_linkers,
        removed_clusters,
        caps,
    })
}

mod variants;
pub use variants::{MixtureOptions, interpenetrate, mixed_framework, stack_layers};
