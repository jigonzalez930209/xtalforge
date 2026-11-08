//! Amorphous polymer chains: packing, build steps and equilibration protocol
//! (Tasks 6.9.1 to 6.9.4).

use serde::Deserialize;
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

use crate::embed::{EmbedOptions, embed_smiles};
use crate::smiles::parse_smiles;

/// Molar mass of carbon, gram per mole.
const MOLAR_MASS_C: f64 = 12.011;
/// Avogadro constant, per mole.
const AVOGADRO: f64 = 6.022_140_76e23;
/// Cubic Angstrom per cubic centimetre.
const A3_PER_CM3: f64 = 1.0e24;
/// Chains closer than this to an already placed chain are skipped, Angstrom.
const MIN_CHAIN_DISTANCE: f64 = 1.0;
/// Empty space around a replicated chain, Angstrom.
const CHAIN_PADDING: f64 = 10.0;
/// Gap between two replicated repeat units, Angstrom.
const REPEAT_GAP: f64 = 1.5;

/// Target conditions for [`build_polymer`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PolymerOptions {
    /// Target density of the packed model in gram per cubic centimetre.
    pub density: f64,
}

impl Default for PolymerOptions {
    fn default() -> Self {
        Self { density: 1.0 }
    }
}

/// A packed polymer and the steps that produced it.
#[derive(Debug, Clone, PartialEq)]
pub struct PolymerBuild {
    /// Packed P1 structure.
    pub structure: Structure,
    /// Human-readable build steps in order.
    pub steps: Vec<String>,
}

/// Packs whole chains into a cubic P1 cell at the requested density.
///
/// The cell side follows from the total atom count, a carbon mass estimate
/// (`12.011 g/mol` per atom) and the target density. Chains are translated to
/// deterministic cubic-grid offsets and a chain is skipped when any of its
/// atoms comes closer than one Angstrom to an already placed atom, including
/// periodic images. Labels are renumbered as element plus a running index.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for no chains, an empty chain, a
/// non-positive or non-finite density, a non-representable cell, or a site
/// that cannot be rebuilt.
pub fn pack_chains(chains: &[Structure], density: f64) -> Result<Structure, Error> {
    if chains.is_empty() {
        return Err(Error::InvalidInput(
            "packing needs at least one chain".to_owned(),
        ));
    }
    if !density.is_finite() || density <= 0.0 {
        return Err(Error::InvalidInput(
            "the packing density must be positive and finite".to_owned(),
        ));
    }
    let mut atom_count = 0usize;
    for chain in chains {
        if chain.sites().is_empty() {
            return Err(Error::InvalidInput(
                "a packed chain must contain at least one atom".to_owned(),
            ));
        }
        atom_count = atom_count
            .checked_add(chain.sites().len())
            .ok_or_else(|| Error::InvalidInput("the packed atom count overflows".to_owned()))?;
    }
    let side = cubic_side(atom_count, density)?;
    let lattice = Lattice::cubic(side)?;
    let grid = grid_side(chains.len());
    let spacing = side / grid as f64;
    let mut builder = StructureBuilder::new("packed_chains", lattice);
    let mut placed: Vec<[f64; 3]> = Vec::new();
    let mut label = 0usize;
    for (chain_index, chain) in chains.iter().enumerate() {
        let offset = grid_offset(chain_index, grid, spacing);
        let points = chain.cart_sites();
        let mut candidates = Vec::with_capacity(points.len());
        for point in &points {
            let cart = [
                point[0] + offset[0],
                point[1] + offset[1],
                point[2] + offset[2],
            ];
            candidates.push(Lattice::wrap(lattice.to_frac(cart)));
        }
        let clashes = candidates.iter().any(|candidate| {
            placed
                .iter()
                .any(|other| lattice.distance(*candidate, *other) < MIN_CHAIN_DISTANCE)
        });
        if clashes {
            continue;
        }
        for (site, frac) in chain.sites().iter().zip(candidates) {
            label += 1;
            builder = builder.site(Site::new(
                format!("{}{}", site.element(), label),
                site.element(),
                frac,
            )?);
            placed.push(frac);
        }
    }
    builder.build()
}

/// Builds one chain from a SMILES repeat unit and packs it.
///
/// The repeat unit is embedded with [`EmbedOptions::default`], copied
/// `monomers` times along `x` with a small gap, and packed with
/// [`pack_chains`] at [`PolymerOptions::density`]. The returned steps record
/// the parse, embed, replicate and pack operations, in that order.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for zero monomers, a bad density, an empty
/// or unparsable SMILES, or any error of [`pack_chains`].
pub fn build_polymer(
    smiles: &str,
    monomers: u32,
    options: &PolymerOptions,
) -> Result<PolymerBuild, Error> {
    if monomers == 0 {
        return Err(Error::InvalidInput(
            "a polymer needs at least one repeat unit".to_owned(),
        ));
    }
    if !options.density.is_finite() || options.density <= 0.0 {
        return Err(Error::InvalidInput(
            "the polymer density must be positive and finite".to_owned(),
        ));
    }
    let graph = parse_smiles(smiles)?;
    let mut steps = vec![format!(
        "parsed {} atoms and {} bonds from the SMILES",
        graph.atoms.len(),
        graph.bonds.len()
    )];
    let monomer = embed_smiles(&graph, &EmbedOptions::default())?;
    steps.push(format!(
        "embedded one repeat unit with {} sites",
        monomer.sites().len()
    ));
    let chain = repeat_chain(&monomer, monomers)?;
    steps.push(format!(
        "replicated the repeat unit {monomers} times along a"
    ));
    let structure = pack_chains(std::slice::from_ref(&chain), options.density)?;
    steps.push(format!(
        "packed one chain at {} g/cm3 in a cubic P1 cell",
        options.density
    ));
    Ok(PolymerBuild { structure, steps })
}

/// The documented equilibration protocol, one line per stage.
///
/// The steps are data, not a simulation: they describe the usual minimise,
/// anneal, NVT, NPT and production sequence for a packed amorphous model. An
/// external molecular-dynamics engine is free to interpret them.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the embedded protocol is invalid or
/// empty.
pub fn polymer_protocol() -> Result<Vec<String>, Error> {
    let protocol: Protocol = serde_json::from_str(include_str!("../data/polymer_protocol.json"))
        .map_err(|error| Error::InvalidInput(format!("polymer protocol: {error}")))?;
    if protocol.steps.is_empty() {
        return Err(Error::InvalidInput(
            "the polymer protocol has no steps".to_owned(),
        ));
    }
    Ok(protocol
        .steps
        .into_iter()
        .map(|step| format!("{}: {}", step.name, step.details))
        .collect())
}

#[derive(Debug, Deserialize)]
struct Protocol {
    steps: Vec<ProtocolStep>,
}

#[derive(Debug, Deserialize)]
struct ProtocolStep {
    name: String,
    details: String,
}

fn cubic_side(atoms: usize, density: f64) -> Result<f64, Error> {
    let mass = atoms as f64 * MOLAR_MASS_C / AVOGADRO;
    let volume_a3 = mass / density * A3_PER_CM3;
    let side = volume_a3.cbrt();
    if !side.is_finite() || side <= 0.0 {
        return Err(Error::InvalidInput(
            "the packing cell side is not representable".to_owned(),
        ));
    }
    Ok(side)
}

fn grid_side(count: usize) -> usize {
    let mut side = 1usize;
    while side.saturating_mul(side).saturating_mul(side) < count {
        side += 1;
    }
    side
}

fn grid_offset(index: usize, grid: usize, spacing: f64) -> [f64; 3] {
    let plane = grid.saturating_mul(grid);
    [
        (index % grid) as f64 * spacing,
        ((index / grid) % grid) as f64 * spacing,
        (index / plane) as f64 * spacing,
    ]
}

fn repeat_chain(monomer: &Structure, monomers: u32) -> Result<Structure, Error> {
    let points = monomer.cart_sites();
    if points.is_empty() {
        return Err(Error::InvalidInput(
            "a repeat unit needs at least one atom".to_owned(),
        ));
    }
    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    for point in &points {
        min_x = min_x.min(point[0]);
        max_x = max_x.max(point[0]);
    }
    let shift = (max_x - min_x) + REPEAT_GAP;
    let mut atoms = Vec::with_capacity(monomer.sites().len().saturating_mul(monomers as usize));
    for copy in 0..monomers {
        let dx = f64::from(copy) * shift;
        for (site, point) in monomer.sites().iter().zip(&points) {
            atoms.push((
                site.element().to_owned(),
                [point[0] + dx, point[1], point[2]],
            ));
        }
    }
    Structure::molecule("polymer_chain", &atoms, CHAIN_PADDING)
}
