//! Amorphous polymer chains: packing, build steps and equilibration protocol
//! (Tasks 6.9.1 to 6.9.4).

use serde::Deserialize;
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder, element_info};

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
/// Closest approach allowed between a chain and its periodic image, Angstrom.
const MIN_IMAGE_DISTANCE: f64 = 2.5;

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

/// Builds one periodic chain from a SMILES repeat unit.
///
/// The repeat unit is embedded with [`EmbedOptions::default`] and rotated so
/// that its first SMILES atom (the head) and its last atom (the tail) lie on
/// the `a` axis. `monomers` copies are joined head to tail, the tail of each
/// unit bonded to the head of the next at the sum of their covalent radii, and
/// the last unit bonds to the first one through the periodic boundary, so the
/// chain is infinite along `a`. The cell is orthorhombic: `a` is the chain
/// period and the square cross-section is chosen from the real atomic masses
/// so the model reaches [`PolymerOptions::density`]. When that cross-section
/// would bring the chain closer than 2.5 Angstrom to its own periodic image,
/// the cell is widened and the lower density is reported in the steps.
/// Hydrogens are not added; the SMILES atoms are the heavy-atom backbone.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for zero monomers, a bad density, an empty
/// or unparsable SMILES, or a non-representable cell.
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
    let (unit, period) = align_unit(&monomer)?;
    let chain_period = period * f64::from(monomers);
    steps.push(format!(
        "joined {monomers} repeat units head to tail along a ({chain_period:.3} A period)"
    ));
    let mut mass = 0.0;
    let mut radius: f64 = 0.0;
    for (element, point) in &unit {
        mass += element_info(element)?.mass;
        radius = radius.max(point[1].hypot(point[2]));
    }
    let mass_g = mass * f64::from(monomers) / AVOGADRO;
    let target_side = (mass_g / options.density * A3_PER_CM3 / chain_period).sqrt();
    let min_side = 2.0 * radius + MIN_IMAGE_DISTANCE;
    let side = target_side.max(min_side);
    if !side.is_finite() || side <= 0.0 || !chain_period.is_finite() {
        return Err(Error::InvalidInput(
            "the polymer cell is not representable".to_owned(),
        ));
    }
    let lattice = Lattice::new([[chain_period, 0.0, 0.0], [0.0, side, 0.0], [0.0, 0.0, side]])?;
    let mut builder = StructureBuilder::new("polymer_chain", lattice);
    let mut label = 0usize;
    for copy in 0..monomers {
        let dx = f64::from(copy) * period;
        for (element, point) in &unit {
            label += 1;
            let cart = [point[0] + dx, point[1] + side / 2.0, point[2] + side / 2.0];
            builder = builder.site(Site::new(
                format!("{element}{label}"),
                element,
                Lattice::wrap(lattice.to_frac(cart)),
            )?);
        }
    }
    let density = mass_g / (chain_period * side * side / A3_PER_CM3);
    if side > target_side {
        steps.push(format!(
            "widened the cell to {side:.3} A so the chain clears its images; density {density:.3} g/cm3 instead of {}",
            options.density
        ));
    } else {
        steps.push(format!(
            "set the {side:.3} A square cross-section for {density:.3} g/cm3"
        ));
    }
    Ok(PolymerBuild {
        structure: builder.build()?,
        steps,
    })
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

/// An element symbol and its Cartesian position.
type Atom = (String, [f64; 3]);

/// Rotates the repeat unit so its head sits at the origin and its tail on +x.
///
/// Returns the aligned atoms and the period: head-tail distance plus the
/// tail-head bond.
fn align_unit(monomer: &Structure) -> Result<(Vec<Atom>, f64), Error> {
    let points = monomer.cart_sites();
    let (Some(head), Some(tail)) = (points.first(), points.last()) else {
        return Err(Error::InvalidInput(
            "a repeat unit needs at least one atom".to_owned(),
        ));
    };
    let elements: Vec<&str> = monomer.sites().iter().map(|site| site.element()).collect();
    let head_element = elements.first().copied().unwrap_or("C");
    let tail_element = elements.last().copied().unwrap_or("C");
    let bond =
        element_info(head_element)?.covalent_radius + element_info(tail_element)?.covalent_radius;
    let axis = [tail[0] - head[0], tail[1] - head[1], tail[2] - head[2]];
    let length = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    let rotation = if length > 1e-8 {
        rotation_to_x([axis[0] / length, axis[1] / length, axis[2] / length])
    } else {
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    };
    let aligned = elements
        .iter()
        .zip(&points)
        .map(|(element, point)| {
            let d = [point[0] - head[0], point[1] - head[1], point[2] - head[2]];
            let rotated = rotation.map(|row| row[0] * d[0] + row[1] * d[1] + row[2] * d[2]);
            ((*element).to_owned(), rotated)
        })
        .collect();
    Ok((aligned, length + bond))
}

/// Rotation matrix taking the unit vector `u` onto `+x` (Rodrigues).
fn rotation_to_x(u: [f64; 3]) -> [[f64; 3]; 3] {
    // Axis k = u x x = (0, u_z, -u_y), cos = u_x, sin = |k|.
    let cosine = u[0];
    let (ky, kz) = (u[2], -u[1]);
    let sine = ky.hypot(kz);
    if sine < 1e-12 {
        return if cosine > 0.0 {
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        } else {
            [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]
        };
    }
    let (y, z) = (ky / sine, kz / sine);
    let t = 1.0 - cosine;
    [
        [cosine, -z * sine, y * sine],
        [z * sine, cosine + y * y * t, y * z * t],
        [-y * sine, y * z * t, cosine + z * z * t],
    ]
}
