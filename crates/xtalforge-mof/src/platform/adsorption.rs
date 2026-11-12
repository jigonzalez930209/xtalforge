#![allow(clippy::indexing_slicing)] // fixed 3x3 and bounded graph arrays

//! Host-guest force fields, Widom insertion, GCMC and energy grids
//! (Tasks 10.4.1-10.4.5).

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};

/// Boltzmann constant in eV/K.
pub const KB_EV: f64 = 8.617_333_262e-5;
/// Avogadro constant.
pub const AVOGADRO: f64 = 6.022_140_76e23;

const GUEST_JSON: &str = include_str!("../../data/guests.json");

#[derive(Debug, Clone, Deserialize)]
struct GuestFile {
    guests: Vec<GuestEntry>,
}

#[derive(Debug, Clone, Deserialize)]
struct GuestEntry {
    name: String,
    atoms: Vec<GuestAtomEntry>,
}

#[derive(Debug, Clone, Deserialize)]
struct GuestAtomEntry {
    element: String,
    position: [f64; 3],
    sigma: f64,
    epsilon_k: f64,
    charge: f64,
}

/// One atom of a rigid guest.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GuestAtom {
    /// Element symbol.
    pub element: String,
    /// Position relative to the guest centre in Angstrom.
    pub position: [f64; 3],
    /// Lennard-Jones diameter in Angstrom.
    pub sigma: f64,
    /// Lennard-Jones well depth in K.
    pub epsilon_k: f64,
    /// Point charge in elementary units.
    pub charge: f64,
}

/// A rigid guest molecule.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Guest {
    /// TraPPE name (`CO2`, `N2`, `CH4`, `H2`, `H2O`).
    pub name: String,
    /// Atoms with their force-field parameters.
    pub atoms: Vec<GuestAtom>,
}

/// One host atom with its force-field parameters.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HostAtom {
    /// P1 fractional position.
    pub frac: [f64; 3],
    /// Element symbol.
    pub element: String,
    /// Lennard-Jones diameter in Angstrom.
    pub sigma: f64,
    /// Lennard-Jones well depth in K.
    pub epsilon_k: f64,
    /// Framework charge in elementary units.
    pub charge: f64,
}

/// A host-guest force field.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HostGuestForceField {
    /// Framework atoms.
    pub hosts: Vec<HostAtom>,
    /// Guest models by name.
    pub guests: BTreeMap<String, Guest>,
    /// Cell matrix rows.
    pub matrix: [[f64; 3]; 3],
    /// Inverse cell matrix rows.
    pub inverse: [[f64; 3]; 3],
}

/// Lists the embedded rigid guests.
///
/// # Errors
///
/// Returns [`crate::Error::Data`] for an invalid guest table.
pub fn guests() -> Result<Vec<Guest>, Error> {
    let data: GuestFile = serde_json::from_str(GUEST_JSON)
        .map_err(|error| Error::InvalidInput(format!("invalid guest table: {error}")))?;
    Ok(data
        .guests
        .into_iter()
        .map(|entry| Guest {
            name: entry.name,
            atoms: entry
                .atoms
                .into_iter()
                .map(|atom| GuestAtom {
                    element: atom.element,
                    position: atom.position,
                    sigma: atom.sigma,
                    epsilon_k: atom.epsilon_k,
                    charge: atom.charge,
                })
                .collect(),
        })
        .collect())
}

/// Finds one embedded guest.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown name or invalid data.
pub fn find_guest(name: &str) -> Result<Guest, Error> {
    guests()?
        .into_iter()
        .find(|guest| guest.name.eq_ignore_ascii_case(name))
        .ok_or_else(|| Error::InvalidInput(format!("unknown guest {name:?}")))
}

/// Builds the UFF host parameters and EQeq charges plus the TraPPE guests.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an element without UFF parameters or a
/// failed P1 expansion.
pub fn force_field(structure: &Structure) -> Result<HostGuestForceField, Error> {
    let p1 = structure.p1()?;
    let snapshot = crate::platform::util::p1_structure(&p1)?;
    let charges = xtalforge_ff::eqeq_charges(&snapshot, &xtalforge_ff::ChargeOptions::eqeq())
        .map(|charges| charges.per_atom)
        .unwrap_or_else(|_| vec![0.0; p1.atoms.len()]);
    let mut hosts = Vec::with_capacity(p1.atoms.len());
    for (index, atom) in p1.atoms.iter().enumerate() {
        let params = xtalforge_ff::uff::params::types_for(&atom.element)
            .map_err(|error| Error::InvalidInput(error.to_string()))?;
        let Some(params) = params.first() else {
            return Err(Error::InvalidInput(format!(
                "no UFF parameters for host element {}",
                atom.element
            )));
        };
        hosts.push(HostAtom {
            frac: atom.frac,
            element: atom.element.clone(),
            sigma: params.x1,
            epsilon_k: params.d1 * 0.043_364_1 / KB_EV,
            charge: charges.get(index).copied().unwrap_or(0.0),
        });
    }
    let guest_list = guests()?;
    let guests = guest_list
        .into_iter()
        .map(|guest| (guest.name.to_ascii_uppercase(), guest))
        .collect();
    let matrix = p1.lattice.matrix();
    let inverse = p1.lattice.inverse();
    Ok(HostGuestForceField {
        hosts,
        guests,
        matrix,
        inverse,
    })
}

pub(crate) fn guest_by_name<'a>(
    force_field: &'a HostGuestForceField,
    name: &str,
) -> Result<&'a Guest, Error> {
    force_field
        .guests
        .get(&name.to_ascii_uppercase())
        .ok_or_else(|| Error::InvalidInput(format!("unknown guest {name:?}")))
}

fn minimum_image(delta: [f64; 3], matrix: [[f64; 3]; 3]) -> f64 {
    let mut best = f64::INFINITY;
    for i in -1..=1 {
        for j in -1..=1 {
            for k in -1..=1 {
                let shift = [
                    delta[0] + f64::from(i),
                    delta[1] + f64::from(j),
                    delta[2] + f64::from(k),
                ];
                let cart = [
                    shift[0] * matrix[0][0] + shift[1] * matrix[1][0] + shift[2] * matrix[2][0],
                    shift[0] * matrix[0][1] + shift[1] * matrix[1][1] + shift[2] * matrix[2][1],
                    shift[0] * matrix[0][2] + shift[1] * matrix[1][2] + shift[2] * matrix[2][2],
                ];
                let distance = (cart[0] * cart[0] + cart[1] * cart[1] + cart[2] * cart[2]).sqrt();
                best = best.min(distance);
            }
        }
    }
    best
}

fn rotation(angles: [f64; 3]) -> [[f64; 3]; 3] {
    let (sin_x, cos_x) = angles[0].sin_cos();
    let (sin_y, cos_y) = angles[1].sin_cos();
    let (sin_z, cos_z) = angles[2].sin_cos();
    let rx = [[1.0, 0.0, 0.0], [0.0, cos_x, -sin_x], [0.0, sin_x, cos_x]];
    let ry = [[cos_y, 0.0, sin_y], [0.0, 1.0, 0.0], [-sin_y, 0.0, cos_y]];
    let rz = [[cos_z, -sin_z, 0.0], [sin_z, cos_z, 0.0], [0.0, 0.0, 1.0]];
    multiply(rz, multiply(ry, rx))
}

fn multiply(first: [[f64; 3]; 3], second: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut result = [[0.0; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            for k in 0..3 {
                result[row][column] += first[row][k] * second[k][column];
            }
        }
    }
    result
}

fn apply(matrix: [[f64; 3]; 3], point: [f64; 3]) -> [f64; 3] {
    [
        matrix[0][0] * point[0] + matrix[0][1] * point[1] + matrix[0][2] * point[2],
        matrix[1][0] * point[0] + matrix[1][1] * point[1] + matrix[1][2] * point[2],
        matrix[2][0] * point[0] + matrix[2][1] * point[1] + matrix[2][2] * point[2],
    ]
}

/// Interaction energy of one rigid guest at a fractional position and
/// orientation.
#[must_use]
pub fn guest_energy(
    force_field: &HostGuestForceField,
    guest: &Guest,
    frac: [f64; 3],
    angles: [f64; 3],
) -> f64 {
    let rotation = rotation(angles);
    let centre = [
        frac[0] * force_field.matrix[0][0]
            + frac[1] * force_field.matrix[1][0]
            + frac[2] * force_field.matrix[2][0],
        frac[0] * force_field.matrix[0][1]
            + frac[1] * force_field.matrix[1][1]
            + frac[2] * force_field.matrix[2][1],
        frac[0] * force_field.matrix[0][2]
            + frac[1] * force_field.matrix[1][2]
            + frac[2] * force_field.matrix[2][2],
    ];
    let mut energy = 0.0;
    for atom in &guest.atoms {
        let offset = apply(rotation, atom.position);
        let point = [
            centre[0] + offset[0],
            centre[1] + offset[1],
            centre[2] + offset[2],
        ];
        let guest_frac = [
            point[0] * force_field.inverse[0][0]
                + point[1] * force_field.inverse[1][0]
                + point[2] * force_field.inverse[2][0],
            point[0] * force_field.inverse[0][1]
                + point[1] * force_field.inverse[1][1]
                + point[2] * force_field.inverse[2][1],
            point[0] * force_field.inverse[0][2]
                + point[1] * force_field.inverse[1][2]
                + point[2] * force_field.inverse[2][2],
        ];
        for host in &force_field.hosts {
            let delta = [
                guest_frac[0] - host.frac[0],
                guest_frac[1] - host.frac[1],
                guest_frac[2] - host.frac[2],
            ];
            let distance = minimum_image(delta, force_field.matrix);
            if distance > 15.0 {
                continue;
            }
            let sigma = 0.5 * (atom.sigma + host.sigma);
            let epsilon = (atom.epsilon_k * host.epsilon_k).sqrt();
            let ratio = sigma / distance.max(1e-6);
            let ratio6 = ratio.powi(6);
            energy += 4.0 * epsilon * (ratio6 * ratio6 - ratio6) * KB_EV;
            energy += 14.399_645_478_6 * atom.charge * host.charge / distance.max(1e-6);
        }
    }
    energy
}

/// Function field mixing rules (Lorentz-Berthelot).
#[must_use]
pub fn lennard_jones(sigma_a: f64, epsilon_a: f64, sigma_b: f64, epsilon_b: f64) -> (f64, f64) {
    (0.5 * (sigma_a + sigma_b), (epsilon_a * epsilon_b).sqrt())
}

mod mc;
pub use mc::{
    EnergyGrid, GcmcOptions, Isotherm, IsothermPoint, RaspaFiles, WidomOptions, WidomResult,
    energy_grid, gcmc_isotherm, raspa_input, widom_insertion,
};
