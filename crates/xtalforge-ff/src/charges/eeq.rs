//! EEQ2019 charge equilibration (Task 4.3.3, functional audit 2026-11-11).
//!
//! Ported from sci-form `src/charges_eeq/eeq2019.rs` (MIT, same maintainer),
//! itself a port of grimme-lab/multicharge `eeq2019` (Apache-2.0). The model
//! of Caldeweyher et al. (J. Chem. Phys. 150, 154122, 2019) solves
//!
//! ```text
//! A_ii = eta_i + sqrt(2/pi) / rad_i
//! A_ij = sum_images erf(r / sqrt(rad_i^2 + rad_j^2)) / (eps r)
//! x_i  = -chi_i + kcn_i sqrt(CN_i)
//! A q = x  with  sum q = total charge,
//! ```
//!
//! in atomic units. `CN_i` is the erf-counted coordination number with the
//! D3 covalent radii (4/3 of the Pyykko radii) and a smooth cap at 8, summed
//! over the periodic images within 25 Bohr. The Coulomb sum runs over the
//! periodic images inside `options.cutoff`, the same real-space treatment as
//! QEq, and `eps` is [`ChargeOptions::dielectric`] (1 reproduces EEQ2019).
//! The parameters cover H to Lr and live in `data/eeq2019.json`.
//!
//! Before the 2026-11-11 audit this module reused the QEq table and
//! reproduced QEq exactly at the default dielectric.

use std::f64::consts::PI;
use std::sync::OnceLock;

use serde::Deserialize;
use xtalforge_core::{Lattice, Structure, element_info};

use crate::Error;
use crate::charges::periodic::image_range;
use crate::charges::qeq::solve_constrained;
use crate::charges::{ChargeOptions, Charges, MAX_CHARGE_ATOMS};
use crate::math::{add, norm, sub, vec_mat};

/// Bohr per Angstrom.
const BOHR_PER_ANGSTROM: f64 = 1.889_726_124_565_062;
/// Steepness of the erf coordination counting function.
const CN_STEEPNESS: f64 = 7.5;
/// Smooth coordination cap.
const CN_MAX: f64 = 8.0;
/// Coordination cutoff in Bohr.
const CN_CUTOFF_BOHR: f64 = 25.0;

#[derive(Deserialize)]
struct Table {
    chi: Vec<f64>,
    eta: Vec<f64>,
    rad: Vec<f64>,
    kcnchi: Vec<f64>,
    covrad_ang: Vec<f64>,
}

fn table() -> Result<&'static Table, Error> {
    static TABLE: OnceLock<Option<Table>> = OnceLock::new();
    TABLE
        .get_or_init(|| serde_json::from_str(include_str!("../../data/eeq2019.json")).ok())
        .as_ref()
        .ok_or_else(|| Error::Data("the embedded EEQ2019 table is invalid".to_owned()))
}

/// One atom's EEQ2019 parameters, in atomic units.
#[derive(Debug, Clone, Copy)]
struct Params {
    chi: f64,
    eta: f64,
    rad: f64,
    kcnchi: f64,
    covrad_bohr: f64,
}

fn params(element: &str) -> Result<Params, Error> {
    let z = usize::from(element_info(element)?.z);
    let table = table()?;
    let get = |values: &[f64]| {
        z.checked_sub(1)
            .and_then(|index| values.get(index).copied())
            .ok_or_else(|| Error::Data(format!("EEQ2019 has no parameters for {element}")))
    };
    Ok(Params {
        chi: get(&table.chi)?,
        eta: get(&table.eta)?,
        rad: get(&table.rad)?,
        kcnchi: get(&table.kcnchi)?,
        covrad_bohr: get(&table.covrad_ang)? * BOHR_PER_ANGSTROM * 4.0 / 3.0,
    })
}

/// Abramowitz-Stegun 7.1.26 error function (the reference's accuracy class).
fn erf(x: f64) -> f64 {
    let (sign, x) = if x >= 0.0 { (1.0, x) } else { (-1.0, -x) };
    let t = 1.0 / (1.0 + 0.327_591_1 * x);
    let poly = t
        * (0.254_829_592
            + t * (-0.284_496_736
                + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
    sign * (1.0 - poly * (-x * x).exp())
}

/// Calls `visit(i, j, r_bohr)` for every pair and image within `cutoff_bohr`.
fn for_each_pair(
    lattice: &Lattice,
    positions: &[[f64; 3]],
    cutoff_bohr: f64,
    mut visit: impl FnMut(usize, usize, f64),
) {
    let range = image_range(lattice, cutoff_bohr / BOHR_PER_ANGSTROM);
    let base = lattice.matrix();
    for na in -range[0]..=range[0] {
        for nb in -range[1]..=range[1] {
            for nc in -range[2]..=range[2] {
                let offset = vec_mat([f64::from(na), f64::from(nb), f64::from(nc)], base);
                let same_cell = na == 0 && nb == 0 && nc == 0;
                for (i, first) in positions.iter().enumerate() {
                    for (j, second) in positions.iter().enumerate() {
                        if same_cell && i == j {
                            continue;
                        }
                        let r = norm(sub(add(*second, offset), *first)) * BOHR_PER_ANGSTROM;
                        if r > 1e-6 && r <= cutoff_bohr {
                            visit(i, j, r);
                        }
                    }
                }
            }
        }
    }
}

/// Solves EEQ2019 for one structure.
pub(crate) fn solve(structure: &Structure, options: &ChargeOptions) -> Result<Charges, Error> {
    let cell = structure.p1()?;
    let n = cell.n_atoms();
    if n == 0 {
        return Ok(Charges::empty());
    }
    if n > MAX_CHARGE_ATOMS {
        return Err(Error::InvalidInput(format!(
            "EEQ solves a dense {n} x {n} system; above {MAX_CHARGE_ATOMS} atoms set site charges instead"
        )));
    }
    let mut indices = Vec::new();
    let mut atoms = Vec::new();
    let mut positions = Vec::new();
    for (index, atom) in cell.atoms.iter().enumerate() {
        if element_info(&atom.element)?.is_dummy() {
            continue;
        }
        indices.push(index);
        atoms.push(params(&atom.element)?);
        positions.push(cell.lattice.to_cart(atom.frac));
    }
    let mut per_atom = vec![0.0; n];
    if indices.is_empty() {
        return Ok(Charges {
            per_atom,
            total: 0.0,
            iterations: 0,
        });
    }
    let count = atoms.len();
    let mut cn = vec![0.0; count];
    for_each_pair(&cell.lattice, &positions, CN_CUTOFF_BOHR, |i, j, r| {
        if let (Some(first), Some(second), Some(slot)) = (atoms.get(i), atoms.get(j), cn.get_mut(i))
        {
            let rc = first.covrad_bohr + second.covrad_bohr;
            *slot += 0.5 * (1.0 + erf(-CN_STEEPNESS * (r - rc) / rc));
        }
    });
    let capped = cn
        .iter()
        .map(|value| (1.0 + CN_MAX.exp()).ln() - (1.0 + (CN_MAX - value).exp()).ln());
    // QEq's solver minimises sum chi q + 1/2 q A q, so the effective
    // electronegativity is chi - kcn sqrt(CN).
    let chi: Vec<f64> = atoms
        .iter()
        .zip(capped)
        .map(|(atom, value)| atom.chi - atom.kcnchi * value.max(0.0).sqrt())
        .collect();
    let mut matrix = vec![vec![0.0; count]; count];
    for (index, atom) in atoms.iter().enumerate() {
        if let Some(slot) = matrix.get_mut(index).and_then(|row| row.get_mut(index)) {
            *slot = atom.eta + (2.0 / PI).sqrt() / atom.rad;
        }
    }
    let cutoff_bohr = options.cutoff * BOHR_PER_ANGSTROM;
    for_each_pair(&cell.lattice, &positions, cutoff_bohr, |i, j, r| {
        if let (Some(first), Some(second)) = (atoms.get(i), atoms.get(j)) {
            let width = (first.rad * first.rad + second.rad * second.rad).sqrt();
            let value = erf(r / width) / (r * options.dielectric);
            if let Some(slot) = matrix.get_mut(i).and_then(|row| row.get_mut(j)) {
                *slot += value;
            }
        }
    });
    let solution = solve_constrained(&chi, &matrix, options.total_charge)?;
    for (slot, index) in indices.into_iter().enumerate() {
        if let (Some(value), Some(target)) = (solution.get(slot), per_atom.get_mut(index)) {
            *target = *value;
        }
    }
    let total = per_atom.iter().sum();
    Ok(Charges {
        per_atom,
        total,
        iterations: 1,
    })
}
