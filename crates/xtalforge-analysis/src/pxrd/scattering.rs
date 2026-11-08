//! Z-based atomic scattering factors and structure factors for PXRD.
//!
//! The approximation is documented on the parent module: `f(s) = Z exp(-b s^2)`
//! with `b = 0.5 Z^(2/3)` Angstrom squared, plus the isotropic Debye-Waller
//! factor `exp(-8 pi^2 U s^2)` when the site carries a `U_iso`.

use std::f64::consts::{PI, TAU};

use xtalforge_core::{Error, P1Cell, Structure, elements::element_by_symbol};

/// A scattering atom in the P1 cell, with its form-factor parameters.
pub(super) struct Scatterer {
    pub(super) frac: [f64; 3],
    pub(super) occupancy: f64,
    atomic_number: f64,
    /// Gaussian decay coefficient `b` from the Z-based approximation.
    gaussian_b: f64,
    /// Isotropic Debye-Waller exponent `B = 8 pi^2 U`.
    debye_b: f64,
}

/// Expands the structure to P1 and builds one scatterer per atom.
pub(super) fn build_scatterers(
    structure: &Structure,
    cell: &P1Cell,
) -> Result<(Vec<Scatterer>, f64), Error> {
    let mut scatterers = Vec::with_capacity(cell.atoms.len());
    let mut total = 0.0;
    for atom in &cell.atoms {
        let element = element_by_symbol(&atom.element)?;
        let atomic_number = f64::from(element.z);
        let gaussian_b = 0.5 * atomic_number.powf(2.0 / 3.0);
        let u_iso = structure
            .sites()
            .get(atom.site_idx)
            .and_then(xtalforge_core::Site::u_iso)
            .unwrap_or(0.0);
        total += atomic_number * atom.occupancy;
        scatterers.push(Scatterer {
            frac: atom.frac,
            occupancy: atom.occupancy,
            atomic_number,
            gaussian_b,
            debye_b: 8.0 * PI * PI * u_iso,
        });
    }
    Ok((scatterers, total))
}

/// `|F(hkl)|^2` from the P1 atoms.
pub(super) fn structure_factor(scatterers: &[Scatterer], hkl: [i32; 3], d: f64) -> f64 {
    let [h, k, l] = hkl.map(f64::from);
    let s2 = 1.0 / (4.0 * d * d);
    let mut real = 0.0;
    let mut imaginary = 0.0;
    for scatterer in scatterers {
        let [x, y, z] = scatterer.frac;
        let phase = TAU * (h * x + k * y + l * z);
        let attenuation = (-(scatterer.gaussian_b + scatterer.debye_b) * s2).exp();
        let amplitude = scatterer.atomic_number * scatterer.occupancy * attenuation;
        real += amplitude * phase.cos();
        imaginary += amplitude * phase.sin();
    }
    real * real + imaginary * imaginary
}
