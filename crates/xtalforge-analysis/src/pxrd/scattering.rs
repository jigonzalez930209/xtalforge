//! Atomic scattering factors and structure factors for PXRD.
//!
//! Neutral-atom form factors use the International Tables Vol. C (1992)
//! Table 6.1.1.4 four-Gaussian fit, `f(s) = sum a_i exp(-b_i s^2) + c`, stored
//! in `data/it92.json` for H to Cf. Heavier elements fall back to the Z-based
//! Gaussian `f(s) = Z exp(-0.5 Z^(2/3) s^2)`. The isotropic Debye-Waller factor
//! `exp(-8 pi^2 U s^2)` is applied when a site carries a `U_iso`.

use std::collections::BTreeMap;
use std::f64::consts::{PI, TAU};
use std::sync::OnceLock;

use serde::Deserialize;
use xtalforge_core::{Error, P1Cell, Structure, elements::element_by_symbol};

/// The embedded IT92 coefficient table.
const IT92: &str = include_str!("../../data/it92.json");

#[derive(Deserialize)]
struct Table {
    coefficients: BTreeMap<String, [f64; 9]>,
}

fn it92_table() -> Option<&'static BTreeMap<String, [f64; 9]>> {
    static TABLE: OnceLock<Option<BTreeMap<String, [f64; 9]>>> = OnceLock::new();
    TABLE
        .get_or_init(|| {
            serde_json::from_str::<Table>(IT92)
                .ok()
                .map(|table| table.coefficients)
        })
        .as_ref()
}

/// Form-factor model of one element.
#[derive(Debug, Clone, Copy)]
enum FormFactor {
    /// IT92 `a1..a4, b1..b4, c`.
    It92([f64; 9]),
    /// Z-based Gaussian `Z exp(-b s^2)`.
    Gaussian { z: f64, b: f64 },
}

impl FormFactor {
    fn value(self, s2: f64) -> f64 {
        match self {
            Self::It92(c) => {
                c[0] * (-c[4] * s2).exp()
                    + c[1] * (-c[5] * s2).exp()
                    + c[2] * (-c[6] * s2).exp()
                    + c[3] * (-c[7] * s2).exp()
                    + c[8]
            }
            Self::Gaussian { z, b } => z * (-b * s2).exp(),
        }
    }
}

/// A scattering atom in the P1 cell, with its form-factor parameters.
pub(super) struct Scatterer {
    pub(super) frac: [f64; 3],
    pub(super) occupancy: f64,
    form: FormFactor,
    /// Isotropic Debye-Waller exponent `B = 8 pi^2 U`.
    debye_b: f64,
}

/// Expands the structure to P1 and builds one scatterer per atom.
pub(super) fn build_scatterers(
    structure: &Structure,
    cell: &P1Cell,
) -> Result<(Vec<Scatterer>, f64), Error> {
    let table = it92_table();
    let mut scatterers = Vec::with_capacity(cell.atoms.len());
    let mut total = 0.0;
    for atom in &cell.atoms {
        let element = element_by_symbol(&atom.element)?;
        let atomic_number = f64::from(element.z);
        let form = table.and_then(|table| table.get(&element.symbol)).map_or(
            FormFactor::Gaussian {
                z: atomic_number,
                b: 0.5 * atomic_number.powf(2.0 / 3.0),
            },
            |coefficients| FormFactor::It92(*coefficients),
        );
        let u_iso = structure
            .sites()
            .get(atom.site_idx)
            .and_then(xtalforge_core::Site::u_iso)
            .unwrap_or(0.0);
        total += atomic_number * atom.occupancy;
        scatterers.push(Scatterer {
            frac: atom.frac,
            occupancy: atom.occupancy,
            form,
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
        let attenuation = (-scatterer.debye_b * s2).exp();
        let amplitude = scatterer.form.value(s2) * scatterer.occupancy * attenuation;
        real += amplitude * phase.cos();
        imaginary += amplitude * phase.sin();
    }
    real * real + imaginary * imaginary
}
