//! Atomic scattering amplitudes, neutron scattering lengths, electron form
//! factors and Debye-Waller attenuation (Tasks 5.1, 14.1.2 and 14.1.4).
//!
//! # X-rays
//!
//! Neutral atoms use the International Tables Vol. C (1992) Table 6.1.1.4
//! four-Gaussian coefficients stored in `data/it92.json`,
//!
//! ```text
//! f(s) = sum_{i=1..4} a_i exp(-b_i s^2) + c,   s = sin(theta) / lambda.
//! ```
//!
//! Elements beyond Cf fall back to `f(s) = Z exp(-0.5 Z^(2/3) s^2)`. Ionic
//! form factors and anomalous dispersion are not modelled.
//!
//! # Neutrons
//!
//! The coherent bound scattering lengths of `data/neutron.json` are used as a
//! constant amplitude per atom, in Angstrom (the table stores femtometre and
//! the API converts with `1 fm = 1e-5 Angstrom`). A negative length (H, Ti,
//! Mn, V) is kept: it changes the phase of the structure factor and therefore
//! the intensities, which is the physical behaviour.
//!
//! # Electrons
//!
//! The Mott-Bethe relation gives the electron scattering factor from the
//! X-ray one,
//!
//! ```text
//! f_e(s) = 0.023934 (Z - f_x(s)) / s^2   [Angstrom],
//! ```
//!
//! with `s` in `1/Angstrom`. The expression is capped at
//! `s^2 = 1e-6 1/Angstrom^2` so the `s -> 0` limit stays finite. It is a
//! kinematic approximation; dynamical scattering is out of scope.
//!
//! # Debye-Waller
//!
//! An isotropic `U_iso` uses `exp(-8 pi^2 U s^2)`. An anisotropic `U_ij`
//! (`[U11, U22, U33, U12, U13, U23]`, Angstrom squared) uses
//!
//! ```text
//! exp(-2 pi^2 sum_ij h_i h_j G*_ij U_ij),
//! ```
//!
//! with `G*` the reciprocal metric tensor. Anisotropic parameters take
//! precedence over `U_iso` when both are present; for `U11 = U22 = U33 = U`
//! the anisotropic expression reduces exactly to the isotropic one.

use std::collections::BTreeMap;
use std::f64::consts::{PI, TAU};
use std::sync::OnceLock;

use serde::Deserialize;
use xtalforge_core::{Error, Lattice, P1Cell, Site, Structure, elements::element_by_symbol};

use crate::instrument::RadiationKind;

/// The embedded IT92 coefficient table.
const IT92: &str = include_str!("../data/it92.json");
/// The embedded neutron scattering-length table.
const NEUTRON: &str = include_str!("../data/neutron.json");
/// Mott-Bethe prefactor for `s` in `1/Angstrom` and `f_e` in Angstrom.
const MOTT_BETHE: f64 = 0.023_934;
/// Smallest `s^2` used by the electron form factor.
const MIN_S2: f64 = 1e-6;

#[derive(Deserialize)]
struct It92Table {
    coefficients: BTreeMap<String, [f64; 9]>,
}

#[derive(Deserialize)]
struct NeutronTable {
    lengths: BTreeMap<String, f64>,
}

fn it92_table() -> Option<&'static BTreeMap<String, [f64; 9]>> {
    static TABLE: OnceLock<Option<BTreeMap<String, [f64; 9]>>> = OnceLock::new();
    TABLE
        .get_or_init(|| {
            serde_json::from_str::<It92Table>(IT92)
                .ok()
                .map(|table| table.coefficients)
        })
        .as_ref()
}

fn neutron_table() -> Option<&'static BTreeMap<String, f64>> {
    static TABLE: OnceLock<Option<BTreeMap<String, f64>>> = OnceLock::new();
    TABLE
        .get_or_init(|| {
            serde_json::from_str::<NeutronTable>(NEUTRON)
                .ok()
                .map(|table| table.lengths)
        })
        .as_ref()
}

/// Coherent bound neutron scattering length of an element, in Angstrom.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown element or an element that is
/// not in the embedded table.
pub fn neutron_scattering_length(symbol: &str) -> Result<f64, Error> {
    let element = element_by_symbol(symbol)?;
    let table = neutron_table().ok_or_else(|| {
        Error::InvalidInput("the neutron scattering-length table could not be read".to_owned())
    })?;
    let femtometre = table.get(&element.symbol).copied().ok_or_else(|| {
        Error::InvalidInput(format!(
            "no neutron scattering length for element {}",
            element.symbol
        ))
    })?;
    Ok(femtometre * 1e-5)
}

/// Mott-Bethe electron scattering factor of an element at `s = sin(theta) /
/// lambda` in `1/Angstrom`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown element or a non-finite or
/// negative `s`.
pub fn electron_scattering_factor(symbol: &str, s: f64) -> Result<f64, Error> {
    let element = element_by_symbol(symbol)?;
    if !s.is_finite() || s < 0.0 {
        return Err(Error::InvalidInput(format!(
            "s must be non-negative and finite, got {s}"
        )));
    }
    let z = f64::from(element.z);
    let coefficients = it92_table().and_then(|table| table.get(&element.symbol).copied());
    let fallback_b = 0.5 * z.powf(2.0 / 3.0);
    Ok(electron_value(z, coefficients, fallback_b, s * s))
}

/// X-ray form factor of an element at `s = sin(theta) / lambda` in
/// `1/Angstrom`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown element or a non-finite or
/// negative `s`.
pub fn xray_form_factor(symbol: &str, s: f64) -> Result<f64, Error> {
    let element = element_by_symbol(symbol)?;
    if !s.is_finite() || s < 0.0 {
        return Err(Error::InvalidInput(format!(
            "s must be non-negative and finite, got {s}"
        )));
    }
    let z = f64::from(element.z);
    let coefficients = it92_table().and_then(|table| table.get(&element.symbol).copied());
    Ok(coefficients.map_or_else(
        || z * (-0.5 * z.powf(2.0 / 3.0) * s * s).exp(),
        |values| it92_value(values, s * s),
    ))
}

/// Form-factor model of one element.
#[derive(Debug, Clone, Copy)]
enum FormFactor {
    /// IT92 `a1..a4, b1..b4, c`.
    It92([f64; 9]),
    /// Z-based Gaussian `Z exp(-b s^2)`.
    Gaussian { z: f64, b: f64 },
    /// Constant neutron scattering length in Angstrom.
    Constant(f64),
    /// Mott-Bethe electron form factor with an optional IT92 X-ray part.
    Electron {
        z: f64,
        coefficients: Option<[f64; 9]>,
        fallback_b: f64,
    },
}

impl FormFactor {
    fn value(self, s2: f64) -> f64 {
        match self {
            Self::It92(coefficients) => it92_value(coefficients, s2),
            Self::Gaussian { z, b } => z * (-b * s2).exp(),
            Self::Constant(length) => length,
            Self::Electron {
                z,
                coefficients,
                fallback_b,
            } => electron_value(z, coefficients, fallback_b, s2),
        }
    }
}

fn it92_value(coefficients: [f64; 9], s2: f64) -> f64 {
    coefficients[0] * (-coefficients[4] * s2).exp()
        + coefficients[1] * (-coefficients[5] * s2).exp()
        + coefficients[2] * (-coefficients[6] * s2).exp()
        + coefficients[3] * (-coefficients[7] * s2).exp()
        + coefficients[8]
}

fn electron_value(z: f64, coefficients: Option<[f64; 9]>, fallback_b: f64, s2: f64) -> f64 {
    let xray = coefficients.map_or_else(
        || z * (-fallback_b * s2).exp(),
        |values| it92_value(values, s2),
    );
    MOTT_BETHE * (z - xray) / s2.max(MIN_S2)
}

/// Debye-Waller attenuation model of one atom.
#[derive(Debug, Clone, Copy)]
enum DebyeWaller {
    /// Isotropic exponent coefficient `B = 8 pi^2 U`.
    Isotropic { b: f64 },
    /// Anisotropic parameters `U11 U22 U33 U12 U13 U23` in Angstrom squared.
    Anisotropic { u: [f64; 6] },
}

/// A scattering atom in the P1 cell, with its amplitude model.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Scatterer {
    pub(crate) frac: [f64; 3],
    pub(crate) occupancy: f64,
    form: FormFactor,
    debye: DebyeWaller,
}

impl Scatterer {
    /// Debye-Waller attenuation `exp(-W)` of one reflection.
    ///
    /// Used by the magnetic diffraction module, which evaluates the nuclear
    /// amplitude itself but must apply the same displacement model.
    pub(crate) fn attenuation(&self, metric: &[[f64; 3]; 3], hkl: [i32; 3]) -> f64 {
        (-debye_exponent(&self.debye, metric, hkl)).exp()
    }
}

/// One structure expanded to P1, its reciprocal metric and the total
/// scattering power used for the absence threshold.
#[derive(Debug, Clone)]
pub(crate) struct Scattering {
    pub(crate) scatterers: Vec<Scatterer>,
    pub(crate) total: f64,
    metric: [[f64; 3]; 3],
}

/// Reciprocal metric tensor `G*_ij = a*_i . a*_j` of a lattice.
///
/// # Errors
///
/// Returns the error of [`Lattice::reciprocal`].
pub(crate) fn reciprocal_metric(lattice: &Lattice) -> Result<[[f64; 3]; 3], Error> {
    let reciprocal = lattice.reciprocal(false)?.matrix();
    let [first, second, third] = reciprocal;
    let dot = |left: [f64; 3], right: [f64; 3]| {
        left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
    };
    Ok([
        [dot(first, first), dot(first, second), dot(first, third)],
        [dot(second, first), dot(second, second), dot(second, third)],
        [dot(third, first), dot(third, second), dot(third, third)],
    ])
}

/// `|h a* + k b* + l c*|^2` from the reciprocal metric.
pub(crate) fn reciprocal_squared(metric: &[[f64; 3]; 3], hkl: [i32; 3]) -> f64 {
    let [h, k, l] = hkl.map(f64::from);
    let [[g00, g01, g02], [_, g11, g12], [_, _, g22]] = *metric;
    h * h * g00 + k * k * g11 + l * l * g22 + 2.0 * (h * k * g01 + h * l * g02 + k * l * g12)
}

/// Expands the structure to P1 and builds one scatterer per atom.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown element or an element
/// without a neutron scattering length when the probe is neutron radiation.
pub(crate) fn build(
    structure: &Structure,
    cell: &P1Cell,
    kind: RadiationKind,
) -> Result<Scattering, Error> {
    let table = it92_table();
    let metric = reciprocal_metric(structure.lattice())?;
    let mut scatterers = Vec::with_capacity(cell.atoms.len());
    let mut total = 0.0;
    for atom in &cell.atoms {
        let element = element_by_symbol(&atom.element)?;
        let z = f64::from(element.z);
        let coefficients = table.and_then(|values| values.get(&element.symbol).copied());
        let fallback_b = 0.5 * z.powf(2.0 / 3.0);
        let (form, weight) = match kind {
            RadiationKind::XRay => {
                let form = coefficients
                    .map_or(FormFactor::Gaussian { z, b: fallback_b }, FormFactor::It92);
                (form, z)
            }
            RadiationKind::Neutron => {
                let length = neutron_scattering_length(&element.symbol)?;
                (FormFactor::Constant(length), length.abs())
            }
            RadiationKind::Electron => (
                FormFactor::Electron {
                    z,
                    coefficients,
                    fallback_b,
                },
                z,
            ),
        };
        let site = structure.sites().get(atom.site_idx);
        let debye = match site.and_then(Site::u_aniso) {
            Some(u) => DebyeWaller::Anisotropic { u },
            None => DebyeWaller::Isotropic {
                b: 8.0 * PI * PI * site.and_then(Site::u_iso).unwrap_or(0.0),
            },
        };
        total += weight * atom.occupancy;
        scatterers.push(Scatterer {
            frac: atom.frac,
            occupancy: atom.occupancy,
            form,
            debye,
        });
    }
    Ok(Scattering {
        scatterers,
        total,
        metric,
    })
}

/// Complex structure factor `(real, imaginary)` of one reflection.
pub(crate) fn amplitude(scattering: &Scattering, hkl: [i32; 3]) -> (f64, f64) {
    let s2 = 0.25 * reciprocal_squared(&scattering.metric, hkl);
    let [h, k, l] = hkl.map(f64::from);
    let mut real = 0.0;
    let mut imaginary = 0.0;
    for scatterer in &scattering.scatterers {
        let [x, y, z] = scatterer.frac;
        let phase = TAU * (h * x + k * y + l * z);
        let attenuation = scatterer.attenuation(&scattering.metric, hkl);
        let value = scatterer.form.value(s2) * scatterer.occupancy * attenuation;
        real += value * phase.cos();
        imaginary += value * phase.sin();
    }
    (real, imaginary)
}

/// `|F(hkl)|^2` of one reflection.
pub(crate) fn squared_amplitude(scattering: &Scattering, hkl: [i32; 3]) -> f64 {
    let (real, imaginary) = amplitude(scattering, hkl);
    real * real + imaginary * imaginary
}

/// Debye-Waller exponent of one reflection.
fn debye_exponent(debye: &DebyeWaller, metric: &[[f64; 3]; 3], hkl: [i32; 3]) -> f64 {
    match *debye {
        DebyeWaller::Isotropic { b } => {
            let s2 = 0.25 * reciprocal_squared(metric, hkl);
            b * s2
        }
        DebyeWaller::Anisotropic { u } => {
            let [h, k, l] = hkl.map(f64::from);
            let [[g00, g01, g02], [_, g11, g12], [_, _, g22]] = *metric;
            let sum = h * h * g00 * u[0]
                + k * k * g11 * u[1]
                + l * l * g22 * u[2]
                + 2.0 * h * k * g01 * u[3]
                + 2.0 * h * l * g02 * u[4]
                + 2.0 * k * l * g12 * u[5];
            2.0 * PI * PI * sum
        }
    }
}
