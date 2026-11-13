//! Single-crystal structure factors, HKL export and electron diffraction
//! (Task 14.1.4).
//!
//! # Single-crystal structure factors
//!
//! [`structure_factors`] lists every `(h k l)` in a `2 theta` window with the
//! complex structure factor of the P1 expansion,
//!
//! ```text
//! F(hkl) = sum_j f_j exp(-W_j) exp(2 pi i (h x_j + k y_j + l z_j)),
//! ```
//!
//! including Friedel partners (`-h -k -l` is listed separately) and without a
//! multiplicity factor. X-ray, neutron and electron radiation are supported
//! through [`RadiationSource`]; anisotropic displacement parameters enter
//! through the reciprocal metric tensor (Task 14.1.2).
//!
//! # HKL export
//!
//! [`hkl_text`] writes the SHELX `HKLF 4` columns `h k l F^2 sigma(F^2)`.
//! A calculated pattern has no counting statistics, so the standard
//! uncertainty is reported as `sqrt(F^2)` (the Poisson estimate of a line with
//! that intensity). The `write_hkl` helper saves the same text to a file; it
//! is not available on WebAssembly.
//!
//! # Electron diffraction
//!
//! [`electron_diffraction`] builds a kinematic electron pattern from an
//! accelerating voltage: the relativistic wavelength comes from
//! [`electron_wavelength`](crate::electron_wavelength), intensities use the
//! Mott-Bethe electron form factor, and an optional camera length converts
//! every reflection to a ring radius `L tan(2 theta)` in millimetre.
//! Dynamical scattering, inelastic background and the thickness-dependent
//! electron Lorentz factor are out of scope.

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};

use crate::instrument::{RadiationKind, RadiationSource};
use crate::pxrd::enumerate;
use crate::scattering;

/// Number of `(h k l)` candidates above which enumeration is refused.
const MAX_REFLECTION_SAMPLES: f64 = 8_000_000.0;
/// Relative amplitude below which a reflection counts as absent.
const ABSENCE_TOLERANCE: f64 = 1e-8;

/// Options for [`structure_factors`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct StructureFactorOptions {
    /// Radiation source.
    pub radiation: RadiationSource,
    /// Inclusive `2 theta` range in degrees.
    pub two_theta_range: [f64; 2],
    /// Reflections with `|F|` below this value are dropped (`0` keeps every
    /// enumerated reflection).
    pub min_amplitude: f64,
}

impl Default for StructureFactorOptions {
    fn default() -> Self {
        Self {
            radiation: RadiationSource::default(),
            two_theta_range: [5.0, 80.0],
            min_amplitude: 0.0,
        }
    }
}

impl StructureFactorOptions {
    /// Validates the radiation, the range and the amplitude threshold.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an invalid radiation, a non-finite
    /// or reversed range, or a negative or non-finite threshold.
    pub fn validate(&self) -> Result<(), Error> {
        self.radiation.validate()?;
        let [minimum, maximum] = self.two_theta_range;
        if !minimum.is_finite() || !maximum.is_finite() || maximum <= minimum {
            return Err(Error::InvalidInput(format!(
                "2 theta range must be finite and increasing, got {minimum}..{maximum}"
            )));
        }
        if minimum < 0.0 || maximum > 180.0 {
            return Err(Error::InvalidInput(format!(
                "2 theta range must satisfy 0 <= min < max <= 180, got {minimum}..{maximum}"
            )));
        }
        if !self.min_amplitude.is_finite() || self.min_amplitude < 0.0 {
            return Err(Error::InvalidInput(format!(
                "minimum amplitude must be non-negative and finite, got {}",
                self.min_amplitude
            )));
        }
        Ok(())
    }
}

/// One single-crystal reflection with its complex structure factor.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StructureFactor {
    /// Miller indices.
    pub hkl: [i32; 3],
    /// Interplanar spacing in Angstrom.
    pub d: f64,
    /// Bragg angle `2 theta` in degrees.
    pub two_theta: f64,
    /// Real part of the structure factor.
    pub real: f64,
    /// Imaginary part of the structure factor.
    pub imaginary: f64,
    /// `|F|^2`.
    pub amplitude_squared: f64,
}

impl StructureFactor {
    /// `|F|`.
    #[must_use]
    pub fn amplitude(&self) -> f64 {
        self.amplitude_squared.max(0.0).sqrt()
    }

    /// `(h k l)` with a minus sign on negative indices.
    #[must_use]
    pub fn label(&self) -> String {
        let [h, k, l] = self.hkl;
        format!("({h} {k} {l})")
    }
}

/// Enumerates every reflection in the range with its complex structure
/// factor.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, a structure without
/// scattering atoms, an element outside the reference tables, or too many
/// `(h k l)` candidates.
pub fn structure_factors(
    structure: &Structure,
    options: &StructureFactorOptions,
) -> Result<Vec<StructureFactor>, Error> {
    options.validate()?;
    let wavelength = options.radiation.wavelength;
    let cell = structure.p1()?;
    let scattering = scattering::build(structure, &cell, options.radiation.kind)?;
    if scattering.scatterers.is_empty() || scattering.total <= 0.0 {
        return Err(Error::InvalidInput(
            "there are no scattering atoms to calculate structure factors from".to_owned(),
        ));
    }
    let lattice = structure.lattice();
    let reciprocal = lattice.reciprocal(false)?.matrix();
    let [minimum, maximum] = options.two_theta_range;
    let theta_max = (maximum.min(179.9) * 0.5).to_radians();
    let sine = theta_max.sin();
    if sine <= 0.0 {
        return Err(Error::InvalidInput(
            "the 2 theta range has no Bragg solution".to_owned(),
        ));
    }
    let d_min = wavelength / (2.0 * sine);
    let [a, b, c] = lattice.lengths();
    let h_max = (a / d_min).ceil() + 1.0;
    let k_max = (b / d_min).ceil() + 1.0;
    let l_max = (c / d_min).ceil() + 1.0;
    let samples = (2.0 * h_max + 1.0) * (2.0 * k_max + 1.0) * (2.0 * l_max + 1.0);
    if !samples.is_finite() || samples > MAX_REFLECTION_SAMPLES {
        return Err(Error::InvalidInput(format!(
            "the requested range needs {samples} hkl candidates; use a shorter range"
        )));
    }
    let h_max = h_max as i64;
    let k_max = k_max as i64;
    let l_max = l_max as i64;
    let threshold = options
        .min_amplitude
        .max(ABSENCE_TOLERANCE * scattering.total);
    let mut factors = Vec::new();
    for h in -h_max..=h_max {
        for k in -k_max..=k_max {
            for l in -l_max..=l_max {
                if h == 0 && k == 0 && l == 0 {
                    continue;
                }
                let (Ok(h), Ok(k), Ok(l)) = (i32::try_from(h), i32::try_from(k), i32::try_from(l))
                else {
                    continue;
                };
                let hkl = [h, k, l];
                let inverse_d = crate::pxrd::reciprocal_length(hkl, &reciprocal);
                if inverse_d <= 0.0 {
                    continue;
                }
                let d = 1.0 / inverse_d;
                let sine = 0.5 * wavelength * inverse_d;
                if sine >= 1.0 {
                    continue;
                }
                let two_theta = 2.0 * sine.asin().to_degrees();
                if two_theta < minimum || two_theta > maximum {
                    continue;
                }
                let (real, imaginary) = scattering::amplitude(&scattering, hkl);
                let amplitude_squared = real * real + imaginary * imaginary;
                if !amplitude_squared.is_finite() || amplitude_squared < threshold * threshold {
                    continue;
                }
                factors.push(StructureFactor {
                    hkl,
                    d,
                    two_theta,
                    real,
                    imaginary,
                    amplitude_squared,
                });
            }
        }
    }
    let key = |value: f64| (value * 1e9).round();
    factors.sort_by(|left, right| {
        key(left.two_theta)
            .total_cmp(&key(right.two_theta))
            .then_with(|| left.hkl.cmp(&right.hkl))
    });
    Ok(factors)
}

/// SHELX `HKLF 4` text: `h k l F^2 sigma(F^2)`.
///
/// The rows follow the order of `reflections`; the trailing newline is
/// included. The standard uncertainty is `sqrt(F^2)`, the Poisson estimate of
/// a calculated line.
#[must_use]
pub fn hkl_text(reflections: &[StructureFactor]) -> String {
    let mut text = String::with_capacity(reflections.len() * 40);
    for reflection in reflections {
        let sigma = reflection.amplitude_squared.max(0.0).sqrt().max(1e-6);
        text.push_str(&format!(
            "{:>4} {:>4} {:>4} {:>14.6} {:>14.6}\n",
            reflection.hkl[0],
            reflection.hkl[1],
            reflection.hkl[2],
            reflection.amplitude_squared,
            sigma
        ));
    }
    text
}

/// Writes [`hkl_text`] to a file.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the file cannot be written.
#[cfg(not(target_arch = "wasm32"))]
pub fn write_hkl(
    path: impl AsRef<std::path::Path>,
    reflections: &[StructureFactor],
) -> Result<(), Error> {
    let path = path.as_ref();
    std::fs::write(path, hkl_text(reflections)).map_err(|error| {
        Error::InvalidInput(format!("could not write {}: {error}", path.display()))
    })
}

/// Options for [`electron_diffraction`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ElectronDiffractionOptions {
    /// Accelerating voltage in kilovolt.
    pub acceleration_kv: f64,
    /// Largest `2 theta` in degrees.
    pub two_theta_max: f64,
    /// Camera length in millimetre; when present every reflection reports a
    /// ring radius.
    pub camera_length_mm: Option<f64>,
    /// Reflections whose kinematic intensity is below this value are dropped.
    pub min_intensity: f64,
}

impl Default for ElectronDiffractionOptions {
    fn default() -> Self {
        Self {
            acceleration_kv: 200.0,
            two_theta_max: 3.0,
            camera_length_mm: None,
            min_intensity: 0.0,
        }
    }
}

impl ElectronDiffractionOptions {
    /// Validates the voltage, angle, camera length and threshold.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an invalid voltage, a non-positive
    /// or non-finite maximum angle, a non-positive camera length, or a negative
    /// threshold.
    pub fn validate(&self) -> Result<(), Error> {
        crate::instrument::electron_wavelength(self.acceleration_kv)?;
        if !self.two_theta_max.is_finite() || self.two_theta_max <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "electron 2 theta maximum must be positive and finite, got {}",
                self.two_theta_max
            )));
        }
        if let Some(length) = self.camera_length_mm
            && (!length.is_finite() || length <= 0.0)
        {
            return Err(Error::InvalidInput(format!(
                "camera length must be positive and finite, got {length}"
            )));
        }
        if !self.min_intensity.is_finite() || self.min_intensity < 0.0 {
            return Err(Error::InvalidInput(
                "minimum intensity must be non-negative and finite".to_owned(),
            ));
        }
        Ok(())
    }
}

/// One electron diffraction reflection.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ElectronReflection {
    /// Miller indices of the canonical representative.
    pub hkl: [i32; 3],
    /// Interplanar spacing in Angstrom.
    pub d: f64,
    /// Scattering vector length `1 / d` in inverse Angstrom.
    pub g: f64,
    /// Bragg angle in degrees two-theta.
    pub two_theta: f64,
    /// Ring radius `L tan(2 theta)` in millimetre, when a camera length is set.
    pub ring_radius_mm: Option<f64>,
    /// Number of symmetry-equivalent reflections.
    pub multiplicity: usize,
    /// Kinematic intensity `|F|^2 * multiplicity * LP`.
    pub intensity: f64,
}

/// A kinematic electron diffraction pattern.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ElectronPattern {
    /// Accelerating voltage in kilovolt.
    pub acceleration_kv: f64,
    /// Relativistic electron wavelength in Angstrom.
    pub wavelength: f64,
    /// Camera length in millimetre, when set.
    pub camera_length_mm: Option<f64>,
    /// Reflections, ordered by increasing angle.
    pub reflections: Vec<ElectronReflection>,
}

/// Builds a kinematic electron diffraction pattern.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, a structure without
/// scattering atoms, an element outside the reference tables, too many
/// reflections, or a voltage whose wavelength is not finite.
pub fn electron_diffraction(
    structure: &Structure,
    options: &ElectronDiffractionOptions,
) -> Result<ElectronPattern, Error> {
    options.validate()?;
    let radiation = RadiationSource::electron(options.acceleration_kv)?;
    let reflections = enumerate(
        structure,
        radiation.wavelength,
        RadiationKind::Electron,
        [0.0, options.two_theta_max],
        0.5,
        true,
    )?;
    let threshold = options.min_intensity;
    let camera = options.camera_length_mm;
    let mapped = reflections
        .into_iter()
        .filter(|reflection| reflection.intensity >= threshold)
        .map(|reflection| ElectronReflection {
            hkl: reflection.hkl,
            d: reflection.d,
            g: 1.0 / reflection.d,
            two_theta: reflection.two_theta,
            ring_radius_mm: camera.map(|length| length * reflection.two_theta.to_radians().tan()),
            multiplicity: reflection.multiplicity,
            intensity: reflection.intensity,
        })
        .collect();
    Ok(ElectronPattern {
        acceleration_kv: options.acceleration_kv,
        wavelength: radiation.wavelength,
        camera_length_mm: camera,
        reflections: mapped,
    })
}
