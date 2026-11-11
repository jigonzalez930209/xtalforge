//! Calculated powder X-ray diffraction (Task 5.1).
//!
//! # What this module computes
//!
//! Every reflection `(h k l)` whose `2 theta = 2 asin(lambda / (2 d))` lies in
//! the requested range is enumerated from the cell. The structure factor is
//! built from the P1 expansion of the structure:
//!
//! ```text
//! F(hkl) = sum_j f_j(s) exp(2 pi i (h x_j + k y_j + l z_j))
//! ```
//!
//! and the powder intensity is
//!
//! ```text
//! I = |F|^2 * multiplicity * LP(2 theta)
//! LP = (1 + cos^2(2 theta)) / (sin^2(theta) cos(theta)).
//! ```
//!
//! # Approximations (deliberate, documented here)
//!
//! * **Atomic form factor.** Neutral atoms use the International Tables
//!   Vol. C (1992) Table 6.1.1.4 four-Gaussian coefficients,
//!
//!   ```text
//!   f(s) = sum_{i=1..4} a_i exp(-b_i s^2) + c,   s = sin(theta) / lambda,
//!   ```
//!
//!   embedded from `data/it92.json` for H to Cf. Elements beyond Cf fall
//!   back to the Gaussian `f(s) = Z exp(-0.5 Z^(2/3) s^2)`. Ionic form
//!   factors are not used; formal charges do not change `f`.
//! * **Anomalous scattering** is ignored (`f' = f'' = 0`), so Friedel pairs
//!   have equal intensity.
//! * **Debye-Waller.** An isotropic `U_iso` present on a site is applied as
//!   `exp(-8 pi^2 U s^2)`. Anisotropic parameters are ignored.
//! * **Profile.** A Gaussian of the requested FWHM is used, with unit area.
//!   A pseudo-Voigt/Caglioti width is not implemented yet.
//! * **Multiplicity.** When the structure carries explicit symmetry
//!   operations, reflections are grouped into orbits (Friedel pairs merged)
//!   and only the canonical representative is listed, with the orbit size as
//!   its multiplicity. A P1 structure has multiplicity 1 and lists each
//!   `(h k l)` separately.
//!
//! The result is deterministic: the same structure and options always produce
//! bit-identical output.

mod profile;
mod scattering;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};

use profile::{lorentz_polarization, sample_pattern};
use scattering::{build_scatterers, structure_factor};

/// Number of samples above which the requested pattern is refused.
const MAX_PATTERN_POINTS: usize = 2_000_000;
/// Number of `(h k l)` candidates above which enumeration is refused.
const MAX_REFLECTION_SAMPLES: f64 = 8_000_000.0;
/// Relative amplitude below which a structure factor counts as absent.
const ABSENCE_TOLERANCE: f64 = 1e-8;

/// X-ray wavelengths of common laboratory sources, in Angstrom.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum Radiation {
    /// Copper K alpha 1, `1.540598 Angstrom`.
    CuKa1,
    /// Copper K alpha average, `1.54184 Angstrom`.
    CuKa,
    /// Molybdenum K alpha, `0.71073 Angstrom`.
    MoKa,
    /// Cobalt K alpha, `1.78897 Angstrom`.
    CoKa,
    /// Iron K alpha, `1.93604 Angstrom`.
    FeKa,
    /// Chromium K alpha, `2.28970 Angstrom`.
    CrKa,
    /// Silver K alpha, `0.55941 Angstrom`.
    AgKa,
    /// An arbitrary wavelength in Angstrom, validated when used.
    Custom(f64),
}

impl Radiation {
    /// Wavelength in Angstrom.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when a custom wavelength is not
    /// positive and finite.
    pub fn wavelength(self) -> Result<f64, Error> {
        let value = match self {
            Self::CuKa1 => 1.540_598,
            Self::CuKa => 1.541_84,
            Self::MoKa => 0.710_73,
            Self::CoKa => 1.788_97,
            Self::FeKa => 1.936_04,
            Self::CrKa => 2.289_70,
            Self::AgKa => 0.559_41,
            Self::Custom(value) => value,
        };
        if !value.is_finite() || value <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "radiation wavelength must be positive and finite, got {value}"
            )));
        }
        Ok(value)
    }

    /// Short human-readable label, including the wavelength.
    #[must_use]
    pub fn label(self) -> String {
        match self {
            Self::CuKa1 => "Cu Ka1 (1.540598 A)".to_owned(),
            Self::CuKa => "Cu Ka (1.54184 A)".to_owned(),
            Self::MoKa => "Mo Ka (0.71073 A)".to_owned(),
            Self::CoKa => "Co Ka (1.78897 A)".to_owned(),
            Self::FeKa => "Fe Ka (1.93604 A)".to_owned(),
            Self::CrKa => "Cr Ka (2.28970 A)".to_owned(),
            Self::AgKa => "Ag Ka (0.55941 A)".to_owned(),
            Self::Custom(value) => format!("custom ({value:.5} A)"),
        }
    }
}

/// Options for [`calculate_pxrd`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PxrdOptions {
    /// Radiation source.
    pub radiation: Radiation,
    /// Inclusive `2 theta` range in degrees.
    pub two_theta_range: [f64; 2],
    /// Sampling step in degrees.
    pub step: f64,
    /// Gaussian full width at half maximum in degrees.
    pub fwhm: f64,
}

impl Default for PxrdOptions {
    fn default() -> Self {
        Self {
            radiation: Radiation::CuKa,
            two_theta_range: [5.0, 80.0],
            step: 0.02,
            fwhm: 0.1,
        }
    }
}

impl PxrdOptions {
    /// Copy with a different radiation.
    #[must_use]
    pub fn with_radiation(mut self, radiation: Radiation) -> Self {
        self.radiation = radiation;
        self
    }

    /// Copy with a different `2 theta` range.
    #[must_use]
    pub fn with_range(mut self, minimum: f64, maximum: f64) -> Self {
        self.two_theta_range = [minimum, maximum];
        self
    }

    /// Copy with a different sampling step.
    #[must_use]
    pub fn with_step(mut self, step: f64) -> Self {
        self.step = step;
        self
    }

    /// Copy with a different Gaussian FWHM.
    #[must_use]
    pub fn with_fwhm(mut self, fwhm: f64) -> Self {
        self.fwhm = fwhm;
        self
    }

    fn validate(&self) -> Result<f64, Error> {
        let wavelength = self.radiation.wavelength()?;
        let [minimum, maximum] = self.two_theta_range;
        if !minimum.is_finite() || !maximum.is_finite() {
            return Err(Error::InvalidInput(
                "2 theta range must be finite".to_owned(),
            ));
        }
        if minimum < 0.0 || maximum <= minimum || maximum > 180.0 {
            return Err(Error::InvalidInput(format!(
                "2 theta range must satisfy 0 <= min < max <= 180, got {minimum}..{maximum}"
            )));
        }
        if !self.step.is_finite() || self.step <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "pattern step must be positive and finite, got {}",
                self.step
            )));
        }
        if !self.fwhm.is_finite() || self.fwhm <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "peak FWHM must be positive and finite, got {}",
                self.fwhm
            )));
        }
        let points = ((maximum - minimum) / self.step).floor() + 1.0;
        if !points.is_finite() || points > MAX_PATTERN_POINTS as f64 {
            return Err(Error::InvalidInput(format!(
                "pattern would have {points} points; use a larger step or a shorter range"
            )));
        }
        Ok(wavelength)
    }
}

/// One calculated reflection.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Reflection {
    /// Miller indices of the canonical representative.
    pub hkl: [i32; 3],
    /// Interplanar spacing in Angstrom.
    pub d: f64,
    /// Bragg angle `2 theta` in degrees.
    pub two_theta: f64,
    /// Number of symmetry-equivalent reflections at this position.
    pub multiplicity: usize,
    /// Unnormalised powder intensity `|F|^2 * multiplicity * LP`.
    pub intensity: f64,
}

impl Reflection {
    /// `(h k l)` with a minus sign on negative indices.
    #[must_use]
    pub fn label(&self) -> String {
        let [h, k, l] = self.hkl;
        format!("({h} {k} {l})")
    }
}

/// One sample of the continuous pattern.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PatternPoint {
    /// `2 theta` in degrees.
    pub two_theta: f64,
    /// Intensity in arbitrary units.
    pub intensity: f64,
}

/// A calculated powder pattern and its reflection list.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Pattern {
    /// Sampled profile, ordered by increasing `2 theta`.
    pub points: Vec<PatternPoint>,
    /// Allowed reflections, ordered by `2 theta`.
    pub reflections: Vec<Reflection>,
}

impl Pattern {
    /// Largest sampled intensity, or `0` for an empty pattern.
    #[must_use]
    pub fn max_intensity(&self) -> f64 {
        self.points
            .iter()
            .map(|point| point.intensity)
            .fold(0.0, f64::max)
    }

    /// Sampled profile scaled so that its maximum is 100.
    #[must_use]
    pub fn normalised(&self) -> Vec<PatternPoint> {
        let maximum = self.max_intensity();
        if maximum <= 0.0 {
            return self.points.clone();
        }
        self.points
            .iter()
            .map(|point| PatternPoint {
                two_theta: point.two_theta,
                intensity: point.intensity * 100.0 / maximum,
            })
            .collect()
    }
}

/// Calculates the powder pattern of a structure.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, a structure without
/// scattering atoms, a cell that would require too many reflections, or an
/// element outside the reference table.
pub fn calculate_pxrd(structure: &Structure, options: &PxrdOptions) -> Result<Pattern, Error> {
    let wavelength = options.validate()?;
    let lattice = structure.lattice();
    let cell = structure.p1()?;
    let (scatterers, total_scattering) = build_scatterers(structure, &cell)?;
    if scatterers.is_empty() || total_scattering <= 0.0 {
        return Err(Error::InvalidInput(
            "there are no scattering atoms to calculate a pattern from".to_owned(),
        ));
    }

    let [minimum, maximum] = options.two_theta_range;
    let theta_max = (maximum.min(179.9) * 0.5).to_radians();
    let d_min = wavelength / (2.0 * theta_max.sin());
    let reciprocal = lattice.reciprocal(false)?.matrix();
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

    let inverse_rotations = inverse_rotations(structure)?;
    let merges_symmetry = !inverse_rotations.is_empty();
    let mut reflections: Vec<Reflection> = Vec::new();

    for h in -h_max..=h_max {
        for k in -k_max..=k_max {
            for l in -l_max..=l_max {
                if h == 0 && k == 0 && l == 0 {
                    continue;
                }
                let Ok(h) = i32::try_from(h) else {
                    continue;
                };
                let Ok(k) = i32::try_from(k) else {
                    continue;
                };
                let Ok(l) = i32::try_from(l) else {
                    continue;
                };
                let hkl = [h, k, l];
                let inverse_d = reciprocal_length(hkl, reciprocal);
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

                let (multiplicity, representative) = if merges_symmetry {
                    let orbit = reflection_orbit(hkl, &inverse_rotations);
                    let representative = canonical_reflection(&orbit);
                    if representative != hkl {
                        continue;
                    }
                    (orbit.len(), representative)
                } else {
                    (1, hkl)
                };

                let squared = structure_factor(&scatterers, hkl, d);
                let cutoff = ABSENCE_TOLERANCE * total_scattering;
                if !squared.is_finite() || squared <= cutoff * cutoff {
                    continue;
                }
                let intensity = squared * multiplicity as f64 * lorentz_polarization(two_theta);
                if intensity <= 0.0 || !intensity.is_finite() {
                    continue;
                }
                reflections.push(Reflection {
                    hkl: representative,
                    d,
                    two_theta,
                    multiplicity,
                    intensity,
                });
            }
        }
    }

    // Equivalent reflections share a two-theta up to the last bits, which
    // differ between native and WebAssembly maths; compare on a 1e-9 degree
    // grid so the hkl tie-break gives the same order on every platform.
    let key = |value: f64| (value * 1e9).round();
    reflections.sort_by(|left, right| {
        key(left.two_theta)
            .total_cmp(&key(right.two_theta))
            .then_with(|| left.hkl.cmp(&right.hkl))
    });

    let points = sample_pattern(&reflections, options);
    Ok(Pattern {
        points,
        reflections,
    })
}

fn reciprocal_length(hkl: [i32; 3], reciprocal: [[f64; 3]; 3]) -> f64 {
    let [h, k, l] = hkl.map(f64::from);
    let vector = [
        h * reciprocal[0][0] + k * reciprocal[1][0] + l * reciprocal[2][0],
        h * reciprocal[0][1] + k * reciprocal[1][1] + l * reciprocal[2][1],
        h * reciprocal[0][2] + k * reciprocal[1][2] + l * reciprocal[2][2],
    ];
    (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]).sqrt()
}

fn inverse_rotations(structure: &Structure) -> Result<Vec<[[i32; 3]; 3]>, Error> {
    if structure.operations().is_empty() {
        return Ok(Vec::new());
    }
    structure
        .operations()
        .iter()
        .map(|operation| operation.inverse().map(|inverse| inverse.rotation()))
        .collect()
}

/// Symmetry-equivalent reflections, with Friedel partners, sorted and unique.
fn reflection_orbit(hkl: [i32; 3], inverse_rotations: &[[[i32; 3]; 3]]) -> Vec<[i32; 3]> {
    let mut orbit = Vec::with_capacity(2 * inverse_rotations.len() + 2);
    orbit.push(hkl);
    orbit.push(hkl.map(|value| -value));
    for rotation in inverse_rotations {
        let [h, k, l] = hkl;
        let transformed = [
            h * rotation[0][0] + k * rotation[1][0] + l * rotation[2][0],
            h * rotation[0][1] + k * rotation[1][1] + l * rotation[2][1],
            h * rotation[0][2] + k * rotation[1][2] + l * rotation[2][2],
        ];
        orbit.push(transformed);
        orbit.push(transformed.map(|value| -value));
    }
    orbit.sort_unstable();
    orbit.dedup();
    orbit
}

/// Canonical member of a reflection orbit: fewest minus signs, then the
/// negative indices sorted, which gives `h >= k >= l >= 0` where the group
/// permits it.
fn canonical_reflection(orbit: &[[i32; 3]]) -> [i32; 3] {
    orbit
        .iter()
        .copied()
        .min_by_key(|hkl| {
            (
                hkl.iter().filter(|value| **value < 0).count(),
                hkl.map(|value| -value),
            )
        })
        .unwrap_or([0, 0, 0])
}
