//! Neutron powder diffraction (Task 19.3.1).
//!
//! # Nuclear intensities
//!
//! The coherent bound scattering lengths of
//! [`neutron_scattering_length`](crate::neutron_scattering_length) are
//! constant amplitudes, so the nuclear structure factor is
//!
//! ```text
//! F_N(hkl) = sum_j b_j exp(-W_j) exp(2 pi i (h x_j + k y_j + l z_j)),
//! ```
//!
//! with the same Debye-Waller model as the X-ray engine (isotropic `U_iso` or
//! anisotropic `U_ij`). A negative length (H, Ti, Mn, V) keeps its sign and
//! therefore changes the phase, which is the physical behaviour. The lengths
//! come from the embedded `data/neutron.json` NIST NCNR table
//! (public-domain factual values); no additional data set is bundled.
//!
//! The powder intensity is `|F_N|^2 * multiplicity * L` with the neutron
//! Lorentz factor `L = 1 / (sin^2 theta cos theta)`. The magnetic module
//! reuses the same reflection enumeration and profile sampler.
//!
//! # Axes
//!
//! * **Constant wavelength.** Positions are `2 theta` in degrees from Bragg's
//!   law and the Phase 7 [`InstrumentProfile`](crate::InstrumentProfile)
//!   supplies the width law, line shape, asymmetry and (ignored for neutrons)
//!   polarisation.
//! * **Time of flight.** One detector bank sits at a fixed `2 theta` and the
//!   wavelength follows from the flight time,
//!
//!   ```text
//!   t[us] = 252.777 L[m] lambda[Angstrom]
//!         = 505.554 L[m] d[Angstrom] sin(theta).
//!   ```
//!
//!   The line shape is the configured Phase 7 shape with a resolution
//!   `FWHM_t = (Delta d / d) * t`; a constant relative `d` resolution is the
//!   documented approximation. Caglioti/TCH widths and the FCJ asymmetry are
//!   functions of `2 theta` and are not used on the time axis. The integrated
//!   intensity is `|F_N|^2 * multiplicity * d^4`, the standard time-of-flight
//!   weighting for a constant effective incident flux.
//!
//! # Approximations (deliberate, documented here)
//!
//! * Coherent elastic nuclear scattering only: no incoherent background, no
//!   absorption, extinction or multiple scattering, no incident-spectrum
//!   weighting and no magnetic contribution (use
//!   [`magnetic_diffraction`](crate::magnetic_diffraction)).
//! * Multiplicities are the orbit sizes under the explicit operations of the
//!   structure. A P1 description lists every `(h k l)` with multiplicity one.
//! * The sampled profile is evaluated within twelve FWHM of each centre, like
//!   the X-ray profile entry point.
//!
//! The result is deterministic: the same structure and options give the same
//! output independent of platform and thread count.

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};

use crate::instrument::RadiationKind;
use crate::scattering;

mod options;
mod reflections;
mod sample;

pub use options::{NeutronAxis, NeutronGeometry, NeutronOptions};
pub(crate) use reflections::{bragg_angle, d_range_from_two_theta, reflections_in_d_range};
pub(crate) use sample::{point_count, sample_cw_peaks};
use sample::{sample_tof_peaks, sort_reflections};

/// Time of flight per unit path length and wavelength in
/// `us / (m Angstrom)`.
pub(crate) const TIME_OF_FLIGHT: f64 = 252.777;
/// Relative amplitude below which a reflection counts as absent.
const ABSENCE_TOLERANCE: f64 = 1e-8;

/// One sample of a neutron pattern; `x` is in the unit of [`NeutronAxis`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct NeutronPoint {
    /// Axis coordinate.
    pub x: f64,
    /// Intensity in arbitrary units.
    pub intensity: f64,
}

/// One allowed neutron reflection.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NeutronReflection {
    /// Miller indices of the canonical representative.
    pub hkl: [i32; 3],
    /// Interplanar spacing in Angstrom.
    pub d: f64,
    /// Bragg angle `2 theta` in degrees; for time of flight this is the fixed
    /// bank angle.
    pub two_theta: f64,
    /// Time of flight in microseconds, on the time-of-flight axis.
    pub tof_us: Option<f64>,
    /// Number of symmetry-equivalent reflections at this position.
    pub multiplicity: usize,
    /// `|F_N|^2` in squared Angstrom.
    pub nuclear_amplitude_squared: f64,
    /// Unnormalised powder intensity, including the axis weighting.
    pub intensity: f64,
}

impl NeutronReflection {
    /// `(h k l)` with a minus sign on negative indices.
    #[must_use]
    pub fn label(&self) -> String {
        let [h, k, l] = self.hkl;
        format!("({h} {k} {l})")
    }
}

/// A calculated neutron powder pattern.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NeutronPattern {
    /// Sampling axis.
    pub axis: NeutronAxis,
    /// Inclusive sampling range in the axis unit.
    pub range: [f64; 2],
    /// Sampling step in the axis unit.
    pub step: f64,
    /// Sampled profile, ordered by increasing `x`.
    pub points: Vec<NeutronPoint>,
    /// Allowed reflections, ordered by increasing position.
    pub reflections: Vec<NeutronReflection>,
}

impl NeutronPattern {
    /// Largest sampled intensity, or `0` for an empty pattern.
    #[must_use]
    pub fn max_intensity(&self) -> f64 {
        self.points
            .iter()
            .map(|point| point.intensity)
            .fold(0.0, f64::max)
    }

    /// Copy of the points scaled so that the maximum is 100.
    #[must_use]
    pub fn normalised(&self) -> Vec<NeutronPoint> {
        let maximum = self.max_intensity();
        if maximum <= 0.0 {
            return self.points.clone();
        }
        self.points
            .iter()
            .map(|point| NeutronPoint {
                x: point.x,
                intensity: point.intensity * 100.0 / maximum,
            })
            .collect()
    }
}

/// Calculates a nuclear neutron powder pattern.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, a structure without
/// scattering atoms, an element without a coherent scattering length, a cell
/// that would need too many reflections, or a pattern that would need too many
/// points.
pub fn neutron_powder(
    structure: &Structure,
    options: &NeutronOptions,
) -> Result<NeutronPattern, Error> {
    options.validate()?;
    let cell = structure.p1()?;
    let scattering = scattering::build(structure, &cell, RadiationKind::Neutron)?;
    if scattering.scatterers.is_empty() || scattering.total <= 0.0 {
        return Err(Error::InvalidInput(
            "there are no scattering atoms to calculate a neutron pattern from".to_owned(),
        ));
    }
    let profile = &options.profile;
    let absence = (ABSENCE_TOLERANCE * scattering.total).powi(2);
    match options.geometry {
        NeutronGeometry::ConstantWavelength { wavelength } => {
            let [d_min, d_max] = d_range_from_two_theta(wavelength, options.range)?;
            let data = reflections_in_d_range(structure, d_min, d_max, true)?;
            let mut reflections = Vec::with_capacity(data.len());
            for reflection in &data {
                let two_theta = bragg_angle(wavelength, reflection.d)?;
                let amplitude_squared = scattering::squared_amplitude(&scattering, reflection.hkl);
                if !amplitude_squared.is_finite() || amplitude_squared <= absence {
                    continue;
                }
                let intensity = amplitude_squared
                    * reflection.multiplicity as f64
                    * profile.lorentz_polarization(two_theta);
                if !intensity.is_finite() || intensity <= 0.0 {
                    continue;
                }
                reflections.push(NeutronReflection {
                    hkl: reflection.hkl,
                    d: reflection.d,
                    two_theta,
                    tof_us: None,
                    multiplicity: reflection.multiplicity,
                    nuclear_amplitude_squared: amplitude_squared,
                    intensity,
                });
            }
            sort_reflections(&mut reflections);
            let peaks: Vec<(f64, f64)> = reflections
                .iter()
                .map(|reflection| (reflection.two_theta, reflection.intensity))
                .collect();
            let points = sample_cw_peaks(&peaks, profile, options.range, options.step)?;
            Ok(NeutronPattern {
                axis: NeutronAxis::TwoTheta { wavelength },
                range: options.range,
                step: options.step,
                points,
                reflections,
            })
        }
        NeutronGeometry::TimeOfFlight {
            total_path_m,
            two_theta_deg,
            relative_d_fwhm,
        } => {
            let sine = (0.5 * two_theta_deg.to_radians()).sin();
            let scale = 2.0 * TIME_OF_FLIGHT * total_path_m * sine;
            if scale <= 0.0 || !scale.is_finite() {
                return Err(Error::InvalidInput(
                    "the time-of-flight calibration is not finite".to_owned(),
                ));
            }
            let d_min = options.range[0] / scale;
            let d_max = options.range[1] / scale;
            let data = reflections_in_d_range(structure, d_min, d_max, true)?;
            let mut reflections = Vec::with_capacity(data.len());
            for reflection in &data {
                let tof = scale * reflection.d;
                let amplitude_squared = scattering::squared_amplitude(&scattering, reflection.hkl);
                if !amplitude_squared.is_finite() || amplitude_squared <= absence {
                    continue;
                }
                let intensity =
                    amplitude_squared * reflection.multiplicity as f64 * reflection.d.powi(4);
                if !intensity.is_finite() || intensity <= 0.0 {
                    continue;
                }
                reflections.push(NeutronReflection {
                    hkl: reflection.hkl,
                    d: reflection.d,
                    two_theta: two_theta_deg,
                    tof_us: Some(tof),
                    multiplicity: reflection.multiplicity,
                    nuclear_amplitude_squared: amplitude_squared,
                    intensity,
                });
            }
            sort_reflections(&mut reflections);
            let points = sample_tof_peaks(
                &reflections,
                profile,
                relative_d_fwhm,
                options.range,
                options.step,
            )?;
            Ok(NeutronPattern {
                axis: NeutronAxis::TimeOfFlight {
                    total_path_m,
                    two_theta_deg,
                },
                range: options.range,
                step: options.step,
                points,
                reflections,
            })
        }
    }
}
