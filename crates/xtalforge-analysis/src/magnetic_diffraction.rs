//! Magnetic neutron powder diffraction (Task 19.3.2).
//!
//! # Magnetic structure factor
//!
//! A moment `m_j` in Bohr magnetons contributes the magnetic scattering length
//!
//! ```text
//! p_j = (gamma r_e / 2) m_j f_m(s) = 2.695e-5 A/ mu_B * m_j * f_m(s),
//! ```
//!
//! with the vector interaction amplitude `p = m - (m . ghat) ghat`. The
//! structure factor is the vector sum
//!
//! ```text
//! F_M(hkl) = sum_j m_j f_m,j(s) exp(-W_j) exp(2 pi i (h x_j + k y_j + l z_j)).
//! ```
//!
//! With the default powder domain average the intensity uses
//! `|F_M|^2 = (2/3) (|F_x|^2 + |F_y|^2 + |F_z|^2)`; a fixed magnetic domain
//! projects the vector `F_M` perpendicular to `ghat` instead.
//!
//! The powder intensity is `|F_M|^2 * multiplicity * L`, the same Lorentz
//! factor and Phase 7 profile as the nuclear neutron engine, so a magnetic
//! pattern can be added to a [`NeutronPattern`] with
//! [`MagneticPattern::with_nuclear`].
//!
//! # Magnetic form factor
//!
//! No third-party magnetic coefficient table is bundled. The module uses the
//! documented approximation
//!
//! ```text
//! f_m(s) = f_x(s) / f_x(0),
//! ```
//!
//! the normalized free-atom X-ray form factor of the embedded International
//! Tables Vol. C data. It is the spin-only `<j0>` shape under the assumption
//! that the spin density follows the total electron density; the orbital
//! `<j2>` term, the anisotropic tensor and the exact 3d radial density are
//! not modelled. It reproduces `f_m(0) = 1` and decays monotonically, which is
//! enough to compare peak positions and relative intensities, but it is not a
//! quantitative replacement for measured ion form factors.
//!
//! # Moments
//!
//! Moments come from the typed site field
//! [`Site::magnetic_moment`](xtalforge_core::Site::magnetic_moment)
//! (`MagneticMoment` in Bohr magnetons, magCIF crystal-axis components). The
//! components are projections onto the unit vectors along `a`, `b` and `c`;
//! [`magnetic_moment_cartesian`] converts them to the Cartesian frame used by
//! the scattering vector. A site without a moment is a missing moment: the
//! default is a typed error that names the site;
//! [`MagneticOptions::allow_missing_moments`] turns it into a zero moment and
//! lists the labels in [`MagneticPattern::missing_moments`].
//!
//! The historical `magmom` site property is not read: `MagneticMoment` is the
//! first-class field introduced by Task 19.1.1.
//!
//! # Domain average
//!
//! With [`MagneticOptions::domain_average`] (the default) the moments are
//! treated as an isotropic powder of magnetic domains, so the intensity carries
//! the standard `2/3` factor. With it disabled the actual perpendicular
//! component of the vector structure factor is used for the fixed magnetic
//! domain of the cell.
//!
//! # Multiplicity and scale
//!
//! Reflection multiplicities are the orbit sizes under the explicit
//! `Structure::operations()`. A magnetic powder should be described by its
//! magnetic cell and magnetic operations (Task 19.1.2); a P1 description lists
//! every `(h k l)` with multiplicity one. Intensities are unnormalised: there
//! is no absorption, extinction, incident-spectrum or absolute calibration
//! factor.

use std::f64::consts::TAU;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Lattice, MagneticMoment, Site, Structure};

use crate::instrument::{InstrumentProfile, RadiationKind, RadiationSource};
use crate::neutron::{
    NeutronAxis, NeutronPattern, NeutronPoint, d_range_from_two_theta, point_count,
    reflections_in_d_range, sample_cw_peaks,
};
use crate::scattering;

/// Magnetic scattering length per Bohr magneton in Angstrom.
const MAGNETIC_PREFACTOR: f64 = 2.695e-5;
/// Relative squared-amplitude threshold below which a line is absent.
const ABSENCE_TOLERANCE: f64 = 1e-12;

/// Converts a magnetic moment to Cartesian components in Bohr magnetons.
///
/// [`MagneticMoment::components`] are projections onto the unit vectors along
/// the direct lattice vectors `a`, `b` and `c`; this function multiplies them
/// by those unit vectors.
#[must_use]
pub fn magnetic_moment_cartesian(lattice: &Lattice, moment: &MagneticMoment) -> [f64; 3] {
    let [ma, mb, mc] = moment.components();
    let [a, b, c] = lattice.matrix();
    let [length_a, length_b, length_c] = lattice.lengths();
    let [ax, ay, az] = a;
    let [bx, by, bz] = b;
    let [cx, cy, cz] = c;
    [
        ma * ax / length_a + mb * bx / length_b + mc * cx / length_c,
        ma * ay / length_a + mb * by / length_b + mc * cy / length_c,
        ma * az / length_a + mb * bz / length_b + mc * cz / length_c,
    ]
}

/// Returns the moment of a site, when present.
///
/// This is a thin re-export of [`Site::magnetic_moment`] that keeps the
/// diffraction entry points independent of the `Site` construction API.
#[must_use]
pub fn site_moment(site: &Site) -> Option<MagneticMoment> {
    site.magnetic_moment()
}

/// Magnetic form factor at `s = sin(theta) / lambda` in inverse Angstrom.
///
/// The approximation is the normalized free-atom X-ray form factor,
/// `f_x(s) / f_x(0)`, clipped to `[0, 1]`; see the module documentation.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown element, a non-finite or
/// negative `s`, or a form factor that cannot be normalized.
pub fn magnetic_form_factor(symbol: &str, s: f64) -> Result<f64, Error> {
    let zero = crate::scattering::xray_form_factor(symbol, 0.0)?;
    let value = crate::scattering::xray_form_factor(symbol, s)?;
    if !zero.is_finite() || zero <= 0.0 || !value.is_finite() {
        return Err(Error::InvalidInput(format!(
            "the magnetic form factor of {symbol} is not finite"
        )));
    }
    Ok((value / zero).clamp(0.0, 1.0))
}

/// Options for [`magnetic_powder`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MagneticOptions {
    /// Phase 7 instrument profile with neutron radiation.
    pub profile: InstrumentProfile,
    /// Inclusive `2 theta` range in degrees.
    pub two_theta_range: [f64; 2],
    /// Sampling step in degrees.
    pub step: f64,
    /// Apply the `2/3` powder average over magnetic domains.
    pub domain_average: bool,
    /// Treat a site without a moment as zero instead of failing.
    pub allow_missing_moments: bool,
}

impl Default for MagneticOptions {
    fn default() -> Self {
        Self {
            profile: InstrumentProfile {
                radiation: RadiationSource {
                    kind: RadiationKind::Neutron,
                    wavelength: 1.54,
                },
                fwhm: 0.2,
                ..InstrumentProfile::default()
            },
            two_theta_range: [5.0, 120.0],
            step: 0.05,
            domain_average: true,
            allow_missing_moments: false,
        }
    }
}

impl MagneticOptions {
    /// Validates the profile, angle range and step.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-neutron profile, an invalid
    /// wavelength, a reversed or out-of-range window, a non-positive step or a
    /// pattern with more than two million points.
    pub fn validate(&self) -> Result<(), Error> {
        self.profile.validate()?;
        if self.profile.radiation.kind != RadiationKind::Neutron {
            return Err(Error::InvalidInput(
                "magnetic neutron diffraction needs a neutron instrument profile".to_owned(),
            ));
        }
        d_range_from_two_theta(self.profile.radiation.wavelength, self.two_theta_range)?;
        let count = point_count(self.two_theta_range, self.step)?;
        if count > 2_000_000 {
            return Err(Error::InvalidInput(format!(
                "magnetic pattern would have {count} points; use a larger step"
            )));
        }
        Ok(())
    }
}

/// One magnetic site and its moment.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MagneticSite {
    /// Site label.
    pub label: String,
    /// Element symbol.
    pub element: String,
    /// Moment in Bohr magnetons as magCIF crystal-axis components along
    /// `a`, `b` and `c`.
    pub moment: [f64; 3],
}

/// One magnetic reflection.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MagneticReflection {
    /// Miller indices of the canonical representative.
    pub hkl: [i32; 3],
    /// Interplanar spacing in Angstrom.
    pub d: f64,
    /// Bragg angle `2 theta` in degrees.
    pub two_theta: f64,
    /// Number of symmetry-equivalent reflections at this position.
    pub multiplicity: usize,
    /// `|F_M|^2` in squared Angstrom, domain factor included.
    pub magnetic_amplitude_squared: f64,
    /// Unnormalised powder intensity, including the Lorentz factor.
    pub intensity: f64,
}

impl MagneticReflection {
    /// `(h k l)` with a minus sign on negative indices.
    #[must_use]
    pub fn label(&self) -> String {
        let [h, k, l] = self.hkl;
        format!("({h} {k} {l})")
    }
}

/// A calculated magnetic neutron powder pattern.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MagneticPattern {
    /// Wavelength in Angstrom.
    pub wavelength: f64,
    /// Inclusive `2 theta` range in degrees.
    pub two_theta_range: [f64; 2],
    /// Sampling step in degrees.
    pub step: f64,
    /// Sampled magnetic profile, ordered by increasing `2 theta`.
    pub points: Vec<NeutronPoint>,
    /// Magnetic reflections with non-zero intensity, ordered by angle.
    pub reflections: Vec<MagneticReflection>,
    /// Sites that carry a moment, including explicit zero moments.
    pub magnetic_sites: Vec<MagneticSite>,
    /// Labels of sites treated as zero because a moment was missing and
    /// [`MagneticOptions::allow_missing_moments`] was set.
    pub missing_moments: Vec<String>,
}

impl MagneticPattern {
    /// Largest sampled intensity, or `0` for an empty pattern.
    #[must_use]
    pub fn max_intensity(&self) -> f64 {
        self.points
            .iter()
            .map(|point| point.intensity)
            .fold(0.0, f64::max)
    }

    /// Adds this magnetic pattern to a constant-wavelength nuclear pattern.
    ///
    /// The two patterns must share the wavelength, range and sampling grid;
    /// the result is the point-by-point sum, ready to plot as the total
    /// scattering.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the nuclear pattern is not a
    /// constant-wavelength pattern at the same wavelength, or its grid differs.
    pub fn with_nuclear(&self, nuclear: &NeutronPattern) -> Result<Vec<NeutronPoint>, Error> {
        let NeutronAxis::TwoTheta { wavelength } = nuclear.axis else {
            return Err(Error::InvalidInput(
                "a magnetic pattern can only be added to a constant-wavelength pattern".to_owned(),
            ));
        };
        if (wavelength - self.wavelength).abs() > 1e-9 {
            return Err(Error::InvalidInput(format!(
                "wavelength mismatch: magnetic {:.6} A, nuclear {wavelength:.6} A",
                self.wavelength
            )));
        }
        if nuclear.points.len() != self.points.len() {
            return Err(Error::InvalidInput(
                "the nuclear and magnetic grids have different lengths".to_owned(),
            ));
        }
        let mut total = Vec::with_capacity(self.points.len());
        for (nuclear_point, magnetic_point) in nuclear.points.iter().zip(&self.points) {
            if (nuclear_point.x - magnetic_point.x).abs() > 1e-9 {
                return Err(Error::InvalidInput(
                    "the nuclear and magnetic grids do not share their x values".to_owned(),
                ));
            }
            total.push(NeutronPoint {
                x: nuclear_point.x,
                intensity: nuclear_point.intensity + magnetic_point.intensity,
            });
        }
        Ok(total)
    }
}

/// Calculates a magnetic neutron powder pattern from the site moments.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, a structure without
/// scattering atoms, a site without a moment when missing moments are not
/// allowed, malformed moment properties, an element outside the reference
/// tables, or a reflection window that is too large.
pub fn magnetic_powder(
    structure: &Structure,
    options: &MagneticOptions,
) -> Result<MagneticPattern, Error> {
    options.validate()?;
    let wavelength = options.profile.radiation.wavelength;
    let cell = structure.p1()?;
    let scattering = scattering::build(structure, &cell, RadiationKind::XRay)?;
    if scattering.scatterers.is_empty() {
        return Err(Error::InvalidInput(
            "there are no atoms to calculate a magnetic pattern from".to_owned(),
        ));
    }
    let metric = scattering::reciprocal_metric(structure.lattice())?;
    let mut site_moments = Vec::with_capacity(structure.sites().len());
    let mut magnetic_sites = Vec::new();
    let mut missing_moments = Vec::new();
    for site in structure.sites() {
        match site.magnetic_moment() {
            Some(moment) => {
                site_moments.push(Some(magnetic_moment_cartesian(
                    structure.lattice(),
                    &moment,
                )));
                magnetic_sites.push(MagneticSite {
                    label: site.label().to_owned(),
                    element: site.element().to_owned(),
                    moment: moment.components(),
                });
            }
            None => {
                if options.allow_missing_moments {
                    missing_moments.push(site.label().to_owned());
                    site_moments.push(None);
                } else {
                    return Err(Error::InvalidInput(format!(
                        "site {} has no magnetic moment; set one with \
                         Site::with_magnetic_moment, or allow missing moments",
                        site.label()
                    )));
                }
            }
        }
    }
    let total_moment: f64 = cell
        .atoms
        .iter()
        .map(|atom| {
            site_moments
                .get(atom.site_idx)
                .and_then(|moment| *moment)
                .map_or(0.0, norm)
                * atom.occupancy
        })
        .sum();
    let threshold = (MAGNETIC_PREFACTOR * total_moment).powi(2) * ABSENCE_TOLERANCE;
    let [d_min, d_max] = d_range_from_two_theta(wavelength, options.two_theta_range)?;
    let data = reflections_in_d_range(structure, d_min, d_max, true)?;
    let profile = &options.profile;
    let mut reflections = Vec::new();
    for reflection in &data {
        let hkl = reflection.hkl;
        let normal = structure.lattice().plane_normal(hkl.map(f64::from))?;
        let [h, k, l] = hkl.map(f64::from);
        let mut real_x = 0.0;
        let mut real_y = 0.0;
        let mut real_z = 0.0;
        let mut imaginary_x = 0.0;
        let mut imaginary_y = 0.0;
        let mut imaginary_z = 0.0;
        for (atom, scatterer) in cell.atoms.iter().zip(&scattering.scatterers) {
            let moment = site_moments
                .get(atom.site_idx)
                .and_then(|moment| *moment)
                .unwrap_or([0.0; 3]);
            if norm(moment) <= 0.0 {
                continue;
            }
            let form = magnetic_form_factor(&atom.element, 0.5 / reflection.d)?;
            let attenuation = scatterer.attenuation(&metric, hkl);
            let phase = TAU * (h * atom.frac[0] + k * atom.frac[1] + l * atom.frac[2]);
            let value = MAGNETIC_PREFACTOR * form * atom.occupancy * attenuation;
            let [mx, my, mz] = moment;
            let (sine, cosine) = phase.sin_cos();
            real_x += value * mx * cosine;
            real_y += value * my * cosine;
            real_z += value * mz * cosine;
            imaginary_x += value * mx * sine;
            imaginary_y += value * my * sine;
            imaginary_z += value * mz * sine;
        }
        let amplitude_squared = if options.domain_average {
            (real_x * real_x
                + real_y * real_y
                + real_z * real_z
                + imaginary_x * imaginary_x
                + imaginary_y * imaginary_y
                + imaginary_z * imaginary_z)
                * 2.0
                / 3.0
        } else {
            let [nx, ny, nz] = normal;
            let along_real = real_x * nx + real_y * ny + real_z * nz;
            let along_imaginary = imaginary_x * nx + imaginary_y * ny + imaginary_z * nz;
            let perp_real_x = real_x - along_real * nx;
            let perp_real_y = real_y - along_real * ny;
            let perp_real_z = real_z - along_real * nz;
            let perp_imaginary_x = imaginary_x - along_imaginary * nx;
            let perp_imaginary_y = imaginary_y - along_imaginary * ny;
            let perp_imaginary_z = imaginary_z - along_imaginary * nz;
            perp_real_x * perp_real_x
                + perp_real_y * perp_real_y
                + perp_real_z * perp_real_z
                + perp_imaginary_x * perp_imaginary_x
                + perp_imaginary_y * perp_imaginary_y
                + perp_imaginary_z * perp_imaginary_z
        };
        if !amplitude_squared.is_finite() || amplitude_squared <= threshold {
            continue;
        }
        let two_theta = crate::neutron::bragg_angle(wavelength, reflection.d)?;
        let intensity = amplitude_squared
            * reflection.multiplicity as f64
            * profile.lorentz_polarization(two_theta);
        if !intensity.is_finite() || intensity <= 0.0 {
            continue;
        }
        reflections.push(MagneticReflection {
            hkl,
            d: reflection.d,
            two_theta,
            multiplicity: reflection.multiplicity,
            magnetic_amplitude_squared: amplitude_squared,
            intensity,
        });
    }
    reflections.sort_by(|left, right| {
        (left.two_theta * 1e9)
            .round()
            .total_cmp(&(right.two_theta * 1e9).round())
            .then_with(|| left.hkl.cmp(&right.hkl))
    });
    let peaks: Vec<(f64, f64)> = reflections
        .iter()
        .map(|reflection| (reflection.two_theta, reflection.intensity))
        .collect();
    let points = sample_cw_peaks(&peaks, profile, options.two_theta_range, options.step)?;
    Ok(MagneticPattern {
        wavelength,
        two_theta_range: options.two_theta_range,
        step: options.step,
        points,
        reflections,
        magnetic_sites,
        missing_moments,
    })
}

fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

fn norm(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}
