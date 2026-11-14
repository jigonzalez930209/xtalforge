//! Kinematic selected-area electron diffraction along a zone axis
//! (Task 19.3.3).
//!
//! # Geometry
//!
//! The electron beam is parallel to the direct lattice vector
//! `z = u a + v b + w c` of the requested zone axis. Only the zero-order Laue
//! zone (`h u + k v + l w = 0`, or a larger `max_laue_index`) is kept. The
//! reciprocal vector `g = h a* + k b* + l c*` is projected on a detector frame
//! perpendicular to `z`,
//!
//! ```text
//! x = L lambda (g . ex),   y = L lambda (g . ey),
//! ```
//!
//! with `L` the camera length in millimetre, `lambda` the relativistic
//! electron wavelength and `(ex, ey)` an orthonormal frame with `ez = z`.
//! The excitation error is reported as `g . ez + lambda g^2 / 2` in inverse
//! Angstrom.
//!
//! # Intensities
//!
//! Spot intensities are the kinematic `|F_e(hkl)|^2` of the electron structure
//! factor, using the Mott-Bethe form factor of
//! [`electron_scattering_factor`](crate::electron_scattering_factor) (the
//! embedded International Tables Vol. C X-ray coefficients) and the same
//! Debye-Waller model as the X-ray engine. There is no Lorentz
//! factor, absorption, thickness or dynamical correction; kinematic
//! intensities are only a first approximation for thin crystals.
//!
//! # Image buffer
//!
//! [`SaedPattern::image`] renders the spots as Gaussian discs in a square
//! `Vec<f32>` buffer normalized to a maximum of one. The direct beam is drawn
//! at the centre with an intensity relative to the brightest reflection
//! ([`SaedOptions::central_intensity`]). [`SaedImage::to_u8`] gives the same
//! image on the 0 to 255 scale for a web canvas.
//!
//! The result is deterministic: fixed hkl enumeration, fixed orientation frame
//! and fixed summation order.

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};

use crate::instrument::{RadiationKind, electron_wavelength};
use crate::scattering;

mod render;

pub use render::SaedImage;
use render::{MAX_IMAGE_SIZE, MIN_IMAGE_SIZE, add_spot};

/// Number of reflection samples above which enumeration is refused.
const MAX_REFLECTION_SAMPLES: f64 = 8_000_000.0;
/// Relative amplitude below which a reflection counts as absent.
const ABSENCE_TOLERANCE: f64 = 1e-8;
/// Largest accepted Laue-zone index.
const MAX_LAUE_INDEX: u32 = 100;

/// Options for [`zone_axis_pattern`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SaedOptions {
    /// Accelerating voltage in kilovolt.
    pub acceleration_kv: f64,
    /// Zone axis `[u v w]` in direct lattice indices.
    pub zone_axis: [i32; 3],
    /// Camera length in millimetre.
    pub camera_length_mm: f64,
    /// Detector width in millimetre.
    pub field_of_view_mm: f64,
    /// Square image edge in pixels, from 16 to 2048.
    pub image_size: usize,
    /// Standard deviation of the Gaussian spot in pixels.
    pub spot_sigma_pixels: f64,
    /// Reflections with a kinematic intensity below this value are dropped.
    pub min_intensity: f64,
    /// Direct-beam intensity relative to the brightest reflection; `0` hides it.
    pub central_intensity: f64,
    /// Largest included `|h u + k v + l w|`; `0` is the zero-order Laue zone.
    pub max_laue_index: u32,
}

impl Default for SaedOptions {
    fn default() -> Self {
        Self {
            acceleration_kv: 200.0,
            zone_axis: [0, 0, 1],
            camera_length_mm: 500.0,
            field_of_view_mm: 20.0,
            image_size: 512,
            spot_sigma_pixels: 2.0,
            min_intensity: 0.0,
            central_intensity: 4.0,
            max_laue_index: 0,
        }
    }
}

impl SaedOptions {
    /// Validates the voltage, zone axis, geometry, image and thresholds.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive voltage, a zero zone
    /// axis, a non-positive camera length or field of view, an image edge
    /// outside `[16, 2048]`, a non-positive spot width, a negative threshold or
    /// central intensity, or a Laue index above 100.
    pub fn validate(&self) -> Result<(), Error> {
        electron_wavelength(self.acceleration_kv)?;
        if self.zone_axis == [0, 0, 0] {
            return Err(Error::InvalidInput(
                "the zone axis must not be the zero vector".to_owned(),
            ));
        }
        for (name, value) in [
            ("camera length", self.camera_length_mm),
            ("field of view", self.field_of_view_mm),
            ("spot sigma", self.spot_sigma_pixels),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(Error::InvalidInput(format!(
                    "{name} must be positive and finite, got {value}"
                )));
            }
        }
        if !(MIN_IMAGE_SIZE..=MAX_IMAGE_SIZE).contains(&self.image_size) {
            return Err(Error::InvalidInput(format!(
                "image size must be in [{MIN_IMAGE_SIZE}, {MAX_IMAGE_SIZE}] pixels, got {}",
                self.image_size
            )));
        }
        if self.spot_sigma_pixels > self.image_size as f64 {
            return Err(Error::InvalidInput(
                "the spot sigma is larger than the image".to_owned(),
            ));
        }
        for (name, value) in [
            ("minimum intensity", self.min_intensity),
            ("central intensity", self.central_intensity),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(Error::InvalidInput(format!(
                    "{name} must be non-negative and finite, got {value}"
                )));
            }
        }
        if self.max_laue_index > MAX_LAUE_INDEX {
            return Err(Error::InvalidInput(format!(
                "the Laue index must not exceed {MAX_LAUE_INDEX}, got {}",
                self.max_laue_index
            )));
        }
        Ok(())
    }

    fn pixels_per_mm(&self) -> f64 {
        self.image_size as f64 / self.field_of_view_mm
    }
}

/// One electron diffraction spot.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SaedSpot {
    /// Miller indices.
    pub hkl: [i32; 3],
    /// Interplanar spacing in Angstrom.
    pub d: f64,
    /// `|g|` in inverse Angstrom.
    pub g: f64,
    /// Detector coordinate along the reference axis, in millimetre.
    pub x_mm: f64,
    /// Detector coordinate along the second axis, in millimetre.
    pub y_mm: f64,
    /// Continuous image column, with the direct beam at `image_size / 2`.
    pub x_px: f64,
    /// Continuous image row, with the direct beam at `image_size / 2`.
    pub y_px: f64,
    /// Excitation error `g . ez + lambda g^2 / 2` in inverse Angstrom.
    pub excitation_error: f64,
    /// Kinematic intensity `|F_e|^2`.
    pub intensity: f64,
}

impl SaedSpot {
    /// Distance of the spot from the direct beam in millimetre.
    #[must_use]
    pub fn radius_mm(&self) -> f64 {
        (self.x_mm * self.x_mm + self.y_mm * self.y_mm).sqrt()
    }

    /// `(h k l)` with a minus sign on negative indices.
    #[must_use]
    pub fn label(&self) -> String {
        let [h, k, l] = self.hkl;
        format!("({h} {k} {l})")
    }
}

/// A kinematic zone-axis electron diffraction pattern.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SaedPattern {
    /// Accelerating voltage in kilovolt.
    pub acceleration_kv: f64,
    /// Relativistic electron wavelength in Angstrom.
    pub wavelength: f64,
    /// Zone axis `[u v w]`.
    pub zone_axis: [i32; 3],
    /// Camera length in millimetre.
    pub camera_length_mm: f64,
    /// Detector width in millimetre.
    pub field_of_view_mm: f64,
    /// Square image edge in pixels.
    pub image_size: usize,
    /// Pixels per millimetre of the image.
    pub pixels_per_mm: f64,
    /// Unit Cartesian zone normal.
    pub zone_normal: [f64; 3],
    /// Unit Cartesian detector x axis.
    pub detector_x: [f64; 3],
    /// Unit Cartesian detector y axis.
    pub detector_y: [f64; 3],
    /// Direct-beam intensity relative to the brightest reflection.
    pub central_intensity: f64,
    /// Standard deviation of the Gaussian spot in pixels.
    pub spot_sigma_pixels: f64,
    /// Spots, ordered by increasing radius and then `hkl`.
    pub spots: Vec<SaedSpot>,
}

impl SaedPattern {
    /// Largest reflection intensity, or `0` without reflections.
    #[must_use]
    pub fn max_intensity(&self) -> f64 {
        self.spots
            .iter()
            .map(|spot| spot.intensity)
            .fold(0.0, f64::max)
    }

    /// Renders the pattern into a normalized square image.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the image size or spot width is
    /// outside the supported range.
    pub fn image(&self) -> Result<SaedImage, Error> {
        if !(MIN_IMAGE_SIZE..=MAX_IMAGE_SIZE).contains(&self.image_size) {
            return Err(Error::InvalidInput(format!(
                "image size must be in [{MIN_IMAGE_SIZE}, {MAX_IMAGE_SIZE}] pixels, got {}",
                self.image_size
            )));
        }
        if !self.spot_sigma_pixels.is_finite() || self.spot_sigma_pixels <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "spot sigma must be positive and finite, got {}",
                self.spot_sigma_pixels
            )));
        }
        let size = self.image_size;
        let center = size as f64 / 2.0;
        let sigma = self.spot_sigma_pixels;
        let brightest = self.max_intensity();
        let central = if brightest > 0.0 {
            self.central_intensity * brightest
        } else {
            self.central_intensity.max(1.0)
        };
        let mut pixels = vec![0.0_f32; size * size];
        if central > 0.0 {
            add_spot(&mut pixels, size, center, center, sigma, central);
        }
        for spot in &self.spots {
            add_spot(
                &mut pixels,
                size,
                spot.x_px,
                spot.y_px,
                sigma,
                spot.intensity,
            );
        }
        let maximum = pixels.iter().copied().fold(0.0_f32, f32::max);
        if maximum > 0.0 {
            for pixel in &mut pixels {
                *pixel /= maximum;
            }
        }
        Ok(SaedImage {
            width: size,
            height: size,
            pixels,
        })
    }
}

/// Builds a kinematic zone-axis electron diffraction pattern.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, a structure without
/// scattering atoms, an element outside the reference tables, a singular
/// lattice frame, or too many `(h k l)` candidates.
pub fn zone_axis_pattern(
    structure: &Structure,
    options: &SaedOptions,
) -> Result<SaedPattern, Error> {
    options.validate()?;
    let wavelength = electron_wavelength(options.acceleration_kv)?;
    let cell = structure.p1()?;
    let scattering = scattering::build(structure, &cell, RadiationKind::Electron)?;
    if scattering.scatterers.is_empty() {
        return Err(Error::InvalidInput(
            "there are no atoms to calculate an electron pattern from".to_owned(),
        ));
    }
    let lattice = structure.lattice();
    let direct = lattice.matrix();
    let [u, v, w] = options.zone_axis.map(f64::from);
    let zone = [
        u * direct[0][0] + v * direct[1][0] + w * direct[2][0],
        u * direct[0][1] + v * direct[1][1] + w * direct[2][1],
        u * direct[0][2] + v * direct[1][2] + w * direct[2][2],
    ];
    let zone_norm = norm(zone);
    if zone_norm <= 1e-12 || !zone_norm.is_finite() {
        return Err(Error::InvalidInput(
            "the zone axis has no finite direct-lattice direction".to_owned(),
        ));
    }
    let zone_normal = scale(zone, 1.0 / zone_norm);
    let reciprocal = lattice.reciprocal(false)?.matrix();
    let detector_x = reference_axis(reciprocal, zone_normal)?;
    let detector_y = cross(zone_normal, detector_x);
    let detector_y_norm = norm(detector_y);
    if detector_y_norm <= 1e-12 {
        return Err(Error::InvalidInput(
            "the detector frame is singular".to_owned(),
        ));
    }
    let detector_y = scale(detector_y, 1.0 / detector_y_norm);

    let pixels_per_mm = options.pixels_per_mm();
    let radius_mm = 0.5 * options.field_of_view_mm;
    let absence = (ABSENCE_TOLERANCE * scattering.total).powi(2);
    let g_max = radius_mm / (options.camera_length_mm * wavelength);
    let [a, b, c] = lattice.lengths();
    let h_max = (a * g_max).ceil() + 1.0;
    let k_max = (b * g_max).ceil() + 1.0;
    let l_max = (c * g_max).ceil() + 1.0;
    let samples = (2.0 * h_max + 1.0) * (2.0 * k_max + 1.0) * (2.0 * l_max + 1.0);
    if !samples.is_finite() || samples > MAX_REFLECTION_SAMPLES {
        return Err(Error::InvalidInput(format!(
            "the requested detector needs {samples} hkl candidates; use a smaller field of view"
        )));
    }
    let h_max = h_max as i64;
    let k_max = k_max as i64;
    let l_max = l_max as i64;
    let center = options.image_size as f64 / 2.0;
    let laue_limit = u64::from(options.max_laue_index);
    let mut spots = Vec::new();
    for h in -h_max..=h_max {
        for k in -k_max..=k_max {
            for l in -l_max..=l_max {
                if h == 0 && k == 0 && l == 0 {
                    continue;
                }
                let zone_index = h * i64::from(options.zone_axis[0])
                    + k * i64::from(options.zone_axis[1])
                    + l * i64::from(options.zone_axis[2]);
                if zone_index.unsigned_abs() > laue_limit {
                    continue;
                }
                let (Ok(h), Ok(k), Ok(l)) = (i32::try_from(h), i32::try_from(k), i32::try_from(l))
                else {
                    continue;
                };
                let hkl = [h, k, l];
                let [hf, kf, lf] = hkl.map(f64::from);
                let g = [
                    hf * reciprocal[0][0] + kf * reciprocal[1][0] + lf * reciprocal[2][0],
                    hf * reciprocal[0][1] + kf * reciprocal[1][1] + lf * reciprocal[2][1],
                    hf * reciprocal[0][2] + kf * reciprocal[1][2] + lf * reciprocal[2][2],
                ];
                let length = norm(g);
                if length <= 0.0 || !length.is_finite() {
                    continue;
                }
                let along = dot(g, zone_normal);
                let perpendicular = [
                    g[0] - along * zone_normal[0],
                    g[1] - along * zone_normal[1],
                    g[2] - along * zone_normal[2],
                ];
                let x_mm = options.camera_length_mm * wavelength * dot(perpendicular, detector_x);
                let y_mm = options.camera_length_mm * wavelength * dot(perpendicular, detector_y);
                if x_mm * x_mm + y_mm * y_mm > radius_mm * radius_mm {
                    continue;
                }
                let intensity = scattering::squared_amplitude(&scattering, hkl);
                if !intensity.is_finite() || intensity <= absence.max(options.min_intensity) {
                    continue;
                }
                spots.push(SaedSpot {
                    hkl,
                    d: 1.0 / length,
                    g: length,
                    x_mm,
                    y_mm,
                    x_px: center + x_mm * pixels_per_mm,
                    y_px: center - y_mm * pixels_per_mm,
                    excitation_error: along + 0.5 * wavelength * length * length,
                    intensity,
                });
            }
        }
    }
    spots.sort_by(|left, right| {
        let key = |spot: &SaedSpot| (spot.radius_mm() * 1e9).round();
        key(left)
            .total_cmp(&key(right))
            .then_with(|| left.hkl.cmp(&right.hkl))
    });
    Ok(SaedPattern {
        acceleration_kv: options.acceleration_kv,
        wavelength,
        zone_axis: options.zone_axis,
        camera_length_mm: options.camera_length_mm,
        field_of_view_mm: options.field_of_view_mm,
        image_size: options.image_size,
        pixels_per_mm,
        zone_normal,
        detector_x,
        detector_y,
        central_intensity: options.central_intensity,
        spot_sigma_pixels: options.spot_sigma_pixels,
        spots,
    })
}

fn reference_axis(reciprocal: [[f64; 3]; 3], zone_normal: [f64; 3]) -> Result<[f64; 3], Error> {
    let mut best = [0.0; 3];
    let mut best_norm = 0.0;
    for row in reciprocal {
        let along = dot(row, zone_normal);
        let projection = [
            row[0] - along * zone_normal[0],
            row[1] - along * zone_normal[1],
            row[2] - along * zone_normal[2],
        ];
        let length = norm(projection);
        if length > best_norm {
            best = projection;
            best_norm = length;
        }
    }
    if best_norm <= 1e-12 || !best_norm.is_finite() {
        return Err(Error::InvalidInput(
            "the zone axis has no perpendicular reciprocal reference".to_owned(),
        ));
    }
    Ok(scale(best, 1.0 / best_norm))
}

fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

fn norm(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}

fn scale(vector: [f64; 3], factor: f64) -> [f64; 3] {
    [vector[0] * factor, vector[1] * factor, vector[2] * factor]
}

fn cross(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}
