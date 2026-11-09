//! Calculated-curve evaluation with zero shift, width scale, texture and the
//! copper K alpha doublet (Task 7.8.1).

use std::f64::consts::{LN_2, PI};

use xtalforge_core::{Error, Structure};
use xtalforge_io::xy::PatternPoint;

use crate::instrument::doublet;
use crate::pxrd::{PxrdOptions, Radiation, calculate_pxrd};

/// Copper K alpha 1 wavelength in Angstrom.
const CU_KA1: f64 = 1.540_598;
/// Copper K alpha 2 wavelength in Angstrom.
const CU_KA2: f64 = 1.544_426;
/// Copper K alpha 2 / K alpha 1 intensity ratio.
const CU_KA2_RATIO: f64 = 0.5;

/// Parameters that move the calculated curve without changing the structure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CurveOptions {
    /// Zero shift in degrees two-theta.
    pub(crate) zero: f64,
    /// Multiplier on the PXRD FWHM.
    pub(crate) width_scale: f64,
    /// March-Dollase coefficient (`1` is random orientation).
    pub(crate) texture: f64,
    /// Preferred-orientation axis in fractional direct coordinates.
    pub(crate) texture_axis: [f64; 3],
    /// Place the copper K alpha 1/K alpha 2 doublet instead of one line.
    pub(crate) doublet: bool,
}

impl Default for CurveOptions {
    fn default() -> Self {
        Self {
            zero: 0.0,
            width_scale: 1.0,
            texture: 1.0,
            texture_axis: [0.0, 0.0, 1.0],
            doublet: false,
        }
    }
}

/// Samples the phase at the observed points.
pub(crate) fn evaluate_curve(
    structure: &Structure,
    observed: &[PatternPoint],
    pxrd: &PxrdOptions,
    options: &CurveOptions,
) -> Result<Vec<f64>, Error> {
    // The doublet is placed from the K alpha 1 position, so the pattern is
    // calculated at that wavelength and split into the two lines below.
    let effective = if options.doublet && pxrd.radiation == Radiation::CuKa {
        PxrdOptions {
            radiation: Radiation::CuKa1,
            ..*pxrd
        }
    } else {
        *pxrd
    };
    let pattern = calculate_pxrd(structure, &effective)?;
    if pattern.reflections.is_empty() {
        return Err(Error::InvalidInput(
            "the structure has no reflections in the observed range".to_owned(),
        ));
    }
    let width_scale = options.width_scale.clamp(0.05, 20.0);
    let fwhm = (pxrd.fwhm * width_scale).max(1e-3);
    let normalisation = 2.0 / fwhm * (LN_2 / PI).sqrt();
    let reciprocal = structure.lattice().reciprocal(false)?.matrix();
    let matrix = structure.lattice().matrix();
    let axis = axis_cartesian(options.texture_axis, matrix);
    let texture = options.texture.clamp(0.1, 10.0);
    let mut curve = vec![0.0; observed.len()];
    for reflection in &pattern.reflections {
        let factor = march_dollase(reflection.hkl, &reciprocal, axis, texture);
        let lines: Vec<(f64, f64)> = if options.doublet && pxrd.radiation == Radiation::CuKa {
            doublet(reflection.two_theta, CU_KA1, CU_KA2, CU_KA2_RATIO)
        } else {
            vec![(reflection.two_theta, 1.0)]
        };
        for (position, weight) in lines {
            let intensity = reflection.intensity * factor * weight;
            let center = position + options.zero;
            for (slot, point) in curve.iter_mut().zip(observed.iter()) {
                let delta = (point.two_theta - center) / fwhm;
                if delta.abs() > 8.0 {
                    continue;
                }
                *slot += intensity * normalisation * (-4.0 * LN_2 * delta * delta).exp();
            }
        }
    }
    Ok(curve)
}

/// March-Dollase preferred-orientation factor for one reflection.
fn march_dollase(
    hkl: [i32; 3],
    reciprocal: &[[f64; 3]; 3],
    axis: Option<[f64; 3]>,
    ratio: f64,
) -> f64 {
    let Some(axis) = axis else {
        return 1.0;
    };
    let [h, k, l] = hkl.map(f64::from);
    let normal = [
        h * reciprocal[0][0] + k * reciprocal[1][0] + l * reciprocal[2][0],
        h * reciprocal[0][1] + k * reciprocal[1][1] + l * reciprocal[2][1],
        h * reciprocal[0][2] + k * reciprocal[1][2] + l * reciprocal[2][2],
    ];
    let norm = normal.iter().map(|value| value * value).sum::<f64>().sqrt();
    if norm <= 1e-12 {
        return 1.0;
    }
    let cosine = (normal[0] * axis[0] + normal[1] * axis[1] + normal[2] * axis[2]) / norm;
    let squared = cosine * cosine;
    let sine = (1.0 - squared).max(0.0);
    let base = ratio * ratio * squared + sine / ratio;
    if base <= 0.0 {
        return 1.0;
    }
    base.powf(-1.5)
}

fn axis_cartesian(axis: [f64; 3], matrix: [[f64; 3]; 3]) -> Option<[f64; 3]> {
    let norm = axis.iter().map(|value| value * value).sum::<f64>().sqrt();
    if !norm.is_finite() || norm <= 1e-12 {
        return None;
    }
    let scaled = axis.map(|value| value / norm);
    let cartesian = [
        scaled[0] * matrix[0][0] + scaled[1] * matrix[1][0] + scaled[2] * matrix[2][0],
        scaled[0] * matrix[0][1] + scaled[1] * matrix[1][1] + scaled[2] * matrix[2][1],
        scaled[0] * matrix[0][2] + scaled[1] * matrix[1][2] + scaled[2] * matrix[2][2],
    ];
    let length = cartesian
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    if length <= 1e-12 {
        None
    } else {
        Some(cartesian.map(|value| value / length))
    }
}
