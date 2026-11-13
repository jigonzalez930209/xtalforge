//! Gaussian peak profile sampling for the historical PXRD entry point.

use std::f64::consts::{LN_2, PI};

use super::{PatternPoint, PxrdOptions, Reflection};

/// How many FWHM from a peak centre the Gaussian is evaluated.
const PROFILE_CUTOFF_FWHM: f64 = 12.0;

/// Samples every reflection profile on the `2 theta` grid described by the
/// options. The Gaussian has unit area and the requested FWHM.
pub(super) fn sample_pattern(
    reflections: &[Reflection],
    options: &PxrdOptions,
) -> Vec<PatternPoint> {
    let [minimum, maximum] = options.two_theta_range;
    let count = ((maximum - minimum) / options.step).floor() as usize + 1;
    let normalisation = 2.0 / options.fwhm * (LN_2 / PI).sqrt();
    let cutoff = PROFILE_CUTOFF_FWHM * options.fwhm;
    let mut points = Vec::with_capacity(count);
    for index in 0..count {
        let two_theta = minimum + index as f64 * options.step;
        let mut intensity = 0.0;
        for reflection in reflections {
            let offset = two_theta - reflection.two_theta;
            if offset.abs() > cutoff {
                continue;
            }
            let scaled = offset / options.fwhm;
            intensity +=
                reflection.intensity * normalisation * (-4.0 * LN_2 * scaled * scaled).exp();
        }
        points.push(PatternPoint {
            two_theta,
            intensity,
        });
    }
    points
}
