//! Profile sampling of the neutron powder engine (Task 19.3.1).

use xtalforge_core::Error;

use crate::instrument::{InstrumentProfile, LineShape};

use super::{NeutronPoint, NeutronReflection};

/// Number of pattern points above which sampling is refused.
pub(crate) const MAX_PATTERN_POINTS: usize = 2_000_000;
/// How many FWHM from a line centre the profile is evaluated.
const PROFILE_CUTOFF_FWHM: f64 = 12.0;

/// Number of points of a sampled range; also rejects non-finite counts.
pub(crate) fn point_count(range: [f64; 2], step: f64) -> Result<usize, Error> {
    let count = ((range[1] - range[0]) / step).floor() + 1.0;
    if !count.is_finite() || count <= 0.0 {
        return Err(Error::InvalidInput(
            "the sampling range and step do not produce a finite pattern".to_owned(),
        ));
    }
    Ok(count as usize)
}

/// Samples a constant-wavelength peak list on `2 theta`.
pub(crate) fn sample_cw_peaks(
    peaks: &[(f64, f64)],
    profile: &InstrumentProfile,
    range: [f64; 2],
    step: f64,
) -> Result<Vec<NeutronPoint>, Error> {
    let count = point_count(range, step)?;
    let mut points = Vec::with_capacity(count);
    for index in 0..count {
        let x = range[0] + index as f64 * step;
        let mut intensity = 0.0;
        for &(center, height) in peaks {
            let width = profile.fwhm_at(center);
            if (x - center).abs() > PROFILE_CUTOFF_FWHM * width {
                continue;
            }
            intensity += height * profile.area_normalised(center, x);
        }
        points.push(NeutronPoint { x, intensity });
    }
    Ok(points)
}

pub(crate) fn sample_tof_peaks(
    reflections: &[NeutronReflection],
    profile: &InstrumentProfile,
    relative_d_fwhm: f64,
    range: [f64; 2],
    step: f64,
) -> Result<Vec<NeutronPoint>, Error> {
    let count = point_count(range, step)?;
    let mut line = *profile;
    line.caglioti = None;
    line.tch = None;
    line.asymmetry = None;
    if line.shape == LineShape::Tch {
        line.shape = LineShape::PseudoVoigt;
        line.eta = 0.5;
    }
    let mut points = Vec::with_capacity(count);
    for index in 0..count {
        let x = range[0] + index as f64 * step;
        let mut intensity = 0.0;
        for reflection in reflections {
            let Some(center) = reflection.tof_us else {
                continue;
            };
            let width = (relative_d_fwhm * center).max(1e-9);
            if (x - center).abs() > PROFILE_CUTOFF_FWHM * width {
                continue;
            }
            line.fwhm = width;
            intensity += reflection.intensity * line.area_normalised(center, x);
        }
        points.push(NeutronPoint { x, intensity });
    }
    Ok(points)
}

pub(crate) fn sort_reflections(reflections: &mut [NeutronReflection]) {
    let position =
        |reflection: &NeutronReflection| reflection.tof_us.unwrap_or(reflection.two_theta);
    reflections.sort_by(|left, right| {
        (position(left) * 1e9)
            .round()
            .total_cmp(&(position(right) * 1e9).round())
            .then_with(|| left.hkl.cmp(&right.hkl))
    });
}
