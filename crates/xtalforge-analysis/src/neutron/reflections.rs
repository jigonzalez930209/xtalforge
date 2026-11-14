//! `d`-window reflection enumeration shared by the neutron engines
//! (Task 19.3.1).

use xtalforge_core::{Error, Structure};

use crate::scattering;

/// Number of reflection samples above which enumeration is refused.
const MAX_REFLECTION_SAMPLES: f64 = 8_000_000.0;
/// Relative tolerance of the `d` window boundaries.
const D_TOLERANCE: f64 = 1e-12;

/// One reflection of a `d`-window enumeration.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ReflectionData {
    /// Miller indices of the canonical representative.
    pub hkl: [i32; 3],
    /// Interplanar spacing in Angstrom.
    pub d: f64,
    /// Orbit size, or one when the structure has no explicit operations.
    pub multiplicity: usize,
}

/// `d` window that covers a `2 theta` window at a wavelength.
pub(crate) fn d_range_from_two_theta(
    wavelength: f64,
    two_theta_range: [f64; 2],
) -> Result<[f64; 2], Error> {
    if !wavelength.is_finite() || wavelength <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "wavelength must be positive and finite, got {wavelength}"
        )));
    }
    let [minimum, maximum] = two_theta_range;
    if minimum <= 0.0 || maximum > 180.0 || maximum <= minimum {
        return Err(Error::InvalidInput(format!(
            "2 theta range must satisfy 0 < min < max <= 180, got {minimum}..{maximum}"
        )));
    }
    let sine_min = (0.5 * minimum.to_radians()).sin();
    let sine_max = (0.5 * maximum.min(179.9).to_radians()).sin();
    if sine_min <= 0.0 || sine_max <= 0.0 {
        return Err(Error::InvalidInput(
            "the 2 theta range has no Bragg solution".to_owned(),
        ));
    }
    let d_max = wavelength / (2.0 * sine_min);
    let d_min = wavelength / (2.0 * sine_max);
    if !d_min.is_finite() || !d_max.is_finite() || d_min <= 0.0 || d_max <= d_min {
        return Err(Error::InvalidInput(
            "the 2 theta range gives no finite d window".to_owned(),
        ));
    }
    Ok([d_min, d_max])
}

/// Enumerates the reflections whose `d` spacing lies in an inclusive window.
pub(crate) fn reflections_in_d_range(
    structure: &Structure,
    d_min: f64,
    d_max: f64,
    merge_symmetry: bool,
) -> Result<Vec<ReflectionData>, Error> {
    if !d_min.is_finite() || !d_max.is_finite() || d_min <= 0.0 || d_max < d_min {
        return Err(Error::InvalidInput(format!(
            "d window must be positive, finite and increasing, got {d_min}..{d_max}"
        )));
    }
    let lattice = structure.lattice();
    let metric = scattering::reciprocal_metric(lattice)?;
    let [a, b, c] = lattice.lengths();
    let h_max = (a / d_min).ceil() + 1.0;
    let k_max = (b / d_min).ceil() + 1.0;
    let l_max = (c / d_min).ceil() + 1.0;
    let samples = (2.0 * h_max + 1.0) * (2.0 * k_max + 1.0) * (2.0 * l_max + 1.0);
    if !samples.is_finite() || samples > MAX_REFLECTION_SAMPLES {
        return Err(Error::InvalidInput(format!(
            "the requested window needs {samples} hkl candidates; use a narrower range"
        )));
    }
    let h_max = h_max as i64;
    let k_max = k_max as i64;
    let l_max = l_max as i64;
    let inverse_rotations = if merge_symmetry {
        inverse_rotations(structure)?
    } else {
        Vec::new()
    };
    let merges = !inverse_rotations.is_empty();
    let mut reflections = Vec::new();
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
                let inverse_d = scattering::reciprocal_squared(&metric, hkl).sqrt();
                if inverse_d <= 0.0 {
                    continue;
                }
                let d = 1.0 / inverse_d;
                if d < d_min * (1.0 - D_TOLERANCE) || d > d_max * (1.0 + D_TOLERANCE) {
                    continue;
                }
                let (representative, multiplicity) = if merges {
                    let orbit = reflection_orbit(hkl, &inverse_rotations);
                    let representative = canonical_reflection(&orbit);
                    if representative != hkl {
                        continue;
                    }
                    (representative, orbit.len())
                } else {
                    (hkl, 1)
                };
                reflections.push(ReflectionData {
                    hkl: representative,
                    d,
                    multiplicity,
                });
            }
        }
    }
    reflections.sort_by(|left, right| {
        (left.d * 1e9)
            .round()
            .total_cmp(&(right.d * 1e9).round())
            .then_with(|| left.hkl.cmp(&right.hkl))
    });
    Ok(reflections)
}

/// Bragg angle of a `d` spacing in degrees.
pub(crate) fn bragg_angle(wavelength: f64, d: f64) -> Result<f64, Error> {
    let sine = 0.5 * wavelength / d;
    if !sine.is_finite() || sine <= 0.0 || sine >= 1.0 {
        return Err(Error::InvalidInput(format!(
            "d = {d} A has no Bragg solution at {wavelength} A"
        )));
    }
    Ok(2.0 * sine.asin().to_degrees())
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
