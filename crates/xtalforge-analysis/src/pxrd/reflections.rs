//! Reflection enumeration shared by the PXRD, single-crystal and electron
//! engines (Tasks 5.1 and 14.1.4).

use xtalforge_core::{Error, Structure};

use crate::instrument::{RadiationKind, profile::lorentz_polarization};
use crate::scattering;

use super::Reflection;

/// Number of `(h k l)` candidates above which enumeration is refused.
const MAX_REFLECTION_SAMPLES: f64 = 8_000_000.0;
/// Relative amplitude below which a structure factor counts as absent.
const ABSENCE_TOLERANCE: f64 = 1e-8;

/// Enumerates the allowed reflections of a structure.
///
/// `merge_symmetry` groups each reflection orbit (Friedel partners included)
/// and keeps the canonical representative with the orbit size as its
/// multiplicity; without it every `(h k l)` is listed with multiplicity one.
/// The intensity is `|F|^2 * multiplicity * LP`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an invalid wavelength or range, a
/// structure without scattering atoms, an element outside the reference
/// tables, or too many `(h k l)` candidates.
pub(crate) fn enumerate(
    structure: &Structure,
    wavelength: f64,
    kind: RadiationKind,
    two_theta_range: [f64; 2],
    polarization: f64,
    merge_symmetry: bool,
) -> Result<Vec<Reflection>, Error> {
    if !wavelength.is_finite() || wavelength <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "radiation wavelength must be positive and finite, got {wavelength}"
        )));
    }
    let [minimum, maximum] = two_theta_range;
    if !minimum.is_finite() || !maximum.is_finite() || maximum <= minimum {
        return Err(Error::InvalidInput(format!(
            "2 theta range must be finite and increasing, got {minimum}..{maximum}"
        )));
    }
    let lattice = structure.lattice();
    let cell = structure.p1()?;
    let scattering = scattering::build(structure, &cell, kind)?;
    if scattering.scatterers.is_empty() || scattering.total <= 0.0 {
        return Err(Error::InvalidInput(
            "there are no scattering atoms to calculate a pattern from".to_owned(),
        ));
    }
    let reciprocal = lattice.reciprocal(false)?.matrix();
    let [a, b, c] = lattice.lengths();
    let theta_max = (maximum.min(179.9) * 0.5).to_radians();
    let sine = theta_max.sin();
    if sine <= 0.0 {
        return Err(Error::InvalidInput(
            "the 2 theta range has no Bragg solution".to_owned(),
        ));
    }
    let d_min = wavelength / (2.0 * sine);
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

    let inverse_rotations = if merge_symmetry {
        inverse_rotations(structure)?
    } else {
        Vec::new()
    };
    let merges_symmetry = !inverse_rotations.is_empty();
    let cutoff = ABSENCE_TOLERANCE * scattering.total;
    let mut reflections: Vec<Reflection> = Vec::new();

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
                let inverse_d = reciprocal_length(hkl, &reciprocal);
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

                let squared = scattering::squared_amplitude(&scattering, hkl);
                if !squared.is_finite() || squared <= cutoff * cutoff {
                    continue;
                }
                let intensity = squared
                    * multiplicity as f64
                    * lorentz_polarization(kind, polarization, two_theta);
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
    Ok(reflections)
}

/// Length of the reciprocal vector of one reflection, in inverse Angstrom.
pub(crate) fn reciprocal_length(hkl: [i32; 3], reciprocal: &[[f64; 3]; 3]) -> f64 {
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
