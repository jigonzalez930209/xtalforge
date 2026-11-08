//! Golden-data comparison helpers for xtalforge tests (Task 0.9.4).
//!
//! The helpers return `Result<(), String>` with a readable diff; the
//! `assert_*` wrappers turn that into a test failure with the same message.
//! Comparisons are deterministic and independent of atom or bond order.

fn fail(message: String) {
    // Every comparator returns a non-empty message; this assertion is the
    // test failure and avoids clippy's panic lint in a test-support library.
    let empty = message.is_empty();
    assert!(empty, "{message}");
}

/// One atom in fractional coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    /// Element symbol, case-sensitive as stored by the reader.
    pub element: String,
    /// Fractional coordinates.
    pub frac: [f64; 3],
}

impl Atom {
    /// Creates an atom.
    #[must_use]
    pub fn new(element: impl Into<String>, frac: [f64; 3]) -> Self {
        Self {
            element: element.into(),
            frac,
        }
    }
}

/// One bond between two atoms, with a lattice image for the second atom.
#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    /// Index of the first atom.
    pub i: usize,
    /// Index of the second atom.
    pub j: usize,
    /// Lattice image of `j`, in lattice translations.
    pub image: [i32; 3],
    /// Bond order (1.0 single, 1.5 aromatic, ...).
    pub order: f64,
}

impl Bond {
    /// Creates a bond.
    #[must_use]
    pub fn new(i: usize, j: usize, image: [i32; 3], order: f64) -> Self {
        Self { i, j, image, order }
    }
}

/// One diffraction peak.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Peak {
    /// Position in degrees two-theta.
    pub position: f64,
    /// Normalized intensity.
    pub intensity: f64,
}

/// Compares two numbers with absolute and relative tolerances.
///
/// # Errors
///
/// Returns a message with both values when they differ by more than
/// `atol + rtol * |expected|`.
pub fn compare_close(actual: f64, expected: f64, atol: f64, rtol: f64) -> Result<(), String> {
    if !actual.is_finite() || !expected.is_finite() {
        return Err(format!(
            "non-finite value: actual={actual}, expected={expected}"
        ));
    }
    let difference = (actual - expected).abs();
    let allowed = atol + rtol * expected.abs();
    if difference > allowed {
        return Err(format!(
            "values differ: actual={actual}, expected={expected}, \
             difference={difference}, allowed={allowed}"
        ));
    }
    Ok(())
}

/// Asserts [`compare_close`] with a label.
pub fn assert_close(actual: f64, expected: f64, atol: f64, rtol: f64) {
    if let Err(message) = compare_close(actual, expected, atol, rtol) {
        fail(message);
    }
}

/// Compares atom sets up to a common fractional translation and permutation.
///
/// The translation is taken from the first expected atom against every actual
/// atom, so the comparison is exact up to `atol` and does not need a global
/// origin. The function is intended for test-sized sets (a few hundred atoms).
///
/// # Errors
///
/// Returns a message when the lengths, elements or positions differ.
pub fn compare_atom_sets(expected: &[Atom], actual: &[Atom], atol: f64) -> Result<(), String> {
    if expected.len() != actual.len() {
        return Err(format!(
            "atom count differs: expected={}, actual={}",
            expected.len(),
            actual.len()
        ));
    }
    let Some(first_expected) = expected.first() else {
        return Ok(());
    };
    let [first_x, first_y, first_z] = first_expected.frac;

    for candidate in actual {
        let [candidate_x, candidate_y, candidate_z] = candidate.frac;
        let [shift_x, shift_y, shift_z] = [
            candidate_x - first_x,
            candidate_y - first_y,
            candidate_z - first_z,
        ];
        let mut used = vec![false; actual.len()];
        let mut matched = true;
        for atom in expected {
            let [x, y, z] = atom.frac;
            let target = [x + shift_x, y + shift_y, z + shift_z];
            let found = actual.iter().enumerate().find(|(index, other)| {
                !used.get(*index).copied().unwrap_or(true)
                    && other.element == atom.element
                    && frac_close(other.frac, target, atol)
            });
            match found {
                Some((index, _)) => {
                    if let Some(slot) = used.get_mut(index) {
                        *slot = true;
                    }
                }
                None => {
                    matched = false;
                    break;
                }
            }
        }
        if matched {
            return Ok(());
        }
    }

    Err("atom sets differ after every candidate translation".to_owned())
}

/// Asserts [`compare_atom_sets`].
pub fn assert_atom_sets_equal(expected: &[Atom], actual: &[Atom], atol: f64) {
    if let Err(message) = compare_atom_sets(expected, actual, atol) {
        fail(message);
    }
}

/// Compares bond sets independent of order.
///
/// Atom indices must refer to the same atom order; use atom-set comparison
/// first when the reference order differs. Bond order is compared with
/// `order_atol`.
///
/// # Errors
///
/// Returns a message naming the first missing or extra bond.
pub fn compare_bond_sets(
    expected: &[Bond],
    actual: &[Bond],
    order_atol: f64,
) -> Result<(), String> {
    if expected.len() != actual.len() {
        return Err(format!(
            "bond count differs: expected={}, actual={}",
            expected.len(),
            actual.len()
        ));
    }
    let key = |bond: &Bond| {
        let (i, j) = if bond.i <= bond.j {
            (bond.i, bond.j)
        } else {
            (bond.j, bond.i)
        };
        (i, j, bond.image)
    };
    let mut used = vec![false; actual.len()];
    for bond in expected {
        let expected_key = key(bond);
        let found = actual.iter().enumerate().find(|(index, other)| {
            !used.get(*index).copied().unwrap_or(true)
                && key(other) == expected_key
                && (other.order - bond.order).abs() <= order_atol
        });
        match found {
            Some((index, _)) => {
                if let Some(slot) = used.get_mut(index) {
                    *slot = true;
                }
            }
            None => {
                return Err(format!(
                    "bond missing or different: {bond:?} (key {expected_key:?})"
                ));
            }
        }
    }
    Ok(())
}

/// Asserts [`compare_bond_sets`].
pub fn assert_bond_sets_equal(expected: &[Bond], actual: &[Bond], order_atol: f64) {
    if let Err(message) = compare_bond_sets(expected, actual, order_atol) {
        fail(message);
    }
}

/// Compares two normalized PXRD peak lists by position, independent of order.
///
/// Peaks are matched greedily within `position_atol` degrees; intensities are
/// compared within `intensity_atol` after normalization by the strongest peak.
///
/// # Errors
///
/// Returns a message when the peak counts or any matched pair differ.
pub fn compare_peaks(
    expected: &[Peak],
    actual: &[Peak],
    position_atol: f64,
    intensity_atol: f64,
) -> Result<(), String> {
    if expected.len() != actual.len() {
        return Err(format!(
            "peak count differs: expected={}, actual={}",
            expected.len(),
            actual.len()
        ));
    }
    let normalize = |peaks: &[Peak]| {
        let max = peaks.iter().map(|peak| peak.intensity).fold(0.0, f64::max);
        if max <= 0.0 {
            return peaks.to_vec();
        }
        peaks
            .iter()
            .map(|peak| Peak {
                position: peak.position,
                intensity: peak.intensity / max,
            })
            .collect::<Vec<_>>()
    };
    let expected = normalize(expected);
    let actual = normalize(actual);
    let mut used = vec![false; actual.len()];
    for peak in &expected {
        let found = actual.iter().enumerate().find(|(index, other)| {
            !used.get(*index).copied().unwrap_or(true)
                && (other.position - peak.position).abs() <= position_atol
                && (other.intensity - peak.intensity).abs() <= intensity_atol
        });
        match found {
            Some((index, _)) => {
                if let Some(slot) = used.get_mut(index) {
                    *slot = true;
                }
            }
            None => {
                return Err(format!(
                    "peak missing or different: position={}, intensity={}",
                    peak.position, peak.intensity
                ));
            }
        }
    }
    Ok(())
}

/// Asserts [`compare_peaks`].
pub fn assert_peaks_equal(
    expected: &[Peak],
    actual: &[Peak],
    position_atol: f64,
    intensity_atol: f64,
) {
    if let Err(message) = compare_peaks(expected, actual, position_atol, intensity_atol) {
        fail(message);
    }
}

fn frac_close(first: [f64; 3], second: [f64; 3], atol: f64) -> bool {
    first.iter().zip(second.iter()).all(|(first, second)| {
        let difference = (first - second).rem_euclid(1.0);
        difference.min(1.0 - difference) <= atol
    })
}
