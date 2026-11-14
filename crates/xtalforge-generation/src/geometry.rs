#![allow(clippy::indexing_slicing)] // fixed three-component vectors

//! Geometry primitives shared by the generators: seeded lattices, orbit
//! deduplication, periodic minimum distances and rigid rotations.

use xtalforge_core::{Error, Lattice, Structure, SymOp};

use crate::rng::SplitMix64;

/// Cartesian distance below which two images are the same point.
pub(crate) const SAME_POINT: f64 = 0.05;

/// Builds a conventional lattice whose parameters match the crystal system.
///
/// `volume` is the target cell volume in cubic Angstrom; the axial ratios are
/// drawn from `rng` inside a documented range so the cell is not accidentally
/// special. The angles are fixed by the crystal system except for the
/// monoclinic `beta` and the triclinic triple, which are drawn.
pub(crate) fn lattice_for_system(
    system: &str,
    volume: f64,
    rng: &mut SplitMix64,
) -> Result<Lattice, Error> {
    if !volume.is_finite() || volume <= 0.0 {
        return Err(Error::InvalidInput(
            "the target cell volume must be positive and finite".to_owned(),
        ));
    }
    let ratio_a = rng.range(0.85, 1.18);
    let ratio_b = rng.range(0.85, 1.18);
    let ratio_c = rng.range(0.85, 1.18);
    let cubic = volume.cbrt();
    let (a, b, c, alpha, beta, gamma) = match system {
        "cubic" => (cubic, cubic, cubic, 90.0, 90.0, 90.0),
        "tetragonal" => {
            let c_over_a = rng.range(0.75, 1.45);
            let a = (volume / c_over_a).cbrt();
            (a, a, a * c_over_a, 90.0, 90.0, 90.0)
        }
        "hexagonal" => {
            let c_over_a = rng.range(0.9, 1.7);
            let a = (volume / (3.0_f64.sqrt() / 2.0 * c_over_a)).cbrt();
            (a, a, a * c_over_a, 90.0, 90.0, 120.0)
        }
        "trigonal" => {
            let c_over_a = rng.range(0.9, 1.9);
            let a = (volume / (3.0_f64.sqrt() / 2.0 * c_over_a)).cbrt();
            (a, a, a * c_over_a, 90.0, 90.0, 120.0)
        }
        "orthorhombic" => {
            let abc = volume.cbrt();
            (
                abc * ratio_a,
                abc * ratio_b,
                abc * ratio_c,
                90.0,
                90.0,
                90.0,
            )
        }
        "monoclinic" => {
            let abc = volume.cbrt();
            (
                abc * ratio_a,
                abc * ratio_b,
                abc * ratio_c,
                90.0,
                rng.range(96.0, 124.0),
                90.0,
            )
        }
        _ => {
            let abc = volume.cbrt();
            (
                abc * ratio_a,
                abc * ratio_b,
                abc * ratio_c,
                rng.range(72.0, 108.0),
                rng.range(72.0, 108.0),
                rng.range(72.0, 108.0),
            )
        }
    };
    Lattice::from_parameters(a, b, c, alpha, beta, gamma)
}

/// Wraps fractional coordinates into `[0, 1)`.
pub(crate) fn wrap(mut frac: [f64; 3]) -> [f64; 3] {
    for value in &mut frac {
        let wrapped = value.rem_euclid(1.0);
        // `rem_euclid` can return exactly 1.0 after rounding for values just
        // below zero; folding it keeps the half-open contract.
        *value = if wrapped >= 1.0 { 0.0 } else { wrapped };
    }
    frac
}

/// Appends `point` unless it is within [`SAME_POINT`] of a stored point.
pub(crate) fn push_unique(points: &mut Vec<[f64; 3]>, lattice: &Lattice, point: [f64; 3]) -> bool {
    let point = wrap(point);
    for existing in points.iter() {
        let (_, distance) = lattice.minimum_image([
            point[0] - existing[0],
            point[1] - existing[1],
            point[2] - existing[2],
        ]);
        if distance < SAME_POINT {
            return false;
        }
    }
    points.push(point);
    true
}

/// Expands points through operations, keeping one image per orbit.
pub(crate) fn orbit(points: &[[f64; 3]], operations: &[SymOp], lattice: &Lattice) -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    for point in points {
        for operation in operations {
            push_unique(&mut out, lattice, operation.apply(*point));
        }
    }
    out
}

/// Smallest Cartesian distance between any two structure atoms, minimum-image.
///
/// Returns `f64::INFINITY` for fewer than two atoms.
pub(crate) fn minimum_distance(structure: &Structure) -> Result<f64, Error> {
    let p1 = structure.p1()?;
    let lattice = &p1.lattice;
    let mut best = f64::INFINITY;
    for (index, first) in p1.atoms.iter().enumerate() {
        for second in p1.atoms.iter().skip(index + 1) {
            let (_, distance) = lattice.minimum_image([
                first.frac[0] - second.frac[0],
                first.frac[1] - second.frac[1],
                first.frac[2] - second.frac[2],
            ]);
            if distance < best && distance > 1e-6 {
                best = distance;
            }
        }
    }
    Ok(best)
}

/// Smallest distance between atoms belonging to different groups.
///
/// `group[i]` labels atom `i`; atoms with the same label are excluded. This is
/// how an intermolecular clash test ignores the short bonds inside a molecule.
pub(crate) fn minimum_intergroup_distance(
    positions: &[[f64; 3]],
    groups: &[usize],
    lattice: &Lattice,
) -> f64 {
    let mut best = f64::INFINITY;
    for (index, first) in positions.iter().enumerate() {
        for (other, second) in positions.iter().enumerate().skip(index + 1) {
            if groups.get(index) == groups.get(other) {
                continue;
            }
            let (_, distance) = lattice.minimum_image([
                first[0] - second[0],
                first[1] - second[1],
                first[2] - second[2],
            ]);
            if distance < best && distance > 1e-6 {
                best = distance;
            }
        }
    }
    best
}

/// Whether every operation maps the expanded atoms onto themselves.
pub(crate) fn operations_preserve(
    structure: &Structure,
    operations: &[SymOp],
) -> Result<bool, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Ok(false);
    }
    for operation in operations {
        for atom in &p1.atoms {
            let image = wrap(operation.apply(atom.frac));
            let mut found = false;
            for other in &p1.atoms {
                if other.element != atom.element {
                    continue;
                }
                let (_, distance) = p1.lattice.minimum_image([
                    image[0] - other.frac[0],
                    image[1] - other.frac[1],
                    image[2] - other.frac[2],
                ]);
                if distance < 0.05 {
                    found = true;
                    break;
                }
            }
            if !found {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// Rotation matrix from a unit quaternion `(w, x, y, z)`.
pub(crate) fn quaternion_matrix([w, x, y, z]: [f64; 4]) -> [[f64; 3]; 3] {
    let norm = (w * w + x * x + y * y + z * z).sqrt();
    let (w, x, y, z) = if norm > 0.0 {
        (w / norm, x / norm, y / norm, z / norm)
    } else {
        (1.0, 0.0, 0.0, 0.0)
    };
    [
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y - z * w),
            2.0 * (x * z + y * w),
        ],
        [
            2.0 * (x * y + z * w),
            1.0 - 2.0 * (x * x + z * z),
            2.0 * (y * z - x * w),
        ],
        [
            2.0 * (x * z - y * w),
            2.0 * (y * z + x * w),
            1.0 - 2.0 * (x * x + y * y),
        ],
    ]
}

/// Applies a 3x3 matrix to a Cartesian vector.
pub(crate) fn apply_matrix(matrix: [[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    matrix.map(|row| row[0] * vector[0] + row[1] * vector[1] + row[2] * vector[2])
}

/// Evaluates a linear coordinate template.
///
/// Moyo stores Wyckoff templates as comma-separated components such as `x`,
/// `1/2-x` or `-y+1/4`. The parser accepts a sum of terms, each a rational or
/// decimal constant, an optional integer coefficient and at most one variable
/// (`x`, `y` or `z`). Unknown syntax returns an error instead of guessing.
pub(crate) fn evaluate_template(template: &str, free: [f64; 3]) -> Result<[f64; 3], Error> {
    let parts: Vec<&str> = template.split(',').collect();
    if parts.len() != 3 {
        return Err(Error::InvalidInput(format!(
            "Wyckoff template `{template}` does not have three components"
        )));
    }
    let mut out = [0.0; 3];
    for (axis, part) in parts.iter().enumerate() {
        out[axis] = evaluate_component(part, free).map_err(|reason| {
            Error::InvalidInput(format!("Wyckoff template `{template}`: {reason}"))
        })?;
    }
    Ok(out)
}

fn evaluate_component(text: &str, free: [f64; 3]) -> Result<f64, String> {
    let compact: String = text.chars().filter(|c| !c.is_whitespace()).collect();
    if compact.is_empty() {
        return Err("empty component".to_owned());
    }
    let mut total = 0.0;
    let mut current = String::new();
    let mut sign = 1.0;
    let mut terms: Vec<(f64, String)> = Vec::new();
    for (index, character) in compact.char_indices() {
        if (character == '+' || character == '-') && index > 0 {
            terms.push((sign, std::mem::take(&mut current)));
            sign = if character == '-' { -1.0 } else { 1.0 };
        } else if index == 0 && (character == '+' || character == '-') {
            sign = if character == '-' { -1.0 } else { 1.0 };
        } else {
            current.push(character);
        }
    }
    terms.push((sign, current));
    for (sign, term) in terms {
        if term.is_empty() {
            return Err("empty term".to_owned());
        }
        if let Some(position) = term.find(['x', 'y', 'z']) {
            let (coefficient, variable) = term.split_at(position);
            if variable.len() != 1 {
                return Err(format!("`{term}` has more than one variable"));
            }
            let coefficient = if coefficient.is_empty() {
                1.0
            } else {
                parse_number(coefficient)?
            };
            let value = match variable {
                "x" => free[0],
                "y" => free[1],
                _ => free[2],
            };
            total += sign * coefficient * value;
        } else {
            total += sign * parse_number(&term)?;
        }
    }
    Ok(total)
}

fn parse_number(text: &str) -> Result<f64, String> {
    if let Some((numerator, denominator)) = text.split_once('/') {
        let numerator = numerator
            .parse::<f64>()
            .map_err(|_| format!("`{text}` is not a number"))?;
        let denominator = denominator
            .parse::<f64>()
            .map_err(|_| format!("`{text}` is not a number"))?;
        if denominator == 0.0 {
            return Err(format!("`{text}` divides by zero"));
        }
        Ok(numerator / denominator)
    } else {
        text.parse::<f64>()
            .map_err(|_| format!("`{text}` is not a number"))
    }
}
