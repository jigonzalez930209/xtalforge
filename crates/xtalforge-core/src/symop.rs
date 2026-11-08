//! Symmetry operations with exact integer rotations and rational translations.

use std::fmt;

use serde::{Deserialize, Serialize};

use crate::Error;

/// One symmetry operation: `x' = R x + t` in fractional coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SymOp {
    rotation: [[i32; 3]; 3],
    translation: [f64; 3],
}

impl SymOp {
    /// Creates an operation and wraps the translation into `[0, 1)`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a singular rotation or non-finite
    /// translation.
    pub fn new(rotation: [[i32; 3]; 3], translation: [f64; 3]) -> Result<Self, Error> {
        if !translation.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "symmetry translation must be finite".to_owned(),
            ));
        }
        if determinant(rotation) == 0 {
            return Err(Error::InvalidInput(
                "symmetry rotation is singular".to_owned(),
            ));
        }
        Ok(Self {
            rotation,
            translation: translation.map(|value| value.rem_euclid(1.0)),
        })
    }

    /// Identity operation.
    #[must_use]
    pub const fn identity() -> Self {
        Self {
            rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            translation: [0.0, 0.0, 0.0],
        }
    }

    /// Integer rotation matrix.
    #[must_use]
    pub const fn rotation(&self) -> [[i32; 3]; 3] {
        self.rotation
    }

    /// Fractional translation in `[0, 1)`.
    #[must_use]
    pub const fn translation(&self) -> [f64; 3] {
        self.translation
    }

    /// Applies the operation to one fractional point.
    #[must_use]
    pub fn apply(&self, frac: [f64; 3]) -> [f64; 3] {
        let [x, y, z] = frac;
        let [[a, b, c], [d, e, f], [g, h, i]] = self.rotation;
        let [tx, ty, tz] = self.translation;
        [
            a as f64 * x + b as f64 * y + c as f64 * z + tx,
            d as f64 * x + e as f64 * y + f as f64 * z + ty,
            g as f64 * x + h as f64 * y + i as f64 * z + tz,
        ]
    }

    /// Applies the operation to a block of points.
    #[must_use]
    pub fn apply_many(&self, points: &[[f64; 3]]) -> Vec<[f64; 3]> {
        points.iter().map(|point| self.apply(*point)).collect()
    }

    /// Composition `self ∘ other`: apply `other`, then `self`.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`SymOp::new`].
    pub fn compose(&self, other: &Self) -> Result<Self, Error> {
        let rotation = multiply_int(self.rotation, other.rotation);
        let shifted = apply_int(self.rotation, other.translation);
        let translation = [
            shifted[0] + self.translation[0],
            shifted[1] + self.translation[1],
            shifted[2] + self.translation[2],
        ];
        Self::new(rotation, translation)
    }

    /// Inverse operation.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the rotation determinant is not
    /// `±1`.
    pub fn inverse(&self) -> Result<Self, Error> {
        let det = determinant(self.rotation);
        if det.abs() != 1 {
            return Err(Error::InvalidInput(
                "symmetry rotation determinant must be ±1".to_owned(),
            ));
        }
        let inverse = inverse_int(self.rotation, det)?;
        let shifted = apply_int(inverse, self.translation);
        Self::new(inverse, [-shifted[0], -shifted[1], -shifted[2]])
    }

    /// Whether this is the identity operation.
    #[must_use]
    pub fn is_identity(&self) -> bool {
        self.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            && self.translation.iter().all(|value| *value == 0.0)
    }

    /// Parses a triplet such as `-x+1/2,y,z`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for the wrong number of components or a
    /// malformed term.
    pub fn parse(text: &str) -> Result<Self, Error> {
        let components: Vec<&str> = text.split(',').map(str::trim).collect();
        if components.len() != 3 {
            return Err(Error::InvalidInput(format!(
                "a symmetry triplet needs three comma-separated components: {text:?}"
            )));
        }
        let mut rotation = [[0i32; 3]; 3];
        let mut translation = [0.0; 3];
        for (row, component) in components.iter().enumerate() {
            let (line, constant) = parse_component(component)?;
            if let (Some(line_slot), Some(constant_slot)) =
                (rotation.get_mut(row), translation.get_mut(row))
            {
                *line_slot = line;
                *constant_slot = constant;
            }
        }
        Self::new(rotation, translation)
    }

    /// Formats the operation as a triplet.
    #[must_use]
    pub fn triplet(&self) -> String {
        (0..3)
            .map(|row| format_component(&self.rotation, self.translation, row))
            .collect::<Vec<_>>()
            .join(",")
    }
}

impl fmt::Display for SymOp {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.triplet())
    }
}

fn parse_component(text: &str) -> Result<([i32; 3], f64), Error> {
    let compact: String = text
        .chars()
        .filter(|character| !character.is_whitespace())
        .collect();
    if compact.is_empty() {
        return Err(Error::InvalidInput("empty symmetry component".to_owned()));
    }
    let mut rotation = [0i32; 3];
    let mut translation = 0.0;
    let mut sign = 1.0;
    let mut term = String::new();
    for character in compact.chars() {
        if character == '+' || character == '-' {
            if !term.is_empty() {
                add_term(&term, sign, &mut rotation, &mut translation)?;
                term.clear();
            }
            sign = if character == '-' { -1.0 } else { 1.0 };
        } else {
            term.push(character);
        }
    }
    if term.is_empty() {
        return Err(Error::InvalidInput(format!(
            "malformed symmetry component: {text:?}"
        )));
    }
    add_term(&term, sign, &mut rotation, &mut translation)?;
    Ok((rotation, translation))
}

fn add_term(
    term: &str,
    sign: f64,
    rotation: &mut [i32; 3],
    translation: &mut f64,
) -> Result<(), Error> {
    if let Some(index) = term.find(['x', 'y', 'z']) {
        let axis = match term.get(index..index + 1) {
            Some("x") => 0,
            Some("y") => 1,
            Some("z") => 2,
            _ => 0,
        };
        let coefficient = term.get(..index).unwrap_or("");
        let after = term.get(index + 1..).unwrap_or("");
        let mut value = if coefficient.is_empty() || coefficient == "+" {
            1.0
        } else if coefficient == "-" {
            -1.0
        } else {
            coefficient.parse::<f64>().map_err(|_| {
                Error::InvalidInput(format!("malformed symmetry coefficient: {term:?}"))
            })?
        };
        if let Some(denominator) = after.strip_prefix('/') {
            let denominator: f64 = denominator.parse().map_err(|_| {
                Error::InvalidInput(format!("malformed symmetry denominator: {term:?}"))
            })?;
            if denominator == 0.0 {
                return Err(Error::InvalidInput(
                    "symmetry denominator must not be zero".to_owned(),
                ));
            }
            value /= denominator;
        }
        let value = (sign * value).round() as i32;
        if let Some(slot) = rotation.get_mut(axis) {
            *slot += value;
        }
        return Ok(());
    }
    let value = parse_rational(term)?;
    *translation += sign * value;
    Ok(())
}

fn parse_rational(text: &str) -> Result<f64, Error> {
    if let Some((numerator, denominator)) = text.split_once('/') {
        let numerator = numerator
            .parse::<f64>()
            .map_err(|_| Error::InvalidInput(format!("malformed rational constant: {text:?}")))?;
        let denominator = denominator.parse::<f64>().map_err(|_| {
            Error::InvalidInput(format!("malformed rational denominator: {text:?}"))
        })?;
        if denominator == 0.0 {
            return Err(Error::InvalidInput(
                "symmetry denominator must not be zero".to_owned(),
            ));
        }
        Ok(numerator / denominator)
    } else {
        text.parse::<f64>()
            .map_err(|_| Error::InvalidInput(format!("malformed constant: {text:?}")))
    }
}

fn format_component(rotation: &[[i32; 3]; 3], translation: [f64; 3], row: usize) -> String {
    let mut out = String::new();
    let line = rotation.get(row).copied().unwrap_or([0; 3]);
    for (axis, coefficient) in line.iter().enumerate() {
        let coefficient = *coefficient;
        if coefficient == 0 {
            continue;
        }
        let negative = coefficient < 0;
        if out.is_empty() {
            if negative {
                out.push('-');
            }
        } else {
            out.push(if negative { '-' } else { '+' });
        }
        let magnitude = coefficient.unsigned_abs();
        if magnitude != 1 {
            out.push_str(&magnitude.to_string());
        }
        out.push(match axis {
            0 => 'x',
            1 => 'y',
            _ => 'z',
        });
    }
    let constant = translation.get(row).copied().unwrap_or(0.0);
    if constant != 0.0 {
        let text = fraction(constant);
        if out.is_empty() {
            out.push_str(&text);
        } else {
            out.push('+');
            out.push_str(&text);
        }
    }
    if out.is_empty() {
        out.push('0');
    }
    out
}

fn fraction(value: f64) -> String {
    let value = value.rem_euclid(1.0);
    for denominator in 1..=12u32 {
        let numerator = (value * f64::from(denominator)).round();
        let candidate = numerator / f64::from(denominator);
        if (candidate - value).abs() < 1e-9 {
            if denominator == 1 {
                return format!("{}", numerator as i64);
            }
            let divisor = gcd(numerator as i64, i64::from(denominator));
            return format!(
                "{}/{}",
                numerator as i64 / divisor,
                i64::from(denominator) / divisor
            );
        }
    }
    format!("{value:.6}")
}

fn gcd(mut first: i64, mut second: i64) -> i64 {
    while second != 0 {
        let remainder = first % second;
        first = second;
        second = remainder;
    }
    first.abs().max(1)
}

fn determinant(matrix: [[i32; 3]; 3]) -> i32 {
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
}

fn multiply_int(first: [[i32; 3]; 3], second: [[i32; 3]; 3]) -> [[i32; 3]; 3] {
    let [[a, b, c], [d, e, f], [g, h, i]] = first;
    let [[j, k, l], [m, n, o], [p, q, r]] = second;
    [
        [
            a * j + b * m + c * p,
            a * k + b * n + c * q,
            a * l + b * o + c * r,
        ],
        [
            d * j + e * m + f * p,
            d * k + e * n + f * q,
            d * l + e * o + f * r,
        ],
        [
            g * j + h * m + i * p,
            g * k + h * n + i * q,
            g * l + h * o + i * r,
        ],
    ]
}

fn apply_int(matrix: [[i32; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    let [x, y, z] = vector;
    [
        a as f64 * x + b as f64 * y + c as f64 * z,
        d as f64 * x + e as f64 * y + f as f64 * z,
        g as f64 * x + h as f64 * y + i as f64 * z,
    ]
}

fn inverse_int(matrix: [[i32; 3]; 3], det: i32) -> Result<[[i32; 3]; 3], Error> {
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    let adjugate = [
        [e * i - f * h, c * h - b * i, b * f - c * e],
        [f * g - d * i, a * i - c * g, c * d - a * f],
        [d * h - e * g, b * g - a * h, a * e - b * d],
    ];
    let mut inverse = [[0i32; 3]; 3];
    for (output_row, adjudicate_row) in inverse.iter_mut().zip(adjugate) {
        for (output, value) in output_row.iter_mut().zip(adjudicate_row) {
            if value % det != 0 {
                return Err(Error::InvalidInput(
                    "symmetry rotation is not unimodular".to_owned(),
                ));
            }
            *output = value / det;
        }
    }
    Ok(inverse)
}
