//! Wyckoff letters in the input origin (functional audit 2026-11-11).
//!
//! moyo reports the Wyckoff letter of each atom in its standardized cell,
//! whose origin may differ from the input by an element of the Euclidean
//! normalizer. Silicon in origin choice 1 (`F 4d 2 3 -1d`) with Si at
//! `0, 0, 0` then comes out as `8b` although the input is literally the
//! `8a` site, which is what spglib, the CIF and International Tables say.
//!
//! When the input basis already is the conventional basis (`std_linear` is
//! the identity, or a signed axis permutation in a cubic group), the letter is taken from the input coordinates instead:
//! every Wyckoff position of the detected Hall setting with the same
//! multiplicity is tried in table order, and a position matches when some
//! operation of the group maps the atom onto its representative coordinate
//! template (free parameters allowed, modulo lattice translations). When no
//! position matches, or the basis differs, moyo's letter is kept.

#![allow(clippy::indexing_slicing)] // fixed 3 x 3 arrays indexed by 0..3

use xtalforge_core::SymOp;

use crate::spacegroup::WyckoffPosition;

/// Linear template `A v + c` of one coordinate triplet.
#[derive(Debug, Clone, Copy)]
struct Template {
    matrix: [[f64; 3]; 3],
    constant: [f64; 3],
}

fn parse_component(text: &str) -> Option<([f64; 3], f64)> {
    let mut coefficients = [0.0; 3];
    let mut constant = 0.0;
    let mut term = String::new();
    let mut flush = |term: &mut String| -> Option<()> {
        if term.is_empty() || term == "+" || term == "-" {
            term.clear();
            return Some(());
        }
        let (sign, body) = match term.strip_prefix('-') {
            Some(rest) => (-1.0, rest),
            None => (1.0, term.trim_start_matches('+')),
        };
        let variable = body.chars().last().and_then(|last| match last {
            'x' => Some(0),
            'y' => Some(1),
            'z' => Some(2),
            _ => None,
        });
        match variable {
            Some(axis) => {
                let factor = &body[..body.len() - 1];
                let factor = if factor.is_empty() {
                    1.0
                } else {
                    factor.parse::<f64>().ok()?
                };
                *coefficients.get_mut(axis)? += sign * factor;
            }
            None => {
                let value = match body.split_once('/') {
                    Some((top, bottom)) => top.parse::<f64>().ok()? / bottom.parse::<f64>().ok()?,
                    None => body.parse::<f64>().ok()?,
                };
                constant += sign * value;
            }
        }
        term.clear();
        Some(())
    };
    for character in text.chars().filter(|character| !character.is_whitespace()) {
        if (character == '+' || character == '-') && !term.is_empty() {
            flush(&mut term)?;
        }
        term.push(character);
    }
    flush(&mut term)?;
    Some((coefficients, constant))
}

fn parse_template(text: &str) -> Option<Template> {
    let parts: Vec<&str> = text.split(',').collect();
    if parts.len() != 3 {
        return None;
    }
    let mut matrix = [[0.0; 3]; 3];
    let mut constant = [0.0; 3];
    for (row, part) in parts.iter().enumerate() {
        let (coefficients, value) = parse_component(part)?;
        *matrix.get_mut(row)? = coefficients;
        *constant.get_mut(row)? = value;
    }
    Some(Template { matrix, constant })
}

fn solve3(matrix: [[f64; 3]; 3], rhs: [f64; 3]) -> Option<[f64; 3]> {
    let det = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let base = det(matrix);
    if base.abs() < 1e-30 {
        return None;
    }
    let mut out = [0.0; 3];
    for (column, slot) in out.iter_mut().enumerate() {
        let mut replaced = matrix;
        for row in 0..3 {
            replaced[row][column] = rhs[row];
        }
        *slot = det(replaced) / base;
    }
    Some(out)
}

/// Whether `point` lies on the template modulo lattice translations.
fn on_template(template: &Template, point: [f64; 3], tolerance: f64) -> bool {
    let a = template.matrix;
    // Normal equations with a small ridge handle rank-deficient templates.
    let mut normal = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            normal[i][j] =
                (0..3).map(|k| a[k][i] * a[k][j]).sum::<f64>() + if i == j { 1e-6 } else { 0.0 };
        }
    }
    for nx in -2..=2 {
        for ny in -2..=2 {
            for nz in -2..=2 {
                let shift = [f64::from(nx), f64::from(ny), f64::from(nz)];
                let target: [f64; 3] =
                    [0, 1, 2].map(|k| point[k] + shift[k] - template.constant[k]);
                let rhs = [0, 1, 2].map(|i| (0..3).map(|k| a[k][i] * target[k]).sum::<f64>());
                let Some(v) = solve3(normal, rhs) else {
                    continue;
                };
                let fits = (0..3).all(|k| {
                    let value: f64 = (0..3).map(|j| a[k][j] * v[j]).sum();
                    (value - target[k]).abs() < tolerance
                });
                if fits {
                    return true;
                }
            }
        }
    }
    false
}

/// Letters in the input origin, falling back to `moyo_letters`.
pub(crate) fn input_origin_letters(
    operations: &[SymOp],
    hall_number: u32,
    std_linear: [[f64; 3]; 3],
    positions: &[[f64; 3]],
    moyo_letters: &[char],
) -> Vec<char> {
    let identity = (0..3)
        .all(|i| (0..3).all(|j| (std_linear[i][j] - if i == j { 1.0 } else { 0.0 }).abs() < 1e-6));
    // A cubic group contains the cyclic axis permutations, so its templates
    // are the same in any signed permutation of the conventional axes.
    let signed_permutation = std_linear.iter().all(|row| {
        row.iter()
            .filter(|value| (value.abs() - 1.0).abs() < 1e-6)
            .count()
            == 1
            && row.iter().filter(|value| value.abs() < 1e-6).count() == 2
    });
    let cubic = crate::spacegroup::SpaceGroup::from_hall_number(hall_number)
        .is_ok_and(|group| (195..=230).contains(&group.number));
    let identity = identity || (signed_permutation && cubic);
    let Ok(table) = WyckoffPosition::list(hall_number) else {
        return moyo_letters.to_vec();
    };
    if !identity || operations.is_empty() {
        return moyo_letters.to_vec();
    }
    let parsed: Vec<(char, usize, Option<Template>)> = table
        .iter()
        .map(|position| {
            (
                position.letter,
                position.multiplicity,
                parse_template(&position.coordinates),
            )
        })
        .collect();
    let multiplicity_of = |letter: char| {
        parsed
            .iter()
            .find(|entry| entry.0 == letter)
            .map(|entry| entry.1)
    };
    // Table order runs from the general position down to `a`; scan from `a`.
    let mut ordered = parsed.clone();
    ordered.sort_by_key(|entry| entry.0);
    positions
        .iter()
        .zip(moyo_letters)
        .map(|(point, letter)| {
            let Some(multiplicity) = multiplicity_of(*letter) else {
                return *letter;
            };
            ordered
                .iter()
                .filter(|(_, size, template)| *size == multiplicity && template.is_some())
                .find(|(_, _, template)| {
                    template.as_ref().is_some_and(|template| {
                        operations
                            .iter()
                            .any(|operation| on_template(template, operation.apply(*point), 2e-3))
                    })
                })
                .map_or(*letter, |entry| entry.0)
        })
        .collect()
}
