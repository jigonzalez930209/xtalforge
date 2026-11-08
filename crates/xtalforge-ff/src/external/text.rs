//! Shared text helpers for the external-engine writers and parsers.

use crate::Error;

/// Appends one line, converting the value with `Display`.
pub(crate) fn push_line(output: &mut String, line: impl std::fmt::Display) {
    output.push_str(&line.to_string());
    output.push('\n');
}

/// Rejects empty tokens and characters that could break a one-line format.
pub(crate) fn validate_token(label: &str, token: &str) -> Result<(), Error> {
    if token.trim().is_empty() {
        return Err(Error::InvalidInput(format!("{label} must not be empty")));
    }
    if !token.chars().all(|character| {
        character.is_ascii_alphanumeric()
            || matches!(character, '-' | '_' | '+' | '.' | '/' | '(' | ')' | '*')
    }) {
        return Err(Error::InvalidInput(format!(
            "{label} contains a character that is not safe in an input line"
        )));
    }
    Ok(())
}

/// Rejects a non-finite parsed quantity.
pub(crate) fn finite(label: &str, value: f64) -> Result<(), Error> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(Error::InvalidInput(format!("{label} is not finite")))
    }
}

/// Numbers in one line, after trimming punctuation the engines put around
/// them.
pub(crate) fn numbers(text: &str) -> Vec<f64> {
    text.split_whitespace()
        .filter_map(|token| {
            let trimmed = token.trim_matches(|character: char| {
                !(character.is_ascii_digit() || matches!(character, '-' | '+' | '.' | 'e' | 'E'))
            });
            trimmed.parse::<f64>().ok()
        })
        .collect()
}

/// First line containing `keyword`, case-insensitively.
pub(crate) fn find_line(lines: &[&str], keyword: &str) -> Option<usize> {
    let needle = keyword.to_ascii_lowercase();
    lines
        .iter()
        .position(|line| line.to_ascii_lowercase().contains(&needle))
}

/// First number after `keyword` on the first line that contains it.
pub(crate) fn find_number(lines: &[&str], keyword: &str) -> Option<f64> {
    let index = find_line(lines, keyword)?;
    let line = lines.get(index)?;
    let lower = line.to_ascii_lowercase();
    let needle = keyword.to_ascii_lowercase();
    let start = lower.find(&needle)? + needle.len();
    numbers(line.get(start..)?).into_iter().next()
}

/// Reads indexed vector rows (`index x y z`) after `start`, in the given unit
/// scale.  Stops at a blank line or at the first line that is not a row.
pub(crate) fn parse_vectors(lines: &[&str], start: usize, scale: f64) -> Vec<[f64; 3]> {
    let mut vectors = Vec::new();
    for line in lines.iter().skip(start) {
        if line.trim().is_empty() {
            break;
        }
        let values = numbers(line);
        let components = if values.len() >= 4 {
            [
                values.get(1).copied().unwrap_or(0.0),
                values.get(2).copied().unwrap_or(0.0),
                values.get(3).copied().unwrap_or(0.0),
            ]
        } else if values.len() == 3 {
            [
                values.first().copied().unwrap_or(0.0),
                values.get(1).copied().unwrap_or(0.0),
                values.get(2).copied().unwrap_or(0.0),
            ]
        } else {
            if vectors.is_empty() {
                continue;
            }
            break;
        };
        let vector = [
            components[0] * scale,
            components[1] * scale,
            components[2] * scale,
        ];
        if vector.iter().any(|value| !value.is_finite()) {
            break;
        }
        vectors.push(vector);
    }
    vectors
}
