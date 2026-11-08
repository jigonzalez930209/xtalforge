//! Whitespace-delimited XY patterns (Task 3.5.3).

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked above
use xtalforge_core::Error;

/// One row of a two-theta/intensity pattern.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct PatternPoint {
    /// Two-theta in degrees.
    pub two_theta: f64,
    /// Intensity (arbitrary units).
    pub intensity: f64,
    /// Optional uncertainty on the intensity.
    pub uncertainty: Option<f64>,
}

/// Reads an XY pattern, ignoring `#` comments and blank lines.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a row without two numeric columns.
pub fn read_xy(text: &str) -> Result<Vec<PatternPoint>, Error> {
    let mut points = Vec::new();
    for (line_index, line) in text.lines().enumerate() {
        let line = line.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        let values: Vec<f64> = line
            .split_whitespace()
            .filter_map(|value| value.parse().ok())
            .collect();
        let (Some(two_theta), Some(intensity)) = (values.first(), values.get(1)) else {
            return Err(Error::InvalidInput(format!(
                "line {}: an XY row needs two numbers",
                line_index + 1
            )));
        };
        points.push(PatternPoint {
            two_theta: *two_theta,
            intensity: *intensity,
            uncertainty: values.get(2).copied(),
        });
    }
    Ok(points)
}

/// Writes an XY pattern with a comment header.
#[must_use]
pub fn write_xy(points: &[PatternPoint]) -> String {
    let mut out = String::from("# 2theta intensity uncertainty\n");
    for point in points {
        match point.uncertainty {
            Some(uncertainty) => out.push_str(&format!(
                "{:.6} {:.6} {:.6}\n",
                point.two_theta, point.intensity, uncertainty
            )),
            None => out.push_str(&format!("{:.6} {:.6}\n", point.two_theta, point.intensity)),
        }
    }
    out
}
