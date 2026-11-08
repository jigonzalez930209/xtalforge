//! Fullerene geometry from the built-in C60 coordinates (Task 6.7).

use serde::Deserialize;
use xtalforge_core::{Error, Structure};

#[derive(Debug, Deserialize)]
struct Data {
    name: String,
    points: Vec<[f64; 3]>,
}

const C60: &str = include_str!("../data/c60.json");

/// Builds a C60 buckyball in a padded P1 box.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a size other than 60 or invalid data.
pub fn fullerene(size: u32) -> Result<Structure, Error> {
    if size != 60 {
        return Err(Error::InvalidInput(format!(
            "only C60 is available, got C{size}"
        )));
    }
    let data: Data = serde_json::from_str(C60)
        .map_err(|error| Error::InvalidInput(format!("C60 data: {error}")))?;
    if data.name != "C60" || data.points.len() != 60 {
        return Err(Error::InvalidInput("C60 data is incomplete".to_owned()));
    }
    let points: Vec<(String, [f64; 3])> = data
        .points
        .iter()
        .map(|point| ("C".to_owned(), *point))
        .collect();
    Structure::molecule("C60", &points, 5.0)
}
