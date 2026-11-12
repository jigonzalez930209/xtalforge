#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Radial distribution functions, bond-angle distributions and fingerprints
//! (Task 9.3.4).

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::{Error, Structure, neighbors};

/// Options for [`radial_distribution`].
#[derive(Debug, Clone, PartialEq)]
pub struct RdfOptions {
    /// Maximum distance in Angstrom.
    pub r_max: f64,
    /// Histogram bin width in Angstrom.
    pub bin_width: f64,
    /// Moving-average window applied to every curve; zero disables it.
    pub smoothing: usize,
}

impl Default for RdfOptions {
    fn default() -> Self {
        Self {
            r_max: 10.0,
            bin_width: 0.05,
            smoothing: 0,
        }
    }
}

/// A total and partial radial distribution function.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RadialDistribution {
    /// Bin centres in Angstrom.
    pub r: Vec<f64>,
    /// Total `g(r)`.
    pub total: Vec<f64>,
    /// Partial `g_ab(r)` per unordered element pair, key `A-B`.
    pub partial: BTreeMap<String, Vec<f64>>,
    /// Running coordination number from the total `g(r)`.
    pub coordination: Vec<f64>,
    /// Number density in atoms per cubic Angstrom.
    pub density: f64,
}

/// Computes the total and partial RDF of a structure.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-positive range or bin width, an
/// empty structure, or the errors of [`neighbors()`].
pub fn radial_distribution(
    structure: &Structure,
    options: &RdfOptions,
) -> Result<RadialDistribution, Error> {
    if !options.r_max.is_finite() || options.r_max <= 0.0 {
        return Err(Error::InvalidInput(
            "RDF range must be positive and finite".to_owned(),
        ));
    }
    if !options.bin_width.is_finite() || options.bin_width <= 0.0 {
        return Err(Error::InvalidInput(
            "RDF bin width must be positive and finite".to_owned(),
        ));
    }
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "the RDF needs at least one atom".to_owned(),
        ));
    }
    let bins = (options.r_max / options.bin_width).ceil() as usize;
    let volume = p1.lattice.volume();
    let density = p1.atoms.len() as f64 / volume;
    let mut total = vec![0.0; bins];
    let mut partial: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    let mut element_counts: BTreeMap<String, usize> = BTreeMap::new();
    for atom in &p1.atoms {
        *element_counts.entry(atom.element.clone()).or_default() += 1;
    }
    for (first, second) in element_counts.keys().enumerate() {
        for other in element_counts.keys().skip(first) {
            let key = pair_key(second, other);
            partial.insert(key, vec![0.0; bins]);
        }
    }
    let snapshot = p1_structure(&p1)?;
    for pair in neighbors(&snapshot, options.r_max)? {
        let Some(bin) = bin_of(pair.distance, options.bin_width, bins) else {
            continue;
        };
        let (Some(a), Some(b)) = (p1.atoms.get(pair.i), p1.atoms.get(pair.j)) else {
            continue;
        };
        if let Some(entry) = total.get_mut(bin) {
            *entry += 2.0;
        }
        let key = pair_key(&a.element, &b.element);
        if a.element == b.element {
            if let Some(entry) = partial.get_mut(&key).and_then(|data| data.get_mut(bin)) {
                *entry += 2.0;
            }
        } else {
            if let Some(entry) = partial.get_mut(&key).and_then(|data| data.get_mut(bin)) {
                *entry += 1.0;
            }
            let swapped = pair_key(&b.element, &a.element);
            if let Some(entry) = partial.get_mut(&swapped).and_then(|data| data.get_mut(bin)) {
                *entry += 1.0;
            }
        }
    }
    let mut r = Vec::with_capacity(bins);
    let mut coordination = Vec::with_capacity(bins);
    let mut running = 0.0;
    for index in 0..bins {
        let radius = (index as f64 + 0.5) * options.bin_width;
        r.push(radius);
        let shell = 4.0 * std::f64::consts::PI * radius * radius * options.bin_width;
        if let Some(entry) = total.get_mut(index) {
            *entry /= p1.atoms.len() as f64 * density * shell;
            running += *entry * density * shell;
        }
        coordination.push(running);
    }
    for (key, data) in &mut partial {
        let (Some(left), Some(right)) = (split_key(key, 0), split_key(key, 1)) else {
            continue;
        };
        let count_left = element_counts.get(left).copied().unwrap_or(0) as f64;
        let count_right = element_counts.get(right).copied().unwrap_or(0) as f64;
        if count_left == 0.0 || count_right == 0.0 {
            continue;
        }
        let rho_right = count_right / volume;
        for (index, value) in data.iter_mut().enumerate() {
            let radius = (index as f64 + 0.5) * options.bin_width;
            let shell = 4.0 * std::f64::consts::PI * radius * radius * options.bin_width;
            *value /= count_left * rho_right * shell;
        }
    }
    if options.smoothing > 0 {
        total = smooth(&total, options.smoothing);
        for data in partial.values_mut() {
            *data = smooth(data, options.smoothing);
        }
    }
    Ok(RadialDistribution {
        r,
        total,
        partial,
        coordination,
        density,
    })
}

/// Bond-angle histogram in degrees.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AngleDistribution {
    /// Bin centres in degrees.
    pub angle: Vec<f64>,
    /// Counts per bin.
    pub counts: Vec<f64>,
}

/// Computes the distribution of angles between neighbour pairs.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-positive cutoff or bin width, or
/// the errors of [`neighbors()`].
pub fn bond_angle_distribution(
    structure: &Structure,
    cutoff: f64,
    bin_width: f64,
) -> Result<AngleDistribution, Error> {
    if !cutoff.is_finite() || cutoff <= 0.0 {
        return Err(Error::InvalidInput(
            "bond-angle cutoff must be positive and finite".to_owned(),
        ));
    }
    if !bin_width.is_finite() || bin_width <= 0.0 {
        return Err(Error::InvalidInput(
            "bond-angle bin width must be positive and finite".to_owned(),
        ));
    }
    let p1 = structure.p1()?;
    let snapshot = p1_structure(&p1)?;
    let pairs = neighbors(&snapshot, cutoff)?;
    let matrix = p1.lattice.matrix();
    let mut shells: Vec<Vec<[f64; 3]>> = vec![Vec::new(); p1.atoms.len()];
    for pair in pairs {
        let (Some(a), Some(b)) = (p1.atoms.get(pair.i), p1.atoms.get(pair.j)) else {
            continue;
        };
        let base = p1.lattice.to_cart(b.frac);
        let image = [
            base[0]
                + f64::from(pair.image[0]) * matrix[0][0]
                + f64::from(pair.image[1]) * matrix[1][0]
                + f64::from(pair.image[2]) * matrix[2][0],
            base[1]
                + f64::from(pair.image[0]) * matrix[0][1]
                + f64::from(pair.image[1]) * matrix[1][1]
                + f64::from(pair.image[2]) * matrix[2][1],
            base[2]
                + f64::from(pair.image[0]) * matrix[0][2]
                + f64::from(pair.image[1]) * matrix[1][2]
                + f64::from(pair.image[2]) * matrix[2][2],
        ];
        let centre_a = p1.lattice.to_cart(a.frac);
        if let Some(entry) = shells.get_mut(pair.i) {
            entry.push([
                image[0] - centre_a[0],
                image[1] - centre_a[1],
                image[2] - centre_a[2],
            ]);
        }
        let centre_b = p1.lattice.to_cart(b.frac);
        let opposite = [
            centre_a[0]
                - f64::from(pair.image[0]) * matrix[0][0]
                - f64::from(pair.image[1]) * matrix[1][0]
                - f64::from(pair.image[2]) * matrix[2][0]
                - centre_b[0],
            centre_a[1]
                - f64::from(pair.image[0]) * matrix[0][1]
                - f64::from(pair.image[1]) * matrix[1][1]
                - f64::from(pair.image[2]) * matrix[2][1]
                - centre_b[1],
            centre_a[2]
                - f64::from(pair.image[0]) * matrix[0][2]
                - f64::from(pair.image[1]) * matrix[1][2]
                - f64::from(pair.image[2]) * matrix[2][2]
                - centre_b[2],
        ];
        if let Some(entry) = shells.get_mut(pair.j) {
            entry.push(opposite);
        }
    }
    let bins = (180.0 / bin_width).ceil() as usize;
    let mut counts = vec![0.0; bins];
    for vectors in &shells {
        for first in 0..vectors.len() {
            for second in (first + 1)..vectors.len() {
                let (Some(a), Some(b)) = (vectors.get(first), vectors.get(second)) else {
                    continue;
                };
                let norm_a = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
                let norm_b = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();
                if norm_a < 1e-12 || norm_b < 1e-12 {
                    continue;
                }
                let cosine = ((a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (norm_a * norm_b))
                    .clamp(-1.0, 1.0);
                let angle = cosine.acos().to_degrees();
                let bin = (angle / bin_width).floor() as usize;
                if let Some(entry) = counts.get_mut(bin.min(bins - 1)) {
                    *entry += 1.0;
                }
            }
        }
    }
    let angle = (0..bins)
        .map(|index| (index as f64 + 0.5) * bin_width)
        .collect();
    Ok(AngleDistribution { angle, counts })
}

/// A compact composition, cell and RDF fingerprint for quick similarity.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StructureFingerprint {
    /// Element fractions in the cell.
    pub composition: BTreeMap<String, f64>,
    /// Cell parameters `a b c alpha beta gamma`.
    pub cell: [f64; 6],
    /// Smoothed total `g(r)` sampled at a fixed grid.
    pub rdf: Vec<f64>,
}

impl StructureFingerprint {
    /// Cosine similarity: `1` for identical fingerprints, `0` for unrelated.
    #[must_use]
    pub fn similarity(&self, other: &Self) -> f64 {
        let composition = composition_similarity(&self.composition, &other.composition);
        let cell = cell_similarity(self.cell, other.cell);
        let rdf = cosine(&self.rdf, &other.rdf);
        0.4 * composition + 0.2 * cell + 0.4 * rdf
    }
}

/// Builds the fingerprint of a structure.
///
/// # Errors
///
/// Returns the errors of [`radial_distribution`].
pub fn fingerprint(
    structure: &Structure,
    options: &RdfOptions,
) -> Result<StructureFingerprint, Error> {
    let distribution = radial_distribution(structure, options)?;
    let p1 = structure.p1()?;
    let total = p1.atoms.len().max(1) as f64;
    let mut counts: BTreeMap<String, f64> = BTreeMap::new();
    for atom in &p1.atoms {
        *counts.entry(atom.element.clone()).or_default() += 1.0;
    }
    let composition = counts
        .into_iter()
        .map(|(element, count)| (element, count / total))
        .collect();
    Ok(StructureFingerprint {
        composition,
        cell: structure.lattice().parameters(),
        rdf: distribution.total,
    })
}

fn composition_similarity(first: &BTreeMap<String, f64>, second: &BTreeMap<String, f64>) -> f64 {
    let mut keys: Vec<&String> = first.keys().chain(second.keys()).collect();
    keys.sort();
    keys.dedup();
    let mut left = Vec::new();
    let mut right = Vec::new();
    for key in keys {
        left.push(first.get(key).copied().unwrap_or(0.0));
        right.push(second.get(key).copied().unwrap_or(0.0));
    }
    cosine(&left, &right)
}

fn cell_similarity(first: [f64; 6], second: [f64; 6]) -> f64 {
    let mut difference = 0.0;
    for (a, b) in first.iter().zip(second.iter()) {
        let scale = a.abs().max(b.abs()).max(1e-6);
        difference += ((a - b) / scale).abs();
    }
    (1.0 - difference / 6.0).clamp(0.0, 1.0)
}

fn cosine(first: &[f64], second: &[f64]) -> f64 {
    let length = first.len().min(second.len());
    let mut dot = 0.0;
    let mut norm_first = 0.0;
    let mut norm_second = 0.0;
    for index in 0..length {
        let a = first[index];
        let b = second[index];
        dot += a * b;
        norm_first += a * a;
        norm_second += b * b;
    }
    if norm_first < 1e-18 || norm_second < 1e-18 {
        return 0.0;
    }
    (dot / (norm_first.sqrt() * norm_second.sqrt())).clamp(0.0, 1.0)
}

fn smooth(values: &[f64], window: usize) -> Vec<f64> {
    if window == 0 || values.is_empty() {
        return values.to_vec();
    }
    let radius = window / 2;
    let mut smoothed = Vec::with_capacity(values.len());
    for index in 0..values.len() {
        let start = index.saturating_sub(radius);
        let end = (index + radius + 1).min(values.len());
        let slice = values.get(start..end).unwrap_or(&[]);
        let mean = if slice.is_empty() {
            0.0
        } else {
            slice.iter().sum::<f64>() / slice.len() as f64
        };
        smoothed.push(mean);
    }
    smoothed
}

fn pair_key(first: &str, second: &str) -> String {
    if first <= second {
        format!("{first}-{second}")
    } else {
        format!("{second}-{first}")
    }
}

fn split_key(key: &str, index: usize) -> Option<&str> {
    key.split('-').nth(index)
}

fn bin_of(distance: f64, width: f64, bins: usize) -> Option<usize> {
    if !distance.is_finite() || distance <= 0.0 {
        return None;
    }
    let index = (distance / width).floor() as usize;
    if index < bins { Some(index) } else { None }
}

fn p1_structure(p1: &crate::P1Cell) -> Result<Structure, Error> {
    let mut builder = crate::StructureBuilder::new("p1", p1.lattice);
    for atom in &p1.atoms {
        builder = builder.site(
            crate::Site::new(atom.label.clone(), &atom.element, atom.frac)?
                .with_occupancy(atom.occupancy)?,
        );
    }
    builder.build()
}
