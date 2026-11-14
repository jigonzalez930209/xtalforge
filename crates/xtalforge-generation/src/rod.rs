//! Rod groups and their generator (Task 20.1.3).
//!
//! moyo has no rod-group tables. This module ships a documented catalogue of
//! fourteen common rod groups as explicit operation sets and generates a
//! chain in a wide cell with an explicit repeat along `c`.

#![allow(clippy::indexing_slicing)] // fixed three-component vectors and bounded orbit arrays

use serde::{Deserialize, Serialize};
use serde_json::json;
use xtalforge_core::{Error, Lattice, Site, StructureBuilder, SymOp};

use crate::geometry::{minimum_distance, operations_preserve, push_unique};
use crate::lowdim::{LowDimCrystal, validate_common};
use crate::rng::SplitMix64;

/// One rod group.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct RodGroup {
    /// Catalogue number (1 to 14 in this subset).
    pub number: u32,
    /// Symbol used by this catalogue.
    pub symbol: String,
    /// Operations about the chain axis.
    pub operations: Vec<SymOp>,
}

/// The documented fourteen-rod-group subset.
///
/// # Errors
///
/// Returns [`Error::Symmetry`] when an operation cannot be built.
pub fn rod_groups() -> Result<Vec<RodGroup>, Error> {
    let mut groups = Vec::new();
    let entries: [(&str, Vec<[[i32; 3]; 3]>); 14] = [
        ("p1", vec![identity()]),
        ("p-1", vec![negate()]),
        ("p2", vec![rot2()]),
        ("p2/m", vec![rot2(), mirror_z()]),
        ("p222", vec![rot2(), rot2x(), rot2y()]),
        ("pmmm", vec![mirror_x(), mirror_y(), mirror_z()]),
        ("p4", vec![rot4()]),
        ("p4/m", vec![rot4(), mirror_z()]),
        ("p4mm", vec![rot4(), mirror_x(), mirror_y()]),
        ("p3", vec![rot3()]),
        ("p-3", vec![rot3(), negate()]),
        ("p6", vec![rot6()]),
        ("p6/m", vec![rot6(), mirror_z()]),
        ("p6mm", vec![rot6(), mirror_x(), mirror_y()]),
    ];
    for (index, (symbol, rotations)) in entries.into_iter().enumerate() {
        let mut operations = vec![SymOp::new(identity(), [0.0; 3])?];
        for rotation in rotations {
            let operation = SymOp::new(rotation, [0.0; 3])?;
            if !operations.contains(&operation) {
                operations.push(operation);
            }
        }
        groups.push(RodGroup {
            number: index as u32 + 1,
            symbol: symbol.to_owned(),
            operations,
        });
    }
    Ok(groups)
}

/// Options for [`generate_rod_crystal`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RodCrystalOptions {
    /// Rod-group number in the [`rod_groups`] catalogue.
    pub rod_group: u32,
    /// Element symbols and counts per repeat unit.
    pub composition: Vec<(String, u32)>,
    /// Repeat length along the chain in Angstrom.
    pub repeat: f64,
    /// Smallest allowed periodic distance in Angstrom.
    pub min_distance: f64,
    /// Number of attempts before the generator reports failure.
    pub max_attempts: usize,
    /// Explicit seed.
    pub seed: u64,
}

impl RodCrystalOptions {
    /// Options with a 5 A repeat and a 1.0 A minimum distance.
    #[must_use]
    pub fn new(rod_group: u32, composition: Vec<(String, u32)>, seed: u64) -> Self {
        Self {
            rod_group,
            composition,
            repeat: 5.0,
            min_distance: 1.0,
            max_attempts: 64,
            seed,
        }
    }
}

/// Generates one random rod crystal.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options or exhausted attempts.
pub fn generate_rod_crystal(options: &RodCrystalOptions) -> Result<LowDimCrystal, Error> {
    if !options.repeat.is_finite() || options.repeat <= 0.0 {
        return Err(Error::InvalidInput(
            "the rod repeat must be positive and finite".to_owned(),
        ));
    }
    validate_common(
        &options.composition,
        options.min_distance,
        options.max_attempts,
    )?;
    let catalogue = rod_groups()?;
    let group = catalogue
        .iter()
        .find(|group| group.number == options.rod_group)
        .ok_or_else(|| {
            Error::InvalidInput(format!(
                "rod-group number must be in 1..={}, got {}",
                catalogue.len(),
                options.rod_group
            ))
        })?;
    let mut rng = SplitMix64::new(options.seed);
    let mut last_reason = String::from("no candidate was tried");
    for attempt in 1..=options.max_attempts {
        let lattice = Lattice::from_parameters(
            20.0,
            20.0,
            options.repeat,
            90.0,
            90.0,
            if matches!(options.rod_group, 10..=14) {
                120.0
            } else {
                90.0
            },
        )?;
        let mut builder = StructureBuilder::new(format!("rod-{}", group.symbol), lattice)
            .hall("P 1")
            .metadata("generator", json!("rod"))
            .metadata("rod_group", json!(group.number))
            .metadata("seed", json!(options.seed));
        for operation in &group.operations {
            builder = builder.operation(*operation);
        }
        let mut points: Vec<[f64; 3]> = Vec::new();
        let mut labels = 0usize;
        let mut ok = true;
        for (element, count) in &options.composition {
            let mut remaining = *count as usize;
            let mut guard = 0;
            while remaining > 0 {
                guard += 1;
                if guard > 64 {
                    ok = false;
                    break;
                }
                let free = [rng.uniform(), rng.uniform(), rng.uniform()];
                let point = [0.05 + 0.15 * free[0], 0.05 + 0.15 * free[1], free[2]];
                let mut images = Vec::new();
                for operation in &group.operations {
                    push_unique(&mut images, &lattice, operation.apply(point));
                }
                if images.len() > remaining {
                    continue;
                }
                let count = images.len();
                labels += 1;
                builder = builder.site(Site::new(format!("{element}{labels}"), element, point)?);
                points.push(point);
                remaining -= count;
            }
            if !ok {
                break;
            }
        }
        let _ = points;
        if !ok {
            last_reason = "the composition cannot be expressed with the rod orbits".to_owned();
            continue;
        }
        let structure = builder.build()?;
        if minimum_distance(&structure)? < options.min_distance {
            last_reason = "a candidate violated the minimum distance".to_owned();
            continue;
        }
        if operations_preserve(&structure, &group.operations)? {
            return Ok(LowDimCrystal {
                structure,
                layer_group: group.number,
                hm_symbol: group.symbol.clone(),
                wyckoffs: Vec::new(),
                attempt,
                seed: options.seed,
            });
        }
    }
    Err(Error::InvalidInput(format!(
        "could not generate rod group {} after {} attempts: {last_reason}",
        options.rod_group, options.max_attempts
    )))
}

// Integer rotation matrices used by the rod catalogue.

fn identity() -> [[i32; 3]; 3] {
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
}

fn negate() -> [[i32; 3]; 3] {
    [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]
}

fn rot2() -> [[i32; 3]; 3] {
    [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
}

fn rot2x() -> [[i32; 3]; 3] {
    [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
}

fn rot2y() -> [[i32; 3]; 3] {
    [[-1, 0, 0], [0, 1, 0], [0, 0, -1]]
}

fn mirror_x() -> [[i32; 3]; 3] {
    [[-1, 0, 0], [0, 1, 0], [0, 0, 1]]
}

fn mirror_y() -> [[i32; 3]; 3] {
    [[1, 0, 0], [0, -1, 0], [0, 0, 1]]
}

fn mirror_z() -> [[i32; 3]; 3] {
    [[1, 0, 0], [0, 1, 0], [0, 0, -1]]
}

fn rot4() -> [[i32; 3]; 3] {
    [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
}

fn rot3() -> [[i32; 3]; 3] {
    [[0, -1, 0], [1, -1, 0], [0, 0, 1]]
}

fn rot6() -> [[i32; 3]; 3] {
    [[1, -1, 0], [1, 0, 0], [0, 0, 1]]
}
