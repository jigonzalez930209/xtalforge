//! Interpenetrated, layered and mixed frameworks.

use std::collections::BTreeMap;

use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};
use xtalforge_symmetry::{SymmetryOptions, detect_symmetry};

use crate::platform::rng::SplitMix64;

/// Builds an `folds`-fold interpenetrated framework.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for fewer than two folds, or the errors of
/// the merge helper.
pub fn interpenetrate(structure: &Structure, folds: u32) -> Result<Structure, Error> {
    if folds < 2 {
        return Err(Error::InvalidInput(
            "interpenetration needs at least two folds".to_owned(),
        ));
    }
    let mut result = structure.clone();
    for fold in 1..folds {
        let fraction = f64::from(fold) / f64::from(folds);
        let shift = [fraction, fraction, fraction];
        result = crate::advanced::merge_frameworks(&result, structure, shift)
            .map_err(|error| Error::InvalidInput(error.to_string()))?;
    }
    Ok(result.with_name(format!("{} x{folds}", structure.name())))
}

/// Stacks 2D layers along `c` in AA or AB fashion.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for zero layers, or the errors of the
/// supercell builder.
pub fn stack_layers(
    structure: &Structure,
    layers: u32,
    offset: [f64; 2],
) -> Result<Structure, Error> {
    if layers == 0 {
        return Err(Error::InvalidInput(
            "a stack needs at least one layer".to_owned(),
        ));
    }
    let p1 = structure.as_p1()?;
    let matrix = p1.lattice().matrix();
    let lattice = Lattice::new([
        matrix[0],
        matrix[1],
        [
            matrix[2][0] * f64::from(layers),
            matrix[2][1] * f64::from(layers),
            matrix[2][2] * f64::from(layers),
        ],
    ])?;
    let mut builder = StructureBuilder::new(structure.name(), lattice);
    for layer in 0..layers {
        let shift = if layer % 2 == 0 { [0.0, 0.0] } else { offset };
        for (index, site) in p1.sites().iter().enumerate() {
            let frac = [
                (site.frac()[0] + shift[0]).rem_euclid(1.0) / 1.0,
                (site.frac()[1] + shift[1]).rem_euclid(1.0) / 1.0,
                (site.frac()[2] + f64::from(layer)) / f64::from(layers),
            ];
            builder = builder.site(
                Site::new(
                    format!("{}{}_{}", site.label(), layer, index),
                    site.element(),
                    frac,
                )?
                .with_occupancy(site.occupancy())?,
            );
        }
    }
    builder.build()
}

/// Options for mixed-metal or mixed-linker frameworks.
#[derive(Debug, Clone, PartialEq)]
pub struct MixtureOptions {
    /// Replacements by element (`Zn` -> `Co`).
    pub map: BTreeMap<String, String>,
    /// Replace one whole symmetry orbit per element rather than a random
    /// fraction.
    pub ordered: bool,
    /// Fraction of P1 atoms replaced when `ordered` is false.
    pub fraction: f64,
    /// Deterministic seed.
    pub seed: u64,
}

impl Default for MixtureOptions {
    fn default() -> Self {
        Self {
            map: BTreeMap::new(),
            ordered: true,
            fraction: 1.0,
            seed: 0,
        }
    }
}

/// Builds a mixed-metal or mixed-linker framework.
///
/// # Errors
///
/// Returns the errors of the symmetry search, species replacement or builder.
pub fn mixed_framework(
    structure: &Structure,
    options: &MixtureOptions,
) -> Result<Structure, Error> {
    let p1 = structure.as_p1()?;
    if options.ordered {
        let dataset = detect_symmetry(&p1, SymmetryOptions::default())?;
        let mut result = p1.clone();
        let mut rng = SplitMix64::new(options.seed);
        let mut elements: Vec<&String> = options.map.keys().collect();
        elements.sort();
        for element in elements {
            let Some(target) = options.map.get(element) else {
                continue;
            };
            let orbits: Vec<usize> = dataset
                .orbits
                .iter()
                .enumerate()
                .filter(|(index, representative)| {
                    **representative == *index
                        && result
                            .sites()
                            .get(*index)
                            .is_some_and(|site| site.element() == element)
                })
                .map(|(index, _)| index)
                .collect();
            if orbits.is_empty() {
                continue;
            }
            let chosen = orbits
                .get((rng.next_u64() as usize) % orbits.len())
                .copied()
                .unwrap_or(*orbits.first().unwrap_or(&0));
            let members = dataset.equivalents(chosen);
            for atom in members {
                result = xtalforge_core::replace_site(&result, atom, target)?;
            }
        }
        Ok(result)
    } else {
        let map = options.map.clone();
        xtalforge_core::replace_species(
            &p1,
            &xtalforge_core::SpeciesReplacement {
                map,
                fraction: options.fraction,
                seed: options.seed,
            },
        )
    }
}
