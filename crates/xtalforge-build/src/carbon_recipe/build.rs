//! Recipe-driven disordered carbon builder (Task 6.8.4).

#![allow(clippy::indexing_slicing)] // bounded graph indices

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};

use super::graph;
use super::surface::solve_surface;
use super::terminate::{Ratios, terminate};
use super::{ATOMS_PER_GCC, Rng};

/// A carbon recipe. The defaults are the reference's example ZTC:
/// diamond net at 2x2x2, 0.42 g/cm3, 38 % of the sheet, one layer.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CarbonRecipe {
    /// Net name (`dia` or `pcu`).
    pub net: String,
    /// Repetition of the net cell.
    pub repeat: [usize; 3],
    /// Carbon density in g/cm3.
    pub density: f64,
    /// Strut radius over the mean net edge.
    pub radius_ratio: f64,
    /// Fraction of the closed sheet kept.
    pub coverage: f64,
    /// Stacked layers.
    pub layers: usize,
    /// Interlayer spacing in Angstrom.
    pub interlayer: f64,
    /// Vertex distortion in Angstrom.
    pub sigma_vertex: f64,
    /// Edge bow in Angstrom.
    pub sigma_edge: f64,
    /// Stone-Wales pairs per 100 carbon rings (accepted, not yet placed).
    pub stone_wales: f64,
    /// Hydrogen per carbon.
    pub hydrogen: f64,
    /// Fluorine per carbon.
    pub fluorine: f64,
    /// Oxygen per carbon.
    pub oxygen: f64,
    /// Fraction of the oxygen that is a ring ether.
    pub ether: f64,
    /// Fraction of the oxygen that is a hydroxyl.
    pub hydroxyl: f64,
    /// Fraction of the oxygen that is a carbonyl.
    pub carbonyl: f64,
    /// Deterministic seed.
    pub seed: u64,
}

impl Default for CarbonRecipe {
    fn default() -> Self {
        Self {
            net: "dia".to_owned(),
            repeat: [2, 2, 2],
            density: 0.42,
            radius_ratio: 0.25,
            coverage: 0.38,
            layers: 1,
            interlayer: 3.35,
            sigma_vertex: 0.6,
            sigma_edge: 0.4,
            stone_wales: 1.0,
            hydrogen: 0.07,
            fluorine: 0.29,
            oxygen: 0.044,
            ether: 0.8,
            hydroxyl: 0.18,
            carbonyl: 0.02,
            seed: 0,
        }
    }
}

impl CarbonRecipe {
    fn validate(&self) -> Result<(), Error> {
        if self.net != "dia" && self.net != "pcu" {
            return Err(Error::InvalidInput(format!(
                "carbon recipes support dia and pcu, not {:?}",
                self.net
            )));
        }
        if self.repeat.contains(&0) || self.repeat.iter().any(|value| *value > 8) {
            return Err(Error::InvalidInput(
                "carbon repeat must be between 1 and 8".to_owned(),
            ));
        }
        if !self.density.is_finite() || self.density <= 0.0 {
            return Err(Error::InvalidInput(
                "carbon density must be positive and finite".to_owned(),
            ));
        }
        if !self.radius_ratio.is_finite() || self.radius_ratio <= 0.0 {
            return Err(Error::InvalidInput(
                "carbon radius ratio must be positive and finite".to_owned(),
            ));
        }
        if !self.coverage.is_finite()
            || !(0.0..=1.0).contains(&self.coverage)
            || self.coverage == 0.0
        {
            return Err(Error::InvalidInput(
                "carbon coverage must be above zero and at most one".to_owned(),
            ));
        }
        if self.layers != 1 {
            return Err(Error::InvalidInput(
                "the carbon recipe supports one layer in this build".to_owned(),
            ));
        }
        for value in [
            self.sigma_vertex,
            self.sigma_edge,
            self.stone_wales,
            self.hydrogen,
            self.fluorine,
            self.oxygen,
            self.ether,
            self.hydroxyl,
            self.carbonyl,
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(Error::InvalidInput(
                    "carbon recipe numbers must be finite and non-negative".to_owned(),
                ));
            }
        }
        Ok(())
    }

    fn ratios(&self) -> Ratios {
        Ratios {
            hydrogen: self.hydrogen,
            fluorine: self.fluorine,
            oxygen: self.oxygen,
            ether: self.ether,
            hydroxyl: self.hydroxyl,
            carbonyl: self.carbonyl,
        }
    }
}

/// A built disordered carbon and its report.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CarbonBuild {
    /// The terminated structure.
    pub structure: Structure,
    /// The recipe that produced it.
    pub recipe: CarbonRecipe,
    /// Ring census up to ten-membered rings.
    pub census: Vec<(usize, usize)>,
    /// Carbon density in g/cm3.
    pub density: f64,
    /// H/C, F/C and O/C.
    pub ratios: [f64; 3],
    /// Fraction of the graph that was an edge.
    pub edge_fraction: f64,
    /// Connected pieces.
    pub pieces: usize,
    /// Directions in which the graph percolates.
    pub periodicity: usize,
    /// Mean net edge in Angstrom.
    pub edge_length: f64,
    /// Carbons on the closed sheet before the cut.
    pub closed_carbons: usize,
    /// Edge carbons left without a termination.
    pub bare_edges: usize,
    /// Mean C-C bond length in Angstrom.
    pub mean_bond: f64,
    /// Human-readable notes.
    pub notes: Vec<String>,
}

/// Builds the carbon a recipe describes.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unsupported net, a bad repeat,
/// density, radius or coverage, a multi-layer recipe, or a sheet that
/// cannot be meshed or terminated.
pub fn build_carbon(recipe: &CarbonRecipe) -> Result<CarbonBuild, Error> {
    recipe.validate()?;
    let mut rng = Rng::new(recipe.seed);
    let solved = solve_surface(
        &recipe.net,
        recipe.repeat,
        recipe.density,
        recipe.radius_ratio,
        recipe.coverage,
        recipe.sigma_vertex,
        recipe.sigma_edge,
        &mut rng,
    )?;
    let closed = graph::dualise(&solved.mesh)?;
    let closed_carbons = closed.len();
    // The ethers take carbons' places, so the sheet keeps that many more
    // for the carbon left to be the density's.
    let target = solved.carbons + recipe.ratios().wanted(solved.carbons).ether;
    let cut = graph::cut(&closed, &solved.segments, &solved.vertices, target)?;
    let periodicity = cut.periodicity();
    let cut = super::settle::settle(&cut);
    let terminated = terminate(&cut, recipe.ratios())?;
    let mut notes = Vec::new();
    let census = match xtalforge_core::rings::ring_counts_explicit(&terminated.carbon_structure, 10)
    {
        Ok(census) => census,
        Err(error) => {
            notes.push(format!("ring census skipped: {}", error.message()));
            Vec::new()
        }
    };
    let volume = cut.lattice.volume();
    let density = terminated.carbon_count as f64 / volume / ATOMS_PER_GCC;
    let carbons = terminated.carbon_count.max(1) as f64;
    if terminated.bare_edges > 0 {
        notes.push(format!(
            "{} edge carbons left bare: raise the coverage or the termination ratios",
            terminated.bare_edges
        ));
    }
    if !terminated.short.is_empty() {
        let missing = terminated
            .short
            .iter()
            .map(|(kind, count)| format!("{count} {kind}"))
            .collect::<Vec<_>>()
            .join(", ");
        notes.push(format!(
            "too few edge carbons for {missing}: raise the coverage or lower the ratios"
        ));
    }
    if recipe.stone_wales > 0.0 {
        notes.push(
            "Stone-Wales pairs are accepted by the recipe but not placed in this build".to_owned(),
        );
    }
    Ok(CarbonBuild {
        structure: terminated.structure,
        recipe: recipe.clone(),
        census,
        density,
        ratios: [
            terminated.hydrogen as f64 / carbons,
            terminated.fluorine as f64 / carbons,
            terminated.oxygen as f64 / carbons,
        ],
        edge_fraction: terminated.edge_carbons as f64 / cut.len().max(1) as f64,
        pieces: terminated.pieces,
        periodicity,
        edge_length: solved.edge_length,
        closed_carbons,
        bare_edges: terminated.bare_edges,
        mean_bond: {
            let matrix = cut.lattice.matrix();
            let mut total = 0.0;
            let mut count = 0usize;
            for bond in &cut.bonds {
                let start = cut.centers[bond.a];
                let end = [
                    cut.centers[bond.b][0] + bond.image[0] as f64,
                    cut.centers[bond.b][1] + bond.image[1] as f64,
                    cut.centers[bond.b][2] + bond.image[2] as f64,
                ];
                let vector = [
                    (end[0] - start[0]) * matrix[0][0]
                        + (end[1] - start[1]) * matrix[1][0]
                        + (end[2] - start[2]) * matrix[2][0],
                    (end[0] - start[0]) * matrix[0][1]
                        + (end[1] - start[1]) * matrix[1][1]
                        + (end[2] - start[2]) * matrix[2][1],
                    (end[0] - start[0]) * matrix[0][2]
                        + (end[1] - start[1]) * matrix[1][2]
                        + (end[2] - start[2]) * matrix[2][2],
                ];
                total +=
                    (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]).sqrt();
                count += 1;
            }
            if count == 0 {
                0.0
            } else {
                total / count as f64
            }
        },
        notes,
    })
}
