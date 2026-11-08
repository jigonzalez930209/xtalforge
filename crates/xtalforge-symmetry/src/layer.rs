//! Layer-group detection for slabs (Task 2.8.1).

use moyo::MoyoLayerDataset;
use moyo::base::AngleTolerance;
use moyo::data::LayerSetting;
use serde::Serialize;
use xtalforge_core::{Error, Structure, SymOp};

use crate::convert::{operations_from_moyo, to_moyo};

/// Layer-group symmetry of a slab.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LayerDataset {
    /// Layer-group number (1 to 80).
    pub number: u32,
    /// Layer Hall number (1 to 116).
    pub hall_number: u32,
    /// Hermann-Mauguin symbol.
    pub hm_symbol: String,
    /// Two-dimensional Bravais type and atom count.
    pub pearson_symbol: String,
    /// Operations in the input cell.
    pub operations: Vec<SymOp>,
    /// Representative atom index for each input atom.
    pub orbits: Vec<usize>,
    /// Wyckoff letter for each input atom.
    pub wyckoffs: Vec<char>,
    /// Site-symmetry symbol for each input atom.
    pub site_symmetry_symbols: Vec<String>,
}

/// Detects the layer group of a slab.
///
/// The lattice must satisfy the layer contract: `c` is perpendicular to `a`
/// and `b` within the tolerance, because moyo treats the third axis as
/// aperiodic.
///
/// # Errors
///
/// Returns [`Error::Symmetry`] when the contract is violated or the search
/// fails.
pub fn detect_layer_symmetry(structure: &Structure, symprec: f64) -> Result<LayerDataset, Error> {
    if !symprec.is_finite() || symprec <= 0.0 {
        return Err(Error::Symmetry(format!(
            "symprec must be positive and finite, got {symprec}"
        )));
    }
    let p1 = structure.p1()?;
    let moyo_cell = to_moyo(&p1)?;
    let dataset = MoyoLayerDataset::new(
        &moyo_cell.cell,
        symprec,
        AngleTolerance::Default,
        LayerSetting::Standard,
        false,
    )
    .map_err(|error| {
        Error::Symmetry(format!(
            "layer-group search failed at symprec={symprec} A: {error}"
        ))
    })?;
    Ok(LayerDataset {
        number: dataset.number as u32,
        hall_number: dataset.hall_number as u32,
        hm_symbol: dataset.hm_symbol.clone(),
        pearson_symbol: dataset.pearson_symbol.clone(),
        operations: operations_from_moyo(&dataset.operations)?,
        orbits: dataset.orbits.clone(),
        wyckoffs: dataset.wyckoffs.clone(),
        site_symmetry_symbols: dataset.site_symmetry_symbols.clone(),
    })
}
