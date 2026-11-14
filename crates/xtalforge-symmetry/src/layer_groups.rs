//! Layer groups and their Wyckoff positions (Phase 20 support).
//!
//! moyo stores the 80 layer groups, the 116 layer Hall settings and the layer
//! Wyckoff tables. This module exposes them through xtalforge's own
//! [`SymOp`] type so generators and analysis code never depend on moyo
//! directly.
//!
//! A layer operation acts on the first two coordinates; the third is the
//! aperiodic stacking direction and is kept as the identity. Generation uses a
//! slab-shaped cell with a large `c`, which is what the layer detector of Task
//! 2.8.1 also expects.

use moyo::data::{
    LayerHallNumber, LayerHallSymbol, LayerNumber, iter_layer_hall_symbol_entry,
    iter_layer_wyckoff_positions_from_hall_number, layer_hall_symbol_entry,
};
use serde::Serialize;
use xtalforge_core::{Error, SymOp};

/// One Hall setting of a layer group.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LayerGroup {
    /// Layer-group number (1 to 80).
    pub number: u32,
    /// Layer Hall number (1 to 116).
    pub hall_number: u32,
    /// Hermann-Mauguin symbol.
    pub hm_symbol: String,
    /// Setting code from the layer-group paper; empty for the default.
    pub setting: String,
    /// `p` or `c` centering.
    pub centering: char,
    /// Two-dimensional Bravais class (`mp`, `op`, `oc`, `tp` or `hp`).
    pub bravais: String,
    /// Operations acting on the layer cell.
    pub operations: Vec<SymOp>,
}

impl LayerGroup {
    /// Builds one layer Hall setting.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Symmetry`] for a Hall number outside 1 to 116 or an
    /// operation moyo cannot express in the core type.
    pub fn from_hall_number(hall_number: u32) -> Result<Self, Error> {
        let entry = layer_hall_symbol_entry(hall_number as LayerHallNumber)
            .ok_or_else(|| Error::Symmetry(format!("unknown layer Hall number {hall_number}")))?;
        let symbol =
            LayerHallSymbol::from_hall_number(hall_number as LayerHallNumber).ok_or_else(|| {
                Error::Symmetry(format!("unknown layer Hall symbol for {hall_number}"))
            })?;
        let base = symbol.traverse();
        let mut expanded = base.clone();
        for operation in &base {
            for translation in symbol.centering_translations() {
                expanded.push(
                    operation.clone()
                        * moyo::base::Operation::new(nalgebra::Matrix3::identity(), *translation),
                );
            }
        }
        let operations = crate::convert::operations_from_moyo(&expanded)?;
        let arithmetic = moyo::data::layer_arithmetic_crystal_class_entry(entry.arithmetic_number)
            .ok_or_else(|| {
                Error::Symmetry(format!(
                    "unknown layer arithmetic class for Hall {hall_number}"
                ))
            })?;
        Ok(Self {
            number: entry.number as u32,
            hall_number,
            hm_symbol: entry.hm_short.to_owned(),
            setting: entry.setting.to_owned(),
            centering: match entry.centering {
                moyo::data::LayerCentering::P => 'p',
                moyo::data::LayerCentering::C => 'c',
            },
            bravais: arithmetic.layer_bravais_class.to_string(),
            operations,
        })
    }

    /// Builds the default setting of a layer-group number.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Symmetry`] for a number outside 1 to 80.
    pub fn from_number(number: u32) -> Result<Self, Error> {
        if !(1..=80).contains(&number) {
            return Err(Error::Symmetry(format!(
                "layer-group number must be in 1..=80, got {number}"
            )));
        }
        let hall = iter_layer_hall_symbol_entry()
            .find(|entry| entry.number as u32 == number)
            .map(|entry| entry.hall_number as u32)
            .ok_or_else(|| Error::Symmetry(format!("no layer Hall setting for {number}")))?;
        Self::from_hall_number(hall)
    }
}

/// The 80 layer groups, one default setting each.
///
/// # Errors
///
/// Returns the conversion errors of [`LayerGroup::from_number`].
pub fn layer_groups() -> Result<Vec<LayerGroup>, Error> {
    let mut groups = Vec::with_capacity(80);
    for number in 1..=80u32 {
        groups.push(LayerGroup::from_number(number)?);
    }
    Ok(groups)
}

/// One layer Wyckoff position.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LayerWyckoffPosition {
    /// Multiplicity in the conventional layer cell.
    pub multiplicity: usize,
    /// Wyckoff letter.
    pub letter: char,
    /// Site-symmetry symbol.
    pub site_symmetry: String,
    /// Coordinate triplet template.
    pub coordinates: String,
}

impl LayerWyckoffPosition {
    /// Lists the positions of a layer Hall setting.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Symmetry`] for a Hall number outside 1 to 116.
    pub fn list(hall_number: u32) -> Result<Vec<Self>, Error> {
        if hall_number == 0 || hall_number > 116 {
            return Err(Error::Symmetry(format!(
                "layer Hall number must be in 1..=116, got {hall_number}"
            )));
        }
        Ok(
            iter_layer_wyckoff_positions_from_hall_number(hall_number as LayerHallNumber)
                .map(|position| Self {
                    multiplicity: position.multiplicity,
                    letter: position.letter,
                    site_symmetry: position.site_symmetry.to_owned(),
                    coordinates: position.coordinates.to_owned(),
                })
                .collect(),
        )
    }
}

/// Layer number for a Hall setting, or `None`.
#[must_use]
pub fn layer_number(hall_number: u32) -> Option<u32> {
    layer_hall_symbol_entry(hall_number as LayerHallNumber)
        .map(|entry| entry.number as LayerNumber as u32)
}
