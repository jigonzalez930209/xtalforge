//! Plain overlay data for view clients (Task 5.9.1).
//!
//! An overlay is something a run draws over a crystal without being part of
//! it: a pore network, atomic charges or an orbital lobe. This module only
//! carries the data; it has no dependency on `xtalforge-view`, and the
//! serialized form is the interchange shape for every surface.

use serde::{Deserialize, Serialize};

use crate::nets::QuotientGraph;

/// How an overlay is drawn.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct OverlayStyle {
    /// RGB colour, one byte per channel.
    pub color: [u8; 3],
    /// Ball radius in Angstrom.
    pub radius: f64,
}

/// One drawable overlay.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Overlay {
    /// Display name.
    pub name: String,
    /// Unit cell whose rows are the lattice vectors in Angstrom.
    pub cell: [[f64; 3]; 3],
    /// Cartesian sites in Angstrom, each with its label.
    pub sites: Vec<([f64; 3], String)>,
    /// Drawing style.
    pub style: OverlayStyle,
}

/// A collection of overlays, serialized as JSON for every surface.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct OverlaySet {
    /// Overlays in insertion order.
    pub overlays: Vec<Overlay>,
}

impl OverlaySet {
    /// Empty set.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds an overlay.
    pub fn push(&mut self, overlay: Overlay) {
        self.overlays.push(overlay);
    }

    /// Number of overlays.
    #[must_use]
    pub fn len(&self) -> usize {
        self.overlays.len()
    }

    /// Whether the set has no overlays.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.overlays.is_empty()
    }
}

/// Turns a quotient graph into an overlay of its nodes.
///
/// Fractional node coordinates are converted to Cartesian Angstrom with the
/// supplied cell, whose rows are the lattice vectors. The overlay is named
/// `net` and keeps the graph's labels.
#[must_use]
pub fn overlay_net(graph: &QuotientGraph, cell: [[f64; 3]; 3], style: OverlayStyle) -> Overlay {
    let sites = graph
        .nodes
        .iter()
        .enumerate()
        .map(|(index, frac)| {
            let label = graph
                .labels
                .get(index)
                .cloned()
                .unwrap_or_else(|| format!("node {}", index + 1));
            (frac_to_cart(*frac, cell), label)
        })
        .collect();
    Overlay {
        name: "net".to_owned(),
        cell,
        sites,
        style,
    }
}

fn frac_to_cart(frac: [f64; 3], cell: [[f64; 3]; 3]) -> [f64; 3] {
    let mut cart = [0.0; 3];
    for (component, row) in frac.iter().zip(cell.iter()) {
        for (value, basis) in cart.iter_mut().zip(row.iter()) {
            *value += component * basis;
        }
    }
    cart
}
