//! View settings and draw styles (Task 3.14.1/3.14.2).

use serde::{Deserialize, Serialize};

/// How atoms and bonds are drawn.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum DrawStyle {
    /// Spheres and cylinders.
    #[default]
    BallStick,
    /// Van der Waals spheres.
    SpaceFill,
    /// Lines only.
    Wireframe,
    /// Thick sticks without atom spheres.
    Sticks,
}

/// View settings a front end stores with the document.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ViewSettings {
    /// Draw style.
    pub style: DrawStyle,
    /// Atom radius multiplier.
    pub atom_scale: f64,
    /// Bond radius in Angstrom.
    pub bond_radius: f64,
    /// Whether to draw the cell edges.
    pub show_cell: bool,
    /// Whether to draw the axes.
    pub show_axes: bool,
    /// Background colour as 0-255 RGB.
    pub background: [u8; 3],
    /// Colour map for scalar colouring, when used.
    pub colormap: crate::Colormap,
    /// Whether to hide hydrogen atoms.
    pub hide_hydrogens: bool,
}

impl Default for ViewSettings {
    fn default() -> Self {
        Self {
            style: DrawStyle::BallStick,
            atom_scale: 1.0,
            bond_radius: 0.15,
            show_cell: true,
            show_axes: true,
            background: [255, 255, 255],
            colormap: crate::Colormap::Viridis,
            hide_hydrogens: false,
        }
    }
}

impl ViewSettings {
    /// Clamps the numeric fields to sensible ranges.
    pub fn sanitize(&mut self) {
        self.atom_scale = self.atom_scale.clamp(0.05, 5.0);
        self.bond_radius = self.bond_radius.clamp(0.01, 1.0);
    }
}
