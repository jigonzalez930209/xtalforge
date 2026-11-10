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

/// A named view-settings preset (Task 8.3.1).
#[derive(Debug, Clone, PartialEq)]
pub struct StylePreset {
    /// Stable preset name.
    pub name: String,
    /// Settings the front end applies.
    pub settings: ViewSettings,
}

/// A registry of named draw-style presets.
///
/// [`DrawStyle`] is a closed enum because the geometry a renderer draws is a
/// contract; the settings around it (radii, colours, cell, hydrogens) are
/// open. A plugin registers a preset and a viewer offers it by name.
#[derive(Debug, Default)]
pub struct StyleRegistry {
    presets: std::collections::BTreeMap<String, ViewSettings>,
}

impl StyleRegistry {
    /// A registry with ball-and-stick, space-fill, wireframe and sticks.
    #[must_use]
    pub fn new() -> Self {
        let mut registry = Self::default();
        for (name, style) in [
            ("ball-and-stick", DrawStyle::BallStick),
            ("space-fill", DrawStyle::SpaceFill),
            ("wireframe", DrawStyle::Wireframe),
            ("sticks", DrawStyle::Sticks),
        ] {
            registry.presets.insert(
                name.to_owned(),
                ViewSettings {
                    style,
                    ..ViewSettings::default()
                },
            );
        }
        registry
    }

    /// An empty registry.
    #[must_use]
    pub fn empty() -> Self {
        Self::default()
    }

    /// Registers a preset.
    ///
    /// # Errors
    ///
    /// Returns [`xtalforge_core::Error::InvalidInput`] for an empty or
    /// duplicate name.
    pub fn register(
        &mut self,
        name: &str,
        settings: ViewSettings,
    ) -> Result<(), xtalforge_core::Error> {
        if name.trim().is_empty() {
            return Err(xtalforge_core::Error::InvalidInput(
                "a style preset needs a non-empty name".to_owned(),
            ));
        }
        if self.presets.contains_key(name) {
            return Err(xtalforge_core::Error::InvalidInput(format!(
                "style preset {name} is already registered"
            )));
        }
        self.presets.insert(name.to_owned(), settings);
        Ok(())
    }

    /// One preset's settings.
    #[must_use]
    pub fn get(&self, name: &str) -> Option<&ViewSettings> {
        self.presets.get(name)
    }

    /// Registered names, sorted.
    #[must_use]
    pub fn names(&self) -> Vec<&str> {
        self.presets.keys().map(String::as_str).collect()
    }

    /// Named presets, sorted.
    #[must_use]
    pub fn presets(&self) -> Vec<StylePreset> {
        self.presets
            .iter()
            .map(|(name, settings)| StylePreset {
                name: name.clone(),
                settings: settings.clone(),
            })
            .collect()
    }
}
