//! Style preset registry (Task 8.3.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_view::{DrawStyle, StyleRegistry, ViewSettings};

#[test]
fn the_builtin_presets_cover_the_draw_styles() {
    let registry = StyleRegistry::new();
    assert_eq!(
        registry.names(),
        vec!["ball-and-stick", "space-fill", "sticks", "wireframe"]
    );
    assert_eq!(
        registry.get("ball-and-stick").unwrap().style,
        DrawStyle::BallStick
    );
    assert_eq!(
        registry.get("space-fill").unwrap().style,
        DrawStyle::SpaceFill
    );
    assert_eq!(
        registry.get("wireframe").unwrap().style,
        DrawStyle::Wireframe
    );
    assert_eq!(registry.get("sticks").unwrap().style, DrawStyle::Sticks);
    assert_eq!(registry.presets().len(), 4);
}

#[test]
fn a_plugin_preset_registers() {
    let mut registry = StyleRegistry::empty();
    let settings = ViewSettings {
        atom_scale: 2.5,
        bond_radius: 0.2,
        hide_hydrogens: true,
        ..ViewSettings::default()
    };
    registry.register("large-atoms", settings.clone()).unwrap();
    assert_eq!(registry.get("large-atoms"), Some(&settings));
    assert_eq!(registry.names(), vec!["large-atoms"]);
}

#[test]
fn duplicates_and_empty_names_are_refused() {
    let mut registry = StyleRegistry::new();
    assert!(
        registry
            .register("ball-and-stick", ViewSettings::default())
            .is_err()
    );
    assert!(registry.register(" ", ViewSettings::default()).is_err());
}
