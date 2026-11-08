//! Integration tests for view settings, colormaps and scene data (Task 3.14).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_view::{Colormap, DrawStyle, ViewSettings, build_scene};

fn nacl() -> xtalforge_core::Structure {
    StructureBuilder::new("diamond", Lattice::cubic(3.566_986).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("C2", "C", [0.25, 0.25, 0.25]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn scene_has_atoms_bonds_and_cell_edges() {
    let settings = ViewSettings::default();
    let scene = build_scene(&nacl(), &settings).unwrap();
    assert_eq!(scene.atoms.len(), 2);
    assert_eq!(scene.bonds.len(), 1);
    assert_eq!(scene.cell_edges.len(), 12);
    assert!(scene.atoms[0].radius > 0.0);
    assert!(scene.atoms[0].position[0].abs() < 1e-12);
    assert_eq!(scene.bonds[0].order, 1.0);
}

#[test]
fn styles_change_radii_and_hidden_atoms() {
    let mut settings = ViewSettings {
        style: DrawStyle::SpaceFill,
        ..Default::default()
    };
    settings.sanitize();
    let scene = build_scene(&nacl(), &settings).unwrap();
    assert!(scene.atoms[0].radius > 0.35);
    let wire = ViewSettings {
        style: DrawStyle::Wireframe,
        ..Default::default()
    };
    let scene = build_scene(&nacl(), &wire).unwrap();
    assert_eq!(scene.atoms[0].radius, 0.0);
    assert_eq!(scene.bonds[0].radius, 0.02);
}

#[test]
fn colormaps_are_clamped_and_capped() {
    assert_eq!(Colormap::Grayscale.sample(0.0), [0, 0, 0]);
    assert_eq!(Colormap::Grayscale.sample(1.0), [255, 255, 255]);
    assert_eq!(
        Colormap::Viridis.sample(-1.0),
        Colormap::Viridis.sample(0.0)
    );
    assert_eq!(Colormap::Viridis.sample(2.0), Colormap::Viridis.sample(1.0));
}

#[test]
fn printable_mesh_is_ascii_stl() {
    use xtalforge_view::{build_mesh, write_stl};
    let scene = build_scene(&nacl(), &ViewSettings::default()).unwrap();
    let triangles = build_mesh(&scene);
    assert!(triangles.len() >= 8);
    let stl = write_stl(&triangles, "test");
    assert!(stl.starts_with("solid test"));
    assert!(stl.ends_with("endsolid test\n"));
    assert_eq!(stl.matches("facet normal").count(), triangles.len());
}
