//! View settings, draw styles, colormaps and scene data (Phase 3).
//!
//! The crate has no windowing dependency: a viewer consumes [`SceneData`] and
//! is free to use wgpu, VTK, three.js or anything else.

pub mod colormap;
pub mod export;
pub mod mesh;
pub mod normals;
pub mod scene;
pub mod settings;

pub use colormap::Colormap;
pub use export::{write_glb, write_gltf, write_obj};
pub use mesh::{Triangle, build_mesh, write_stl};
pub use normals::{IndexedMesh, MeshVertex, build_indexed_mesh, vertex_normals};
pub use scene::{SceneAtom, SceneBond, SceneData, build_scene};
pub use settings::{DrawStyle, StylePreset, StyleRegistry, ViewSettings};
