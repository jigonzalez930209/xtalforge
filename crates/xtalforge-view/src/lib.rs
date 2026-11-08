//! View settings, draw styles, colormaps and scene data (Phase 3).
//!
//! The crate has no windowing dependency: a viewer consumes [`SceneData`] and
//! is free to use wgpu, VTK, three.js or anything else.

pub mod colormap;
pub mod mesh;
pub mod scene;
pub mod settings;

pub use colormap::Colormap;
pub use mesh::{Triangle, build_mesh, write_stl};
pub use scene::{SceneAtom, SceneBond, SceneData, build_scene};
pub use settings::{DrawStyle, ViewSettings};
