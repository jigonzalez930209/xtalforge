//! MOF building blocks, topologies and framework assembly (Phase 6).
//!
//! The builder follows the PORMAKE idea (MIT; see `THIRD_PARTY.md`): a net
//! supplies the topology, a node building block sits on every vertex and a
//! linker spans every edge. The data here is a curated subset, not the full
//! PORMAKE catalogue.
//!
//! This crate implements Tasks 6.11 to 6.16 as a documented subset:
//!
//! - [`topologies`] reads `data/topologies.json` (six nets: `pcu`, `dia`,
//!   `srs`, `sql`, `hcb`, `bcu`).
//! - [`blocks`] reads `data/blocks.json` (four nodes and four linkers).
//! - [`build_mof`] places the node block on every vertex and the linker block
//!   on every edge, producing a P1 [`xtalforge_core::Structure`].
//!
//! It deliberately does not implement the full PORMAKE catalogue, the RMSD
//! locator, rod or merged nets, or the supertetrahedral decoration. See the
//! crate documentation of [`build`] for the orientation and scale rules.
//!
//! Coordinates in `data/blocks.json` are illustrative geometries, not
//! quantum-chemistry-optimized structures. The builder does not relax them.

pub mod block;
pub mod build;
pub mod error;
pub mod orient;
pub mod topology;

pub use block::{Block, BlockKind, blocks, find_block};
pub use build::{BuiltMof, MofRequest, build_mof, build_request_json};
pub use error::Error;
pub use orient::node_oriented;
pub use topology::{Topology, find_topology, topologies};
