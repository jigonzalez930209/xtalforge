//! Builders: graphene, nanotubes and the carbon family (Phase 6).

pub mod carbon;
pub mod embed;
pub mod fragments;
pub mod fullerene;
pub mod guest;
pub mod hydrogens;
pub mod net;
pub mod sketch;
pub mod smiles;
pub mod substitute;

pub use carbon::{graphene, nanotube};
pub use embed::{
    CoordinationGeometry, EmbedOptions, apply_coordination, coordination_geometry, embed_smiles,
};
pub use fragments::{Fragment, find_fragment, fragments};
pub use fullerene::fullerene;
pub use guest::place_guest;
pub use hydrogens::{HydrogenOptions, add_hydrogens};
pub use net::build_net;
pub use sketch::{RDKIT_BEHAVIOUR_INVENTORY, Sketch, SketchAtom, SketchBond, from_smiles};
pub use smiles::{SmilesAtom, SmilesBond, SmilesMolecule, parse_smiles};
pub use substitute::{clearance, substitute};
