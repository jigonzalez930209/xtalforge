//! Builders: graphene, nanotubes and the carbon family (Phase 6).

pub mod carbon;
pub mod carbon_recipe;
pub use carbon_recipe::{CarbonBuild, CarbonRecipe, build_carbon, carbon_surface};
pub mod embed;
pub mod fragments;
pub mod fullerene;
pub mod guest;
pub mod hydrogens;
pub mod net;
pub mod polymer;
pub mod ribbon;
pub mod sketch;
pub mod smiles;
pub mod spiral;
pub mod substitute;
pub mod twist;
pub mod ztc;

pub use carbon::{graphene, nanotube};
pub use embed::{
    CoordinationGeometry, EmbedOptions, apply_coordination, coordination_geometry, embed_smiles,
    embed_smiles_with_hydrogens,
};
pub use fragments::{Fragment, find_fragment, fragments};
pub use fullerene::fullerene;
pub use guest::place_guest;
pub use hydrogens::{HydrogenOptions, add_hydrogens};
pub use net::build_net;
pub use polymer::{PolymerBuild, PolymerOptions, build_polymer, pack_chains, polymer_protocol};
pub use ribbon::nanoribbon;
pub use sketch::{RDKIT_BEHAVIOUR_INVENTORY, Sketch, SketchAtom, SketchBond, from_smiles};
pub use smiles::{SmilesAtom, SmilesBond, SmilesMolecule, parse_smiles};
pub use spiral::{face_spiral, spiral_sequence};
pub use substitute::{clearance, substitute};
pub use twist::commensurate_twist;
pub use ztc::{ZtcOptions, build_schwarzite, build_ztc, schwarzite_surface_mesh, ztc_surface_mesh};
