//! Builders: graphene, nanotubes and the carbon family (Phase 6).

pub mod carbon;
pub mod fragments;
pub mod fullerene;
pub mod guest;
pub mod net;
pub mod smiles;
pub mod substitute;

pub use carbon::{graphene, nanotube};
pub use fragments::{Fragment, find_fragment, fragments};
pub use fullerene::fullerene;
pub use guest::place_guest;
pub use net::build_net;
pub use smiles::{SmilesAtom, SmilesBond, SmilesMolecule, parse_smiles};
pub use substitute::{clearance, substitute};
