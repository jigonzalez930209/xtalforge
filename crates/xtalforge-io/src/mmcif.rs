//! mmCIF/PDBx reader (Task 13.2.3).
//!
//! mmCIF is the macromolecular dialect of CIF: tags use dots
//! (`_atom_site.Cartn_x`) and atoms are usually Cartesian. The shared CIF
//! parser matches dots and underscores, and the structure extractor accepts
//! `_atom_site.Cartn_x` as well as `_atom_site.fract_x`, reads the first
//! `pdbx_PDB_model_num`, converts B factors to `U_iso`, and keeps the chain,
//! residue and sequence identifiers as site properties.

use xtalforge_core::{Error, Structure};

use crate::cif::{CifDocument, structure_from_block};

/// Reads the first mmCIF data block that holds an atom-site loop.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unparsable document, no atom-site
/// loop, a missing cell or a malformed coordinate.
pub fn read_mmcif(text: &str) -> Result<Structure, Error> {
    let document = CifDocument::parse(text)?;
    for block in &document.blocks {
        let has_atoms = block.loop_with("_atom_site_Cartn_x").is_some()
            || block.loop_with("_atom_site_fract_x").is_some();
        if has_atoms {
            return structure_from_block(block);
        }
    }
    Err(Error::InvalidInput(
        "mmCIF has no _atom_site loop with coordinates".to_owned(),
    ))
}

/// Reads a named mmCIF data block.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the block is missing or cannot become
/// a structure.
pub fn read_mmcif_block(text: &str, name: &str) -> Result<Structure, Error> {
    let document = CifDocument::parse(text)?;
    let block = document
        .blocks
        .iter()
        .find(|block| block.name == name)
        .ok_or_else(|| Error::InvalidInput(format!("mmCIF block not found: {name}")))?;
    structure_from_block(block)
}
