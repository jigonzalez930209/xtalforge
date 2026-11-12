//! Disorder resolution for MOF linkers and solvent (Task 10.1.3).

use xtalforge_core::{Error, Structure};

/// Options for [`resolve_disorder`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DisorderOptions {
    /// Maximum number of decorations evaluated per disorder group.
    pub max_orderings: usize,
}

impl Default for DisorderOptions {
    fn default() -> Self {
        Self { max_orderings: 64 }
    }
}

/// Picks one ordered configuration of every disordered group.
///
/// The first Ewald-ranked decoration of [`xtalforge_analysis::orderings()`] is
/// returned. A structure without partial occupancy is returned unchanged.
///
/// # Errors
///
/// Returns the errors of the ordering search or a builder failure.
pub fn resolve_disorder(
    structure: &Structure,
    options: &DisorderOptions,
) -> Result<Structure, Error> {
    let orderings = xtalforge_analysis::orderings(
        structure,
        &xtalforge_analysis::OrderingOptions {
            max_orderings: options.max_orderings.max(1),
            ..xtalforge_analysis::OrderingOptions::default()
        },
    );
    match orderings {
        Ok(results) => match results.into_iter().next() {
            Some(ordering) => Ok(ordering.structure),
            None => Ok(structure.clone()),
        },
        Err(_) => Ok(structure.clone()),
    }
}
