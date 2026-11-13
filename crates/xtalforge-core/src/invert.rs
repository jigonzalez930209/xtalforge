//! Inverting a structure to its other hand (Task 14.3.3 support).
//!
//! Coordinates are negated about the origin. A centrosymmetric crystal already
//! contains the inversion, so its other hand is the same crystal and it is
//! returned unchanged. A setting whose operation list is invariant under
//! negating the coordinates keeps its Hall symbol and operations; every other
//! setting is expanded to P1, because a chiral crystal cannot obey a group that
//! does not describe it.

use crate::{Error, Lattice, Structure, StructureBuilder, SymOp};

/// Whether the structure's operations contain an inversion.
///
/// A structure with no operations describes the identity only and is not
/// centrosymmetric.
#[must_use]
pub fn is_centrosymmetric(structure: &Structure) -> bool {
    structure.operations().iter().any(|operation| {
        let rotation = operation.rotation();
        rotation.iter().enumerate().all(|(row, values)| {
            values.iter().enumerate().all(|(column, value)| {
                let expected = if row == column { -1 } else { 0 };
                *value == expected
            })
        })
    })
}

/// Returns the structure in the other hand, or an error for no atoms.
///
/// A centrosymmetric structure already contains the inversion, so its other
/// hand is the same crystal and it is returned unchanged with every bond.
/// P1 and settings whose operation list is invariant under negating the
/// coordinates keep their operations. Every other setting is expanded to P1:
/// negating the coordinates of a chiral crystal cannot satisfy a group that
/// does not contain an inversion.
///
/// Explicit and topology bonds keep their atom indices and have their lattice
/// images negated. A record that names an asymmetric site cannot name an image
/// once a chiral structure is expanded, so it is dropped and perception
/// provides the bonds of the inverted cell.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the structure has no sites, or the
/// errors of the P1 expansion and the rebuild.
pub fn invert(structure: &Structure) -> Result<Structure, Error> {
    if structure.sites().is_empty() {
        return Err(Error::InvalidInput(
            "there are no atoms to invert".to_owned(),
        ));
    }
    if is_centrosymmetric(structure) {
        return Ok(structure.clone());
    }
    if structure.is_p1() || operations_are_inversion_invariant(structure.operations()) {
        return negate_keeping_group(structure);
    }
    let expanded = structure.as_p1()?;
    let inverted = negate_keeping_group(&expanded)?;
    let mut bonds: Vec<crate::ExplicitBond> = Vec::new();
    for bond in structure.explicit_bonds() {
        if bond.i < expanded.sites().len() && bond.j < expanded.sites().len() {
            bonds.push(crate::ExplicitBond {
                image: bond.image.map(|value| -value),
                ..*bond
            });
        }
    }
    let mut topology: Vec<crate::TopologyBond> = Vec::new();
    for bond in structure.topology_bonds() {
        if bond.i < expanded.sites().len() && bond.j < expanded.sites().len() {
            topology.push(crate::TopologyBond {
                image: bond.image.map(|value| -value),
                ..*bond
            });
        }
    }
    rebuild_bonds(&inverted, &bonds, &topology)
}

fn operations_are_inversion_invariant(operations: &[SymOp]) -> bool {
    if operations.is_empty() {
        return false;
    }
    operations.iter().all(|operation| {
        let rotation = operation.rotation().map(|row| row.map(|value| -value));
        let translation = operation.translation().map(|value| -value);
        operations.iter().any(|candidate| {
            candidate.rotation() == rotation
                && candidate.translation().iter().zip(translation.iter()).all(
                    |(actual, expected)| {
                        let difference = (actual - expected).rem_euclid(1.0);
                        difference < 1e-9 || (1.0 - difference) < 1e-9
                    },
                )
        })
    })
}

fn negate_keeping_group(structure: &Structure) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new(structure.name(), *structure.lattice())
        .hall(structure.hall().to_owned());
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for site in structure.sites() {
        let frac = negate(site.frac());
        builder = builder.site(site.clone().with_frac(frac)?);
    }
    for bond in structure.explicit_bonds() {
        builder = builder.explicit_bond(crate::ExplicitBond {
            image: bond.image.map(|value| -value),
            ..*bond
        });
    }
    for bond in structure.topology_bonds() {
        builder = builder.topology_bond(crate::TopologyBond {
            image: bond.image.map(|value| -value),
            ..*bond
        });
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

fn rebuild_bonds(
    source: &Structure,
    explicit: &[crate::ExplicitBond],
    topology: &[crate::TopologyBond],
) -> Result<Structure, Error> {
    let mut builder =
        StructureBuilder::new(source.name(), *source.lattice()).hall(source.hall().to_owned());
    for operation in source.operations() {
        builder = builder.operation(*operation);
    }
    for site in source.sites() {
        builder = builder.site(site.clone());
    }
    for bond in explicit {
        builder = builder.explicit_bond(*bond);
    }
    for bond in topology {
        builder = builder.topology_bond(*bond);
    }
    for (key, value) in source.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

fn negate(frac: [f64; 3]) -> [f64; 3] {
    Lattice::wrap([-frac[0], -frac[1], -frac[2]])
}
