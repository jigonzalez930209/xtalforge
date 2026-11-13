//! Undoable bond edits and perception rules (Task 3.10.2).

use xtalforge_core::{BondRules, Error, ExplicitBond, Structure};

use crate::command::Command;
use crate::edit::{check_index, metadata_with, rebuild_explicit, rebuild_metadata};

/// Adds an explicit bond between two P1 atoms.
#[derive(Debug, Clone, Copy)]
pub struct AddBond {
    bond: ExplicitBond,
}

impl AddBond {
    /// Creates the command from an explicit bond.
    #[must_use]
    pub fn new(bond: ExplicitBond) -> Self {
        Self {
            bond: ExplicitBond {
                remove: false,
                ..bond
            },
        }
    }

    /// Creates a single bond with the default order.
    #[must_use]
    pub fn between(i: usize, j: usize) -> Self {
        Self::new(ExplicitBond {
            i,
            j,
            image: [0, 0, 0],
            order: 1.0,
            remove: false,
        })
    }

    /// The bond the command adds.
    #[must_use]
    pub fn bond(&self) -> &ExplicitBond {
        &self.bond
    }
}

impl Command for AddBond {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        validate_bond(structure, &self.bond)?;
        let same = |candidate: &ExplicitBond| {
            candidate.i == self.bond.i
                && candidate.j == self.bond.j
                && candidate.image == self.bond.image
        };
        let mut bonds: Vec<ExplicitBond> = structure
            .explicit_bonds()
            .iter()
            .filter(|bond| !same(bond))
            .copied()
            .collect();
        bonds.push(self.bond);
        rebuild_explicit(structure, bonds)
    }

    fn description(&self) -> String {
        "Add bond".to_owned()
    }
}

/// Removes a perceived or explicit bond between two P1 atoms.
#[derive(Debug, Clone, Copy)]
pub struct RemoveBond {
    i: usize,
    j: usize,
    image: [i32; 3],
}

impl RemoveBond {
    /// Creates the command.
    #[must_use]
    pub fn new(i: usize, j: usize, image: [i32; 3]) -> Self {
        Self { i, j, image }
    }

    /// First atom.
    #[must_use]
    pub fn i(&self) -> usize {
        self.i
    }

    /// Second atom.
    #[must_use]
    pub fn j(&self) -> usize {
        self.j
    }
}

impl Command for RemoveBond {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        let probe = ExplicitBond {
            i: self.i,
            j: self.j,
            image: self.image,
            order: 1.0,
            remove: false,
        };
        validate_bond(structure, &probe)?;
        let same = |candidate: &ExplicitBond| {
            candidate.i == self.i && candidate.j == self.j && candidate.image == self.image
        };
        if let Some(position) = structure
            .explicit_bonds()
            .iter()
            .position(|bond| same(bond) && !bond.remove)
        {
            let bonds: Vec<ExplicitBond> = structure
                .explicit_bonds()
                .iter()
                .enumerate()
                .filter(|(index, _)| *index != position)
                .map(|(_, bond)| *bond)
                .collect();
            return rebuild_explicit(structure, bonds);
        }
        if structure
            .explicit_bonds()
            .iter()
            .any(|bond| same(bond) && bond.remove)
        {
            return Ok(structure.clone());
        }
        let mut bonds: Vec<ExplicitBond> = structure.explicit_bonds().to_vec();
        bonds.push(ExplicitBond {
            remove: true,
            ..probe
        });
        rebuild_explicit(structure, bonds)
    }

    fn description(&self) -> String {
        "Remove bond".to_owned()
    }
}

/// Changes the distance criteria stored with the structure.
#[derive(Debug, Clone, PartialEq)]
pub struct SetBondRules {
    rules: BondRules,
}

impl SetBondRules {
    /// Creates the command from the global scale and minimum distance.
    #[must_use]
    pub fn new(scale: f64, minimum_distance: f64) -> Self {
        Self {
            rules: BondRules::new(scale, minimum_distance),
        }
    }

    /// Creates the command from a complete rule set, including pair-specific
    /// ranges and the metal-metal switch.
    #[must_use]
    pub fn from_rules(rules: BondRules) -> Self {
        Self { rules }
    }

    /// Rules the command installs.
    #[must_use]
    pub fn rules(&self) -> &BondRules {
        &self.rules
    }
}

impl Command for SetBondRules {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        self.rules.validate()?;
        let metadata = metadata_with(structure, "bond_rules", self.rules.to_json());
        rebuild_metadata(structure, structure.sites().to_vec(), &metadata)
    }

    fn description(&self) -> String {
        format!(
            "Set bond rules (scale {}, minimum {})",
            self.rules.scale, self.rules.minimum_distance
        )
    }
}

fn validate_bond(structure: &Structure, bond: &ExplicitBond) -> Result<(), Error> {
    if !bond.order.is_finite() || bond.order <= 0.0 {
        return Err(Error::InvalidInput(
            "bond order must be positive and finite".to_owned(),
        ));
    }
    let n_atoms = structure.p1()?.n_atoms();
    check_index(bond.i, n_atoms, "atom")?;
    check_index(bond.j, n_atoms, "atom")?;
    if bond.i == bond.j && bond.image == [0, 0, 0] {
        return Err(Error::InvalidInput(
            "a bond needs two distinct atoms".to_owned(),
        ));
    }
    Ok(())
}
