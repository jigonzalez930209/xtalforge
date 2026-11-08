//! Functional-group catalogue and matcher (Task 1.13).

use crate::{BondGraph, Error, Structure};

/// A functional group the matcher knows.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum GroupKind {
    /// O-H.
    Hydroxyl,
    /// C with two oxygen neighbours.
    Carboxylate,
    /// N with one to three neighbours.
    Amine,
    /// C-F.
    Fluoride,
}

impl GroupKind {
    /// Stable name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Hydroxyl => "hydroxyl",
            Self::Carboxylate => "carboxylate",
            Self::Amine => "amine",
            Self::Fluoride => "fluoride",
        }
    }
}

/// One matched group.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GroupMatch {
    /// Group type.
    pub kind: GroupKind,
    /// Atom indices in the P1 cell.
    pub atoms: Vec<usize>,
}

/// Finds every occurrence of the catalogue on the perceived bond graph.
///
/// # Errors
///
/// Returns the errors of the P1 expansion or bond perception.
pub fn find_groups(structure: &Structure) -> Result<Vec<GroupMatch>, Error> {
    let graph = structure.bond_graph()?;
    let p1 = structure.p1()?;
    let element_of = |index: usize| {
        p1.atoms
            .get(index)
            .map(|atom| atom.element.clone())
            .unwrap_or_default()
    };
    let mut matches = Vec::new();
    for index in 0..graph.n_atoms() {
        let element = element_of(index);
        let neighbours = neighbours(&graph, index);
        match element.as_str() {
            "O" if neighbours
                .iter()
                .any(|neighbour| element_of(*neighbour) == "H") =>
            {
                if let Some(hydrogen) = neighbours
                    .iter()
                    .find(|neighbour| element_of(**neighbour) == "H")
                {
                    matches.push(GroupMatch {
                        kind: GroupKind::Hydroxyl,
                        atoms: sorted(vec![index, *hydrogen]),
                    });
                }
            }
            "C" => {
                let oxygens: Vec<usize> = neighbours
                    .iter()
                    .copied()
                    .filter(|neighbour| element_of(*neighbour) == "O")
                    .collect();
                if oxygens.len() >= 2
                    && let (Some(first), Some(second)) = (oxygens.first(), oxygens.get(1))
                {
                    matches.push(GroupMatch {
                        kind: GroupKind::Carboxylate,
                        atoms: sorted(vec![index, *first, *second]),
                    });
                }
            }
            "N" if (1..=3).contains(&neighbours.len()) => {
                matches.push(GroupMatch {
                    kind: GroupKind::Amine,
                    atoms: vec![index],
                });
            }
            "F" if neighbours.len() == 1 => {
                if let Some(carbon) = neighbours.first() {
                    matches.push(GroupMatch {
                        kind: GroupKind::Fluoride,
                        atoms: sorted(vec![index, *carbon]),
                    });
                }
            }
            _ => {}
        }
    }
    matches.sort_by(|first, second| (first.kind, &first.atoms).cmp(&(second.kind, &second.atoms)));
    matches.dedup();
    Ok(matches)
}

/// Counts per group kind.
///
/// # Errors
///
/// Returns the errors of [`find_groups`].
pub fn group_counts(structure: &Structure) -> Result<Vec<(GroupKind, usize)>, Error> {
    let mut counts: Vec<(GroupKind, usize)> = Vec::new();
    for group in find_groups(structure)? {
        if let Some(entry) = counts.iter_mut().find(|(kind, _)| *kind == group.kind) {
            entry.1 += 1;
        } else {
            counts.push((group.kind, 1));
        }
    }
    Ok(counts)
}

fn neighbours(graph: &BondGraph, index: usize) -> Vec<usize> {
    let mut out = Vec::new();
    for bond_index in graph.bonds_of(index) {
        if let Some(bond) = graph.bonds.get(*bond_index) {
            let other = if bond.i == index { bond.j } else { bond.i };
            if !out.contains(&other) {
                out.push(other);
            }
        }
    }
    out
}

fn sorted(mut atoms: Vec<usize>) -> Vec<usize> {
    atoms.sort_unstable();
    atoms.dedup();
    atoms
}
