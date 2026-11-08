//! Bond graph plus geometry in, UFF atom types out (Task 4.2.2).
//!
//! Ported from Crystal-Builder `xtal/ff/uff/typer.py`, split into
//! `rules` (element rules), `coordination` (the general table rule) and
//! `perception` (aromatic rings,
//! terminal refinement and bond orders). Every assignment carries a reason
//! and a confidence, because a mistyped carbon gives a plausible-looking
//! number that is simply wrong.
//!
//! UFF4MOF's fitted rows are only reachable for a framework node, detected
//! from the chemistry around the metal; a type set by hand through the site
//! property `uff_type` always wins.

use std::collections::BTreeMap;

use crate::Error;
use crate::uff::coordination::by_coordination;
use crate::uff::geometry::{AtomView, Geometry, LocalBond};
use crate::uff::params::{ParameterSet, get, has_element};
use crate::uff::perception::{aromatic_atoms, assign_bond_orders, refine_terminal};
use crate::uff::rules::{rule_for, simple_type};

/// How much a typing decision should be trusted.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Confidence {
    /// The rule is unambiguous.
    Certain,
    /// The rule used a geometric heuristic.
    Likely,
    /// Nothing fitted; the user should look at the atom.
    Uncertain,
}

/// One atom's UFF type, and why it was chosen.
#[derive(Debug, Clone, PartialEq)]
pub struct AtomType {
    /// UFF type name.
    pub name: String,
    /// Confidence in the assignment.
    pub confidence: Confidence,
    /// Human-readable reason, for example `"four neighbours"`.
    pub reason: String,
    /// Whether the type came from a site property rather than perception.
    pub overridden: bool,
}

impl AtomType {
    pub(crate) fn new(name: &str, confidence: Confidence, reason: impl Into<String>) -> Self {
        Self {
            name: name.to_owned(),
            confidence,
            reason: reason.into(),
            overridden: false,
        }
    }

    /// Whether the assignment is certain or set by hand.
    #[must_use]
    pub fn is_sure(&self) -> bool {
        self.overridden || self.confidence == Confidence::Certain
    }
}

/// The typing of a whole P1 cell.
#[derive(Debug, Clone, PartialEq)]
pub struct Typing {
    /// One assignment per real (non-dummy) atom.
    pub types: Vec<AtomType>,
    /// One bond order per perceived bond.
    pub bond_orders: Vec<f64>,
    /// Aromatic rings as atom indices.
    pub aromatic_rings: Vec<Vec<usize>>,
    /// Decisions worth showing before the numbers are trusted.
    pub warnings: Vec<String>,
}

impl Typing {
    /// Type names in atom order.
    #[must_use]
    pub fn names(&self) -> Vec<&str> {
        self.types.iter().map(|atom| atom.name.as_str()).collect()
    }

    /// Atoms whose assignment is genuinely uncertain.
    #[must_use]
    pub fn unsure(&self) -> Vec<usize> {
        self.types
            .iter()
            .enumerate()
            .filter(|(_, atom)| atom.confidence == Confidence::Uncertain && !atom.overridden)
            .map(|(index, _)| index)
            .collect()
    }

    /// One-line summary used by the calculator and run records.
    #[must_use]
    pub fn summary(&self) -> String {
        let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
        for atom in &self.types {
            *counts.entry(atom.name.as_str()).or_default() += 1;
        }
        let listed = counts
            .iter()
            .map(|(name, count)| format!("{count}x{name}"))
            .collect::<Vec<_>>()
            .join(", ");
        let unsure = self.unsure().len();
        if unsure == 0 {
            format!("{} atoms typed: {listed}", self.types.len())
        } else {
            format!(
                "{} atoms typed: {listed} -- {unsure} uncertain",
                self.types.len()
            )
        }
    }
}

/// What a hand-written rule's UFF4MOF answer becomes under plain UFF.
const UFF_EQUIVALENT: &[(&str, &str)] = &[("O_3_f", "O_3"), ("O_2_z", "O_2"), ("S_3_f", "S_3+2")];

pub(crate) fn assign(
    view: &AtomView,
    bonds: &[LocalBond],
    geometry: &Geometry,
    rings: &[Vec<usize>],
    parameter_set: ParameterSet,
) -> Result<Typing, Error> {
    for element in &view.elements {
        if !has_element(element)? {
            return Err(Error::UnknownElement(element.clone()));
        }
    }
    let aromatic = aromatic_atoms(rings, view);
    let mut warnings = Vec::new();
    let mut types = Vec::with_capacity(view.n_atoms());
    for index in 0..view.n_atoms() {
        let mut assigned = None;
        if let Some(name) = view.overrides.get(index).and_then(Option::as_deref) {
            match get(name) {
                Ok(_) => {
                    assigned = Some(AtomType {
                        name: name.to_owned(),
                        confidence: Confidence::Certain,
                        reason: "set by hand".to_owned(),
                        overridden: true,
                    });
                }
                Err(_) => {
                    warnings.push(format!("ignored unknown UFF type {name:?} at atom {index}"));
                }
            }
        }
        let atom = if let Some(assigned) = assigned {
            assigned
        } else {
            let element = view
                .elements
                .get(index)
                .ok_or_else(|| Error::Typing(format!("atom index {index} is out of range")))?;
            let mut atom = if let Some(rule) = rule_for(element) {
                rule(
                    index,
                    view,
                    geometry,
                    aromatic.get(index).copied().unwrap_or(false),
                )?
            } else if let Some(name) = simple_type(element) {
                AtomType::new(
                    name,
                    Confidence::Certain,
                    "the only type UFF gives this element",
                )
            } else {
                by_coordination(index, view, geometry, parameter_set)?
            };
            if parameter_set == ParameterSet::Uff
                && let Some((_, replacement)) = UFF_EQUIVALENT
                    .iter()
                    .find(|(mof, _)| mof == &atom.name.as_str())
            {
                atom = AtomType::new(
                    replacement,
                    Confidence::Likely,
                    format!(
                        "{}; {} is UFF4MOF's and plain UFF was asked for",
                        atom.reason, atom.name
                    ),
                );
            }
            atom
        };
        types.push(atom);
    }
    refine_terminal(&mut types, view, geometry)?;
    let aromatic_rings: Vec<Vec<usize>> = rings
        .iter()
        .filter(|ring| {
            ring.iter()
                .all(|index| aromatic.get(*index).copied().unwrap_or(false))
        })
        .cloned()
        .collect();
    let orders = assign_bond_orders(bonds, &types, &aromatic, &aromatic_rings);
    let unsure = types
        .iter()
        .filter(|atom| atom.confidence == Confidence::Uncertain && !atom.overridden)
        .count();
    if unsure > 0 {
        warnings.push(format!(
            "{unsure} atom(s) have a type the typer is not sure of; look at them before trusting the numbers"
        ));
    }
    Ok(Typing {
        types,
        bond_orders: orders,
        aromatic_rings,
        warnings,
    })
}
