//! Sketch data model and a deterministic SMILES writer (Tasks 6.1.5, 6.1.6).
//!
//! A sketch is a molecule as drawn: atoms at 2D points, bonds with orders,
//! formal charges. Nothing is perceived until the sketch is embedded or written
//! as SMILES. The editor itself is a GUI and out of scope.
//!
//! The writer is deterministic, not canonical: it walks a depth-first spanning
//! tree in atom order, marks the remaining edges as ring closures with digits
//! `1`-`9`, and writes charged or non-organic atoms in brackets. Aromaticity is
//! carried as bond order `1.5`, because [`SketchAtom`] has no aromatic flag;
//! writing an aromatic ring therefore produces explicit `:` bonds instead of
//! lower-case atoms. More than nine simultaneously open ring closures is a
//! typed error, since the parser accepts only single digits.
//!
//! [`RDKIT_BEHAVIOUR_INVENTORY`] records what the Crystal-Builder RDKit layer
//! asks for and what this crate actually provides.

mod writer;

use std::collections::BTreeSet;

use xtalforge_core::Error;

use crate::embed::{EmbedOptions, embed_smiles};
use crate::smiles::parse_smiles;
use writer::{Tree, bond_symbol, emit};

/// One atom of a sketch.
#[derive(Debug, Clone, PartialEq)]
pub struct SketchAtom {
    /// Element symbol, for example `C`, `N` or `Fe`.
    pub element: String,
    /// Point on the drawing plane, Angstrom.
    pub position: [f64; 2],
    /// Formal charge.
    pub charge: i32,
}

/// One bond of a sketch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SketchBond {
    /// First atom index.
    pub a: usize,
    /// Second atom index.
    pub b: usize,
    /// Bond order (`1`, `1.5`, `2` or `3`).
    pub order: f64,
}

/// A molecule as drawn.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Sketch {
    /// Atoms in drawing order.
    pub atoms: Vec<SketchAtom>,
    /// Bonds in drawing order.
    pub bonds: Vec<SketchBond>,
}

/// Reference behaviours of the Crystal-Builder RDKit layer and their status.
///
/// The first element is the behaviour, the second what this crate does. The
/// list is data, not documentation: it is meant to be consumed by report code
/// and tests.
pub const RDKIT_BEHAVIOUR_INVENTORY: &[(&str, &str)] = &[
    (
        "aromaticity perception",
        "not implemented: the parser keeps lower-case atoms and written bond orders; sketches store aromatic bonds as order 1.5",
    ),
    (
        "sanitization",
        "not implemented: valences and aromaticity are not checked, malformed graphs are accepted until embedding",
    ),
    (
        "conformer generation (ETKDG)",
        "reduced deterministic distance geometry in embed.rs; no stereochemistry, no ETKDG torsion preferences",
    ),
    (
        "stereochemistry",
        "not implemented: '@', '/' and '\\\\' are parsed and then ignored; no stereochemical claim is made",
    ),
    (
        "hydrogen handling",
        "explicit hydrogens only in the graph; hydrogens.rs adds missing hydrogens geometrically",
    ),
    (
        "fragment library and substructure matching",
        "fragments.rs stores named fragments; no RDKit substructure search or attachment-point matching",
    ),
    (
        "canonical SMILES",
        "sketch::to_smiles writes a deterministic depth-first SMILES, not a canonical one",
    ),
    (
        "MMFF/UFF cleanup",
        "not implemented here; xtalforge-ff exposes the UFF calculator separately",
    ),
    (
        "fragments and disconnected components",
        "disconnected sketches are written with '.' separators and read back as one graph",
    ),
];

impl Sketch {
    /// Creates an empty sketch.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Appends an atom and returns its index.
    pub fn add_atom(&mut self, atom: SketchAtom) -> usize {
        self.atoms.push(atom);
        self.atoms.len() - 1
    }

    /// Appends an atom from its components and returns its index.
    pub fn add_atom_element(
        &mut self,
        element: impl Into<String>,
        position: [f64; 2],
        charge: i32,
    ) -> usize {
        self.add_atom(SketchAtom {
            element: element.into(),
            position,
            charge,
        })
    }

    /// Appends a bond and returns its index.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an atom index outside the sketch, a
    /// self bond, a duplicate bond, a non-positive or non-finite order, or an
    /// order the writer cannot express.
    pub fn add_bond(&mut self, bond: SketchBond) -> Result<usize, Error> {
        let count = self.atoms.len();
        if bond.a >= count || bond.b >= count {
            return Err(Error::InvalidInput(format!(
                "bond ({}, {}) is outside a sketch with {count} atoms",
                bond.a, bond.b
            )));
        }
        if bond.a == bond.b {
            return Err(Error::InvalidInput(
                "a sketch bond needs two different atoms".to_owned(),
            ));
        }
        if !bond.order.is_finite() || bond.order <= 0.0 {
            return Err(Error::InvalidInput(
                "bond order must be positive and finite".to_owned(),
            ));
        }
        bond_symbol(bond.order)?;
        let pair = ordered(bond.a, bond.b);
        if self
            .bonds
            .iter()
            .any(|existing| ordered(existing.a, existing.b) == pair)
        {
            return Err(Error::InvalidInput(format!(
                "a bond between atoms {} and {} already exists",
                bond.a, bond.b
            )));
        }
        self.bonds.push(bond);
        Ok(self.bonds.len() - 1)
    }

    /// Appends a bond from its components and returns its index.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`Sketch::add_bond`].
    pub fn connect(&mut self, a: usize, b: usize, order: f64) -> Result<usize, Error> {
        self.add_bond(SketchBond { a, b, order })
    }

    /// Removes an atom and every bond that touched it, renumbering the rest.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when `index` is outside the sketch.
    pub fn remove_atom(&mut self, index: usize) -> Result<(), Error> {
        if index >= self.atoms.len() {
            return Err(Error::InvalidInput(format!(
                "atom {index} is outside a sketch with {} atoms",
                self.atoms.len()
            )));
        }
        self.atoms.remove(index);
        let mut kept = Vec::with_capacity(self.bonds.len());
        for bond in &self.bonds {
            if bond.a == index || bond.b == index {
                continue;
            }
            kept.push(SketchBond {
                a: renumber(bond.a, index),
                b: renumber(bond.b, index),
                order: bond.order,
            });
        }
        self.bonds = kept;
        Ok(())
    }

    /// Writes the sketch as a deterministic depth-first SMILES string.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty sketch, an empty element,
    /// an invalid bond or more than nine simultaneously open ring closures.
    pub fn to_smiles(&self) -> Result<String, Error> {
        validate(self)?;
        let tree = Tree::new(self)?;
        emit(self, &tree)
    }
}

/// Reads a SMILES string into a sketch with 2D positions from the embedder.
///
/// Aromatic bonds are stored as order `1.5`. The 2D point is the embedded
/// Cartesian `x` and `y`; the `z` coordinate is discarded.
///
/// # Errors
///
/// Returns the errors of [`parse_smiles`] and [`embed_smiles`].
pub fn from_smiles(smiles: &str) -> Result<Sketch, Error> {
    let graph = parse_smiles(smiles)?;
    let structure = embed_smiles(&graph, &EmbedOptions::default())?;
    let points = structure.cart_sites();
    let mut sketch = Sketch::new();
    for (atom, point) in graph.atoms.iter().zip(points.iter()) {
        sketch.add_atom(SketchAtom {
            element: atom.element.clone(),
            position: [point[0], point[1]],
            charge: atom.charge,
        });
    }
    for bond in &graph.bonds {
        let aromatic = (bond.order - 1.0).abs() < 1e-9
            && graph.atoms.get(bond.a).is_some_and(|atom| atom.aromatic)
            && graph.atoms.get(bond.b).is_some_and(|atom| atom.aromatic);
        sketch.connect(bond.a, bond.b, if aromatic { 1.5 } else { bond.order })?;
    }
    Ok(sketch)
}

/// Checks the public fields of a sketch, which callers may mutate freely.
fn validate(sketch: &Sketch) -> Result<(), Error> {
    if sketch.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "a sketch needs at least one atom".to_owned(),
        ));
    }
    for (index, atom) in sketch.atoms.iter().enumerate() {
        if atom.element.trim().is_empty() {
            return Err(Error::InvalidInput(format!(
                "sketch atom {index} has an empty element"
            )));
        }
        if !atom.position.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(format!(
                "sketch atom {index} has a non-finite position"
            )));
        }
    }
    let count = sketch.atoms.len();
    let mut pairs = BTreeSet::new();
    for bond in &sketch.bonds {
        if bond.a >= count || bond.b >= count {
            return Err(Error::InvalidInput(format!(
                "bond ({}, {}) is outside a sketch with {count} atoms",
                bond.a, bond.b
            )));
        }
        if bond.a == bond.b {
            return Err(Error::InvalidInput(
                "a sketch bond needs two different atoms".to_owned(),
            ));
        }
        if !bond.order.is_finite() || bond.order <= 0.0 {
            return Err(Error::InvalidInput(
                "bond order must be positive and finite".to_owned(),
            ));
        }
        bond_symbol(bond.order)?;
        if !pairs.insert(ordered(bond.a, bond.b)) {
            return Err(Error::InvalidInput(format!(
                "duplicate bond between atoms {} and {}",
                bond.a, bond.b
            )));
        }
    }
    Ok(())
}

fn ordered(a: usize, b: usize) -> (usize, usize) {
    if a <= b { (a, b) } else { (b, a) }
}

fn renumber(value: usize, removed: usize) -> usize {
    if value > removed { value - 1 } else { value }
}
