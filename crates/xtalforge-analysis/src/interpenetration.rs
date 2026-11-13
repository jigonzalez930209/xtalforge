//! Blatov interpenetration classes from a net (Task 14.4.2).
//!
//! The component count is topological: split the quotient graph into connected
//! pieces, and count the independent copies each piece holds by the index of
//! the lattice its own cycles generate in the lattice of everything parallel
//! to it (the saturation index of its voltages). A pair of copies written into
//! one cell and a single copy written into a doubled cell therefore both count
//! two.
//!
//! The class follows Blatov's classification (V. A. Blatov, L. Carlucci,
//! G. Ciani, D. M. Proserpio, *CrystEngComm* **6** (2004) 377):
//!
//! * **Class Ia** — the copies are related by translations of the array, so
//!   one net's lattice is a proper sublattice of the array's.
//! * **Class II** — the copies are related by a symmetry operation that is not
//!   a translation.
//! * **Class III** — no symmetry operation relates the copies.
//!
//! # Deliberate subset
//!
//! A quotient graph carries no symmetry operations. [`classify`] can therefore
//! establish Class Ia (from the voltages, or from meaningful node coordinates)
//! but cannot tell Class II from Class III; it reports
//! [`BlatovClass::Undetermined`] in that case rather than guessing.
//! [`classify_with_symmetry`] takes the array's operations (for a structure,
//! [`xtalforge_symmetry::detect_symmetry`] returns them) and resolves the
//! class. Both need node coordinates that are not all the same point; a net
//! parsed from a canonical key stores placeholders, so it yields
//! [`BlatovClass::Undetermined`] unless the doubled cell already shows in its
//! voltages.

#![allow(clippy::indexing_slicing)] // bounded 3x3 lattice arithmetic

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, SymOp};

use crate::nets::QuotientGraph;
use crate::topology::{adjacency, periodicity};

/// Blatov's interpenetration classes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BlatovClass {
    /// One net: not interpenetrated.
    None,
    /// Copies related by translations of the array.
    ClassIa,
    /// Copies related by a non-translation symmetry operation.
    ClassII,
    /// No symmetry operation relates the copies.
    ClassIII,
    /// The net does not carry enough information to decide between II and III.
    Undetermined,
}

impl BlatovClass {
    /// The class number (0 for not interpenetrated), or `None` when the class
    /// could not be established.
    #[must_use]
    pub fn number(&self) -> Option<u8> {
        match self {
            Self::None => Some(0),
            Self::ClassIa => Some(1),
            Self::ClassII => Some(2),
            Self::ClassIII => Some(3),
            Self::Undetermined => None,
        }
    }

    /// The label the reports print.
    #[must_use]
    pub fn label(&self) -> &'static str {
        match self {
            Self::None => "not interpenetrated",
            Self::ClassIa => "Class Ia",
            Self::ClassII => "Class II",
            Self::ClassIII => "Class III",
            Self::Undetermined => "class undetermined",
        }
    }
}

/// The component count and class of one net.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct InterpenetrationReport {
    /// Class, or [`BlatovClass::Undetermined`].
    pub class: BlatovClass,
    /// Connected components with at least one edge.
    pub components: usize,
    /// Components with a non-zero periodicity.
    pub frameworks: usize,
    /// Independent periodic copies, summed over the periodic components.
    pub fold: usize,
    /// Non-periodic components (solvent, ions, dust).
    pub other: usize,
}

impl InterpenetrationReport {
    /// The class number, or `None` when undetermined.
    #[must_use]
    pub fn class_number(&self) -> Option<u8> {
        self.class.number()
    }

    /// Whether a symmetry search could change the answer.
    #[must_use]
    pub fn needs_symmetry(&self) -> bool {
        self.class == BlatovClass::Undetermined
    }
}

/// Classifies a net from its topology and node coordinates.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty graph.
pub fn classify(graph: &QuotientGraph) -> Result<InterpenetrationReport, Error> {
    let (components, fold) = components(graph)?;
    let frameworks = components
        .iter()
        .filter(|component| component.periodicity > 0)
        .count();
    let other = components.len() - frameworks;
    let mut report = InterpenetrationReport {
        class: BlatovClass::Undetermined,
        components: components.len(),
        frameworks,
        fold,
        other,
    };
    if report.fold <= 1 && report.frameworks <= 1 {
        report.class = BlatovClass::None;
        return Ok(report);
    }
    if components
        .iter()
        .any(|component| component.multiplicity > 1)
    {
        report.class = BlatovClass::ClassIa;
        return Ok(report);
    }
    if !coordinates_are_meaningful(graph) {
        return Ok(report);
    }
    let periodic: Vec<&QuotientGraph> = components
        .iter()
        .filter(|component| component.periodicity > 0)
        .map(|component| &component.graph)
        .collect();
    if periodic.len() >= 2 && all_translates(&periodic) {
        report.class = BlatovClass::ClassIa;
    }
    Ok(report)
}

/// Resolves the class with the array's symmetry operations.
///
/// A translation that maps one periodic component onto another is Class Ia; the
/// first non-translation operation that does is Class II; when none does, the
/// class is III.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty graph, a non-positive or
/// non-finite tolerance, or a non-translation relation that cannot be checked
/// because the node coordinates are all the same point.
pub fn classify_with_symmetry(
    graph: &QuotientGraph,
    operations: &[SymOp],
    tolerance: f64,
) -> Result<InterpenetrationReport, Error> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(Error::InvalidInput(
            "interpenetration tolerance must be positive and finite".to_owned(),
        ));
    }
    let mut report = classify(graph)?;
    if report.class != BlatovClass::Undetermined {
        return Ok(report);
    }
    if !coordinates_are_meaningful(graph) {
        return Err(Error::InvalidInput(
            "the class needs node coordinates or a doubled cell in the voltages".to_owned(),
        ));
    }
    let (components, _) = components(graph)?;
    let periodic: Vec<&QuotientGraph> = components
        .iter()
        .filter(|component| component.periodicity > 0)
        .map(|component| &component.graph)
        .collect();
    let Some(first) = periodic.first() else {
        report.class = BlatovClass::None;
        return Ok(report);
    };
    for other in periodic.iter().skip(1) {
        if translation_maps(first, other, tolerance).is_some() {
            report.class = BlatovClass::ClassIa;
            return Ok(report);
        }
    }
    for operation in operations {
        if operation.rotation() == [[1, 0, 0], [0, 1, 0], [0, 0, 1]] {
            continue;
        }
        for other in periodic.iter().skip(1) {
            if operation_maps(*operation, first, other, tolerance) {
                report.class = BlatovClass::ClassII;
                return Ok(report);
            }
        }
    }
    report.class = BlatovClass::ClassIII;
    Ok(report)
}

struct Component {
    graph: QuotientGraph,
    periodicity: usize,
    multiplicity: usize,
}

fn components(graph: &QuotientGraph) -> Result<(Vec<Component>, usize), Error> {
    if graph.nodes.is_empty() {
        return Err(Error::InvalidInput(
            "an interpenetration class needs at least one node".to_owned(),
        ));
    }
    let mut found = Vec::new();
    for nodes in graph.components() {
        if nodes.len() < 2 && graph.edges.is_empty() {
            continue;
        }
        let subgraph = graph.subgraph(&nodes)?;
        if subgraph.edges.is_empty() {
            continue;
        }
        let rank = periodicity(&subgraph);
        let multiplicity = saturation_index(&cycle_voltages(&subgraph));
        found.push(Component {
            graph: subgraph,
            periodicity: rank,
            multiplicity,
        });
    }
    if found.is_empty() {
        return Err(Error::InvalidInput(
            "the net has no connected component with an edge".to_owned(),
        ));
    }
    let fold: usize = found
        .iter()
        .filter(|component| component.periodicity > 0)
        .map(|component| component.multiplicity)
        .sum();
    Ok((found, fold))
}

/// The translations the net's own cycles close onto.
#[must_use]
pub fn cycle_voltages(graph: &QuotientGraph) -> Vec<[i64; 3]> {
    let adjacency = adjacency(graph);
    let mut offset: Vec<Option<[i64; 3]>> = vec![None; graph.nodes.len()];
    let mut voltages = Vec::new();
    for start in 0..graph.nodes.len() {
        if offset.get(start).is_some_and(Option::is_some) {
            continue;
        }
        if let Some(slot) = offset.get_mut(start) {
            *slot = Some([0, 0, 0]);
        }
        let mut stack = vec![start];
        while let Some(node) = stack.pop() {
            let base = offset.get(node).and_then(|value| *value).unwrap_or([0; 3]);
            for (other, shift) in adjacency.get(node).map(Vec::as_slice).unwrap_or(&[]) {
                let image = [
                    base[0] + i64::from(shift[0]),
                    base[1] + i64::from(shift[1]),
                    base[2] + i64::from(shift[2]),
                ];
                match offset.get(*other).and_then(|value| *value) {
                    Some(existing) => {
                        if existing != image {
                            voltages.push([
                                image[0] - existing[0],
                                image[1] - existing[1],
                                image[2] - existing[2],
                            ]);
                        }
                    }
                    None => {
                        if let Some(slot) = offset.get_mut(*other) {
                            *slot = Some(image);
                        }
                        stack.push(*other);
                    }
                }
            }
        }
    }
    voltages
}

/// How many independent copies of a net its voltages generate.
///
/// The index of the lattice the cycles generate in the lattice of everything
/// parallel to it: 1 for a net described in its own cell, `n` for one written
/// into an `n`-fold multiple of it. Exact integer arithmetic.
#[must_use]
pub fn saturation_index(vectors: &[[i64; 3]]) -> usize {
    let basis = row_basis(vectors);
    let rank = basis.len();
    let minors: Vec<i128> = match rank {
        0 => return 1,
        1 => basis[0].iter().map(|value| value.abs()).collect(),
        2 => (0..3)
            .map(|skip| {
                let columns: Vec<usize> = (0..3).filter(|column| *column != skip).collect();
                let a = basis[0];
                let b = basis[1];
                (a[columns[0]] * b[columns[1]] - a[columns[1]] * b[columns[0]]).abs()
            })
            .collect(),
        _ => {
            vec![
                (basis[0][0] * (basis[1][1] * basis[2][2] - basis[1][2] * basis[2][1])
                    - basis[0][1] * (basis[1][0] * basis[2][2] - basis[1][2] * basis[2][0])
                    + basis[0][2] * (basis[1][0] * basis[2][1] - basis[1][1] * basis[2][0]))
                    .abs(),
            ]
        }
    };
    let index = minors.into_iter().fold(0i128, gcd);
    usize::try_from(index).unwrap_or(1).max(1)
}

fn gcd(first: i128, second: i128) -> i128 {
    let (mut a, mut b) = (first.abs(), second.abs());
    while b != 0 {
        let remainder = a % b;
        a = b;
        b = remainder;
    }
    a
}

/// A lattice basis in exact integers, by Euclidean row reduction.
fn row_basis(vectors: &[[i64; 3]]) -> Vec<[i128; 3]> {
    let mut rows: Vec<[i128; 3]> = vectors
        .iter()
        .map(|vector| vector.map(i128::from))
        .collect();
    let mut top = 0usize;
    for column in 0..3 {
        loop {
            let mut nonzero: Vec<usize> = (top..rows.len())
                .filter(|row| rows[*row][column] != 0)
                .collect();
            nonzero.sort_by_key(|row| rows[*row][column].abs());
            let (Some(&pivot), Some(&second)) = (nonzero.first(), nonzero.get(1)) else {
                break;
            };
            let factor = rows[second][column] / rows[pivot][column];
            let pivot_row = rows[pivot];
            for (index, value) in rows[second].iter_mut().enumerate() {
                *value -= factor * pivot_row[index];
            }
        }
        if let Some(pivot) = (top..rows.len()).find(|row| rows[*row][column] != 0) {
            rows.swap(top, pivot);
            top += 1;
        }
    }
    rows.truncate(top);
    rows
}

fn coordinates_are_meaningful(graph: &QuotientGraph) -> bool {
    let Some(first) = graph.nodes.first() else {
        return false;
    };
    graph.nodes.iter().skip(1).any(|node| {
        (0..3).any(|axis| {
            let delta = node[axis] - first[axis];
            (delta - delta.round()).abs() > 1e-6
        })
    })
}

fn all_translates(graphs: &[&QuotientGraph]) -> bool {
    let Some(first) = graphs.first() else {
        return false;
    };
    graphs
        .iter()
        .skip(1)
        .all(|other| translation_maps(first, other, 1e-4).is_some())
}

fn translation_maps(
    source: &QuotientGraph,
    target: &QuotientGraph,
    tolerance: f64,
) -> Option<[f64; 3]> {
    if source.nodes.len() != target.nodes.len() || source.nodes.is_empty() {
        return None;
    }
    let anchor = source.nodes.first()?;
    for candidate in &target.nodes {
        let shift = [
            candidate[0] - anchor[0],
            candidate[1] - anchor[1],
            candidate[2] - anchor[2],
        ];
        let mut used = vec![false; target.nodes.len()];
        let mut matched = true;
        for node in &source.nodes {
            let moved = [node[0] + shift[0], node[1] + shift[1], node[2] + shift[2]];
            let Some(index) = target
                .nodes
                .iter()
                .enumerate()
                .position(|(index, other)| !used[index] && close(moved, *other, tolerance))
            else {
                matched = false;
                break;
            };
            used[index] = true;
        }
        if matched {
            return Some(shift);
        }
    }
    None
}

fn operation_maps(
    operation: SymOp,
    source: &QuotientGraph,
    target: &QuotientGraph,
    tolerance: f64,
) -> bool {
    if source.nodes.len() != target.nodes.len() {
        return false;
    }
    let mut used = vec![false; target.nodes.len()];
    for node in &source.nodes {
        let moved = operation.apply(*node);
        let Some(index) = target
            .nodes
            .iter()
            .enumerate()
            .position(|(index, other)| !used[index] && close(moved, *other, tolerance))
        else {
            return false;
        };
        used[index] = true;
    }
    true
}

fn close(first: [f64; 3], second: [f64; 3], tolerance: f64) -> bool {
    (0..3).all(|axis| {
        let delta = first[axis] - second[axis];
        (delta - delta.round()).abs() <= tolerance
    })
}
