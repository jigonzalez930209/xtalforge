//! Depth-first SMILES emission for a sketch.

use std::collections::BTreeMap;

use xtalforge_core::Error;

use super::{Sketch, SketchAtom};

/// Bond symbol accepted by the parser and writer.
pub(super) fn bond_symbol(order: f64) -> Result<&'static str, Error> {
    if (order - 1.0).abs() < 1e-9 {
        Ok("")
    } else if (order - 1.5).abs() < 1e-9 {
        Ok(":")
    } else if (order - 2.0).abs() < 1e-9 {
        Ok("=")
    } else if (order - 3.0).abs() < 1e-9 {
        Ok("#")
    } else {
        Err(Error::InvalidInput(format!(
            "SMILES cannot express bond order {order}"
        )))
    }
}

/// Depth-first spanning tree plus the edges that close rings.
pub(super) struct Tree {
    roots: Vec<usize>,
    children: Vec<Vec<(usize, f64)>>,
    openings: Vec<Vec<(usize, usize)>>,
    closings: Vec<Vec<usize>>,
}

impl Tree {
    pub(super) fn new(sketch: &Sketch) -> Result<Self, Error> {
        let count = sketch.atoms.len();
        let mut adjacency: Vec<Vec<(usize, usize)>> = vec![Vec::new(); count];
        for (index, bond) in sketch.bonds.iter().enumerate() {
            if let Some(list) = adjacency.get_mut(bond.a) {
                list.push((bond.b, index));
            }
            if let Some(list) = adjacency.get_mut(bond.b) {
                list.push((bond.a, index));
            }
        }
        let mut visited = vec![false; count];
        let mut parent = vec![None; count];
        let mut depth = vec![0_usize; count];
        let mut children: Vec<Vec<(usize, f64)>> = vec![Vec::new(); count];
        let mut seen_edge = vec![false; sketch.bonds.len()];
        let mut roots = Vec::new();
        let mut rings = Vec::new();
        for root in 0..count {
            if visited.get(root).copied().unwrap_or(true) {
                continue;
            }
            roots.push(root);
            if let Some(slot) = visited.get_mut(root) {
                *slot = true;
            }
            let mut stack = vec![(root, 0_usize)];
            while let Some((atom, next)) = stack.last_mut() {
                let current = *atom;
                let list = adjacency.get(current).map(Vec::as_slice).unwrap_or(&[]);
                if *next >= list.len() {
                    stack.pop();
                    continue;
                }
                let Some((neighbour, edge)) = list.get(*next).copied() else {
                    stack.pop();
                    continue;
                };
                *next += 1;
                if parent.get(current).copied().flatten() == Some(neighbour) {
                    continue;
                }
                let order = sketch.bonds.get(edge).map(|bond| bond.order).unwrap_or(1.0);
                if visited.get(neighbour).copied().unwrap_or(true) {
                    if seen_edge.get(edge).copied().unwrap_or(true) {
                        continue;
                    }
                    if let Some(slot) = seen_edge.get_mut(edge) {
                        *slot = true;
                    }
                    let current_depth = depth.get(current).copied().unwrap_or(0);
                    let neighbour_depth = depth.get(neighbour).copied().unwrap_or(0);
                    if current_depth < neighbour_depth {
                        rings.push((current, neighbour, edge));
                    } else {
                        rings.push((neighbour, current, edge));
                    }
                    continue;
                }
                if let Some(slot) = seen_edge.get_mut(edge) {
                    *slot = true;
                }
                if let Some(slot) = visited.get_mut(neighbour) {
                    *slot = true;
                }
                let next_depth = depth.get(current).copied().unwrap_or(0) + 1;
                if let Some(slot) = depth.get_mut(neighbour) {
                    *slot = next_depth;
                }
                if let Some(slot) = parent.get_mut(neighbour) {
                    *slot = Some(current);
                }
                if let Some(list) = children.get_mut(current) {
                    list.push((neighbour, order));
                }
                stack.push((neighbour, 0));
            }
        }
        let mut openings: Vec<Vec<(usize, usize)>> = vec![Vec::new(); count];
        let mut closings: Vec<Vec<usize>> = vec![Vec::new(); count];
        for &(ancestor, descendant, edge) in &rings {
            let mut node = descendant;
            let mut guard = 0;
            while let Some(up) = parent.get(node).copied().flatten() {
                if up == ancestor {
                    break;
                }
                node = up;
                guard += 1;
                if guard > count {
                    return Err(Error::InvalidInput(
                        "the sketch graph is inconsistent".to_owned(),
                    ));
                }
            }
            if let Some(list) = openings.get_mut(ancestor) {
                list.push((node, edge));
            }
            if let Some(list) = closings.get_mut(descendant) {
                list.push(edge);
            }
        }
        for list in &mut openings {
            list.sort_unstable();
        }
        for list in &mut closings {
            list.sort_unstable();
        }
        Ok(Self {
            roots,
            children,
            openings,
            closings,
        })
    }
}

enum Action {
    Atom(usize),
    Subtree(usize),
    Open(usize),
    Close(usize),
    Text(&'static str),
}

pub(super) fn emit(sketch: &Sketch, tree: &Tree) -> Result<String, Error> {
    let mut actions = Vec::new();
    let mut components = Vec::new();
    for (index, root) in tree.roots.iter().enumerate() {
        if index > 0 {
            components.push(Action::Text("."));
        }
        components.push(Action::Subtree(*root));
    }
    actions.extend(components.into_iter().rev());
    let mut output = String::new();
    let mut open: BTreeMap<usize, usize> = BTreeMap::new();
    while let Some(action) = actions.pop() {
        match action {
            Action::Text(text) => output.push_str(text),
            Action::Atom(index) => {
                let atom = sketch.atoms.get(index).ok_or_else(|| {
                    Error::InvalidInput(format!("sketch atom {index} is missing"))
                })?;
                output.push_str(&atom_symbol(atom));
            }
            Action::Subtree(index) => {
                let children = tree.children.get(index).map(Vec::as_slice).unwrap_or(&[]);
                for (position, (child, order)) in children.iter().enumerate().rev() {
                    if position > 0 {
                        actions.push(Action::Text(")"));
                    }
                    actions.push(Action::Subtree(*child));
                    actions.push(Action::Text(bond_symbol(*order)?));
                    for edge in openings_for(tree, index, *child).iter().rev() {
                        actions.push(Action::Open(*edge));
                    }
                    if position > 0 {
                        actions.push(Action::Text("("));
                    }
                }
                let closings = tree.closings.get(index).map(Vec::as_slice).unwrap_or(&[]);
                for edge in closings.iter().rev() {
                    actions.push(Action::Close(*edge));
                }
                actions.push(Action::Atom(index));
            }
            Action::Open(edge) => {
                let digit = free_digit(&open).ok_or_else(|| {
                    Error::InvalidInput("more than nine ring closures are open at once".to_owned())
                })?;
                open.insert(edge, digit);
                let order = sketch.bonds.get(edge).map(|bond| bond.order).unwrap_or(1.0);
                output.push_str(bond_symbol(order)?);
                output.push_str(&digit.to_string());
            }
            Action::Close(edge) => {
                let digit = open.remove(&edge).ok_or_else(|| {
                    Error::InvalidInput("a ring closure was closed twice".to_owned())
                })?;
                let order = sketch.bonds.get(edge).map(|bond| bond.order).unwrap_or(1.0);
                output.push_str(bond_symbol(order)?);
                output.push_str(&digit.to_string());
            }
        }
    }
    Ok(output)
}

fn openings_for(tree: &Tree, atom: usize, child: usize) -> Vec<usize> {
    tree.openings
        .get(atom)
        .map(|list| {
            list.iter()
                .filter(|(target, _)| *target == child)
                .map(|(_, edge)| *edge)
                .collect()
        })
        .unwrap_or_default()
}

fn free_digit(open: &BTreeMap<usize, usize>) -> Option<usize> {
    (1..=9).find(|candidate| !open.values().any(|digit| digit == candidate))
}

fn atom_symbol(atom: &SketchAtom) -> String {
    if atom.charge == 0 && organic(&atom.element) {
        return atom.element.clone();
    }
    let charge = match atom.charge {
        0 => String::new(),
        value if value > 0 => {
            if value == 1 {
                "+".to_owned()
            } else {
                format!("+{value}")
            }
        }
        value => {
            if value == -1 {
                "-".to_owned()
            } else {
                format!("{value}")
            }
        }
    };
    format!("[{}{charge}]", atom.element)
}

fn organic(element: &str) -> bool {
    matches!(
        element,
        "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
    )
}
