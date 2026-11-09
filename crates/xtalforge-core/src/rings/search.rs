//! Periodic cycle walk and Franzblau primitive filter (Task 1.12).
//!
//! Ported from Crystal-Builder `xtal/core/rings.py` (MIT): the cycle walk
//! follows the reference `_cycles`, `_translated_home`, `_is_primitive` and
//! `_distances` functions, with the same path budget and depth.

use std::collections::hash_map::Entry;
use std::collections::{BTreeSet, HashMap, VecDeque};

use crate::{BondGraph, Error, Structure};

use super::{Ring, RingStep, translated_home};

const BUDGET_FLOOR: usize = 1_000_000;
const BUDGET_PER_ATOM: usize = 500;

pub(super) fn find(structure: &Structure, max_size: usize) -> Result<Vec<Ring>, Error> {
    let graph = structure.bond_graph()?;
    find_graph(&graph, max_size)
}

pub(super) fn find_graph(graph: &BondGraph, max_size: usize) -> Result<Vec<Ring>, Error> {
    let adjacency = adjacency(graph);
    let n_atoms = graph.n_atoms();
    let _ = graph;
    let budget = BUDGET_FLOOR.saturating_add(BUDGET_PER_ATOM.saturating_mul(n_atoms));
    let mut left = budget;
    let mut seen: BTreeSet<BTreeSet<RingStep>> = BTreeSet::new();
    let mut candidates: Vec<Vec<RingStep>> = Vec::new();

    for start in 0..n_atoms {
        let origin = RingStep {
            atom: start,
            image: [0, 0, 0],
        };
        let mut stack: Vec<Vec<RingStep>> = vec![vec![origin]];
        while let Some(path) = stack.pop() {
            if path.len() > max_size {
                continue;
            }
            let Some(last) = path.last().copied() else {
                continue;
            };
            let Some(neighbours) = adjacency.get(last.atom) else {
                continue;
            };
            for &(atom, step) in neighbours {
                if atom < start {
                    continue;
                }
                let next = RingStep {
                    atom,
                    image: add(last.image, step),
                };
                if next == origin {
                    if path.len() >= 3 && seen.insert(translated_home(&path)) {
                        candidates.push(path.clone());
                    }
                    continue;
                }
                if path.contains(&next) {
                    continue;
                }
                left = left.checked_sub(1).ok_or_else(|| {
                    Error::InvalidInput("bond graph is too dense for the ring search".to_owned())
                })?;
                let mut extended = path.clone();
                extended.push(next);
                stack.push(extended);
            }
        }
    }

    let depth = max_size / 2;
    let mut reach: Vec<Option<HashMap<RingStep, usize>>> = vec![None; n_atoms];
    let mut rings = Vec::new();
    for path in candidates {
        if primitive(&path, &adjacency, depth, &mut reach) {
            let mut atoms: Vec<usize> = path.iter().map(|step| step.atom).collect();
            atoms.sort_unstable();
            atoms.dedup();
            rings.push(Ring { atoms, path });
        }
    }
    rings.sort_by(|first, second| first.path.cmp(&second.path));
    Ok(rings)
}

fn adjacency(graph: &BondGraph) -> Vec<Vec<(usize, [i32; 3])>> {
    let mut adjacency = vec![Vec::new(); graph.n_atoms()];
    for bond in &graph.bonds {
        if let Some(list) = adjacency.get_mut(bond.i) {
            list.push((bond.j, bond.image));
        }
        if bond.i != bond.j {
            if let Some(list) = adjacency.get_mut(bond.j) {
                list.push((bond.i, negate(bond.image)));
            }
        } else if let Some(list) = adjacency.get_mut(bond.i) {
            list.push((bond.i, negate(bond.image)));
        }
    }
    for list in &mut adjacency {
        list.sort_unstable();
        list.dedup();
    }
    adjacency
}

fn primitive(
    path: &[RingStep],
    adjacency: &[Vec<(usize, [i32; 3])>],
    depth: usize,
    reach: &mut [Option<HashMap<RingStep, usize>>],
) -> bool {
    let n = path.len();
    for (k, step_k) in path.iter().enumerate() {
        let Some(slot) = reach.get_mut(step_k.atom) else {
            continue;
        };
        if slot.is_none() {
            *slot = Some(distances(step_k.atom, adjacency, depth));
        }
        let Some(near) = slot.as_ref() else {
            continue;
        };
        for (m, step_m) in path.iter().enumerate().skip(k + 2) {
            let around = (m - k).min(n - (m - k));
            if around < 2 {
                continue;
            }
            let key = RingStep {
                atom: step_m.atom,
                image: subtract(step_m.image, step_k.image),
            };
            if near.get(&key).is_some_and(|through| *through < around) {
                return false;
            }
        }
    }
    true
}

fn distances(
    start: usize,
    adjacency: &[Vec<(usize, [i32; 3])>],
    depth: usize,
) -> HashMap<RingStep, usize> {
    let origin = RingStep {
        atom: start,
        image: [0, 0, 0],
    };
    let mut seen = HashMap::from([(origin, 0)]);
    let mut queue = VecDeque::from([origin]);
    while let Some(node) = queue.pop_front() {
        let Some(&distance) = seen.get(&node) else {
            continue;
        };
        if distance == depth {
            continue;
        }
        let Some(neighbours) = adjacency.get(node.atom) else {
            continue;
        };
        for &(atom, step) in neighbours {
            let next = RingStep {
                atom,
                image: add(node.image, step),
            };
            if let Entry::Vacant(slot) = seen.entry(next) {
                slot.insert(distance + 1);
                queue.push_back(next);
            }
        }
    }
    seen
}

fn add(first: [i32; 3], second: [i32; 3]) -> [i32; 3] {
    [
        first[0] + second[0],
        first[1] + second[1],
        first[2] + second[2],
    ]
}

fn subtract(first: [i32; 3], second: [i32; 3]) -> [i32; 3] {
    [
        first[0] - second[0],
        first[1] - second[1],
        first[2] - second[2],
    ]
}

fn negate(image: [i32; 3]) -> [i32; 3] {
    [-image[0], -image[1], -image[2]]
}
