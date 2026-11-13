//! Two-connected simplification of a contracted net (Task 14.4.2).
//!
//! A node with exactly two edges is dissolved: its two edges become one edge
//! whose lattice shift is the sum of theirs. This is the topological
//! simplification that turns a chain of linker atoms into the edge a chemist
//! sees. Self-loops are never dissolved, so a chain written as one vertex with
//! a self-loop stays one vertex.

#![allow(clippy::indexing_slicing)] // bounded local node and edge arrays

use crate::contraction::{ContractedEdge, make_edge};

pub(crate) fn strip_degree_two(
    nodes: Vec<[f64; 3]>,
    labels: Vec<String>,
    edges: Vec<ContractedEdge>,
) -> (Vec<[f64; 3]>, Vec<String>, Vec<ContractedEdge>, Vec<usize>) {
    let mut edges = edges;
    let mut removed = vec![false; nodes.len()];
    let mut stripped = Vec::new();
    loop {
        let mut contracted = false;
        for vertex in 0..nodes.len() {
            if removed[vertex] {
                continue;
            }
            let incident: Vec<usize> = edges
                .iter()
                .enumerate()
                .filter(|(_, edge)| {
                    !removed[edge.a]
                        && !removed[edge.b]
                        && edge.a != edge.b
                        && (edge.a == vertex || edge.b == vertex)
                })
                .map(|(index, _)| index)
                .collect();
            if incident.len() != 2 {
                continue;
            }
            let (Some(&first), Some(&second)) = (incident.first(), incident.get(1)) else {
                continue;
            };
            let (left, left_shift) = step(&edges[first], vertex);
            let (right, right_shift) = step(&edges[second], vertex);
            let (high, low) = if first > second {
                (first, second)
            } else {
                (second, first)
            };
            edges.swap_remove(high);
            edges.swap_remove(low);
            removed[vertex] = true;
            stripped.push(vertex);
            edges.push(make_edge(
                left,
                right,
                [
                    left_shift[0] + right_shift[0],
                    left_shift[1] + right_shift[1],
                    left_shift[2] + right_shift[2],
                ],
                Vec::new(),
            ));
            contracted = true;
            break;
        }
        if !contracted {
            break;
        }
    }
    let mut renumber = vec![usize::MAX; nodes.len()];
    let mut kept_nodes = Vec::new();
    let mut kept_labels = Vec::new();
    for (index, node) in nodes.iter().enumerate() {
        if removed[index] {
            continue;
        }
        renumber[index] = kept_nodes.len();
        kept_nodes.push(*node);
        kept_labels.push(labels.get(index).cloned().unwrap_or_default());
    }
    let mut kept_edges = Vec::new();
    for mut edge in edges {
        let (Some(&a), Some(&b)) = (renumber.get(edge.a), renumber.get(edge.b)) else {
            continue;
        };
        if a == usize::MAX || b == usize::MAX {
            continue;
        }
        edge.a = a;
        edge.b = b;
        kept_edges.push(edge);
    }
    kept_edges.sort();
    (kept_nodes, kept_labels, kept_edges, stripped)
}

fn step(edge: &ContractedEdge, from: usize) -> (usize, [i32; 3]) {
    if edge.a == from {
        (edge.b, edge.shift)
    } else {
        (edge.a, edge.shift.map(|value| -value))
    }
}
