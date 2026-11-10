//! Cell scale fit for the MOF builder.

use xtalforge_core::Lattice;

use crate::block::Block;
use crate::error::Error;
use crate::orient;
use crate::topology::Topology;

/// Uniform scale of the topology cell that makes every node-to-node gap match
/// the span of the linker, found by golden-section search on the squared
/// relative mismatch. Edges without a connection pair do not take part.
pub(crate) fn fit_cell_scale(
    topology: &Topology,
    lattice: Lattice,
    rotations: &[[[f64; 3]; 3]],
    assignments: &[[Option<usize>; 2]],
    node: &Block,
    linker: &Block,
) -> f64 {
    let Some((local_first, local_second)) = widest_connection_pair(linker) else {
        return 1.0;
    };
    let (Some(start), Some(end)) = (
        linker.connection_points.get(local_first),
        linker.connection_points.get(local_second),
    ) else {
        return 1.0;
    };
    let span = orient::norm(orient::sub(*end, *start));
    if span < 1e-6 {
        return 1.0;
    }
    // Per edge: the unscaled node-to-node vector and the fixed offset between
    // the two connection points that the linker has to bridge.
    let mut edges: Vec<([f64; 3], [f64; 3])> = Vec::new();
    for (index, edge) in topology.edges.iter().enumerate() {
        let [first, second] = *edge;
        let [Some(first_connection), Some(second_connection)] =
            assignments.get(index).copied().unwrap_or([None, None])
        else {
            continue;
        };
        let (Some(first_fractional), Some(second_fractional)) =
            (topology.nodes.get(first), topology.nodes.get(second))
        else {
            continue;
        };
        let (wrapped, _) =
            lattice.minimum_image(orient::sub(*second_fractional, *first_fractional));
        let (Some(first_point), Some(second_point), Some(first_rotation), Some(second_rotation)) = (
            node.connection_points.get(first_connection),
            node.connection_points.get(second_connection),
            rotations.get(first),
            rotations.get(second),
        ) else {
            continue;
        };
        let offset = orient::sub(
            orient::matrix_apply(*second_rotation, *second_point),
            orient::matrix_apply(*first_rotation, *first_point),
        );
        edges.push((lattice.to_cart(wrapped), offset));
    }
    if edges.is_empty() {
        return 1.0;
    }
    let cost = |scale: f64| -> f64 {
        edges
            .iter()
            .map(|(delta, offset)| {
                let gap = orient::norm(orient::add(orient::scale(*delta, scale), *offset));
                let ratio = gap / span - 1.0;
                ratio * ratio
            })
            .sum()
    };
    let golden = 0.5 * (5.0_f64.sqrt() - 1.0);
    let (mut low, mut high) = (0.1_f64, 20.0_f64);
    for _ in 0..80 {
        let first = high - golden * (high - low);
        let second = low + golden * (high - low);
        if cost(first) < cost(second) {
            high = second;
        } else {
            low = first;
        }
    }
    0.5 * (low + high)
}

pub(crate) fn scaled_lattice(lattice: Lattice, scale: f64) -> Result<Lattice, Error> {
    let matrix = lattice.matrix();
    Ok(Lattice::new(
        matrix.map(|row| row.map(|value| value * scale)),
    )?)
}

pub(crate) fn widest_connection_pair(block: &Block) -> Option<(usize, usize)> {
    let mut best: Option<(f64, usize, usize)> = None;
    for (first, point) in block.connection_points.iter().enumerate() {
        for (second, other) in block.connection_points.iter().enumerate().skip(first + 1) {
            let distance = orient::norm(orient::sub(*other, *point));
            if best.is_none_or(|(value, _, _)| distance > value + 1e-12) {
                best = Some((distance, first, second));
            }
        }
    }
    best.map(|(_, first, second)| (first, second))
}
