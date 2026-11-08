//! A periodic net drawn as a structure (Task 6.10.1).

#![allow(clippy::indexing_slicing)] // fixed cell parameter slots

use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder, TopologyBond};
use xtalforge_io::cgd::CgdEntry;

/// Builds a P1 structure whose sites are the net vertices and whose topology
/// bonds are the edges.
///
/// The cell is the CGD `CELL` when present; otherwise the vertices are placed
/// in a unit cubic cell padded around the fractional range.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty net, an unknown element or a
/// bad cell.
pub fn build_net(entry: &CgdEntry, element: &str) -> Result<Structure, Error> {
    if entry.nodes.is_empty() {
        return Err(Error::InvalidInput(
            "a net needs at least one node".to_owned(),
        ));
    }
    let lattice = match &entry.cell {
        Some(parameters) if parameters.len() >= 6 => Lattice::from_parameters(
            parameters[0],
            parameters[1],
            parameters[2],
            parameters[3],
            parameters[4],
            parameters[5],
        )?,
        _ => Lattice::cubic(1.0)?,
    };
    let mut builder = StructureBuilder::new(entry.name.clone(), lattice);
    for (index, node) in entry.nodes.iter().enumerate() {
        let [x, y, z] = [
            node.frac.first().copied().unwrap_or(0.0),
            node.frac.get(1).copied().unwrap_or(0.0),
            node.frac.get(2).copied().unwrap_or(0.0),
        ];
        builder = builder.site(Site::new(
            format!("{element}{}", index + 1),
            element,
            [x, y, z],
        )?);
    }
    for (first, second) in &entry.edges {
        let a = nearest_node(&entry.nodes, first);
        let b = nearest_node(&entry.nodes, second);
        if let (Some(a), Some(b)) = (a, b)
            && a != b
        {
            builder = builder.topology_bond(TopologyBond {
                i: a,
                j: b,
                image: [0, 0, 0],
            });
        }
    }
    builder.build()
}

fn nearest_node(nodes: &[xtalforge_io::cgd::CgdNode], point: &[f64]) -> Option<usize> {
    let mut best: Option<(f64, usize)> = None;
    for (index, node) in nodes.iter().enumerate() {
        let distance = node
            .frac
            .iter()
            .zip(point.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>();
        if best.is_none_or(|(value, _)| distance < value) {
            best = Some((distance, index));
        }
    }
    best.map(|(_, index)| index)
}
