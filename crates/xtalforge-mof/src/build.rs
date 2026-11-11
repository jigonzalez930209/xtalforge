//! Framework assembly (Tasks 6.12, 6.13 and the Phase 6 request surface).
//!
//! # Placement
//!
//! A node block is copied to every topology vertex without scaling. For a
//! vertex `v`, the incident edge directions `d_e` are unit vectors computed
//! with the minimum image of the periodic cell. The node rotation follows the
//! deterministic rule in [`crate::orient`]: the local reference is the mean
//! connection direction (or the first connection point when the mean is
//! zero), the target is the mean incident direction (or the first edge when
//! that mean is zero), and a greedy matching assigns connection points to
//! edges by descending cosine.
//!
//! # Scale rule
//!
//! For an edge `i-j`, let `p_i` and `p_j` be the node connection points chosen
//! for that edge. The linker connection-point pair `(a, b)` with the largest
//! separation is placed by:
//!
//! - rotation `R` mapping `(b - a)` onto the unit vector from `p_i` to `p_j`;
//! - uniform scale `s = |p_j - p_i| / |b - a|`;
//! - translation putting `a` on `p_i`, so every atom is
//!   `p_i + s R (x - a)`.
//!
//! Then `a` lands on `p_i` and `b` lands exactly on `p_j`; their midpoint is
//! the midpoint between the two node connection points, which is the topology
//! edge midpoint when the same node block terminates both ends. The scale is
//! purely geometric: a long or short edge stretches the linker, and no
//! chemical span check is attempted. When one endpoint has no free connection
//! point, the edge is skipped and reported in [`BuiltMof::warnings`].
//!
//! Atoms are wrapped into the P1 cell and near-duplicates within `1e-4` in
//! fractional coordinates are dropped. Fractional coordinates are rounded to
//! 12 decimal places, which is far below double precision for these cells and
//! keeps the serialized structure stable when it is read back. Labels are
//! renumbered as `<element><index>` in placement order.

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};
use serde_json::Value;
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder, TopologyBond};

use crate::block::{Block, BlockKind, find_block};
use crate::error::Error;
use crate::fit::{fit_cell_scale, scaled_lattice, widest_connection_pair};
use crate::orient::{self, Slot};
use crate::topology::{Topology, find_topology};

/// Smallest linker scale that is reported as an acceptable fit.
pub const LINKER_SCALE_MIN: f64 = 0.8;
/// Largest linker scale that is reported as an acceptable fit.
pub const LINKER_SCALE_MAX: f64 = 1.25;

/// Distance in Angstrom below which two atoms of a built framework are reported
/// as overlapping.
pub const CLOSE_CONTACT: f64 = 0.9;

/// A request to build one framework.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct MofRequest {
    /// Topology name, for example `pcu`.
    pub topology: String,
    /// Node block name, for example `Zn4O`.
    pub node: String,
    /// Linker block name, for example `benzene`.
    pub linker: String,
}

impl MofRequest {
    /// Creates a request from the three catalogue names.
    #[must_use]
    pub fn new(
        topology: impl Into<String>,
        node: impl Into<String>,
        linker: impl Into<String>,
    ) -> Self {
        Self {
            topology: topology.into(),
            node: node.into(),
            linker: linker.into(),
        }
    }
}

/// The result of [`build_mof`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BuiltMof {
    /// The assembled P1 structure.
    pub structure: Structure,
    /// Topology name that was requested.
    pub topology: String,
    /// Node block name that was requested.
    pub node: String,
    /// Linker block name that was requested.
    pub linker: String,
    /// Non-fatal problems found while placing blocks.
    pub warnings: Vec<String>,
}

impl BuiltMof {
    /// Serializes the result and its warnings to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Json`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self)
            .map_err(|error| Error::Json(format!("could not serialize built MOF: {error}")))
    }
}

/// Builds a framework from three catalogue names.
///
/// # Errors
///
/// Returns [`Error::UnknownTopology`] or [`Error::UnknownBlock`] for an
/// unknown name, [`Error::InvalidInput`] when the node or linker name has the
/// wrong kind or an embedded entry fails validation, and the errors of the
/// core crate for an invalid cell or site.
pub fn build_mof(request: &MofRequest) -> Result<BuiltMof, Error> {
    let topology = find_topology(&request.topology)?
        .ok_or_else(|| Error::UnknownTopology(request.topology.clone()))?;
    let node =
        find_block(&request.node)?.ok_or_else(|| Error::UnknownBlock(request.node.clone()))?;
    let linker =
        find_block(&request.linker)?.ok_or_else(|| Error::UnknownBlock(request.linker.clone()))?;
    if node.kind != BlockKind::Node {
        return Err(Error::InvalidInput(format!(
            "block `{}` is a linker, not a node",
            node.name
        )));
    }
    if linker.kind != BlockKind::Linker {
        return Err(Error::InvalidInput(format!(
            "block `{}` is a node, not a linker",
            linker.name
        )));
    }
    topology.validate()?;
    node.validate()?;
    linker.validate()?;
    check_compatible(&node, &linker)?;

    let base_lattice = topology.lattice()?;
    let mut warnings = Vec::new();
    let mut degree_mismatch: Vec<(usize, usize)> = Vec::new();
    let mut scale_range = (f64::INFINITY, 0.0_f64);
    let mut rotations = Vec::with_capacity(topology.nodes.len());
    let mut assignments = vec![[None, None]; topology.edges.len()];
    let mut atoms: Vec<(String, [f64; 3])> = Vec::new();
    let mut seen = BTreeSet::new();

    for vertex in 0..topology.nodes.len() {
        let slots = incident_slots(&topology, vertex, base_lattice);
        if slots.len() != node.connection_count() {
            degree_mismatch.push((vertex, slots.len()));
        }
        let (rotation, assignment) =
            if let Some(fit) = crate::align::kabsch_orientation(&node.connection_points, &slots) {
                fit
            } else {
                let rotation = match mean_slot_direction(&slots) {
                    Some(direction) => orient::node_oriented(&node, direction)?,
                    None => orient::identity(),
                };
                let rotated: Vec<[f64; 3]> = node
                    .connection_points
                    .iter()
                    .map(|connection| orient::matrix_apply(rotation, *connection))
                    .collect();
                (rotation, orient::assign_connections(&rotated, &slots))
            };
        rotations.push(rotation);
        for (slot, connection) in slots.iter().zip(assignment.iter()) {
            if let Some(connection) = connection
                && let Some(entry) = assignments.get_mut(slot.edge)
                && let Some(value) = entry.get_mut(slot.side)
            {
                *value = Some(*connection);
            }
        }
    }

    if let Some((vertex, degree)) = degree_mismatch.first().copied() {
        warnings.push(format!(
            "topology `{}`: {} of {} vertices have a degree different from the {} connection points of node block `{}` (for example vertex {vertex} has degree {degree}); edges without a connection point are skipped",
            topology.name,
            degree_mismatch.len(),
            topology.nodes.len(),
            node.connection_count(),
            node.name
        ));
    }

    let cell_scale = fit_cell_scale(
        &topology,
        base_lattice,
        &rotations,
        &assignments,
        &node,
        &linker,
    );
    let lattice = scaled_lattice(base_lattice, cell_scale)?;
    for (vertex, fractional) in topology.nodes.iter().enumerate() {
        let anchor = lattice.to_cart(*fractional);
        let Some(rotation) = rotations.get(vertex).copied() else {
            continue;
        };
        for (point, element) in node.points.iter().zip(node.elements.iter()) {
            let world = orient::add(anchor, orient::matrix_apply(rotation, *point));
            push_atom(element, world, lattice, &mut atoms, &mut seen);
        }
    }

    for (edge_index, edge) in topology.edges.iter().enumerate() {
        let [first, second] = *edge;
        let [first_connection, second_connection] =
            assignments.get(edge_index).copied().unwrap_or([None, None]);
        let (Some(first_connection), Some(second_connection)) =
            (first_connection, second_connection)
        else {
            let mut missing = Vec::new();
            if first_connection.is_none() {
                missing.push(format!("vertex {first}"));
            }
            if second_connection.is_none() {
                missing.push(format!("vertex {second}"));
            }
            warnings.push(format!(
                "topology `{}` edge {first}-{second} was skipped: no node connection point at {}",
                topology.name,
                missing.join(" and ")
            ));
            continue;
        };
        let Some(first_fractional) = topology.nodes.get(first) else {
            continue;
        };
        let Some(second_fractional) = topology.nodes.get(second) else {
            continue;
        };
        let delta = orient::sub(*second_fractional, *first_fractional);
        let (wrapped, _) = lattice.minimum_image(delta);
        let base_first = lattice.to_cart(*first_fractional);
        let base_second = lattice.to_cart(orient::add(*first_fractional, wrapped));
        let first_rotation = rotations
            .get(first)
            .copied()
            .unwrap_or_else(orient::identity);
        let second_rotation = rotations
            .get(second)
            .copied()
            .unwrap_or_else(orient::identity);
        let Some(first_point) = node.connection_points.get(first_connection) else {
            continue;
        };
        let Some(second_point) = node.connection_points.get(second_connection) else {
            continue;
        };
        let meet_first = orient::add(
            base_first,
            orient::matrix_apply(first_rotation, *first_point),
        );
        let meet_second = orient::add(
            base_second,
            orient::matrix_apply(second_rotation, *second_point),
        );
        let gap_vector = orient::sub(meet_second, meet_first);
        let gap = orient::norm(gap_vector);
        if gap < 1e-6 {
            warnings.push(format!(
                "topology `{}` edge {first}-{second} was skipped: the node connection points coincide",
                topology.name
            ));
            continue;
        }
        let Some(direction) = orient::normalize(gap_vector) else {
            continue;
        };
        let Some((local_first, local_second)) = widest_connection_pair(&linker) else {
            continue;
        };
        let Some(local_start) = linker.connection_points.get(local_first) else {
            continue;
        };
        let Some(local_end) = linker.connection_points.get(local_second) else {
            continue;
        };
        let span_vector = orient::sub(*local_end, *local_start);
        let span = orient::norm(span_vector);
        if span < 1e-6 {
            warnings.push(format!(
                "topology `{}` edge {first}-{second} was skipped: linker `{}` has coincident connection points",
                topology.name, linker.name
            ));
            continue;
        }
        let Some(span_direction) = orient::normalize(span_vector) else {
            continue;
        };
        let linker_rotation = orient::rotation_from_to(span_direction, direction);
        let scale = gap / span;
        scale_range = (scale_range.0.min(scale), scale_range.1.max(scale));
        for (point, element) in linker.points.iter().zip(linker.elements.iter()) {
            let offset = orient::sub(*point, *local_start);
            let world = orient::add(
                meet_first,
                orient::scale(orient::matrix_apply(linker_rotation, offset), scale),
            );
            push_atom(element, world, lattice, &mut atoms, &mut seen);
        }
    }

    if scale_range.0.is_finite()
        && (scale_range.0 < LINKER_SCALE_MIN || scale_range.1 > LINKER_SCALE_MAX)
    {
        warnings.push(format!(
            "linker `{}` is scaled by {:.2} to {:.2} to fit the fixed `{}` cell; outside {LINKER_SCALE_MIN}-{LINKER_SCALE_MAX} the linker geometry is distorted and atoms may overlap",
            linker.name, scale_range.0, scale_range.1, topology.name
        ));
    }

    let name = format!("{}-{}-{}", topology.name, node.name, linker.name);
    let mut builder = StructureBuilder::new(name, lattice)
        .metadata("topology", Value::String(topology.name.clone()))
        .metadata("node", Value::String(node.name.clone()))
        .metadata("linker", Value::String(linker.name.clone()))
        .metadata("cell_scale", Value::from(cell_scale));
    for (index, (element, fractional)) in atoms.iter().enumerate() {
        builder = builder.site(Site::new(
            format!("{element}{}", index + 1),
            element,
            *fractional,
        )?);
    }
    for edge in &topology.edges {
        let [first, second] = *edge;
        let image = edge_image(&topology, lattice, first, second);
        builder = builder.topology_bond(TopologyBond {
            i: first,
            j: second,
            image,
        });
    }
    let structure = builder.build()?;
    if linker.connection_count() > 2 {
        warnings.push(format!(
            "linker `{}` has {} connection points but spans each edge with its widest pair only; the other arms are left dangling",
            linker.name,
            linker.connection_count()
        ));
    }
    let close = xtalforge_core::neighbors(&structure, CLOSE_CONTACT)?.len();
    if close > 0 {
        warnings.push(format!(
            "{close} atom pairs are closer than {CLOSE_CONTACT} Angstrom; the blocks may be chemically incompatible (two carboxylate-terminated blocks, for example) or the net may not suit the node, so the geometry needs editing before use"
        ));
    }
    Ok(BuiltMof {
        structure,
        topology: topology.name,
        node: node.name,
        linker: linker.name,
        warnings,
    })
}

/// Decodes a [`MofRequest`] from JSON and builds it.
///
/// # Errors
///
/// Returns [`Error::Json`] for invalid request JSON and the errors of
/// [`build_mof`] for an unknown name or an invalid block.
pub fn build_request_json(text: &str) -> Result<BuiltMof, Error> {
    let request: MofRequest = serde_json::from_str(text)
        .map_err(|error| Error::Json(format!("invalid MOF request: {error}")))?;
    build_mof(&request)
}

fn incident_slots(topology: &Topology, vertex: usize, lattice: Lattice) -> Vec<Slot> {
    let mut slots = Vec::new();
    for (edge_index, edge) in topology.edges.iter().enumerate() {
        let [first, second] = *edge;
        let (side, other) = if first == vertex {
            (0, second)
        } else if second == vertex {
            (1, first)
        } else {
            continue;
        };
        let (Some(start), Some(end)) = (topology.nodes.get(vertex), topology.nodes.get(other))
        else {
            continue;
        };
        let delta = orient::sub(*end, *start);
        let (wrapped, _) = lattice.minimum_image(delta);
        let direction = orient::normalize(lattice.to_cart(wrapped)).unwrap_or([0.0; 3]);
        slots.push(Slot {
            edge: edge_index,
            side,
            direction,
        });
    }
    slots
}

fn mean_slot_direction(slots: &[Slot]) -> Option<[f64; 3]> {
    let mut sum = [0.0; 3];
    let mut first = None;
    for slot in slots {
        if let Some(direction) = orient::normalize(slot.direction) {
            sum = orient::add(sum, direction);
            first.get_or_insert(direction);
        }
    }
    orient::normalize(sum).or(first)
}

fn push_atom(
    element: &str,
    world: [f64; 3],
    lattice: Lattice,
    atoms: &mut Vec<(String, [f64; 3])>,
    seen: &mut BTreeSet<(i64, i64, i64)>,
) {
    let fractional = Lattice::wrap(lattice.to_frac(world)).map(canonical_fraction);
    let [x, y, z] = fractional;
    let key = (
        (snap(x) * 1e4).round() as i64,
        (snap(y) * 1e4).round() as i64,
        (snap(z) * 1e4).round() as i64,
    );
    if seen.insert(key) {
        atoms.push((element.to_owned(), fractional));
    }
}

fn snap(value: f64) -> f64 {
    if 1.0 - value < 1e-4 { 0.0 } else { value }
}

fn canonical_fraction(value: f64) -> f64 {
    let rounded = (value * 1e12).round() / 1e12;
    if rounded >= 1.0 { 0.0 } else { rounded }
}

fn edge_image(topology: &Topology, lattice: Lattice, first: usize, second: usize) -> [i32; 3] {
    let (Some(start), Some(end)) = (topology.nodes.get(first), topology.nodes.get(second)) else {
        return [0, 0, 0];
    };
    let delta = orient::sub(*end, *start);
    let (wrapped, _) = lattice.minimum_image(delta);
    let [dx, dy, dz] = delta;
    let [wx, wy, wz] = wrapped;
    [
        (dx - wx).round() as i32,
        (dy - wy).round() as i32,
        (dz - wz).round() as i32,
    ]
}

/// Refuses two metal-free blocks that both carry oxygen.
///
/// The catalogue's linkers end in carboxylates, and so does the organic
/// `benzene_tri` node; without a metal cluster between them the two
/// carboxylates are placed on top of each other (O-O at 0.6 Angstrom).
fn check_compatible(node: &Block, linker: &Block) -> Result<(), Error> {
    let metal = |block: &Block| {
        block.elements.iter().any(|element| {
            xtalforge_core::elements::element_by_symbol(element).is_ok_and(|info| info.is_metal)
        })
    };
    let oxygen = |block: &Block| block.elements.iter().any(|element| element == "O");
    if !metal(node) && !metal(linker) && oxygen(node) && oxygen(linker) {
        return Err(Error::InvalidInput(format!(
            "node `{}` and linker `{}` are both carboxylate-terminated and neither holds a metal, so their oxygens would overlap; use a metal node (Zn4O, Cu2, Cr3O) or an oxygen-free linker such as triphenylene",
            node.name, linker.name
        )));
    }
    Ok(())
}
