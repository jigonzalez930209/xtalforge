//! Rods, layers, merged nets, decoration and invented linkers
//! (Tasks 6.14-6.16).
//!
//! These are documented geometric subsets of the PORMAKE catalogue: a rod is a
//! metal chain with linkers, a layered framework stacks honeycomb sheets, an
//! invented linker is a linear carbon chain between two connection points, and
//! the pore partition reports the accessible volume of each framework.

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};

use crate::block::{Block, BlockKind, blocks};
use crate::error::Error;

/// A linker name with its connection-point span in Angstrom.
#[derive(Debug, Clone, PartialEq)]
pub struct LinkerSpan {
    /// Block name.
    pub name: String,
    /// Distance between the first and last connection point.
    pub span: f64,
}

/// Lists the linker blocks whose span lies in `[min, max]`.
///
/// # Errors
///
/// Returns the errors of [`blocks`] or [`Error::InvalidInput`] for a reversed
/// range.
pub fn linkers_by_span(min: f64, max: f64) -> Result<Vec<LinkerSpan>, Error> {
    if !min.is_finite() || !max.is_finite() || min < 0.0 || max < min {
        return Err(Error::InvalidInput(
            "linker span range must be finite and ordered".to_owned(),
        ));
    }
    let mut linkers = Vec::new();
    for block in blocks()? {
        if block.kind != BlockKind::Linker {
            continue;
        }
        let Some(first) = block.connection_points.first() else {
            continue;
        };
        let Some(last) = block.connection_points.last() else {
            continue;
        };
        let span = distance(*first, *last);
        if (min..=max).contains(&span) {
            linkers.push(LinkerSpan {
                name: block.name,
                span,
            });
        }
    }
    linkers.sort_by(|left, right| {
        left.span
            .partial_cmp(&right.span)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(linkers)
}

/// Invents a linear carbon linker of the requested span.
///
/// The chain is a straight line with 1.4 Angstrom spacing; `aromatic` is
/// accepted for the API and currently selects the same geometry (documented).
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-finite or too-small span.
pub fn invent_linker(span: f64, aromatic: bool) -> Result<Block, Error> {
    if !span.is_finite() || span < 1.4 {
        return Err(Error::InvalidInput(format!(
            "an invented linker needs a span of at least 1.4 A, got {span}"
        )));
    }
    let _ = aromatic;
    let count = ((span / 1.4).round() as usize).max(1) + 1;
    let points: Vec<[f64; 3]> = (0..count)
        .map(|index| [index as f64 * span / (count - 1) as f64, 0.0, 0.0])
        .collect();
    let block = Block {
        name: format!("invented-{span:.2}"),
        kind: BlockKind::Linker,
        cell: [[span + 2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]],
        connection_points: vec![[0.0, 0.0, 0.0], [span, 0.0, 0.0]],
        points: points.clone(),
        elements: points.iter().map(|_| "C".to_owned()).collect(),
    };
    block.validate()?;
    Ok(block)
}

/// Builds a rod framework: a metal chain along `c` with benzene linkers in
/// the `ab` plane.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for zero repeats or a builder failure.
pub fn build_rod_framework(repeats: u32) -> Result<Structure, Error> {
    if repeats == 0 {
        return Err(Error::InvalidInput(
            "rod repeats must be positive".to_owned(),
        ));
    }
    let lattice = Lattice::new([
        [10.0, 0.0, 0.0],
        [0.0, 10.0, 0.0],
        [0.0, 0.0, 3.0 * f64::from(repeats)],
    ])
    .map_err(Error::Core)?;
    let mut builder = StructureBuilder::new("rod", lattice);
    for index in 0..repeats {
        let z = (f64::from(index) + 0.5) / f64::from(repeats);
        builder = builder
            .site(Site::new(format!("Zn{}", index + 1), "Zn", [0.5, 0.5, z]).map_err(Error::Core)?);
        for (label, x) in [("C1", 0.2), ("C2", 0.3), ("C3", 0.3), ("C4", 0.2)] {
            builder = builder.site(
                Site::new(format!("{label}_{}", index + 1), "C", [x, 0.5, z])
                    .map_err(Error::Core)?,
            );
        }
    }
    builder.build().map_err(Error::Core)
}

/// Stacks honeycomb carbon layers separated by 3.35 Angstrom.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for zero layers or a builder failure.
pub fn build_layered_framework(layers: u32) -> Result<Structure, Error> {
    if layers == 0 {
        return Err(Error::InvalidInput(
            "a layered framework needs at least one layer".to_owned(),
        ));
    }
    let height = 3.35 * f64::from(layers - 1) + 15.0;
    let lattice = Lattice::new([
        [2.46, 0.0, 0.0],
        [1.23, 2.130_422_3, 0.0],
        [0.0, 0.0, height],
    ])
    .map_err(Error::Core)?;
    let mut builder = StructureBuilder::new("layered", lattice);
    for layer in 0..layers {
        let z = (f64::from(layer) * 3.35 + 1.0) / height;
        for (index, [x, y]) in [[0.0, 0.0], [1.0 / 3.0, 1.0 / 3.0], [2.0 / 3.0, 2.0 / 3.0]]
            .iter()
            .enumerate()
        {
            builder = builder.site(
                Site::new(format!("C{layer}_{index}"), "C", [*x, *y, z]).map_err(Error::Core)?,
            );
        }
    }
    builder.build().map_err(Error::Core)
}

/// Merges two structures, shifting the second by a fractional offset.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-finite shift or a builder failure.
pub fn merge_frameworks(
    first: &Structure,
    second: &Structure,
    shift: [f64; 3],
) -> Result<Structure, Error> {
    if !shift.iter().all(|value| value.is_finite()) {
        return Err(Error::InvalidInput("merge shift must be finite".to_owned()));
    }
    let mut builder = StructureBuilder::new(
        format!("{}-{}", first.name(), second.name()),
        *first.lattice(),
    );
    for site in first.sites() {
        builder = builder.site(site.clone());
    }
    for (index, site) in second.sites().iter().enumerate() {
        let frac = site.frac();
        builder = builder.site(
            Site::new(
                format!("{}{}", site.label(), index + 1),
                site.element(),
                [
                    (frac[0] + shift[0]).rem_euclid(1.0),
                    (frac[1] + shift[1]).rem_euclid(1.0),
                    (frac[2] + shift[2]).rem_euclid(1.0),
                ],
            )
            .map_err(Error::Core)?
            .with_occupancy(site.occupancy())
            .map_err(Error::Core)?,
        );
    }
    builder.build().map_err(Error::Core)
}

/// Adds a tetrahedral `O` decoration around one node.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an out-of-range node or a builder
/// failure.
pub fn decorate_super_tetrahedron(
    structure: &Structure,
    node_index: usize,
) -> Result<Structure, Error> {
    let node = structure
        .sites()
        .get(node_index)
        .ok_or_else(|| {
            Error::InvalidInput(format!(
                "node {node_index} is outside the structure ({} sites)",
                structure.sites().len()
            ))
        })?
        .clone();
    let mut builder = StructureBuilder::new(structure.name(), *structure.lattice())
        .hall(structure.hall().to_owned());
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for site in structure.sites() {
        builder = builder.site(site.clone());
    }
    let center = structure.lattice().to_cart(node.frac());
    let directions: [[f64; 3]; 4] = [
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ];
    for (index, direction) in directions.iter().enumerate() {
        let length = (direction[0] * direction[0]
            + direction[1] * direction[1]
            + direction[2] * direction[2])
            .sqrt();
        let point = [
            center[0] + direction[0] / length * 2.0,
            center[1] + direction[1] / length * 2.0,
            center[2] + direction[2] / length * 2.0,
        ];
        let frac = structure.lattice().to_frac(point);
        builder =
            builder.site(Site::new(format!("O_t{}", index + 1), "O", frac).map_err(Error::Core)?);
    }
    builder.build().map_err(Error::Core)
}

/// Reports the accessible volume of each framework at a probe radius.
///
/// # Errors
///
/// Returns the errors of the distance grid or [`Error::InvalidInput`] for an
/// empty framework list.
pub fn partition_pore_space(
    frameworks: &[Structure],
    spacing: f64,
    probe: f64,
) -> Result<Vec<f64>, Error> {
    if frameworks.is_empty() {
        return Err(Error::InvalidInput(
            "at least one framework is required".to_owned(),
        ));
    }
    let mut volumes = Vec::with_capacity(frameworks.len());
    for framework in frameworks {
        let radii = xtalforge_analysis::default_radii(framework).map_err(Error::Core)?;
        let grid =
            xtalforge_analysis::distance_grid(framework, &radii, spacing).map_err(Error::Core)?;
        volumes.push(grid.accessible_volume(probe).map_err(Error::Core)?);
    }
    Ok(volumes)
}

fn distance(first: [f64; 3], second: [f64; 3]) -> f64 {
    ((first[0] - second[0]).powi(2)
        + (first[1] - second[1]).powi(2)
        + (first[2] - second[2]).powi(2))
    .sqrt()
}
