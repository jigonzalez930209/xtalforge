//! magCIF reader and writer (Task 19.1.3).
//!
//! magCIF is the magnetic CIF dialect. It keeps the ordinary cell, symmetry
//! and atom-site fields and adds magnetic operation and moment loops:
//!
//! - `_space_group_symop_magn_operation_xyz` with a trailing `,+1`/`,-1` or a
//!   prime marking time reversal, and `_space_group_symop_magn_centering_xyz`
//!   with the centering translations. Both the CIF (`_a_b`) and DDLm
//!   (`_a.b`) tag spellings are accepted by the shared parser.
//! - `_atom_site_moment_label` linked to `_atom_site_label` and the moment
//!   components in `_atom_site_moment_crystalaxis_*`,
//!   `_atom_site_moment_cartn_*` or plain `_atom_site_moment_*` fields.
//!
//! Magnetic moments become typed [`xtalforge_core::MagneticMoment`] site
//! data. Cartesian moments are projected onto the unit lattice vectors to
//! obtain the magCIF `crystalaxis` components. The magnetic operation loop is
//! kept as JSON in the [`MAGNETIC_OPERATIONS_KEY`] metadata entry, because
//! typed magnetic space groups arrive in Task 19.1.2; [`write_mcif`] reads it
//! back and falls back to the ordinary operations of the structure.
//!
//! The atom-site list is read as written, like the ordinary CIF reader: the
//! magnetic operations are not applied to the sites in this task. A file in
//! the magnetic group P1, which lists the whole cell, round-trips exactly.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Lattice, MagneticMoment, Structure, StructureBuilder, SymOp};

use crate::cif::{CifBlock, CifDocument, structure_from_block, write_cif};

/// Metadata key holding the magnetic operation loop read from a magCIF block.
///
/// The value is a JSON array of [`MagneticSymOp`] objects.
pub const MAGNETIC_OPERATIONS_KEY: &str = "magnetic_operations";

/// One magnetic symmetry operation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MagneticSymOp {
    /// Spatial operation in fractional coordinates.
    pub operation: SymOp,
    /// Whether the operation includes time reversal.
    pub time_reversal: bool,
}

/// Reads the first data block with an atom-site loop as a magnetic structure.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unparsable document, no atom-site
/// loop, a malformed magnetic operation or an invalid moment value.
pub fn read_mcif(text: &str) -> Result<Structure, Error> {
    let document = CifDocument::parse(text)?;
    for block in &document.blocks {
        let has_atoms = block.loop_with("_atom_site_fract_x").is_some()
            || block.loop_with("_atom_site_Cartn_x").is_some();
        if has_atoms {
            return block_to_structure(block);
        }
    }
    Err(Error::InvalidInput(
        "mCIF has no _atom_site loop with coordinates".to_owned(),
    ))
}

/// Reads a named magCIF data block.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the block is missing or invalid.
pub fn read_mcif_block(text: &str, name: &str) -> Result<Structure, Error> {
    let document = CifDocument::parse(text)?;
    let block = document
        .blocks
        .iter()
        .find(|block| block.name == name)
        .ok_or_else(|| Error::InvalidInput(format!("mCIF block not found: {name}")))?;
    block_to_structure(block)
}

/// Writes a structure as a magCIF data block.
///
/// The ordinary CIF fields are written by [`write_cif`]. A magnetic operation
/// loop follows, from the [`MAGNETIC_OPERATIONS_KEY`] metadata when present
/// and from the structure's ordinary operations otherwise. A moment loop is
/// added when at least one site carries a magnetic moment.
///
/// # Errors
///
/// Returns the errors of [`write_cif`] or an invalid magnetic metadata value.
pub fn write_mcif(structure: &Structure) -> Result<String, Error> {
    let mut text = write_cif(structure)?;
    if !text.ends_with('\n') {
        text.push('\n');
    }
    let moments: Vec<_> = structure
        .sites()
        .iter()
        .filter(|site| site.magnetic_moment().is_some())
        .collect();
    let has_operations = structure.metadata().contains_key(MAGNETIC_OPERATIONS_KEY);
    if has_operations || !moments.is_empty() {
        text.push_str("loop_\n_space_group_symop_magn_operation_xyz\n");
        for operation in structure_magnetic_operations(structure)? {
            text.push_str(&format!("'{}'\n", magnetic_triplet(&operation)));
        }
    }
    if !moments.is_empty() {
        text.push_str(
            "loop_\n_atom_site_moment_label\n\
             _atom_site_moment_crystalaxis_x\n\
             _atom_site_moment_crystalaxis_y\n\
             _atom_site_moment_crystalaxis_z\n",
        );
        for site in moments {
            let Some(moment) = site.magnetic_moment() else {
                continue;
            };
            let [x, y, z] = moment.components();
            text.push_str(&format!(
                "{} {} {} {}\n",
                site.label(),
                format_number(x),
                format_number(y),
                format_number(z)
            ));
        }
    }
    Ok(text)
}

fn block_to_structure(block: &CifBlock) -> Result<Structure, Error> {
    let mut structure = structure_from_block(block)?;
    if let Some(moments) = site_moments(block, &structure)? {
        structure = rebuild_with_moments(&structure, &moments)?;
    }
    let operations = parse_magnetic_operations(block)?;
    if !operations.is_empty() {
        let value = serde_json::to_value(&operations).map_err(|error| {
            Error::InvalidInput(format!(
                "could not serialize the magnetic operations: {error}"
            ))
        })?;
        structure = structure.with_metadata(MAGNETIC_OPERATIONS_KEY, value);
    }
    Ok(structure)
}

/// Moment components read from one row, grouped by source spelling.
#[derive(Debug, Clone, Copy)]
enum MomentFrame {
    /// Components are already crystal-axis projections.
    CrystalAxis,
    /// Components are Cartesian and must be projected onto the axes.
    Cartesian,
}

fn moment_loop(block: &CifBlock) -> Option<&crate::cif::CifLoop> {
    block
        .loops
        .iter()
        .find(|loop_| moment_tags(loop_).is_some())
}

fn moment_tags(loop_: &crate::cif::CifLoop) -> Option<(MomentFrame, [&'static str; 3])> {
    let crystal = [
        "_atom_site_moment_crystalaxis_x",
        "_atom_site_moment_crystalaxis_y",
        "_atom_site_moment_crystalaxis_z",
    ];
    if crystal.iter().any(|tag| loop_.column(tag).is_some()) {
        return Some((MomentFrame::CrystalAxis, crystal));
    }
    let cartesian = [
        "_atom_site_moment_cartn_x",
        "_atom_site_moment_cartn_y",
        "_atom_site_moment_cartn_z",
    ];
    if cartesian.iter().any(|tag| loop_.column(tag).is_some()) {
        return Some((MomentFrame::Cartesian, cartesian));
    }
    let plain = [
        "_atom_site_moment_x",
        "_atom_site_moment_y",
        "_atom_site_moment_z",
    ];
    if plain.iter().any(|tag| loop_.column(tag).is_some()) {
        return Some((MomentFrame::CrystalAxis, plain));
    }
    None
}

fn site_moments(
    block: &CifBlock,
    structure: &Structure,
) -> Result<Option<BTreeMap<String, MagneticMoment>>, Error> {
    let Some(loop_) = moment_loop(block) else {
        return Ok(None);
    };
    let Some((frame, [tag_x, tag_y, tag_z])) = moment_tags(loop_) else {
        return Ok(None);
    };
    let mut moments = BTreeMap::new();
    for row in 0..loop_.rows.len() {
        let target = match loop_.value(row, "_atom_site_moment_label") {
            Some(label) => structure
                .sites()
                .iter()
                .position(|site| site.label() == label.trim()),
            None => (row < structure.sites().len()).then_some(row),
        };
        let Some(target) = target else {
            continue;
        };
        let Some(site) = structure.sites().get(target) else {
            continue;
        };
        let components = [
            moment_component(loop_.value(row, tag_x))?,
            moment_component(loop_.value(row, tag_y))?,
            moment_component(loop_.value(row, tag_z))?,
        ];
        let components = match frame {
            MomentFrame::CrystalAxis => components,
            MomentFrame::Cartesian => crystal_axis_components(structure.lattice(), components),
        };
        moments.insert(site.label().to_owned(), MagneticMoment::vector(components)?);
    }
    if moments.is_empty() {
        Ok(None)
    } else {
        Ok(Some(moments))
    }
}

fn moment_component(value: Option<&str>) -> Result<f64, Error> {
    let Some(text) = value.map(str::trim) else {
        return Ok(0.0);
    };
    if text.is_empty() || text == "." || text == "?" {
        return Ok(0.0);
    }
    let text = text.split('(').next().unwrap_or(text);
    let value = text
        .parse::<f64>()
        .map_err(|_| Error::InvalidInput(format!("invalid magnetic moment value: {text:?}")))?;
    if !value.is_finite() {
        return Err(Error::InvalidInput(
            "magnetic moment values must be finite".to_owned(),
        ));
    }
    Ok(value)
}

/// Projects a Cartesian moment onto the unit lattice vectors.
fn crystal_axis_components(lattice: &Lattice, cartesian: [f64; 3]) -> [f64; 3] {
    let [a, b, c] = lattice.matrix();
    let [length_a, length_b, length_c] = lattice.lengths();
    [
        clean(dot(cartesian, a) / length_a),
        clean(dot(cartesian, b) / length_b),
        clean(dot(cartesian, c) / length_c),
    ]
}

/// Drops the rounding noise of the right-angle cosine in cell parameters.
fn clean(value: f64) -> f64 {
    if value.abs() < 1e-12 { 0.0 } else { value }
}

fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    let [x, y, z] = first;
    let [u, v, w] = second;
    x * u + y * v + z * w
}

fn rebuild_with_moments(
    structure: &Structure,
    moments: &BTreeMap<String, MagneticMoment>,
) -> Result<Structure, Error> {
    let mut builder =
        StructureBuilder::new(structure.name(), *structure.lattice()).hall(structure.hall());
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for bond in structure.explicit_bonds() {
        builder = builder.explicit_bond(*bond);
    }
    for bond in structure.topology_bonds() {
        builder = builder.topology_bond(*bond);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    for site in structure.sites() {
        let site = match moments.get(site.label()) {
            Some(moment) => site.clone().with_magnetic_moment(*moment),
            None => site.clone(),
        };
        builder = builder.site(site);
    }
    builder.build()
}

fn parse_magnetic_operations(block: &CifBlock) -> Result<Vec<MagneticSymOp>, Error> {
    let mut base = Vec::new();
    if let Some(loop_) = block.loops.iter().find(|loop_| {
        loop_
            .column("_space_group_symop_magn_operation_xyz")
            .is_some()
    }) {
        for row in 0..loop_.rows.len() {
            if let Some(text) = loop_.value(row, "_space_group_symop_magn_operation_xyz") {
                base.push(parse_magnetic_operation(text)?);
            }
        }
    }
    let mut centering = Vec::new();
    if let Some(loop_) = block.loops.iter().find(|loop_| {
        loop_
            .column("_space_group_symop_magn_centering_xyz")
            .is_some()
    }) {
        for row in 0..loop_.rows.len() {
            if let Some(text) = loop_.value(row, "_space_group_symop_magn_centering_xyz") {
                centering.push(parse_vector(text)?);
            }
        }
    }
    let mut operations = Vec::new();
    for operation in base {
        if centering.is_empty() {
            push_unique(&mut operations, operation);
        } else {
            for shift in &centering {
                push_unique(&mut operations, compose_translation(operation, *shift)?);
            }
        }
    }
    Ok(operations)
}

fn parse_magnetic_operation(text: &str) -> Result<MagneticSymOp, Error> {
    let mut body = text
        .trim()
        .trim_matches(|character: char| character == '\'' || character == '"')
        .trim()
        .to_owned();
    let mut time_reversal = body.contains('\'');
    body = body.replace('\'', "");
    if let Some((head, flag)) = body.rsplit_once(',') {
        match flag.trim() {
            "+1" | "1" => body = head.to_owned(),
            "-1" => {
                time_reversal = true;
                body = head.to_owned();
            }
            _ => {}
        }
    }
    let operation = SymOp::parse(body.trim())?;
    Ok(MagneticSymOp {
        operation,
        time_reversal,
    })
}

fn compose_translation(magnetic: MagneticSymOp, shift: [f64; 3]) -> Result<MagneticSymOp, Error> {
    let [x, y, z] = magnetic.operation.translation();
    let [u, v, w] = shift;
    let translation = [x + u, y + v, z + w].map(|value| value.rem_euclid(1.0));
    Ok(MagneticSymOp {
        operation: SymOp::new(magnetic.operation.rotation(), translation)?,
        time_reversal: magnetic.time_reversal,
    })
}

fn parse_vector(text: &str) -> Result<[f64; 3], Error> {
    let text = text
        .trim()
        .trim_matches(|character: char| character == '\'' || character == '"');
    let parts: Vec<&str> = text.split(',').map(str::trim).collect();
    let [x, y, z] = <[&str; 3]>::try_from(parts)
        .map_err(|_| Error::InvalidInput(format!("invalid magnetic centering vector: {text:?}")))?;
    Ok([rational(x)?, rational(y)?, rational(z)?])
}

fn rational(text: &str) -> Result<f64, Error> {
    let value = if let Some((numerator, denominator)) = text.split_once('/') {
        let numerator = numerator.trim().parse::<f64>().map_err(|_| {
            Error::InvalidInput(format!("invalid magnetic centering value: {text:?}"))
        })?;
        let denominator = denominator.trim().parse::<f64>().map_err(|_| {
            Error::InvalidInput(format!("invalid magnetic centering value: {text:?}"))
        })?;
        if denominator == 0.0 {
            return Err(Error::InvalidInput(
                "magnetic centering denominator must not be zero".to_owned(),
            ));
        }
        numerator / denominator
    } else {
        text.parse::<f64>().map_err(|_| {
            Error::InvalidInput(format!("invalid magnetic centering value: {text:?}"))
        })?
    };
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::InvalidInput(
            "magnetic centering values must be finite".to_owned(),
        ))
    }
}

fn push_unique(operations: &mut Vec<MagneticSymOp>, operation: MagneticSymOp) {
    if !operations.contains(&operation) {
        operations.push(operation);
    }
}

fn structure_magnetic_operations(structure: &Structure) -> Result<Vec<MagneticSymOp>, Error> {
    if let Some(value) = structure.metadata().get(MAGNETIC_OPERATIONS_KEY) {
        let operations: Vec<MagneticSymOp> =
            serde_json::from_value(value.clone()).map_err(|error| {
                Error::InvalidInput(format!(
                    "invalid {MAGNETIC_OPERATIONS_KEY} metadata: {error}"
                ))
            })?;
        if !operations.is_empty() {
            return Ok(operations);
        }
    }
    if structure.operations().is_empty() {
        return Ok(vec![MagneticSymOp {
            operation: SymOp::identity(),
            time_reversal: false,
        }]);
    }
    Ok(structure
        .operations()
        .iter()
        .map(|operation| MagneticSymOp {
            operation: *operation,
            time_reversal: false,
        })
        .collect())
}

fn magnetic_triplet(operation: &MagneticSymOp) -> String {
    let flag = if operation.time_reversal { "-1" } else { "+1" };
    format!("{},{}", operation.operation.triplet(), flag)
}

fn format_number(value: f64) -> String {
    let text = format!("{value:.6}");
    text.trim_end_matches('0').trim_end_matches('.').to_owned()
}
