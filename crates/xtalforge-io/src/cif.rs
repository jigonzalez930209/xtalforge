//! Minimal CIF reader and writer (Task 1.15).
//!
//! Supports data blocks, tags, loops, quoted and semicolon strings, comments,
//! standard-uncertainty numbers, atom sites, symmetry operations and
//! anisotropic displacement parameters. It is a pure Rust reader for the
//! corpus; the full crystal-builder priority rules are documented in the
//! module comments and the API page.

use std::collections::BTreeMap;

use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder, SymOp};

/// One parsed CIF data block.
#[derive(Debug, Clone, PartialEq)]
pub struct CifBlock {
    /// Block name without the `data_` prefix.
    pub name: String,
    /// Single-valued items.
    pub items: BTreeMap<String, String>,
    /// Loop tables.
    pub loops: Vec<CifLoop>,
}

/// One CIF loop.
#[derive(Debug, Clone, PartialEq)]
pub struct CifLoop {
    /// Column tags.
    pub tags: Vec<String>,
    /// Rows of values.
    pub rows: Vec<Vec<String>>,
}

impl CifBlock {
    /// Single item value.
    #[must_use]
    pub fn item(&self, tag: &str) -> Option<&str> {
        self.items.get(tag).map(String::as_str)
    }

    /// First loop containing `tag`.
    #[must_use]
    pub fn loop_with(&self, tag: &str) -> Option<&CifLoop> {
        self.loops
            .iter()
            .find(|loop_| loop_.tags.iter().any(|candidate| candidate == tag))
    }
}

impl CifLoop {
    /// Column index of a tag.
    #[must_use]
    pub fn column(&self, tag: &str) -> Option<usize> {
        self.tags.iter().position(|candidate| candidate == tag)
    }

    /// Value at `row` and `tag`.
    #[must_use]
    pub fn value(&self, row: usize, tag: &str) -> Option<&str> {
        self.column(tag)
            .and_then(|column| self.rows.get(row)?.get(column))
            .map(String::as_str)
    }
}

/// A parsed CIF document.
#[derive(Debug, Clone, PartialEq)]
pub struct CifDocument {
    /// Data blocks in file order.
    pub blocks: Vec<CifBlock>,
}

impl CifDocument {
    /// Parses a CIF text.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an unbalanced loop or an empty
    /// document.
    pub fn parse(text: &str) -> Result<Self, Error> {
        let tokens = tokenize(text);
        let mut blocks = Vec::new();
        let mut index = 0;
        while index < tokens.len() {
            let token = tokens.get(index).cloned().unwrap_or_default();
            if let Some(name) = token.strip_prefix("data_") {
                blocks.push(CifBlock {
                    name: name.to_owned(),
                    items: BTreeMap::new(),
                    loops: Vec::new(),
                });
                index += 1;
                continue;
            }
            if token == "loop_" {
                let mut tags = Vec::new();
                index += 1;
                while let Some(next) = tokens.get(index) {
                    if next.starts_with('_') {
                        tags.push(next.clone());
                        index += 1;
                    } else {
                        break;
                    }
                }
                let mut values = Vec::new();
                while let Some(next) = tokens.get(index) {
                    if next.starts_with('_')
                        || next == "loop_"
                        || next.starts_with("data_")
                        || next.starts_with("save_")
                    {
                        break;
                    }
                    values.push(next.clone());
                    index += 1;
                }
                if tags.is_empty() {
                    return Err(Error::InvalidInput("CIF loop has no tags".to_owned()));
                }
                if values.len() % tags.len() != 0 {
                    return Err(Error::InvalidInput(format!(
                        "CIF loop has {} values for {} tags",
                        values.len(),
                        tags.len()
                    )));
                }
                let rows = values.chunks(tags.len()).map(<[String]>::to_vec).collect();
                if let Some(block) = blocks.last_mut() {
                    block.loops.push(CifLoop { tags, rows });
                }
                continue;
            }
            if token.starts_with('_') {
                let value = tokens.get(index + 1).cloned().unwrap_or_default();
                if let Some(block) = blocks.last_mut() {
                    block.items.insert(token, value);
                }
                index += 2;
                continue;
            }
            index += 1;
        }
        if blocks.is_empty() {
            return Err(Error::InvalidInput("CIF has no data block".to_owned()));
        }
        Ok(Self { blocks })
    }
}

/// Reads the first data block as a structure.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an invalid cell, missing sites or an
/// unsupported atom loop.
pub fn read_cif(text: &str) -> Result<Structure, Error> {
    let document = CifDocument::parse(text)?;
    let block = document
        .blocks
        .first()
        .ok_or_else(|| Error::InvalidInput("CIF has no data block".to_owned()))?;
    structure_from_block(block)
}

/// Reads a named block as a structure.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the block is missing.
pub fn read_cif_block(text: &str, name: &str) -> Result<Structure, Error> {
    let document = CifDocument::parse(text)?;
    let block = document
        .blocks
        .iter()
        .find(|block| block.name == name)
        .ok_or_else(|| Error::InvalidInput(format!("CIF block not found: {name}")))?;
    structure_from_block(block)
}

/// Reads every block as a structure.
///
/// # Errors
///
/// Returns the error of the first block that cannot become a structure.
pub fn read_cif_all(text: &str) -> Result<Vec<Structure>, Error> {
    let document = CifDocument::parse(text)?;
    document.blocks.iter().map(structure_from_block).collect()
}

/// Writes a structure as a CIF data block.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] if a site has no usable label.
pub fn write_cif(structure: &Structure) -> Result<String, Error> {
    let mut out = String::new();
    out.push_str(&format!("data_{}\n", sanitize_name(structure.name())));
    let [a, b, c, alpha, beta, gamma] = structure.lattice().parameters();
    for (tag, value) in [
        ("_cell_length_a", a),
        ("_cell_length_b", b),
        ("_cell_length_c", c),
        ("_cell_angle_alpha", alpha),
        ("_cell_angle_beta", beta),
        ("_cell_angle_gamma", gamma),
    ] {
        out.push_str(&format!("{tag} {value:.6}\n"));
    }
    out.push_str(&format!("_space_group_name_Hall '{}'\n", structure.hall()));
    if !structure.operations().is_empty() {
        out.push_str("loop_\n_symmetry_equiv_pos_as_xyz\n");
        for operation in structure.operations() {
            out.push_str(&format!("'{}'\n", operation.triplet()));
        }
    }
    out.push_str("loop_\n_atom_site_label\n_atom_site_type_symbol\n");
    out.push_str("_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n");
    out.push_str("_atom_site_occupancy\n_atom_site_U_iso_or_equiv\n_atom_site_charge\n");
    for site in structure.sites() {
        let [x, y, z] = site.frac();
        out.push_str(&format!(
            "{} {} {} {} {} {} {} {}\n",
            site.label(),
            site.element(),
            format_number(x),
            format_number(y),
            format_number(z),
            format_number(site.occupancy()),
            site.u_iso().map_or("?".to_owned(), format_number),
            site.charge().map_or("?".to_owned(), format_number)
        ));
    }
    Ok(out)
}

fn structure_from_block(block: &CifBlock) -> Result<Structure, Error> {
    let a = number(block.item("_cell_length_a")).ok_or_else(|| missing("_cell_length_a"))?;
    let b = number(block.item("_cell_length_b")).ok_or_else(|| missing("_cell_length_b"))?;
    let c = number(block.item("_cell_length_c")).ok_or_else(|| missing("_cell_length_c"))?;
    let alpha =
        number(block.item("_cell_angle_alpha")).ok_or_else(|| missing("_cell_angle_alpha"))?;
    let beta = number(block.item("_cell_angle_beta")).ok_or_else(|| missing("_cell_angle_beta"))?;
    let gamma =
        number(block.item("_cell_angle_gamma")).ok_or_else(|| missing("_cell_angle_gamma"))?;
    let lattice = Lattice::from_parameters(a, b, c, alpha, beta, gamma)?;
    let hall = block
        .item("_space_group_name_Hall")
        .or_else(|| block.item("_symmetry_space_group_name_Hall"))
        .or_else(|| block.item("_symmetry_space_group_name_H-M"))
        .or_else(|| block.item("_space_group_name_H-M_alt"))
        .unwrap_or("P 1")
        .trim_matches('\'')
        .to_owned();
    let mut operations = Vec::new();
    for loop_ in &block.loops {
        let column = loop_.tags.iter().position(|tag| {
            tag == "_symmetry_equiv_pos_as_xyz" || tag == "_space_group_symop_operation_xyz"
        });
        if let Some(column) = column {
            for row in &loop_.rows {
                if let Some(text) = row.get(column) {
                    operations.push(SymOp::parse(text)?);
                }
            }
        }
    }
    let atoms = block
        .loop_with("_atom_site_fract_x")
        .ok_or_else(|| Error::InvalidInput("CIF has no _atom_site loop".to_owned()))?;
    let mut sites = Vec::new();
    for row in 0..atoms.rows.len() {
        let label = atoms
            .value(row, "_atom_site_label")
            .or_else(|| atoms.value(row, "_atom_site_type_symbol"))
            .unwrap_or("X");
        let element = atoms.value(row, "_atom_site_type_symbol").unwrap_or(label);
        let x = number(atoms.value(row, "_atom_site_fract_x"))
            .ok_or_else(|| missing("_atom_site_fract_x"))?;
        let y = number(atoms.value(row, "_atom_site_fract_y"))
            .ok_or_else(|| missing("_atom_site_fract_y"))?;
        let z = number(atoms.value(row, "_atom_site_fract_z"))
            .ok_or_else(|| missing("_atom_site_fract_z"))?;
        let occupancy = number(atoms.value(row, "_atom_site_occupancy")).unwrap_or(1.0);
        let mut site = Site::new(label, element, [x, y, z])?.with_occupancy(occupancy)?;
        if let Some(u_iso) = number(atoms.value(row, "_atom_site_U_iso_or_equiv")) {
            site = site.with_u_iso(u_iso)?;
        }
        if let Some(charge) = number(atoms.value(row, "_atom_site_charge")) {
            site = site.with_charge(charge)?;
        }
        sites.push(site);
    }
    if let Some(aniso) = block.loop_with("_atom_site_aniso_U_11") {
        for (row, site) in sites.iter_mut().enumerate() {
            let values = [
                number(aniso.value(row, "_atom_site_aniso_U_11")),
                number(aniso.value(row, "_atom_site_aniso_U_22")),
                number(aniso.value(row, "_atom_site_aniso_U_33")),
                number(aniso.value(row, "_atom_site_aniso_U_12")),
                number(aniso.value(row, "_atom_site_aniso_U_13")),
                number(aniso.value(row, "_atom_site_aniso_U_23")),
            ];
            if let [
                Some(u11),
                Some(u22),
                Some(u33),
                Some(u12),
                Some(u13),
                Some(u23),
            ] = values
            {
                *site = site.clone().with_u_aniso([u11, u22, u33, u12, u13, u23])?;
            }
        }
    }
    let has_operations = !operations.is_empty();
    let mut builder = StructureBuilder::new(block.name.clone(), lattice).hall(hall.clone());
    for operation in operations {
        builder = builder.operation(operation);
    }
    if !has_operations && hall != "P 1" && hall != "P1" {
        builder = builder.metadata(
            "symmetry_warning",
            serde_json::Value::String(
                "no symmetry operation loop; P1 expansion uses the identity".to_owned(),
            ),
        );
    } else if has_operations {
        builder = builder.metadata(
            "symmetry_source",
            serde_json::Value::String("operation loop".to_owned()),
        );
    }
    for site in sites {
        builder = builder.site(site);
    }
    builder = builder.metadata("cif_block", serde_json::Value::String(block.name.clone()));
    builder.build()
}

fn tokenize(text: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    let mut text_block: Option<String> = None;
    for line in text.lines() {
        if let Some(mut block) = text_block.take() {
            if line.starts_with(';') {
                tokens.push(block);
            } else {
                block.push_str(line);
                block.push('\n');
                text_block = Some(block);
            }
            continue;
        }
        if line.starts_with(';') {
            text_block = Some(line.chars().skip(1).collect());
            continue;
        }
        let mut current = String::new();
        let mut quote: Option<char> = None;
        for character in line.chars() {
            if let Some(active) = quote {
                if character == active {
                    tokens.push(std::mem::take(&mut current));
                    quote = None;
                } else {
                    current.push(character);
                }
                continue;
            }
            match character {
                '\'' | '"' => {
                    if !current.is_empty() {
                        tokens.push(std::mem::take(&mut current));
                    }
                    quote = Some(character);
                }
                '#' => break,
                value if value.is_whitespace() => {
                    if !current.is_empty() {
                        tokens.push(std::mem::take(&mut current));
                    }
                }
                value => current.push(value),
            }
        }
        if !current.is_empty() {
            tokens.push(current);
        }
    }
    if let Some(block) = text_block {
        tokens.push(block);
    }
    tokens
}

fn number(text: Option<&str>) -> Option<f64> {
    let text = text?.trim();
    if text.is_empty() || text == "?" || text == "." {
        return None;
    }
    let text = text.split('(').next().unwrap_or(text);
    text.parse().ok()
}

fn missing(tag: &str) -> Error {
    Error::InvalidInput(format!("CIF is missing {tag}"))
}

fn format_number(value: f64) -> String {
    let text = format!("{value:.6}");
    text.trim_end_matches('0').trim_end_matches('.').to_owned()
}

fn sanitize_name(name: &str) -> String {
    name.chars()
        .map(|character| {
            if character.is_ascii_alphanumeric() || character == '_' {
                character
            } else {
                '_'
            }
        })
        .collect()
}
