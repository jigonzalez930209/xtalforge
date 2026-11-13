//! Streaming CIF atom-loop reader (Task 15.2.2).
//!
//! [`read_cif`](super::read_cif) tokenizes the whole file into a
//! [`CifDocument`](super::CifDocument), which duplicates every atom row as
//! strings. This reader keeps one line and one loop row at a time: the cell
//! items, symmetry operations and atom sites are the only growing state, and
//! the sites are the structure the caller asked for.
//!
//! Supported: one data block, cell items anywhere in the block, a symmetry
//! operation loop, one atom-site loop in fractional or Cartesian coordinates,
//! occupancy, `U_iso`/`B_iso`, charge, chain/residue properties, disorder
//! group/assembly, an anisotropic displacement loop and atom-type oxidation
//! numbers.
//!
//! Not supported, with an explicit error: `_xtal_bond_kind` and
//! `_xtal_bond_atom_1` loops (use [`read_cif`](super::read_cif)), and a
//! Cartesian atom loop that appears before the cell items. Per-atom properties
//! outside the atom-site loop are ignored, as in the document reader.

use std::collections::{BTreeMap, VecDeque};
use std::io::BufRead;

use serde_json::Value;
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder, SymOp};

use super::coincident;
use super::stream_atoms::{AtomColumns, apply_aniso, parse_atom};
use super::tag_matches;
use super::values::{missing, number, same_operations};

/// Reads a CIF from a buffered reader without materialising its tokens.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a truncated loop, a missing cell or
/// atom site, a malformed number, an unsupported explicit-bond loop, or a
/// Cartesian atom loop without a preceding cell.
pub fn read_cif_stream<R: BufRead>(reader: R) -> Result<Structure, Error> {
    let mut tokens = TokenStream::new(reader);
    let mut items: BTreeMap<String, String> = BTreeMap::new();
    let mut sites: Vec<Site> = Vec::new();
    let mut operations: Vec<SymOp> = Vec::new();
    let mut oxidation: BTreeMap<String, f64> = BTreeMap::new();
    let mut aniso_tags: Vec<String> = Vec::new();
    let mut aniso_rows: Vec<Vec<String>> = Vec::new();
    let mut block_name: Option<String> = None;
    let mut first_model: Option<String> = None;
    let mut occupancy_clamped = 0_u32;
    let mut atom_seen = false;
    while let Some(token) = tokens.next()? {
        if let Some(name) = token.strip_prefix("data_") {
            if block_name.is_some() {
                break;
            }
            block_name = Some(name.to_owned());
            continue;
        }
        if token == "loop_" {
            let mut tags = Vec::new();
            while let Some(next) = tokens.peek()? {
                if next.starts_with('_') {
                    tags.push(next);
                    tokens.next()?;
                } else {
                    break;
                }
            }
            if tags.is_empty() {
                return Err(Error::InvalidInput("CIF loop has no tags".to_owned()));
            }
            let kind = LoopKind::of(&tags);
            if kind == LoopKind::ExplicitBonds {
                return Err(Error::InvalidInput(
                    "the streaming CIF reader does not support xtal explicit-bond loops; use read_cif"
                        .to_owned(),
                ));
            }
            let accept_atoms = kind == LoopKind::Atoms && !atom_seen;
            let atom_columns = accept_atoms.then(|| AtomColumns::new(&tags));
            let atom_lattice = match &atom_columns {
                Some(columns) if columns.cartesian => {
                    if cell_items_present(&items) {
                        Some(cell_lattice(&items)?)
                    } else {
                        return Err(Error::InvalidInput(
                            "the streaming CIF reader needs the cell items before a Cartesian atom loop; use read_cif"
                                .to_owned(),
                        ));
                    }
                }
                _ => None,
            };
            let mut atom_row = 0_usize;
            read_rows(&mut tokens, tags.len(), |row| {
                match kind {
                    LoopKind::Symmetry => {
                        if let Some(text) = loop_value(&tags, row, "_symmetry_equiv_pos_as_xyz")
                            .or_else(|| loop_value(&tags, row, "_space_group_symop_operation_xyz"))
                        {
                            operations.push(SymOp::parse(text)?);
                        }
                    }
                    LoopKind::Atoms => {
                        let Some(columns) = &atom_columns else {
                            return Ok(());
                        };
                        let site = parse_atom(
                            row,
                            atom_row,
                            columns,
                            atom_lattice.as_ref(),
                            &mut first_model,
                            &mut occupancy_clamped,
                        )?;
                        atom_row += 1;
                        if let Some(site) = site {
                            sites.push(site);
                        }
                    }
                    LoopKind::Oxidation => {
                        if let (Some(symbol), Some(value)) = (
                            loop_value(&tags, row, "_atom_type_symbol"),
                            loop_value(&tags, row, "_atom_type_oxidation_number"),
                        ) && let Some(value) = number(Some(value))
                        {
                            oxidation.insert(symbol.to_owned(), value);
                        }
                    }
                    LoopKind::Aniso => {
                        if aniso_tags.is_empty() {
                            aniso_tags = tags.clone();
                        }
                        aniso_rows.push(row.to_vec());
                    }
                    LoopKind::Other => {}
                    LoopKind::ExplicitBonds => {}
                }
                Ok(())
            })?;
            atom_seen |= accept_atoms;
            continue;
        }
        if token.starts_with('_') {
            let value = tokens.next()?.unwrap_or_default();
            items.insert(token, value);
            continue;
        }
    }
    let Some(name) = block_name else {
        return Err(Error::InvalidInput("CIF has no data block".to_owned()));
    };
    if !atom_seen {
        return Err(Error::InvalidInput("CIF has no _atom_site loop".to_owned()));
    }
    let lattice = cell_lattice(&items)?;
    let hall = item(&items, "_space_group_name_Hall")
        .or_else(|| item(&items, "_symmetry_space_group_name_Hall"))
        .or_else(|| item(&items, "_symmetry_space_group_name_H-M"))
        .or_else(|| item(&items, "_space_group_name_H-M_alt"))
        .unwrap_or("P 1")
        .trim_matches('\'')
        .to_owned();
    for site in &mut sites {
        if site.charge().is_none()
            && let Some(charge) = oxidation.get(site.element())
        {
            *site = site.clone().with_charge(*charge)?;
        }
    }
    apply_aniso(&mut sites, &aniso_tags, &aniso_rows)?;
    let mut source = "operation loop";
    if operations.is_empty() && hall != "P 1" && hall != "P1" {
        match xtalforge_symmetry::SpaceGroup::parse(&hall) {
            Ok(group) => {
                operations = group.operations;
                source = "Hall symbol";
            }
            Err(_) => {
                source = "identity";
            }
        }
    }
    let has_operations = !operations.is_empty();
    let operation_count = operations.len();
    let operation_warning = if has_operations && hall != "P 1" && hall != "P1" {
        match xtalforge_symmetry::SpaceGroup::parse(&hall) {
            Ok(group) if group.order() != operation_count => Some(format!(
                "operation loop has {operation_count} operations but {hall:?} (Hall {}) has {}",
                group.hall_symbol,
                group.order()
            )),
            Ok(group) if !same_operations(&group.operations, &operations) => Some(format!(
                "operation loop does not match {hall:?} (Hall {})",
                group.hall_symbol
            )),
            Ok(_) => None,
            Err(_) => Some(format!(
                "operation loop is present but the space-group symbol {hall:?} was not recognised"
            )),
        }
    } else {
        None
    };
    let mut builder = StructureBuilder::new(name.clone(), lattice).hall(hall.clone());
    for operation in operations {
        builder = builder.operation(operation);
    }
    if !has_operations && hall != "P 1" && hall != "P1" {
        builder = builder.metadata(
            "symmetry_warning",
            Value::String(
                "no symmetry operation loop and the Hall symbol was not recognised; P1 expansion uses the identity".to_owned(),
            ),
        );
    } else if has_operations {
        builder = builder.metadata("symmetry_source", Value::String(source.to_owned()));
        if let Some(warning) = operation_warning {
            builder = builder.metadata("symmetry_warning", Value::String(warning));
        }
    }
    for site in sites {
        builder = builder.site(site);
    }
    builder = builder.metadata("cif_block", Value::String(name.clone()));
    if let Some(formula) = item(&items, "_chemical_formula_sum")
        .or_else(|| item(&items, "_chemical_formula_structural"))
    {
        builder = builder.metadata(
            "chemical_formula",
            Value::String(formula.trim_matches('\'').to_owned()),
        );
    }
    if occupancy_clamped > 0 {
        builder = builder.metadata("occupancy_clamped", Value::from(occupancy_clamped));
    }
    let structure = builder.clone().build()?;
    match coincident::warning(&structure) {
        Some(message) => builder
            .metadata("coincident_atoms_warning", Value::String(message))
            .build(),
        None => Ok(structure),
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum LoopKind {
    Symmetry,
    Atoms,
    Oxidation,
    Aniso,
    ExplicitBonds,
    Other,
}

impl LoopKind {
    fn of(tags: &[String]) -> Self {
        let has = |name: &str| tags.iter().any(|tag| tag_matches(tag, name));
        if has("_symmetry_equiv_pos_as_xyz") || has("_space_group_symop_operation_xyz") {
            return Self::Symmetry;
        }
        if has("_atom_site_fract_x") || has("_atom_site_Cartn_x") {
            return Self::Atoms;
        }
        if has("_atom_type_oxidation_number") {
            return Self::Oxidation;
        }
        if has("_atom_site_aniso_U_11") {
            return Self::Aniso;
        }
        if has("_xtal_bond_kind") || has("_xtal_bond_atom_1") {
            return Self::ExplicitBonds;
        }
        Self::Other
    }
}

fn item<'a>(items: &'a BTreeMap<String, String>, tag: &str) -> Option<&'a str> {
    items
        .iter()
        .find(|(key, _)| tag_matches(key, tag))
        .map(|(_, value)| value.as_str())
}

fn cell_items_present(items: &BTreeMap<String, String>) -> bool {
    [
        "_cell_length_a",
        "_cell_length_b",
        "_cell_length_c",
        "_cell_angle_alpha",
        "_cell_angle_beta",
        "_cell_angle_gamma",
    ]
    .iter()
    .all(|tag| item(items, tag).is_some())
}

fn cell_lattice(items: &BTreeMap<String, String>) -> Result<Lattice, Error> {
    let a = number(item(items, "_cell_length_a")).ok_or_else(|| missing("_cell_length_a"))?;
    let b = number(item(items, "_cell_length_b")).ok_or_else(|| missing("_cell_length_b"))?;
    let c = number(item(items, "_cell_length_c")).ok_or_else(|| missing("_cell_length_c"))?;
    let alpha =
        number(item(items, "_cell_angle_alpha")).ok_or_else(|| missing("_cell_angle_alpha"))?;
    let beta =
        number(item(items, "_cell_angle_beta")).ok_or_else(|| missing("_cell_angle_beta"))?;
    let gamma =
        number(item(items, "_cell_angle_gamma")).ok_or_else(|| missing("_cell_angle_gamma"))?;
    Lattice::from_parameters(a, b, c, alpha, beta, gamma)
}

pub(super) fn loop_value<'a>(tags: &[String], row: &'a [String], tag: &str) -> Option<&'a str> {
    let column = tags
        .iter()
        .position(|candidate| tag_matches(candidate, tag))?;
    row.get(column).map(String::as_str)
}

fn read_rows(
    tokens: &mut TokenStream<impl BufRead>,
    width: usize,
    mut visit: impl FnMut(&[String]) -> Result<(), Error>,
) -> Result<(), Error> {
    loop {
        if let Some(next) = tokens.peek()?
            && is_delimiter(&next)
        {
            return Ok(());
        }
        let mut row = Vec::with_capacity(width);
        while row.len() < width {
            match tokens.next()? {
                Some(value) => row.push(value),
                None => {
                    if row.is_empty() {
                        return Ok(());
                    }
                    return Err(Error::InvalidInput(
                        "CIF loop ended in the middle of a row".to_owned(),
                    ));
                }
            }
        }
        visit(&row)?;
    }
}

fn is_delimiter(token: &str) -> bool {
    token.starts_with('_')
        || token == "loop_"
        || token.starts_with("data_")
        || token.starts_with("save_")
}

struct TokenStream<R: BufRead> {
    reader: R,
    buffer: String,
    queue: VecDeque<String>,
    done: bool,
}

impl<R: BufRead> TokenStream<R> {
    fn new(reader: R) -> Self {
        Self {
            reader,
            buffer: String::new(),
            queue: VecDeque::new(),
            done: false,
        }
    }

    fn next(&mut self) -> Result<Option<String>, Error> {
        self.fill()?;
        Ok(self.queue.pop_front())
    }

    fn peek(&mut self) -> Result<Option<String>, Error> {
        self.fill()?;
        Ok(self.queue.front().cloned())
    }

    fn fill(&mut self) -> Result<(), Error> {
        while self.queue.is_empty() && !self.done {
            let Some(line) = self.read_line()? else {
                self.done = true;
                break;
            };
            if line.starts_with(';') {
                let mut block: String = line.chars().skip(1).collect();
                while let Some(next) = self.read_line()? {
                    if next.starts_with(';') {
                        break;
                    }
                    block.push_str(&next);
                    block.push('\n');
                }
                self.queue.push_back(block);
            } else {
                tokenize_line(&line, &mut self.queue);
            }
        }
        Ok(())
    }

    fn read_line(&mut self) -> Result<Option<String>, Error> {
        self.buffer.clear();
        let read = self
            .reader
            .read_line(&mut self.buffer)
            .map_err(|error| Error::InvalidInput(format!("could not read the CIF: {error}")))?;
        if read == 0 {
            return Ok(None);
        }
        Ok(Some(self.buffer.trim_end_matches(['\n', '\r']).to_owned()))
    }
}

fn tokenize_line(line: &str, queue: &mut VecDeque<String>) {
    let mut current = String::new();
    let mut quote: Option<char> = None;
    for character in line.chars() {
        if let Some(active) = quote {
            if character == active {
                queue.push_back(std::mem::take(&mut current));
                quote = None;
            } else {
                current.push(character);
            }
            continue;
        }
        match character {
            '\'' | '"' => {
                if !current.is_empty() {
                    queue.push_back(std::mem::take(&mut current));
                }
                quote = Some(character);
            }
            '#' => break,
            value if value.is_whitespace() => {
                if !current.is_empty() {
                    queue.push_back(std::mem::take(&mut current));
                }
            }
            value => current.push(value),
        }
    }
    if !current.is_empty() {
        queue.push_back(current);
    }
}
