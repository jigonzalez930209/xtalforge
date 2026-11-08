//! CGD periodic nets, the Systre and RCSR format (Task 3.4.1).

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked above
use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

/// One vertex of a net.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CgdNode {
    /// Node label.
    pub label: String,
    /// Declared coordination number.
    pub coordination: usize,
    /// Fractional coordinates.
    pub frac: Vec<f64>,
}

/// One `CRYSTAL ... END` block.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CgdEntry {
    /// Net name.
    pub name: String,
    /// Space-group symbol as written.
    pub group: String,
    /// Vertices.
    pub nodes: Vec<CgdNode>,
    /// Edges as pairs of fractional points.
    pub edges: Vec<(Vec<f64>, Vec<f64>)>,
    /// Optional cell parameters.
    pub cell: Option<Vec<f64>>,
    /// Source line of the `CRYSTAL` keyword.
    pub line: usize,
}

/// A parsed CGD file and the blocks that failed.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CgdFile {
    /// Successfully read entries.
    pub entries: Vec<CgdEntry>,
    /// One message per skipped or incomplete block.
    pub problems: Vec<String>,
}

/// Reads every `CRYSTAL ... END` block; malformed blocks are reported, not
/// fatal.
#[must_use]
pub fn read_cgd(text: &str) -> CgdFile {
    let mut entries = Vec::new();
    let mut problems = Vec::new();
    let mut current: Option<(usize, Vec<String>)> = None;
    for (index, raw) in text.lines().enumerate() {
        let line = raw.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        let upper = line.to_ascii_uppercase();
        if upper.starts_with("CRYSTAL") {
            if let Some((start, block)) = current.take() {
                match parse_block(start + 1, &block) {
                    Ok(entry) => entries.push(entry),
                    Err(error) => problems.push(error.message()),
                }
            }
            current = Some((index + 1, Vec::new()));
            continue;
        }
        if upper == "END" {
            if let Some((start, block)) = current.take() {
                match parse_block(start + 1, &block) {
                    Ok(entry) => entries.push(entry),
                    Err(error) => problems.push(error.message()),
                }
            }
            continue;
        }
        if let Some((_, block)) = current.as_mut() {
            block.push(line.to_owned());
        }
    }
    if let Some((start, block)) = current {
        match parse_block(start + 1, &block) {
            Ok(entry) => entries.push(entry),
            Err(error) => problems.push(error.message()),
        }
    }
    CgdFile { entries, problems }
}

fn parse_block(line: usize, rows: &[String]) -> Result<CgdEntry, Error> {
    let mut name = String::new();
    let mut group = "P1".to_owned();
    let mut cell = None;
    let mut nodes = Vec::new();
    let mut edges = Vec::new();
    let mut keyword = String::new();
    for row in rows {
        let first = row
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_ascii_uppercase();
        let rest = row
            .split_once(char::is_whitespace)
            .map_or("", |(_, rest)| rest);
        match first.as_str() {
            "NAME" => name = rest.trim().to_owned(),
            "GROUP" => group = rest.trim().to_owned(),
            "CELL" => {
                cell = Some(
                    rest.split_whitespace()
                        .filter_map(|value| value.parse().ok())
                        .collect(),
                );
            }
            "NODE" | "ATOM" | "EDGE" | "EDGE_CENTER" => {
                keyword = first.clone();
                if first == "EDGE" || first == "EDGE_CENTER" {
                    if let Some(edge) = parse_edge(rest) {
                        edges.push(edge);
                    }
                } else if let Some(node) = parse_node(rest) {
                    nodes.push(node);
                }
            }
            _ => {
                if keyword == "EDGE" || keyword == "EDGE_CENTER" {
                    if let Some(edge) = parse_edge(row) {
                        edges.push(edge);
                    }
                } else if (keyword == "NODE" || keyword == "ATOM")
                    && let Some(node) = parse_node(row)
                {
                    nodes.push(node);
                }
            }
        }
    }
    if name.is_empty() {
        return Err(Error::InvalidInput(format!(
            "CGD block at line {line} has no NAME"
        )));
    }
    if nodes.is_empty() && edges.is_empty() {
        return Err(Error::InvalidInput(format!(
            "CGD block {name} has no nodes or edges"
        )));
    }
    Ok(CgdEntry {
        name,
        group,
        nodes,
        edges,
        cell,
        line,
    })
}

fn parse_node(rest: &str) -> Option<CgdNode> {
    let parts: Vec<&str> = rest.split_whitespace().collect();
    if parts.len() < 3 {
        return None;
    }
    let (label, coordination, start) = if parts[0].chars().all(|c| c.is_ascii_digit()) {
        (format!("N{}", parts[0]), parts[1].parse().ok()?, 2)
    } else {
        (parts[0].to_owned(), parts[1].parse().ok()?, 2)
    };
    let frac: Vec<f64> = parts[start..]
        .iter()
        .filter_map(|value| value.parse().ok())
        .collect();
    if frac.len() < 2 {
        return None;
    }
    Some(CgdNode {
        label,
        coordination,
        frac,
    })
}

fn parse_edge(rest: &str) -> Option<(Vec<f64>, Vec<f64>)> {
    let values: Vec<f64> = rest
        .split_whitespace()
        .filter_map(|value| value.parse().ok())
        .collect();
    if values.len() == 6 {
        Some((values[..3].to_vec(), values[3..6].to_vec()))
    } else if values.len() == 4 {
        Some((values[..3].to_vec(), values[3..].to_vec()))
    } else {
        None
    }
}

/// Writes one block in the common dialect, in `P1`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for empty names or non-finite coordinates.
pub fn write_cgd(entries: &[CgdEntry]) -> Result<String, Error> {
    let mut out = String::new();
    for entry in entries {
        if entry.name.trim().is_empty() {
            return Err(Error::InvalidInput("a CGD net needs a name".to_owned()));
        }
        out.push_str(&format!("CRYSTAL\nNAME {}\nGROUP P1\n", entry.name));
        if let Some(cell) = &entry.cell {
            out.push_str("CELL ");
            out.push_str(
                &cell
                    .iter()
                    .map(|value| format!("{value:.6}"))
                    .collect::<Vec<_>>()
                    .join(" "),
            );
            out.push('\n');
        }
        for node in &entry.nodes {
            if !node.frac.iter().all(|value| value.is_finite()) {
                return Err(Error::InvalidInput(format!(
                    "node {} has non-finite coordinates",
                    node.label
                )));
            }
            out.push_str(&format!(
                "NODE {} {} {}\n",
                node.label,
                node.coordination,
                node.frac
                    .iter()
                    .map(|value| format!("{value:.6}"))
                    .collect::<Vec<_>>()
                    .join(" ")
            ));
        }
        for (first, second) in &entry.edges {
            out.push_str(&format!(
                "EDGE {} {}\n",
                first
                    .iter()
                    .map(|value| format!("{value:.6}"))
                    .collect::<Vec<_>>()
                    .join(" "),
                second
                    .iter()
                    .map(|value| format!("{value:.6}"))
                    .collect::<Vec<_>>()
                    .join(" ")
            ));
        }
        out.push_str("END\n");
    }
    Ok(out)
}
