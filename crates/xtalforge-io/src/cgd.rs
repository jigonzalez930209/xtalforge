//! CGD periodic nets, the Systre and RCSR format (Task 3.4.1 and 5.6.1).
//!
//! Ported from Crystal-Builder `xtal/io/cgd.py` (MIT). The reader accepts
//! both dialects the RCSR file uses: an `EDGE` as a pair of points, or as a
//! node label followed by one neighbour position. Keywords may stand alone
//! with their rows underneath, and `# EDGE_CENTER` lines are comments.
#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked above
use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

/// Every keyword a `.cgd` file uses.
const KEYWORDS: [&str; 9] = [
    "crystal",
    "end",
    "name",
    "group",
    "cell",
    "node",
    "atom",
    "edge",
    "edge_center",
];

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

/// One edge as its two endpoint points.
type EdgePoints = (Vec<f64>, Vec<f64>);

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

impl CgdEntry {
    /// Coordinate dimension of the entry (2 or 3).
    #[must_use]
    pub fn dimension(&self) -> usize {
        if let Some(node) = self.nodes.first() {
            node.frac.len()
        } else if self.cell.as_ref().is_some_and(|cell| cell.len() == 3) {
            2
        } else {
            3
        }
    }
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
        let head = line
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_ascii_lowercase();
        if head == "crystal" {
            if let Some((start, block)) = current.take() {
                match parse_block(start + 1, &block) {
                    Ok(entry) => entries.push(entry),
                    Err(error) => problems.push(error.message()),
                }
            }
            current = Some((index + 1, Vec::new()));
            continue;
        }
        if head == "end" {
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
    let mut cell: Option<Vec<f64>> = None;
    let mut nodes = Vec::new();
    let mut raw_edges: Vec<Vec<f64>> = Vec::new();
    let mut pending = String::new();
    for row in rows {
        let tokens: Vec<&str> = row.split_whitespace().collect();
        let head = tokens.first().copied().unwrap_or("").to_ascii_lowercase();
        let rest = if KEYWORDS.contains(&head.as_str()) {
            pending = head;
            tokens.get(1..).unwrap_or_default()
        } else {
            tokens.as_slice()
        };
        if rest.is_empty() {
            continue;
        }
        match pending.as_str() {
            "name" => name = rest[0].to_owned(),
            "group" => group = rest[0].to_owned(),
            "cell" => {
                cell = Some(rest.iter().filter_map(|value| value.parse().ok()).collect());
            }
            "node" | "atom" => nodes.push(parse_node(rest, line)?),
            "edge" => raw_edges.push(rest.iter().filter_map(|value| value.parse().ok()).collect()),
            "edge_center" => {}
            _ => {}
        }
    }
    if name.is_empty() {
        return Err(Error::InvalidInput(format!(
            "CGD block at line {line} has no NAME"
        )));
    }
    let dimension = if let Some(node) = nodes.first() {
        node.frac.len()
    } else if cell.as_ref().is_some_and(|value| value.len() == 3) {
        2
    } else {
        3
    };
    let edges = resolve_edges(&raw_edges, &nodes, dimension, &name, line)?;
    if edges.is_empty() {
        return Err(Error::InvalidInput(format!(
            "CGD block {name} has no edges"
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

fn parse_node(rest: &[&str], line: usize) -> Result<CgdNode, Error> {
    if rest.len() < 3 {
        return Err(Error::InvalidInput(format!(
            "CGD node at line {line} is incomplete"
        )));
    }
    let label = rest[0].to_owned();
    let coordination = rest[1]
        .parse::<usize>()
        .map_err(|_| Error::InvalidInput(format!("CGD node {label} has no coordination number")))?;
    let frac: Vec<f64> = rest[2..]
        .iter()
        .filter_map(|value| value.parse().ok())
        .collect();
    if frac.len() < 2 {
        return Err(Error::InvalidInput(format!(
            "CGD node {label} has no coordinates"
        )));
    }
    Ok(CgdNode {
        label,
        coordination,
        frac,
    })
}

fn resolve_edges(
    raw: &[Vec<f64>],
    nodes: &[CgdNode],
    dimension: usize,
    name: &str,
    line: usize,
) -> Result<Vec<EdgePoints>, Error> {
    let mut edges = Vec::with_capacity(raw.len());
    for values in raw {
        if values.len() == 2 * dimension {
            edges.push((values[..dimension].to_vec(), values[dimension..].to_vec()));
        } else if values.len() == dimension + 1 {
            let label = format!("{}", values[0] as i64);
            let node = nodes
                .iter()
                .find(|node| node.label == label)
                .ok_or_else(|| {
                    Error::InvalidInput(format!(
                        "CGD block {name} at line {line} has an edge at node {label}, which it never declares"
                    ))
                })?;
            edges.push((node.frac.clone(), values[1..].to_vec()));
        } else {
            return Err(Error::InvalidInput(format!(
                "CGD block {name} at line {line} has an edge of {} numbers, which is neither a pair of {dimension}D points nor a node and a neighbour",
                values.len()
            )));
        }
    }
    Ok(edges)
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
