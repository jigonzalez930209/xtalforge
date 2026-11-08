//! The `mof` command family: catalogues and framework assembly.

use clap::Subcommand;
use serde_json::json;

use crate::format::StructureFormat;
use crate::{Cli, CliError, output};

/// Subcommands of `xtalforge mof`.
#[derive(Debug, Subcommand)]
pub(crate) enum MofAction {
    /// List the topology catalogue.
    Topologies,
    /// List the building blocks.
    Blocks,
    /// Assemble a framework from a topology, node and linker.
    Build {
        /// Topology name.
        topology: String,
        /// Node block name.
        node: String,
        /// Linker block name.
        linker: String,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format; defaults to the output extension or CIF.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
}

/// Runs one `mof` action.
///
/// # Errors
///
/// Returns [`CliError`] for catalogue or builder failures.
pub(crate) fn run(cli: &Cli, action: &MofAction) -> Result<(), CliError> {
    match action {
        MofAction::Topologies => {
            let topologies = xtalforge::mof::topologies().map_err(CliError::mof)?;
            let items: Vec<serde_json::Value> = topologies
                .iter()
                .map(|topology| {
                    json!({
                        "name": topology.name,
                        "nodes": topology.nodes.len(),
                        "edges": topology.edges.len(),
                    })
                })
                .collect();
            if cli.json {
                println!("{}", json!(items));
            } else {
                for topology in topologies {
                    println!(
                        "{:<8} {} nodes, {} edges",
                        topology.name,
                        topology.nodes.len(),
                        topology.edges.len()
                    );
                }
            }
            Ok(())
        }
        MofAction::Blocks => {
            let blocks = xtalforge::mof::blocks().map_err(CliError::mof)?;
            let items: Vec<serde_json::Value> = blocks
                .iter()
                .map(|block| {
                    json!({
                        "name": block.name,
                        "kind": format!("{:?}", block.kind),
                        "atoms": block.elements.len(),
                        "connections": block.connection_points.len(),
                    })
                })
                .collect();
            if cli.json {
                println!("{}", json!(items));
            } else {
                for block in blocks {
                    println!(
                        "{:<12} {:?} {} atoms, {} connections",
                        block.name,
                        block.kind,
                        block.elements.len(),
                        block.connection_points.len()
                    );
                }
            }
            Ok(())
        }
        MofAction::Build {
            topology,
            node,
            linker,
            out,
            format,
        } => {
            let request =
                xtalforge::mof::MofRequest::new(topology.clone(), node.clone(), linker.clone());
            let built = xtalforge::mof::build_mof(&request).map_err(CliError::mof)?;
            if let Some(destination) = out {
                output::emit_structure(cli, &built.structure, Some(destination), *format)?;
                if !cli.quiet && !cli.json {
                    println!(
                        "built {} ({}/{}) with {} warnings",
                        built.structure.name(),
                        built.topology,
                        built.linker,
                        built.warnings.len()
                    );
                }
                return Ok(());
            }
            if cli.json {
                println!(
                    "{}",
                    json!({
                        "structure": built.structure,
                        "topology": built.topology,
                        "node": built.node,
                        "linker": built.linker,
                        "warnings": built.warnings,
                    })
                );
                Ok(())
            } else {
                output::emit_structure(cli, &built.structure, None, *format)
            }
        }
    }
}
