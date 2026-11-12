//! The `mof` command family: catalogues and framework assembly.

use clap::Subcommand;
use serde_json::json;

use crate::format::StructureFormat;
use crate::{Cli, CliError, output};

/// Adsorption actions under `mof adsorption`.
#[derive(Debug, Subcommand)]
pub(crate) enum AdsorptionAction {
    /// List the embedded TraPPE guests.
    Guests,
    /// Host UFF/EQeq parameters and the guest names.
    ForceField { path: String },
    /// Widom insertion and Henry coefficient.
    Widom {
        path: String,
        #[arg(long, default_value = "CO2")]
        guest: String,
        #[arg(long, default_value_t = 298.0)]
        temperature: f64,
        #[arg(long, default_value_t = 1000)]
        samples: usize,
        #[arg(long, default_value_t = 0)]
        seed: u64,
    },
    /// Single-component GCMC isotherm.
    Gcmc {
        path: String,
        #[arg(long, default_value = "CO2")]
        guest: String,
        #[arg(long, default_value_t = 298.0)]
        temperature: f64,
        /// Comma-separated pressures in bar.
        #[arg(long, default_value = "1,5,10")]
        pressures: String,
        #[arg(long, default_value_t = 1000)]
        steps: usize,
        #[arg(long, default_value_t = 0)]
        seed: u64,
    },
    /// Guest energy grid summary.
    Grid {
        path: String,
        #[arg(long, default_value = "CO2")]
        guest: String,
        #[arg(long, default_value_t = 1.0)]
        spacing: f64,
    },
    /// Write a RASPA input set into a directory.
    Raspa {
        path: String,
        #[arg(long, default_value = "CO2")]
        guest: String,
        #[arg(long, default_value_t = 298.0)]
        temperature: f64,
        #[arg(long, default_value_t = 1.0)]
        pressure: f64,
        #[arg(long)]
        dir: Option<String>,
    },
}

/// Post-synthetic modification actions under `mof modify`.
#[derive(Debug, Subcommand)]
pub(crate) enum ModifyAction {
    /// Add -NH2, -OH, -CH3 or -F groups.
    Functionalize {
        path: String,
        #[arg(long, default_value = "amino")]
        group: String,
        #[arg(long)]
        positions: Vec<usize>,
        #[arg(long, default_value_t = 1.0)]
        fraction: f64,
        #[arg(long, default_value_t = 0)]
        seed: u64,
        #[arg(long)]
        out: Option<String>,
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Missing linker and cluster defects with capping.
    Defects {
        path: String,
        #[arg(long, default_value_t = 0.0)]
        linker_fraction: f64,
        #[arg(long, default_value_t = 0.0)]
        cluster_fraction: f64,
        #[arg(long, default_value = "formate")]
        cap: String,
        #[arg(long, default_value_t = 0)]
        seed: u64,
        #[arg(long)]
        out: Option<String>,
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// n-fold interpenetration.
    Interpenetrate {
        path: String,
        #[arg(long, default_value_t = 2)]
        folds: u32,
        #[arg(long)]
        out: Option<String>,
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// AA/AB layer stacking.
    StackLayers {
        path: String,
        #[arg(long, default_value_t = 2)]
        layers: u32,
        /// Offset `u,v` for odd layers.
        #[arg(long, default_value = "0.5,0.5")]
        offset: String,
        #[arg(long)]
        out: Option<String>,
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Mixed-metal or mixed-linker framework.
    Mixed {
        path: String,
        /// Replacement `old:new`, repeated.
        #[arg(long = "map")]
        map: Vec<String>,
        #[arg(long)]
        ordered: bool,
        #[arg(long, default_value_t = 1.0)]
        fraction: f64,
        #[arg(long, default_value_t = 0)]
        seed: u64,
        #[arg(long)]
        out: Option<String>,
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
}

/// Database actions under `mof database`.
#[derive(Debug, Subcommand)]
pub(crate) enum DatabaseAction {
    /// Read a CoRE MOF or QMOF directory and report its entries.
    Read {
        path: String,
        #[arg(long, default_value = "core")]
        kind: String,
    },
    /// Clean, deconstruct and screen a directory.
    Screen {
        path: String,
        #[arg(long, default_value = "CO2")]
        guest: String,
        #[arg(long, default_value_t = 0.8)]
        spacing: f64,
        #[arg(long, default_value_t = 100)]
        samples: usize,
        #[arg(long, default_value_t = 0)]
        limit: usize,
    },
    /// Deduplicate a list of structure paths.
    Deduplicate { paths: Vec<String> },
}

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
    /// Check a MOF for overlaps, fragments, valence, disorder and charge.
    Check {
        /// Input structure.
        path: String,
    },
    /// Remove free solvent and optionally counter-ions.
    SolventRemove {
        /// Input structure.
        path: String,
        /// Remove free counter-ions too.
        #[arg(long)]
        remove_ions: bool,
        /// Output path.
        #[arg(long)]
        out: Option<String>,
        /// Output format.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Resolve disorder into one ordered configuration.
    DisorderResolve {
        /// Input structure.
        path: String,
        /// Output path.
        #[arg(long)]
        out: Option<String>,
        /// Output format.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
    /// Metal oxidation states and framework charge balance.
    ChargeBalance {
        /// Input structure.
        path: String,
    },
    /// Deconstruct into nodes, linkers, solvent and the underlying net.
    Deconstruct {
        /// Input structure.
        path: String,
    },
    /// MOFid and MOFkey identity.
    Mofid {
        /// Input structure.
        path: String,
    },
    /// Coordinatively unsaturated metal sites.
    OpenMetalSites {
        /// Input structure.
        path: String,
    },
    /// Linker, node and framework descriptors.
    Descriptors {
        /// Input structure.
        path: String,
    },
    /// Native PLD, LCD, ASA, accessible volume and channel dimension.
    Pores {
        /// Input structure.
        path: String,
        /// Grid spacing in Angstrom.
        #[arg(long, default_value_t = 0.8)]
        spacing: f64,
        /// Probe radius in Angstrom.
        #[arg(long, default_value_t = 1.2)]
        probe: f64,
    },
    /// Host-guest force fields, Widom, GCMC, grids and RASPA export.
    Adsorption {
        #[command(subcommand)]
        action: AdsorptionAction,
    },
    /// Post-synthetic modification, defects and framework variants.
    Modify {
        #[command(subcommand)]
        action: ModifyAction,
    },
    /// CoRE MOF/QMOF readers, screening and deduplication.
    Database {
        #[command(subcommand)]
        action: DatabaseAction,
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
            if !cli.quiet {
                for warning in &built.warnings {
                    eprintln!("warning: {warning}");
                }
            }
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
        MofAction::Check { path } => crate::mof_platform_commands::check(cli, path),
        MofAction::SolventRemove {
            path,
            remove_ions,
            out,
            format,
        } => crate::mof_platform_commands::solvent_remove(
            cli,
            path,
            *remove_ions,
            out.as_deref(),
            *format,
        ),
        MofAction::DisorderResolve { path, out, format } => {
            crate::mof_platform_commands::disorder(cli, path, out.as_deref(), *format)
        }
        MofAction::ChargeBalance { path } => {
            crate::mof_platform_commands::charge_balance(cli, path)
        }
        MofAction::Deconstruct { path } => crate::mof_platform_commands::deconstruct(cli, path),
        MofAction::Mofid { path } => crate::mof_platform_commands::mofid(cli, path),
        MofAction::OpenMetalSites { path } => {
            crate::mof_platform_commands::open_metal_sites(cli, path)
        }
        MofAction::Descriptors { path } => crate::mof_platform_commands::descriptors(cli, path),
        MofAction::Pores {
            path,
            spacing,
            probe,
        } => crate::mof_platform_commands::pores(cli, path, *spacing, *probe),
        MofAction::Adsorption { action } => crate::mof_adsorption_commands::run(cli, action),
        MofAction::Modify { action } => crate::mof_platform_commands::modify(cli, action),
        MofAction::Database { action } => crate::mof_platform_commands::database(cli, action),
    }
}
