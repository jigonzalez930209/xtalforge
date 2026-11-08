//! The `analysis` command family: porosity, nets and Zeo++.

use clap::Subcommand;
use serde_json::json;
use xtalforge::analysis::{
    RcsrIndex, ZeoOptions, default_radii, distance_grid, overlay_net, quotient_graph,
};

use crate::format::StructureFormat;
use crate::{Cli, CliError, input, output};

/// Subcommands of `xtalforge analysis`.
#[derive(Debug, Subcommand)]
pub(crate) enum AnalysisAction {
    /// Surface-distance grid and porosity numbers.
    Porosity {
        /// Input path, or `-` for standard input.
        path: String,
        /// Grid spacing in Angstrom.
        #[arg(long, default_value_t = 0.35)]
        spacing: f64,
        /// Probe radius in Angstrom.
        #[arg(long, default_value_t = 1.2)]
        probe: f64,
    },
    /// Read a CGD file and identify its nets.
    Nets {
        /// CGD file path.
        path: String,
        /// Endpoint matching tolerance.
        #[arg(long, default_value_t = 0.01)]
        tolerance: f64,
    },
    /// Run Zeo++ on a structure.
    Zeo {
        /// Input path, or `-` for standard input.
        path: String,
        /// Probe radius in Angstrom.
        #[arg(long, default_value_t = 1.2)]
        probe: f64,
        /// Number of samples.
        #[arg(long, default_value_t = 1000)]
        samples: u32,
    },
    /// Write a net overlay as JSON.
    Overlay {
        /// CGD file path.
        path: String,
        /// Endpoint matching tolerance.
        #[arg(long, default_value_t = 0.01)]
        tolerance: f64,
        /// Output path; defaults to standard output.
        #[arg(long, value_name = "PATH")]
        out: Option<String>,
        /// Output format for a one-net structure.
        #[arg(long, value_enum)]
        format: Option<StructureFormat>,
    },
}

/// Runs one `analysis` action.
///
/// # Errors
///
/// Returns [`CliError`] for read, analysis or output failures.
pub(crate) fn run(cli: &Cli, action: &AnalysisAction) -> Result<(), CliError> {
    match action {
        AnalysisAction::Porosity {
            path,
            spacing,
            probe,
        } => porosity(cli, path, *spacing, *probe),
        AnalysisAction::Nets { path, tolerance } => nets(cli, path, *tolerance),
        AnalysisAction::Zeo {
            path,
            probe,
            samples,
        } => zeo(cli, path, *probe, *samples),
        AnalysisAction::Overlay {
            path,
            tolerance,
            out,
            format,
        } => overlay(cli, path, *tolerance, out.as_deref(), *format),
    }
}

fn porosity(cli: &Cli, path: &str, spacing: f64, probe: f64) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let radii = default_radii(&structure).map_err(CliError::analysis)?;
    let grid = distance_grid(&structure, &radii, spacing).map_err(CliError::analysis)?;
    let accessible = grid.accessible_volume(probe).map_err(CliError::analysis)?;
    let pore = grid.pore_volume();
    let surface = grid.surface_area();
    let largest = grid.largest_included_diameter();
    if cli.json {
        println!(
            "{}",
            json!({
                "spacing": spacing,
                "probe_radius": probe,
                "cell_volume": grid.cell_volume(),
                "accessible_volume": accessible,
                "pore_volume": pore,
                "surface_area": surface,
                "largest_included_diameter": largest,
            })
        );
    } else {
        println!("cell volume {:.3} A^3", grid.cell_volume());
        println!("accessible volume {accessible:.3} A^3 (probe {probe} A)");
        println!("pore volume {pore:.3} A^3");
        println!("surface area {surface:.3} A^2");
        println!("largest included diameter {largest:.3} A");
    }
    Ok(())
}

fn nets(cli: &Cli, path: &str, tolerance: f64) -> Result<(), CliError> {
    let text = input::read_source(path).map_err(CliError::core)?;
    let file = xtalforge::io::read_cgd(&text);
    let index = RcsrIndex::from_subset();
    let mut entries = Vec::new();
    for entry in &file.entries {
        let graph = quotient_graph(entry, tolerance).map_err(CliError::analysis)?;
        let invariants = graph.invariants();
        entries.push(json!({
            "name": entry.name,
            "nodes": invariants.nodes,
            "edges": invariants.edges,
            "cycle_rank": invariants.cycle_rank,
            "coordination": invariants.coordination,
            "identified": index.identify(&graph),
        }));
    }
    if cli.json {
        println!("{}", json!({"entries": entries, "problems": file.problems}));
    } else {
        for item in &entries {
            println!(
                "{}: {} nodes, {} edges, rank {}, RCSR {}",
                item["name"], item["nodes"], item["edges"], item["cycle_rank"], item["identified"]
            );
        }
        for problem in &file.problems {
            println!("problem: {problem}");
        }
    }
    Ok(())
}

fn zeo(cli: &Cli, path: &str, probe: f64, samples: u32) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let result = xtalforge::analysis::run_zeo(
        &structure,
        &ZeoOptions {
            probe_radius: probe,
            samples,
            ..ZeoOptions::default()
        },
    )
    .map_err(CliError::analysis)?;
    if cli.json {
        println!(
            "{}",
            serde_json::to_string_pretty(&result).map_err(json_error)?
        );
    } else {
        println!("pore volume {:.3} A^3", result.pore_volume);
        println!("accessible volume {:.3} A^3", result.accessible_volume);
        println!("surface area {:.3} A^2", result.surface_area);
        println!(
            "largest included sphere {:.3} A, free sphere {:.3} A",
            result.largest_included_sphere, result.largest_free_sphere
        );
    }
    Ok(())
}

fn overlay(
    cli: &Cli,
    path: &str,
    tolerance: f64,
    out: Option<&str>,
    _format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let text = input::read_source(path).map_err(CliError::core)?;
    let file = xtalforge::io::read_cgd(&text);
    let Some(entry) = file.entries.first() else {
        return Err(CliError::input("the CGD file has no readable entry"));
    };
    let graph = quotient_graph(entry, tolerance).map_err(CliError::analysis)?;
    let cell = match entry.cell.as_deref() {
        Some([a, b, c, alpha, beta, gamma, ..]) => {
            let lattice =
                xtalforge::core::Lattice::from_parameters(*a, *b, *c, *alpha, *beta, *gamma)
                    .map_err(CliError::core)?;
            lattice.matrix()
        }
        _ => [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    };
    let overlay = overlay_net(
        &graph,
        cell,
        xtalforge::analysis::OverlayStyle {
            color: [90, 90, 90],
            radius: 0.2,
        },
    );
    let text = serde_json::to_string_pretty(&overlay).map_err(json_error)?;
    if let Some(destination) = out {
        return output::write_destination(destination, &text);
    }
    if cli.json {
        println!("{text}");
    } else {
        println!(
            "{}: {} sites in a {} A cell",
            overlay.name,
            overlay.sites.len(),
            overlay.cell[0][0]
        );
    }
    Ok(())
}

fn json_error(error: serde_json::Error) -> CliError {
    CliError::new("json_error", format!("could not serialize: {error}"))
}
