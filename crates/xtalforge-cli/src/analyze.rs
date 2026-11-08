//! Connectivity and group commands: `bonds`, `rings`, `groups`, `neighbors`.

use serde_json::json;
use xtalforge::core::{BondRules, perceive_bonds};

use crate::{Cli, CliError, input};

/// Lists the perceived bonds of the P1 cell.
///
/// # Errors
///
/// Returns [`CliError`] for a read error, a failed bond perception or a JSON
/// serialization failure.
pub(crate) fn bonds(cli: &Cli, path: &str) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let graph = perceive_bonds(&structure, BondRules::default()).map_err(CliError::core)?;
    if cli.json {
        let bonds: Vec<serde_json::Value> = graph
            .bonds
            .iter()
            .map(|bond| {
                json!({
                    "i": bond.i,
                    "j": bond.j,
                    "image": bond.image,
                    "distance": bond.distance,
                    "order": bond.order,
                    "explicit": bond.explicit,
                })
            })
            .collect();
        println!(
            "{}",
            json!({
                "n_atoms": graph.n_atoms(),
                "n_bonds": bonds.len(),
                "bonds": bonds,
            })
        );
        return Ok(());
    }
    if !cli.quiet {
        println!(
            "# {} atoms, {} bonds (i j, image, distance A, order)",
            graph.n_atoms(),
            graph.bonds.len()
        );
    }
    for bond in &graph.bonds {
        let [x, y, z] = bond.image;
        let explicit = if bond.explicit { " explicit" } else { "" };
        println!(
            "{} {} {x} {y} {z} {:.6} {:.3}{explicit}",
            bond.i, bond.j, bond.distance, bond.order
        );
    }
    Ok(())
}

/// Lists the chordless rings up to `max_size` atoms.
///
/// # Errors
///
/// Returns [`CliError`] for a read error, a failed ring search or a JSON
/// serialization failure.
pub(crate) fn rings(cli: &Cli, path: &str, max_size: usize) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let rings = xtalforge::core::rings(&structure, max_size).map_err(CliError::core)?;
    if cli.json {
        let payload: Vec<serde_json::Value> = rings
            .iter()
            .map(|ring| json!({ "size": ring.size(), "atoms": ring.atoms }))
            .collect();
        println!(
            "{}",
            json!({
                "max_size": max_size,
                "n_rings": payload.len(),
                "rings": payload,
            })
        );
        return Ok(());
    }
    if !cli.quiet {
        println!("# {} rings up to size {max_size}", rings.len());
    }
    for ring in &rings {
        println!("{} {}", ring.size(), join(&ring.atoms));
    }
    Ok(())
}

/// Lists the matched functional groups.
///
/// # Errors
///
/// Returns [`CliError`] for a read error, a failed match or a JSON
/// serialization failure.
pub(crate) fn groups(cli: &Cli, path: &str) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let groups = xtalforge::core::find_groups(&structure).map_err(CliError::core)?;
    if cli.json {
        let payload: Vec<serde_json::Value> = groups
            .iter()
            .map(|group| json!({ "kind": group.kind.name(), "atoms": group.atoms }))
            .collect();
        println!(
            "{}",
            json!({
                "n_groups": payload.len(),
                "groups": payload,
            })
        );
        return Ok(());
    }
    if !cli.quiet {
        println!("# {} groups", groups.len());
    }
    for group in &groups {
        println!("{} {}", group.kind.name(), join(&group.atoms));
    }
    Ok(())
}

/// Lists the periodic neighbors within `cutoff` Angstrom.
///
/// # Errors
///
/// Returns [`CliError`] for a read error, an invalid cutoff or a JSON
/// serialization failure.
pub(crate) fn neighbors(cli: &Cli, path: &str, cutoff: f64) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let neighbors = xtalforge::core::neighbors(&structure, cutoff).map_err(CliError::core)?;
    if cli.json {
        let payload: Vec<serde_json::Value> = neighbors
            .iter()
            .map(|neighbor| {
                json!({
                    "i": neighbor.i,
                    "j": neighbor.j,
                    "image": neighbor.image,
                    "distance": neighbor.distance,
                })
            })
            .collect();
        println!(
            "{}",
            json!({
                "cutoff": cutoff,
                "n_neighbors": payload.len(),
                "neighbors": payload,
            })
        );
        return Ok(());
    }
    if !cli.quiet {
        println!(
            "# {} pairs within {cutoff} A (i j, image, distance A)",
            neighbors.len()
        );
    }
    for neighbor in &neighbors {
        let [x, y, z] = neighbor.image;
        println!(
            "{} {} {x} {y} {z} {:.6}",
            neighbor.i, neighbor.j, neighbor.distance
        );
    }
    Ok(())
}

fn join(atoms: &[usize]) -> String {
    atoms
        .iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join(" ")
}
