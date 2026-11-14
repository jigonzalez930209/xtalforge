//! Molecule, sketch, fragment, substitution and guest commands (Task 11.3).

use serde_json::json;

use crate::feature_commands::{load, parse_numbers};
use crate::format::StructureFormat;
use crate::{Cli, CliError, output};

pub(crate) fn molecule(
    cli: &Cli,
    smiles: &str,
    seed: u64,
    iterations: usize,
    hydrogens: bool,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let graph = xtalforge::build::parse_smiles(smiles).map_err(CliError::core)?;
    let options = xtalforge::build::EmbedOptions {
        seed,
        iterations,
        ..xtalforge::build::EmbedOptions::default()
    };
    let structure = if hydrogens {
        xtalforge::build::embed_smiles_with_hydrogens(&graph, &options)
    } else {
        xtalforge::build::embed_smiles(&graph, &options)
    }
    .map_err(CliError::core)?;
    if cli.json && out.is_none() {
        println!(
            "{}",
            json!({
                "formula": structure.formula(),
                "n_atoms": structure.p1().map(|cell| cell.n_atoms()).unwrap_or(0),
                "smiles": smiles,
            })
        );
    }
    output::emit_structure(cli, &structure, out, format)
}

pub(crate) fn fragment(
    cli: &Cli,
    name: Option<&str>,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let Some(name) = name else {
        let fragments = xtalforge::build::fragments().map_err(CliError::core)?;
        if cli.json {
            println!("{}", serde_json::to_string(&fragments).unwrap_or_default());
        } else {
            for fragment in fragments {
                println!(
                    "{:<20} {:<12} {}",
                    fragment.name, fragment.category, fragment.smiles
                );
            }
        }
        return Ok(());
    };
    let fragment = xtalforge::build::find_fragment(name)
        .map_err(CliError::core)?
        .ok_or_else(|| CliError::input(format!("unknown fragment {name:?}")))?;
    molecule(cli, &fragment.smiles, 7, 1500, false, out, format)
}

pub(crate) fn sketch(cli: &Cli, smiles: &str) -> Result<(), CliError> {
    let sketch = xtalforge::build::from_smiles(smiles).map_err(CliError::core)?;
    let roundtrip = sketch.to_smiles().map_err(CliError::core)?;
    if cli.json {
        println!(
            "{}",
            json!({ "smiles": smiles, "roundtrip": roundtrip, "atoms": sketch.atoms.len() })
        );
    } else {
        println!("{roundtrip}");
    }
    Ok(())
}

pub(crate) fn substitute(
    cli: &Cli,
    path: &str,
    index: usize,
    element: &str,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure =
        xtalforge::build::substitute(&load(path)?, index, element).map_err(CliError::core)?;
    output::emit_structure(cli, &structure, out, format)
}

/// Inputs of the `guest` command.
pub(crate) struct GuestRequest<'a> {
    /// Host structure text.
    pub host: &'a str,
    /// Guest SMILES string.
    pub smiles: &'a str,
    /// Cartesian offset `x,y,z`.
    pub offset: &'a str,
    /// Minimum host-guest clearance in Angstrom.
    pub clearance: f64,
    /// Add the hydrogens implied by the guest bond orders.
    pub hydrogens: bool,
}

pub(crate) fn guest(
    cli: &Cli,
    request: &GuestRequest<'_>,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let host = load(request.host)?;
    let graph = xtalforge::build::parse_smiles(request.smiles).map_err(CliError::core)?;
    let options = xtalforge::build::EmbedOptions::default();
    let guest = if request.hydrogens {
        xtalforge::build::embed_smiles_with_hydrogens(&graph, &options)
    } else {
        xtalforge::build::embed_smiles(&graph, &options)
    }
    .map_err(CliError::core)?;
    let points: Vec<(String, [f64; 3])> = guest
        .p1()
        .map_err(CliError::core)?
        .atoms
        .iter()
        .map(|atom| (atom.element.clone(), guest.lattice().to_cart(atom.frac)))
        .collect();
    let values = parse_numbers(request.offset)?;
    let [x, y, z] = match values.as_slice() {
        [x, y, z] => [*x, *y, *z],
        _ => return Err(CliError::input("--offset needs three numbers")),
    };
    let structure = xtalforge::build::place_guest(&host, &points, [x, y, z], request.clearance)
        .map_err(CliError::core)?;
    output::emit_structure(cli, &structure, out, format)
}
