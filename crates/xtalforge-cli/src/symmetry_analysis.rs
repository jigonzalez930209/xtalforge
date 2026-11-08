//! Symmetry analysis and construction commands (Task 2.5-2.12).

use xtalforge::symmetry::{SpaceGroup, SymmetryOptions};

use crate::format::StructureFormat;
use crate::{Cli, CliError, input, output};

/// Standardizes a structure to the conventional cell.
///
/// # Errors
///
/// Returns [`CliError`] for a read or detection failure.
pub(crate) fn standardize(
    cli: &Cli,
    path: &str,
    symprec: f64,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let standardized = xtalforge::symmetry::standardize(
        &structure,
        SymmetryOptions {
            symprec,
            ..SymmetryOptions::default()
        },
    )
    .map_err(CliError::core)?;
    output::emit_structure(cli, &standardized.structure, out, format)
}

/// Builds the standardized primitive cell.
///
/// # Errors
///
/// Returns [`CliError`] for a read or detection failure.
pub(crate) fn primitive(
    cli: &Cli,
    path: &str,
    symprec: f64,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let primitive = xtalforge::symmetry::primitive_cell(
        &structure,
        SymmetryOptions {
            symprec,
            ..SymmetryOptions::default()
        },
    )
    .map_err(CliError::core)?;
    output::emit_structure(cli, &primitive.structure, out, format)
}

/// Reduces a structure to a Niggli, Delaunay or Minkowski cell.
///
/// # Errors
///
/// Returns [`CliError`] for a read or reduction failure.
pub(crate) fn reduce(
    cli: &Cli,
    path: &str,
    algorithm: &str,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let reduction = match algorithm {
        "niggli" => xtalforge::symmetry::Reduction::Niggli,
        "delaunay" => xtalforge::symmetry::Reduction::Delaunay,
        "minkowski" => xtalforge::symmetry::Reduction::Minkowski,
        other => {
            return Err(CliError::input(format!(
                "unknown reduction {other:?}; use niggli, delaunay or minkowski"
            )));
        }
    };
    let reduced =
        xtalforge::symmetry::reduce_structure(&structure, reduction).map_err(CliError::core)?;
    output::emit_structure(cli, &reduced.structure, out, format)
}

/// Prints Wyckoff assignments for every atom.
///
/// # Errors
///
/// Returns [`CliError`] for a read or detection failure.
pub(crate) fn wyckoff(cli: &Cli, path: &str, symprec: f64) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let assignments = xtalforge::symmetry::wyckoff_assignments(
        &structure,
        SymmetryOptions {
            symprec,
            ..SymmetryOptions::default()
        },
    )
    .map_err(CliError::core)?;
    if cli.json {
        println!(
            "{}",
            serde_json::to_string(&assignments).map_err(|error| {
                CliError::new(
                    "json_error",
                    format!("could not serialize Wyckoff data: {error}"),
                )
            })?
        );
    } else {
        for assignment in assignments {
            println!(
                "atom {} {}: Wyckoff {} ({}), site symmetry {}, template {}",
                assignment.atom,
                assignment.element,
                assignment.letter,
                assignment.multiplicity,
                assignment.site_symmetry,
                assignment.coordinates
            );
        }
    }
    Ok(())
}

/// Prints the reciprocal-space Brillouin zone.
///
/// # Errors
///
/// Returns [`CliError`] for a read or zone failure.
pub(crate) fn brillouin(cli: &Cli, path: &str) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let zone = xtalforge::symmetry::brillouin_zone(structure.lattice()).map_err(CliError::core)?;
    if cli.json {
        println!(
            "{}",
            serde_json::to_string(&zone).map_err(|error| {
                CliError::new(
                    "json_error",
                    format!("could not serialize the zone: {error}"),
                )
            })?
        );
    } else {
        println!(
            "Brillouin zone: {} faces, {} vertices, {} edges",
            zone.n_faces(),
            zone.n_vertices(),
            zone.edges.len()
        );
    }
    Ok(())
}

/// Prints the tolerance ladder of a structure.
///
/// # Errors
///
/// Returns [`CliError`] for a read or detection failure.
pub(crate) fn scan(cli: &Cli, path: &str) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let ladder =
        xtalforge::symmetry::tolerance_ladder(&structure, &xtalforge::symmetry::DEFAULT_LADDER)
            .map_err(CliError::core)?;
    if cli.json {
        println!(
            "{}",
            serde_json::to_string(&ladder).map_err(|error| {
                CliError::new(
                    "json_error",
                    format!("could not serialize the ladder: {error}"),
                )
            })?
        );
    } else {
        for range in ladder {
            println!(
                "{:.5}..{:.5} A: {} ({})",
                range.min, range.max, range.hm_symbol, range.number
            );
        }
    }
    Ok(())
}

/// Prints per-atom deviations from a requested group.
///
/// # Errors
///
/// Returns [`CliError`] for a read, parse or detection failure.
pub(crate) fn diagnose(cli: &Cli, path: &str, symbol: &str) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let group = SpaceGroup::parse(symbol).map_err(CliError::core)?;
    let diagnostics = xtalforge::symmetry::diagnose(&structure, &group).map_err(CliError::core)?;
    if cli.json {
        println!(
            "{}",
            serde_json::to_string(&diagnostics).map_err(|error| {
                CliError::new(
                    "json_error",
                    format!("could not serialize diagnostics: {error}"),
                )
            })?
        );
    } else {
        println!(
            "max {:.6} A, mean {:.6} A",
            diagnostics.max_displacement, diagnostics.mean_displacement
        );
        for atom in diagnostics.atoms.iter().take(20) {
            println!(
                "  atom {} {} {}: {:.6} A",
                atom.atom, atom.element, atom.label, atom.displacement
            );
        }
    }
    Ok(())
}

/// Detects the layer group of a slab.
///
/// # Errors
///
/// Returns [`CliError`] for a read or layer-search failure.
pub(crate) fn layer(cli: &Cli, path: &str, symprec: f64) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let dataset =
        xtalforge::symmetry::detect_layer_symmetry(&structure, symprec).map_err(CliError::core)?;
    if cli.json {
        println!(
            "{}",
            serde_json::to_string(&dataset).map_err(|error| {
                CliError::new(
                    "json_error",
                    format!("could not serialize the layer group: {error}"),
                )
            })?
        );
    } else {
        println!(
            "{} (layer group {}, Hall {}) operations {}",
            dataset.hm_symbol,
            dataset.number,
            dataset.hall_number,
            dataset.operations.len()
        );
    }
    Ok(())
}

/// Detects the collinear magnetic space group.
///
/// # Errors
///
/// Returns [`CliError`] for a read or magnetic-search failure.
pub(crate) fn magnetic(cli: &Cli, path: &str, symprec: f64) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let dataset = xtalforge::symmetry::detect_magnetic_symmetry(&structure, symprec)
        .map_err(CliError::core)?;
    if cli.json {
        println!(
            "{}",
            serde_json::to_string(&dataset).map_err(|error| {
                CliError::new(
                    "json_error",
                    format!("could not serialize the magnetic group: {error}"),
                )
            })?
        );
    } else {
        println!(
            "UNI {} with {} magnetic operations",
            dataset.uni_number,
            dataset.operations.len()
        );
    }
    Ok(())
}

/// Changes the setting of a structure to a Hall number.
///
/// # Errors
///
/// Returns [`CliError`] for a read or setting change failure.
pub(crate) fn setting(
    cli: &Cli,
    path: &str,
    hall_number: u32,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let changed =
        xtalforge::symmetry::change_setting(&structure, hall_number).map_err(CliError::core)?;
    output::emit_structure(cli, &changed.structure, out, format)
}

/// Prints the band path of a structure.
///
/// # Errors
///
/// Returns [`CliError`] for a read or k-path failure.
pub(crate) fn kpath(cli: &Cli, path: &str) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let path = xtalforge::symmetry::band_path_for_structure(
        &structure,
        SymmetryOptions {
            symprec: 1e-3,
            ..SymmetryOptions::default()
        },
    )
    .map_err(CliError::core)?;
    let vertices = path.vertices(structure.lattice()).map_err(CliError::core)?;
    if cli.json {
        println!(
            "{}",
            serde_json::json!({
                "bravais": path.bravais,
                "special_points": path.special_points,
                "path": path.path,
                "vertices": vertices.iter().map(|(label, point, distance)| serde_json::json!({
                    "label": label,
                    "cartesian": point,
                    "distance": distance,
                })).collect::<Vec<_>>(),
            })
        );
    } else {
        println!("Bravais {}: {} segments", path.bravais, path.path.len());
        for (start, end) in &path.path {
            println!("  {start} -> {end}");
        }
    }
    Ok(())
}

/// Detects the molecular point group of a finite structure.
///
/// # Errors
///
/// Returns [`CliError`] for a read or point-group failure.
pub(crate) fn pointgroup(cli: &Cli, path: &str, tolerance: f64) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let p1 = structure.p1().map_err(CliError::core)?;
    let points: Vec<(String, [f64; 3])> = p1
        .atoms
        .iter()
        .map(|atom| (atom.element.clone(), p1.lattice.to_cart(atom.frac)))
        .collect();
    let group =
        xtalforge::symmetry::molecular_point_group(&points, tolerance).map_err(CliError::core)?;
    if cli.json {
        println!(
            "{}",
            serde_json::to_string(&group).map_err(|error| {
                CliError::new(
                    "json_error",
                    format!("could not serialize the point group: {error}"),
                )
            })?
        );
    } else {
        println!(
            "{} (order {}, {} proper, {} improper, {} mirrors)",
            group.schoenflies,
            group.order,
            group.proper_rotations,
            group.improper_operations,
            group.mirrors
        );
    }
    Ok(())
}
