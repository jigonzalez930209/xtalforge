//! The `supercell` and `slab` commands.

use xtalforge::core::limits::check_supercell;
use xtalforge::core::{Limits as CoreLimits, Structure, make_slab};

use crate::format::StructureFormat;
use crate::{Cli, CliError, input, output};

/// Builds an `na x nb x nc` supercell and emits it.
///
/// # Errors
///
/// Returns [`CliError`] for a read error, a refused size, a failed supercell
/// or an output failure.
pub(crate) fn supercell(
    cli: &Cli,
    path: &str,
    repeats: [u32; 3],
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    check_size(cli, &structure, repeats)?;
    let built = structure
        .supercell_with_limits(
            repeats[0],
            repeats[1],
            repeats[2],
            &CoreLimits::new(cli.limits.profile()),
        )
        .map_err(CliError::core)?;
    output::emit_structure(cli, &built, out, format)
}

/// Cuts a slab and emits it.
///
/// # Errors
///
/// Returns [`CliError`] for a read error, a degenerate plane, a failed cut or
/// an output failure.
pub(crate) fn slab(
    cli: &Cli,
    path: &str,
    hkl: [i32; 3],
    layers: u32,
    vacuum: f64,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let slab = make_slab(&structure, hkl, layers, vacuum, 0.0).map_err(CliError::core)?;
    output::emit_structure(cli, &slab.structure, out, format)
}

fn check_size(cli: &Cli, structure: &Structure, repeats: [u32; 3]) -> Result<(), CliError> {
    let [na, nb, nc] = repeats;
    let matrix = [
        [f64::from(na), 0.0, 0.0],
        [0.0, f64::from(nb), 0.0],
        [0.0, 0.0, f64::from(nc)],
    ];
    let verdict = check_supercell(
        &CoreLimits::new(cli.limits.profile()),
        structure.p1().map_err(CliError::core)?.n_atoms() as u64,
        matrix,
    )
    .map_err(CliError::core)?;
    if verdict.refused() {
        return Err(CliError::new("limit_exceeded", verdict.to_string()));
    }
    if verdict.warned() && !cli.quiet {
        eprintln!("warning: {verdict}");
    }
    Ok(())
}
