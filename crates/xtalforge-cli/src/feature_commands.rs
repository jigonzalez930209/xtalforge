//! Phase 9 commands: comparison, deformation and the feature dispatcher.

use xtalforge::analysis;
use xtalforge::core::{SiteSort, SpeciesReplacement, Structure, replace_species, sort_sites};

use crate::feature_args::FeatureAction;
use crate::output;
use crate::{Cli, CliError, input};

/// Reads a structure from a path or standard input.
pub(crate) fn load(path: &str) -> Result<Structure, CliError> {
    let text = input::read_source(path).map_err(CliError::core)?;
    input::structure_from_text(&text, path).map_err(CliError::core)
}

/// Maps a symmetry or analysis (core) error to the CLI.
pub(crate) fn core_error(error: xtalforge::core::Error) -> CliError {
    CliError::core(error)
}

/// Maps an analysis error to the CLI.
pub(crate) fn analysis_error(error: analysis::Error) -> CliError {
    CliError::core(error)
}

/// Maps a force-field error to the CLI.
pub(crate) fn ff_error(error: xtalforge::ff::Error) -> CliError {
    CliError::new(error.code(), error.message())
}

/// Writes one structure and optionally prints its summary.
pub(crate) fn emit_structure(
    cli: &Cli,
    structure: &Structure,
    out: Option<&str>,
    format: Option<crate::format::StructureFormat>,
) -> Result<(), CliError> {
    output::emit_structure(cli, structure, out, format)
}

/// Writes all frames as one XYZ trajectory.
pub(crate) fn write_frames(
    cli: &Cli,
    frames: &[Structure],
    out: Option<&str>,
) -> Result<(), CliError> {
    if let Some(path) = out.filter(|path| path.ends_with(".xyz") || path.ends_with(".xdatcar")) {
        let mut text = String::new();
        for frame in frames {
            text.push_str(&xtalforge::io::write_structure(frame, "xyz").map_err(CliError::core)?);
        }
        output::write_destination(path, &text)?;
        if !cli.quiet {
            println!("wrote {} frames to {path}", frames.len());
        }
        return Ok(());
    }
    if let Some(first) = frames.first() {
        emit_structure(cli, first, out, None)?;
    }
    Ok(())
}

/// Decides whether two structures are the same.
pub(crate) fn compare(
    cli: &Cli,
    first: &str,
    second: &str,
    ltol: f64,
    stol: f64,
    angle_tol: f64,
    allow_supercell: bool,
) -> Result<(), CliError> {
    let first = load(first)?;
    let second = load(second)?;
    let options = xtalforge::symmetry::MatcherOptions {
        ltol,
        stol,
        angle_tol,
        allow_supercell,
        ..xtalforge::symmetry::MatcherOptions::default()
    };
    let outcome =
        xtalforge::symmetry::match_structures(&first, &second, &options).map_err(core_error)?;
    if cli.json {
        println!(
            "{}",
            serde_json::json!({
                "matched": outcome.matched,
                "rms": outcome.rms,
                "n_atoms": outcome.n_atoms,
                "mapping": outcome.mapping,
                "translation": outcome.translation,
                "lattice_mapping": outcome.lattice_mapping,
                "supercell": outcome.supercell,
            })
        );
    } else if outcome.matched {
        println!(
            "match: yes rms={:.6} A atoms={}",
            outcome.rms, outcome.n_atoms
        );
    } else {
        println!("match: no");
    }
    Ok(())
}

/// Groups structures into equality classes.
pub(crate) fn group(cli: &Cli, paths: &[String], ltol: f64, stol: f64) -> Result<(), CliError> {
    let mut structures = Vec::with_capacity(paths.len());
    for path in paths {
        structures.push(load(path)?);
    }
    let options = xtalforge::symmetry::MatcherOptions {
        ltol,
        stol,
        allow_supercell: true,
        ..xtalforge::symmetry::MatcherOptions::default()
    };
    let groups =
        xtalforge::symmetry::group_structures(&structures, &options).map_err(core_error)?;
    if cli.json {
        println!(
            "{}",
            serde_json::json!({
                "n_structures": structures.len(),
                "n_groups": groups.len(),
                "groups": groups,
            })
        );
    } else {
        println!("{} structures in {} groups", structures.len(), groups.len());
        for (index, members) in groups.iter().enumerate() {
            let names: Vec<&str> = members
                .iter()
                .filter_map(|member| paths.get(*member).map(String::as_str))
                .collect();
            println!("  {index}: {}", names.join(", "));
        }
    }
    Ok(())
}

/// Applies a strain, deformation, volume or density change.
#[allow(clippy::too_many_arguments)] // one flag per deformation mode
pub(crate) fn deform(
    cli: &Cli,
    path: &str,
    voigt: Option<&str>,
    gradient: Option<&str>,
    volume: Option<f64>,
    density: Option<f64>,
    out: Option<&str>,
    format: Option<crate::format::StructureFormat>,
) -> Result<(), CliError> {
    let structure = load(path)?;
    let mut result = structure.clone();
    if let Some(text) = voigt {
        let values = parse_numbers(text)?;
        let [xx, yy, zz, yz, xz, xy] = match values.as_slice() {
            [xx, yy, zz, yz, xz, xy] => [*xx, *yy, *zz, *yz, *xz, *xy],
            _ => {
                return Err(CliError::input(
                    "--voigt needs six components: xx,yy,zz,yz,xz,xy",
                ));
            }
        };
        let tensor =
            xtalforge::core::voigt_to_tensor([xx, yy, zz, yz, xz, xy]).map_err(core_error)?;
        result = xtalforge::core::strain(&result, tensor).map_err(core_error)?;
    }
    if let Some(text) = gradient {
        let values = parse_numbers(text)?;
        let [f11, f12, f13, f21, f22, f23, f31, f32, f33] = match values.as_slice() {
            [f11, f12, f13, f21, f22, f23, f31, f32, f33] => {
                [*f11, *f12, *f13, *f21, *f22, *f23, *f31, *f32, *f33]
            }
            _ => {
                return Err(CliError::input(
                    "--gradient needs nine components: f11,...,f33",
                ));
            }
        };
        let matrix = [[f11, f12, f13], [f21, f22, f23], [f31, f32, f33]];
        result = xtalforge::core::deform(&result, matrix).map_err(core_error)?;
    }
    if let Some(volume) = volume {
        result = xtalforge::core::scale_to_volume(&result, volume).map_err(core_error)?;
    }
    if let Some(density) = density {
        result = xtalforge::core::scale_to_density(&result, density).map_err(core_error)?;
    }
    if voigt.is_none() && gradient.is_none() && volume.is_none() && density.is_none() {
        return Err(CliError::input(
            "deform needs --voigt, --gradient, --volume or --density",
        ));
    }
    emit_structure(cli, &result, out, format)
}

/// Runs one Phase 9 feature action.
pub(crate) fn features(cli: &Cli, action: &FeatureAction) -> Result<(), CliError> {
    match action {
        FeatureAction::Rattle { .. }
        | FeatureAction::Edit { .. }
        | FeatureAction::Interpolate { .. }
        | FeatureAction::OptimalSupercell { .. } => crate::feature_edits::run(cli, action),
        FeatureAction::Voronoi { .. }
        | FeatureAction::Coordination { .. }
        | FeatureAction::Bvs { .. }
        | FeatureAction::Oxidation { .. }
        | FeatureAction::Rdf { .. }
        | FeatureAction::Fingerprint { .. } => crate::feature_local::run(cli, action),
        FeatureAction::Miller { .. }
        | FeatureAction::Terminations { .. }
        | FeatureAction::Wulff { .. }
        | FeatureAction::Adsorption { .. } => crate::feature_surfaces::run(cli, action),
        FeatureAction::Vacancies { .. }
        | FeatureAction::Substitutions { .. }
        | FeatureAction::Interstitials { .. }
        | FeatureAction::Orderings { .. }
        | FeatureAction::Sqs { .. } => crate::feature_defects::run(cli, action),
        FeatureAction::LatticeMatch { .. }
        | FeatureAction::Stack { .. }
        | FeatureAction::Elastic { .. }
        | FeatureAction::Eos { .. }
        | FeatureAction::Hull { .. } => crate::feature_mechanics::run(cli, action),
        FeatureAction::MagneticOrderings { .. }
        | FeatureAction::Disorder { .. }
        | FeatureAction::Derivatives { .. } => crate::phase19_commands::feature(cli, action),
    }
}

/// Applies the site-editing flags in a deterministic order.
#[allow(clippy::too_many_arguments)] // one flag per editing operation
pub(crate) fn edit(
    cli: &Cli,
    path: &str,
    add: &[String],
    remove: &[usize],
    replace: &[String],
    move_site: &[String],
    species: &[String],
    fraction: f64,
    seed: u64,
    sort: Option<&str>,
    wrap: bool,
    out: Option<&str>,
    format: Option<crate::format::StructureFormat>,
) -> Result<(), CliError> {
    let mut structure = load(path)?;
    if !species.is_empty() {
        let mut map = std::collections::BTreeMap::new();
        for entry in species {
            let (from, to) = entry
                .split_once(':')
                .ok_or_else(|| CliError::input("--species needs old:new pairs"))?;
            map.insert(from.trim().to_owned(), to.trim().to_owned());
        }
        structure = replace_species(
            &structure,
            &SpeciesReplacement {
                map,
                fraction,
                seed,
            },
        )
        .map_err(core_error)?;
    }
    if !add.is_empty() {
        let mut sites = Vec::new();
        for entry in add {
            let (element, coordinates) = entry
                .split_once(':')
                .ok_or_else(|| CliError::input("--add needs element:x,y,z"))?;
            let values = parse_triple(coordinates)?;
            sites.push((element.trim().to_owned(), values));
        }
        structure = xtalforge::core::add_sites(&structure, &sites).map_err(core_error)?;
    }
    if !remove.is_empty() {
        structure = xtalforge::core::remove_sites(&structure, remove).map_err(core_error)?;
    }
    for entry in replace {
        let (index, element) = entry
            .split_once(':')
            .ok_or_else(|| CliError::input("--replace needs index:element"))?;
        let index: usize = index
            .trim()
            .parse()
            .map_err(|_| CliError::input("--replace index must be an integer"))?;
        structure =
            xtalforge::core::replace_site(&structure, index, element.trim()).map_err(core_error)?;
    }
    for entry in move_site {
        let (index, coordinates) = entry
            .split_once(':')
            .ok_or_else(|| CliError::input("--move needs index:x,y,z"))?;
        let index: usize = index
            .trim()
            .parse()
            .map_err(|_| CliError::input("--move index must be an integer"))?;
        let frac = parse_triple(coordinates)?;
        structure = xtalforge::core::move_site(&structure, index, frac).map_err(core_error)?;
    }
    if let Some(sort) = sort {
        let key = match sort {
            "element" => SiteSort::Element,
            "electronegativity" | "en" => SiteSort::Electronegativity,
            "mass" => SiteSort::Mass,
            "label" => SiteSort::Label,
            _ => {
                return Err(CliError::input(
                    "--sort must be element, electronegativity, mass or label",
                ));
            }
        };
        structure = sort_sites(&structure, key).map_err(core_error)?;
    }
    if wrap {
        structure = xtalforge::core::wrap_sites(&structure).map_err(core_error)?;
    }
    emit_structure(cli, &structure, out, format)
}

/// Parses comma-separated numbers.
pub(crate) fn parse_numbers(text: &str) -> Result<Vec<f64>, CliError> {
    text.split(',')
        .map(|value| {
            value
                .trim()
                .parse::<f64>()
                .map_err(|_| CliError::input(format!("invalid number '{value}'")))
        })
        .collect()
}

/// Parses three comma-separated numbers.
pub(crate) fn parse_triple(text: &str) -> Result<[f64; 3], CliError> {
    let values: Vec<f64> = text
        .split(',')
        .map(|value| {
            value
                .trim()
                .parse::<f64>()
                .map_err(|_| CliError::input(format!("invalid number '{value}'")))
        })
        .collect::<Result<_, _>>()?;
    match values.as_slice() {
        [x, y, z] => Ok([*x, *y, *z]),
        _ => Err(CliError::input("expected three comma-separated numbers")),
    }
}
