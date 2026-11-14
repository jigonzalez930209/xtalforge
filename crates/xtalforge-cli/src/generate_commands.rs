//! Phase 20 `generate` commands (Tasks 20.1 to 20.4).

use xtalforge::generation::{
    CrystalOptions, LayerCrystalOptions, MolecularCrystalOptions, PrototypeOptions,
    RodCrystalOptions, SubstitutionOptions, build_prototype, generate_crystal,
    generate_crystal_range, generate_layer_crystal, generate_molecular_crystal,
    generate_rod_crystal, identify_prototype, prototype_library, substitution_report,
};

use xtalforge::mof::FunctionalGroup;

use crate::generate_args::{
    CrystalArgs, CrystalRangeArgs, IdentifyArgs, LayerArgs, MoleculeArgs, PrototypeArgs, RodArgs,
    SubstituteArgs,
};
use crate::{Cli, CliError, format::StructureFormat, input, output};

fn generation_error(error: xtalforge::core::Error) -> CliError {
    CliError::core(error)
}

/// Parses `Element:count,Element:count`.
pub(crate) fn composition(text: &str) -> Result<Vec<(String, u32)>, CliError> {
    let mut out = Vec::new();
    for entry in text.split(',') {
        let entry = entry.trim();
        if entry.is_empty() {
            continue;
        }
        let Some((element, count)) = entry.split_once(':') else {
            return Err(CliError::input(format!(
                "composition entry `{entry}` must be Element:count"
            )));
        };
        let count = count.trim().parse::<u32>().map_err(|error| {
            CliError::input(format!(
                "composition count `{count}` is not an integer: {error}"
            ))
        })?;
        out.push((element.trim().to_owned(), count));
    }
    if out.is_empty() {
        return Err(CliError::input("the composition is empty"));
    }
    Ok(out)
}

/// Parses `1-10` or `225,227`.
pub(crate) fn groups(text: &str) -> Result<Vec<u32>, CliError> {
    let mut out = Vec::new();
    for part in text.split(',') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        if let Some((start, end)) = part.split_once('-') {
            let start = start.trim().parse::<u32>().map_err(|error| {
                CliError::input(format!("space-group range `{part}` is invalid: {error}"))
            })?;
            let end = end.trim().parse::<u32>().map_err(|error| {
                CliError::input(format!("space-group range `{part}` is invalid: {error}"))
            })?;
            if start > end {
                return Err(CliError::input(format!(
                    "space-group range `{part}` starts after it ends"
                )));
            }
            out.extend(start..=end);
        } else {
            out.push(part.parse::<u32>().map_err(|error| {
                CliError::input(format!("space-group `{part}` is not a number: {error}"))
            })?);
        }
    }
    if out.is_empty() {
        return Err(CliError::input("no space groups were given"));
    }
    Ok(out)
}

pub(crate) fn csv(text: &str) -> Vec<String> {
    text.split(',')
        .map(str::trim)
        .filter(|entry| !entry.is_empty())
        .map(str::to_owned)
        .collect()
}

pub(crate) fn fractions(text: &str) -> Result<Vec<f64>, CliError> {
    let values = csv(text);
    let mut out = Vec::with_capacity(values.len());
    for value in values {
        out.push(value.parse::<f64>().map_err(|error| {
            CliError::input(format!("fraction `{value}` is not a number: {error}"))
        })?);
    }
    if out.is_empty() {
        return Err(CliError::input("no fractions were given"));
    }
    Ok(out)
}

pub(crate) fn group_name(text: &str) -> Result<FunctionalGroup, CliError> {
    match text.trim().to_ascii_lowercase().as_str() {
        "amino" | "nh2" => Ok(FunctionalGroup::Amino),
        "hydroxyl" | "oh" => Ok(FunctionalGroup::Hydroxyl),
        "methyl" | "ch3" => Ok(FunctionalGroup::Methyl),
        "fluorine" | "f" => Ok(FunctionalGroup::Fluorine),
        other => Err(CliError::input(format!(
            "unknown functional group `{other}`; use amino, hydroxyl, methyl or fluorine"
        ))),
    }
}

pub(crate) fn read_structure(path: &str) -> Result<xtalforge::core::Structure, CliError> {
    let text = input::read_source(path).map_err(CliError::core)?;
    input::structure_from_text(&text, path).map_err(CliError::core)
}

fn emit_crystal(
    cli: &Cli,
    structure: &xtalforge::core::Structure,
    out: Option<&str>,
    format: Option<StructureFormat>,
    summary: &str,
) -> Result<(), CliError> {
    if let Some(path) = out {
        output::emit_structure(cli, structure, Some(path), format)?;
        if !cli.quiet {
            println!("{summary}");
        }
        return Ok(());
    }
    if cli.json {
        println!(
            "{}",
            structure.to_document().to_json().map_err(CliError::core)?
        );
    } else {
        println!("{summary}");
    }
    Ok(())
}

pub(crate) fn run(
    cli: &Cli,
    action: &crate::generate_args::GenerateAction,
) -> Result<(), CliError> {
    use crate::generate_args::GenerateAction;
    match action {
        GenerateAction::Crystal(args) => crystal(cli, args),
        GenerateAction::CrystalRange(args) => crystal_range(cli, args),
        GenerateAction::Molecule(args) => molecule(cli, args),
        GenerateAction::LayerGroups => layer_group_list(cli),
        GenerateAction::Layer(args) => layer(cli, args),
        GenerateAction::RodGroups => rod_group_list(cli),
        GenerateAction::Rod(args) => rod(cli, args),
        GenerateAction::Prototypes => prototypes(cli),
        GenerateAction::Prototype(args) => prototype(cli, args),
        GenerateAction::Identify(args) => identify(cli, args),
        GenerateAction::Substitute(args) => substitute(cli, args),
        GenerateAction::Search(args) => crate::generate_search_commands::search(cli, args),
        GenerateAction::Evolve(args) => crate::generate_search_commands::evolve(cli, args),
        GenerateAction::Workflow(args) => crate::generate_search_commands::workflow(cli, args),
        GenerateAction::MofDataset(args) => crate::generate_search_commands::mof_dataset(cli, args),
        GenerateAction::Functionalize(args) => {
            crate::generate_search_commands::functionalize(cli, args)
        }
    }
}

fn crystal(cli: &Cli, args: &CrystalArgs) -> Result<(), CliError> {
    let options = CrystalOptions {
        space_group: args.space_group,
        composition: composition(&args.composition)?,
        volume_per_atom: args.volume_per_atom,
        min_distance: args.min_distance,
        max_attempts: args.attempts,
        seed: args.seed,
    };
    let crystal = generate_crystal(&options).map_err(generation_error)?;
    if cli.json && args.out.is_none() {
        println!("{}", crystal.to_json().map_err(CliError::core)?);
        return Ok(());
    }
    emit_crystal(
        cli,
        &crystal.structure,
        args.out.as_deref(),
        args.format,
        &crystal.summary(),
    )
}

fn crystal_range(_cli: &Cli, args: &CrystalRangeArgs) -> Result<(), CliError> {
    let options = CrystalOptions {
        space_group: 1,
        composition: composition(&args.composition)?,
        volume_per_atom: args.volume_per_atom,
        min_distance: args.min_distance,
        max_attempts: args.attempts,
        seed: args.seed,
    };
    let groups = groups(&args.groups)?;
    let results = generate_crystal_range(&groups, &options).map_err(generation_error)?;
    let mut reports = Vec::new();
    for (group, result) in groups.iter().zip(results) {
        match result {
            Ok(crystal) => reports.push(serde_json::json!({
                "space_group": group,
                "status": "ok",
                "attempt": crystal.attempt,
                "wyckoffs": crystal.wyckoffs,
                "structure": crystal.structure,
            })),
            Err(error) => reports.push(serde_json::json!({
                "space_group": group,
                "status": "failed",
                "error": error.message(),
            })),
        }
    }
    println!(
        "{}",
        serde_json::to_string_pretty(&reports)
            .map_err(|error| CliError::new("json", error.to_string()))?
    );
    Ok(())
}

fn molecule(cli: &Cli, args: &MoleculeArgs) -> Result<(), CliError> {
    let molecule = read_structure(&args.molecule)?;
    let options = MolecularCrystalOptions {
        space_group: args.space_group,
        molecule,
        copies: args.copies,
        volume_per_atom: args.volume_per_atom,
        min_distance: args.min_distance,
        tolerance: 0.25,
        max_attempts: args.attempts,
        seed: args.seed,
    };
    let crystal = generate_molecular_crystal(&options).map_err(generation_error)?;
    let summary = format!(
        "molecular crystal in {} with {} molecules per cell (seed {}, attempt {})",
        crystal.space_group, crystal.molecules_per_cell, crystal.seed, crystal.attempt
    );
    emit_crystal(
        cli,
        &crystal.structure,
        args.out.as_deref(),
        args.format,
        &summary,
    )
}

fn layer(cli: &Cli, args: &LayerArgs) -> Result<(), CliError> {
    let options = LayerCrystalOptions {
        layer_group: args.layer_group,
        composition: composition(&args.composition)?,
        area_per_atom: args.area_per_atom,
        min_distance: args.min_distance,
        max_attempts: args.attempts,
        seed: args.seed,
    };
    let crystal = generate_layer_crystal(&options).map_err(generation_error)?;
    let summary = format!(
        "layer group {} ({}) with Wyckoff {} (seed {}, attempt {})",
        crystal.layer_group,
        crystal.hm_symbol,
        crystal.wyckoffs.join(" "),
        crystal.seed,
        crystal.attempt
    );
    emit_crystal(
        cli,
        &crystal.structure,
        args.out.as_deref(),
        args.format,
        &summary,
    )
}

fn rod(cli: &Cli, args: &RodArgs) -> Result<(), CliError> {
    let options = RodCrystalOptions {
        rod_group: args.rod_group,
        composition: composition(&args.composition)?,
        repeat: args.repeat,
        min_distance: args.min_distance,
        max_attempts: args.attempts,
        seed: args.seed,
    };
    let crystal = generate_rod_crystal(&options).map_err(generation_error)?;
    let summary = format!(
        "rod group {} with repeat {} A (seed {}, attempt {})",
        crystal.hm_symbol, args.repeat, crystal.seed, crystal.attempt
    );
    emit_crystal(
        cli,
        &crystal.structure,
        args.out.as_deref(),
        args.format,
        &summary,
    )
}

fn layer_group_list(cli: &Cli) -> Result<(), CliError> {
    let groups = xtalforge::generation::layer_groups().map_err(CliError::core)?;
    if cli.json {
        println!(
            "{}",
            serde_json::to_string_pretty(&groups)
                .map_err(|error| CliError::new("json", error.to_string()))?
        );
    } else {
        for group in groups {
            println!("{} {} ({})", group.number, group.hm_symbol, group.bravais);
        }
    }
    Ok(())
}

fn rod_group_list(cli: &Cli) -> Result<(), CliError> {
    let groups = xtalforge::generation::rod_groups().map_err(CliError::core)?;
    if cli.json {
        println!(
            "{}",
            serde_json::to_string_pretty(&groups)
                .map_err(|error| CliError::new("json", error.to_string()))?
        );
    } else {
        for group in groups {
            println!("{} {}", group.number, group.symbol);
        }
    }
    Ok(())
}

fn prototypes(cli: &Cli) -> Result<(), CliError> {
    let library = prototype_library().map_err(CliError::core)?;
    if cli.json {
        println!(
            "{}",
            serde_json::to_string_pretty(&library)
                .map_err(|error| CliError::new("json", error.to_string()))?
        );
    } else {
        for prototype in library {
            println!(
                "{} {} (space group {}, {} sites)",
                prototype.label,
                prototype.name,
                prototype.space_group,
                prototype.sites.len()
            );
        }
    }
    Ok(())
}

fn prototype(cli: &Cli, args: &PrototypeArgs) -> Result<(), CliError> {
    let structure = build_prototype(
        &args.label,
        &PrototypeOptions {
            scale: args.scale,
            ..PrototypeOptions::default()
        },
    )
    .map_err(CliError::core)?;
    let summary = format!("built prototype {}", args.label);
    emit_crystal(cli, &structure, args.out.as_deref(), args.format, &summary)
}

fn identify(cli: &Cli, args: &IdentifyArgs) -> Result<(), CliError> {
    let structure = read_structure(&args.path)?;
    let identification = identify_prototype(
        &structure,
        &PrototypeOptions {
            matcher: xtalforge::symmetry::MatcherOptions {
                stol: args.stol,
                ..Default::default()
            },
            ..PrototypeOptions::default()
        },
    )
    .map_err(CliError::core)?;
    if cli.json {
        println!(
            "{}",
            serde_json::to_string_pretty(&identification)
                .map_err(|error| CliError::new("json", error.to_string()))?
        );
    } else if let Some(label) = &identification.label {
        println!(
            "{} {} (space group {}, signature {}, RMS {:.4} A)",
            label,
            identification.name.as_deref().unwrap_or(""),
            identification.space_group,
            identification.signature,
            identification.rms.unwrap_or(f64::NAN)
        );
    } else {
        println!(
            "no prototype match (space group {}, signature {})",
            identification.space_group, identification.signature
        );
    }
    Ok(())
}

fn substitute(_cli: &Cli, args: &SubstituteArgs) -> Result<(), CliError> {
    let structure = read_structure(&args.path)?;
    let report = substitution_report(
        &structure,
        &SubstitutionOptions {
            radius_tolerance: args.radius_tolerance,
            same_charge: !args.allow_charge_change,
            limit: args.limit,
        },
    )
    .map_err(CliError::core)?;
    println!("{}", report.to_json().map_err(CliError::core)?);
    Ok(())
}
