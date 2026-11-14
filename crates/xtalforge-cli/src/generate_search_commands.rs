//! Phase 20 search, MOF dataset and functionalisation commands.

use xtalforge::generation::{
    EvolutionaryOptions, MofDatasetOptions, RandomSearchOptions, SearchKind, evolutionary_search,
    functionalization_sweep, generate_mof_dataset, random_search, search_workflow_json,
};
use xtalforge::mof::FunctionalGroup;

use crate::generate_args::{
    EvolveArgs, FunctionalizeArgs, MofDatasetArgs, SearchArgs, WorkflowArgs,
};
use crate::generate_commands::{composition, csv, fractions, group_name, groups, read_structure};
use crate::{Cli, CliError, output};
pub(crate) fn search(cli: &Cli, args: &SearchArgs) -> Result<(), CliError> {
    let options = RandomSearchOptions {
        composition: composition(&args.composition)?,
        space_groups: groups(&args.space_groups)?,
        n_structures: args.structures,
        seed: args.seed,
        volume_per_atom: 20.0,
        min_distance: 1.0,
        max_attempts: 32,
        relax: !args.no_relax,
        relax_cell: args.cell,
        max_steps: args.steps,
        deduplicate: true,
        matcher: Default::default(),
        minimizer: Default::default(),
    };
    let report = random_search(&options).map_err(CliError::core)?;
    let text = report.to_json().map_err(CliError::core)?;
    if let Some(path) = &args.report {
        output::write_destination(path, &text)?;
    } else if cli.json {
        println!("{text}");
    } else if let Some(best) = report.best() {
        println!(
            "{} hits, best {:.6} eV/atom in space group {} (candidate {})",
            report.hits.len(),
            best.energy_per_atom,
            best.space_group,
            best.candidate
        );
    } else {
        println!("no hits");
    }
    if let Some(best) = report.best()
        && args.out.is_some()
    {
        output::emit_structure(cli, &best.structure, args.out.as_deref(), args.format)?;
    }
    Ok(())
}

pub(crate) fn evolve(cli: &Cli, args: &EvolveArgs) -> Result<(), CliError> {
    let options = EvolutionaryOptions {
        population: args.population,
        generations: args.generations,
        relax_cell: args.cell,
        max_steps: args.steps,
        ..EvolutionaryOptions::new(
            composition(&args.composition)?,
            groups(&args.space_groups)?,
            args.seed,
        )
    };
    let report = evolutionary_search(&options).map_err(CliError::core)?;
    if cli.json {
        println!("{}", report.to_json().map_err(CliError::core)?);
    } else if let Some(best) = report.best() {
        println!(
            "{} generations, best {:.6} eV/atom ({} relaxations, converged {})",
            report.generations.len(),
            best.energy_per_atom,
            report.relaxations,
            report.converged
        );
    }
    if let Some(best) = report.best()
        && args.out.is_some()
    {
        output::emit_structure(cli, &best.structure, args.out.as_deref(), args.format)?;
    }
    Ok(())
}

pub(crate) fn workflow(cli: &Cli, args: &WorkflowArgs) -> Result<(), CliError> {
    let kind = match args.kind.as_str() {
        "random" => SearchKind::Random,
        "evolutionary" => SearchKind::Evolutionary,
        other => {
            return Err(CliError::input(format!(
                "search kind `{other}` is not random or evolutionary"
            )));
        }
    };
    let options: serde_json::Value = serde_json::from_str(&args.options)
        .map_err(|error| CliError::input(format!("search options are not JSON: {error}")))?;
    let text = search_workflow_json(
        &args.name,
        kind,
        &args.seed_path,
        options,
        args.checkpoint.as_deref(),
        &args.output,
    )
    .map_err(CliError::core)?;
    if cli.json {
        println!(
            "{}",
            serde_json::from_str::<serde_json::Value>(&text)
                .map_err(|error| CliError::new("json", error.to_string()))?
        );
    } else {
        print!("{text}");
    }
    Ok(())
}

pub(crate) fn mof_dataset(cli: &Cli, args: &MofDatasetArgs) -> Result<(), CliError> {
    let options = MofDatasetOptions {
        topologies: csv(&args.topologies),
        nodes: csv(&args.nodes),
        linkers: csv(&args.linkers),
        limit: args.limit,
        relax: !args.no_relax,
        max_steps: 100,
        descriptors: !args.no_descriptors,
        pore_spacing: 1.0,
        seed: args.seed,
    };
    let report = generate_mof_dataset(&options).map_err(CliError::core)?;
    let text = report.to_json().map_err(CliError::core)?;
    if let Some(path) = &args.out {
        output::write_destination(path, &text)?;
        if !cli.quiet {
            println!(
                "wrote {} records ({} built, {} failed, {} duplicates) to {path}",
                report.records.len(),
                report.built,
                report.failed,
                report.duplicates
            );
        }
    } else if cli.json {
        println!("{text}");
    } else {
        println!(
            "{} records: {} built, {} failed, {} duplicates (full bundle {})",
            report.records.len(),
            report.built,
            report.failed,
            report.duplicates,
            report.full_bundle
        );
    }
    Ok(())
}

pub(crate) fn functionalize(cli: &Cli, args: &FunctionalizeArgs) -> Result<(), CliError> {
    let structure = read_structure(&args.path)?;
    let groups: Vec<FunctionalGroup> = csv(&args.groups)
        .iter()
        .map(|name| group_name(name))
        .collect::<Result<_, _>>()?;
    let fractions = fractions(&args.fractions)?;
    let sweep = functionalization_sweep(
        &structure,
        &groups,
        &fractions,
        args.seed,
        !args.no_descriptors,
    )
    .map_err(CliError::core)?;
    let text = sweep.to_json().map_err(CliError::core)?;
    if let Some(path) = &args.out {
        output::write_destination(path, &text)?;
        if !cli.quiet {
            println!(
                "wrote {} variants ({} built, {} failed) to {path}",
                sweep.records.len(),
                sweep.built,
                sweep.failed
            );
        }
    } else if cli.json {
        println!("{text}");
    } else {
        println!(
            "{} variants: {} built, {} failed",
            sweep.records.len(),
            sweep.built,
            sweep.failed
        );
    }
    Ok(())
}
