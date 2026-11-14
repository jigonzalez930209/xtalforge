//! Phase 21 `ml` commands.

use xtalforge::analysis::dataset::Dataset;
use xtalforge::ml::{
    AcsfOptions, GraphOptions, MbtrOptions, MofFeaturesOptions, Predictor, PropertyModel,
    SelectionOptions, SelectionStrategy, SplitOptions, SplitStrategy, acsf,
    active_learning_workflow, append_predictions, crystal_graph, import_dataset, mbtr,
    mof_features, predict_dataset, select_batch, soap_power_spectrum, split_dataset,
};

use crate::ml_args::*;
use crate::{Cli, CliError, input, output};

fn structure(path: &str) -> Result<xtalforge::core::Structure, CliError> {
    let text = input::read_source(path).map_err(CliError::core)?;
    input::structure_from_text(&text, path).map_err(CliError::core)
}

fn emit_json(cli: &Cli, value: &impl serde::Serialize, out: Option<&str>) -> Result<(), CliError> {
    let text = serde_json::to_string_pretty(value)
        .map_err(|error| CliError::new("json", error.to_string()))?;
    if let Some(path) = out {
        output::write_destination(path, &text)?;
        if !cli.quiet {
            println!("wrote {path}");
        }
    } else {
        println!("{text}");
    }
    Ok(())
}

pub(crate) fn run(cli: &Cli, action: &MlAction) -> Result<(), CliError> {
    match action {
        MlAction::Graph(args) => graph(cli, args),
        MlAction::Soap(args) => soap(cli, args),
        MlAction::Acsf(args) => acsf_command(cli, args),
        MlAction::Mbtr(args) => mbtr_command(cli, args),
        MlAction::MofFeatures(args) => mof_features_command(cli, args),
        MlAction::Split(args) => split(cli, args),
        MlAction::Import(args) => import(cli, args),
        MlAction::Predict(args) => predict(cli, args),
        MlAction::Select(args) => select(cli, args),
        MlAction::ActiveWorkflow(args) => active_workflow(cli, args),
    }
}

fn graph(cli: &Cli, args: &GraphArgs) -> Result<(), CliError> {
    let graph = crystal_graph(
        &structure(&args.path)?,
        &GraphOptions {
            cutoff: args.cutoff,
            max_neighbors: args.k,
        },
    )
    .map_err(CliError::core)?;
    emit_json(cli, &graph, args.out.as_deref())
}

fn soap(cli: &Cli, args: &SoapArgs) -> Result<(), CliError> {
    let descriptor = soap_power_spectrum(
        &structure(&args.path)?,
        &xtalforge::ml::SoapOptions {
            cutoff: args.cutoff,
            n_max: args.n_max,
            l_max: args.l_max,
            sigma: args.sigma,
            average: args.average,
        },
    )
    .map_err(CliError::core)?;
    emit_json(cli, &descriptor, args.out.as_deref())
}

fn acsf_command(cli: &Cli, args: &AcsfArgs) -> Result<(), CliError> {
    let descriptor = acsf(
        &structure(&args.path)?,
        &AcsfOptions {
            cutoff: args.cutoff,
            average: args.average,
            ..AcsfOptions::default()
        },
    )
    .map_err(CliError::core)?;
    emit_json(cli, &descriptor, args.out.as_deref())
}

fn mbtr_command(cli: &Cli, args: &MbtrArgs) -> Result<(), CliError> {
    let descriptor = mbtr(
        &structure(&args.path)?,
        &MbtrOptions {
            cutoff: args.cutoff,
            n_grid: args.grid,
            sigma: args.sigma,
            angular: !args.no_angles,
            ..MbtrOptions::default()
        },
    )
    .map_err(CliError::core)?;
    emit_json(cli, &descriptor, args.out.as_deref())
}

fn mof_features_command(cli: &Cli, args: &MofFeaturesArgs) -> Result<(), CliError> {
    let features = mof_features(
        &structure(&args.path)?,
        &MofFeaturesOptions {
            max_path: args.max_path,
            pores: !args.no_pores,
            ..MofFeaturesOptions::default()
        },
    )
    .map_err(CliError::core)?;
    emit_json(cli, &features, args.out.as_deref())
}

fn split(cli: &Cli, args: &SplitArgs) -> Result<(), CliError> {
    let dataset = Dataset::load(&args.path).map_err(CliError::core)?;
    let strategy = match args.strategy.as_str() {
        "random" => SplitStrategy::Random,
        "composition" => SplitStrategy::Composition,
        "cluster" => SplitStrategy::Cluster {
            matcher: Default::default(),
        },
        "prototype" => SplitStrategy::Prototype,
        other => {
            return Err(CliError::input(format!("unknown split strategy `{other}`")));
        }
    };
    let assignment = split_dataset(
        &dataset,
        &SplitOptions {
            train: args.train,
            validation: args.validation,
            test: args.test,
            seed: args.seed,
            strategy,
        },
    )
    .map_err(CliError::core)?;
    emit_json(cli, &assignment, args.out.as_deref())
}

fn import(cli: &Cli, args: &ImportArgs) -> Result<(), CliError> {
    let text = std::fs::read_to_string(&args.path)
        .map_err(|error| CliError::input(format!("could not read {}: {error}", args.path)))?;
    let format = match args.format.as_str() {
        "extxyz" => xtalforge::ml::ImportFormat::Extxyz,
        "json" => xtalforge::ml::ImportFormat::Json,
        "optimade" => xtalforge::ml::ImportFormat::Optimade,
        other => {
            return Err(CliError::input(format!("unknown import format `{other}`")));
        }
    };
    let report = import_dataset(
        &text,
        &xtalforge::ml::ImportOptions::new(format, &args.name),
    )
    .map_err(CliError::core)?;
    if let Some(path) = &args.out {
        report.dataset.save(path).map_err(CliError::core)?;
        if !cli.quiet {
            println!(
                "wrote {} records ({} skipped) to {path}",
                report.n_records, report.skipped
            );
        }
        return Ok(());
    }
    println!("{} records, {} skipped", report.n_records, report.skipped);
    Ok(())
}

fn feature_vector(
    structure: &xtalforge::core::Structure,
    kind: &str,
) -> Result<Vec<f64>, CliError> {
    match kind {
        "soap" => {
            let descriptor =
                soap_power_spectrum(structure, &Default::default()).map_err(CliError::core)?;
            descriptor
                .rows
                .first()
                .cloned()
                .ok_or_else(|| CliError::input("SOAP produced no row"))
        }
        "acsf" => {
            let descriptor = acsf(
                structure,
                &AcsfOptions {
                    average: true,
                    ..AcsfOptions::default()
                },
            )
            .map_err(CliError::core)?;
            descriptor
                .rows
                .first()
                .cloned()
                .ok_or_else(|| CliError::input("ACSF produced no row"))
        }
        "mbtr" => {
            let descriptor = mbtr(structure, &MbtrOptions::default()).map_err(CliError::core)?;
            Ok(descriptor.features)
        }
        other => Err(CliError::input(format!(
            "unknown feature kind `{other}`; use soap, acsf or mbtr"
        ))),
    }
}

fn predict(cli: &Cli, args: &PredictArgs) -> Result<(), CliError> {
    let dataset = Dataset::load(&args.path).map_err(CliError::core)?;
    let model_text = std::fs::read_to_string(&args.model)
        .map_err(|error| CliError::input(format!("could not read {}: {error}", args.model)))?;
    let model = PropertyModel::from_json(&model_text).map_err(CliError::core)?;
    let kind = args.features.clone();
    let prediction = predict_dataset(&dataset, &model, |structure| {
        feature_vector(structure, &kind)
            .map_err(|error| xtalforge::core::Error::InvalidInput(error.message))
    })
    .map_err(CliError::core)?;
    let updated =
        append_predictions(&dataset, &prediction, &args.target).map_err(CliError::core)?;
    if let Some(path) = &args.out {
        updated.save(path).map_err(CliError::core)?;
        if !cli.quiet {
            println!("wrote {} records to {path}", updated.len());
        }
    } else {
        emit_json(cli, &prediction, None)?;
    }
    Ok(())
}

fn select(cli: &Cli, args: &SelectArgs) -> Result<(), CliError> {
    let text = std::fs::read_to_string(&args.path)
        .map_err(|error| CliError::input(format!("could not read {}: {error}", args.path)))?;
    let value: serde_json::Value = serde_json::from_str(&text)
        .map_err(|error| CliError::input(format!("feature table: {error}")))?;
    let features: Vec<Vec<f64>> = serde_json::from_value(
        value
            .get("features")
            .cloned()
            .ok_or_else(|| CliError::input("feature table has no `features`"))?,
    )
    .map_err(|error| CliError::input(format!("features: {error}")))?;
    let uncertainties: Vec<f64> = serde_json::from_value(
        value
            .get("uncertainties")
            .cloned()
            .unwrap_or_else(|| serde_json::json!(vec![0.0; features.len()])),
    )
    .map_err(|error| CliError::input(format!("uncertainties: {error}")))?;
    let strategy = match args.strategy.as_str() {
        "uncertainty" => SelectionStrategy::Uncertainty,
        "farthest_point" => SelectionStrategy::FarthestPoint,
        "diversity" => SelectionStrategy::Diversity,
        other => {
            return Err(CliError::input(format!(
                "unknown selection strategy `{other}`"
            )));
        }
    };
    let selected = select_batch(
        &features,
        &uncertainties,
        &SelectionOptions {
            batch: args.batch,
            strategy,
            diversity_radius: args.diversity_radius,
        },
    )
    .map_err(CliError::core)?;
    emit_json(cli, &selected, args.out.as_deref())
}

fn active_workflow(cli: &Cli, args: &ActiveWorkflowArgs) -> Result<(), CliError> {
    let document = active_learning_workflow(
        &args.name,
        &args.dataset,
        &args.model,
        args.batch,
        args.rounds,
        args.checkpoint.as_deref(),
    );
    emit_json(cli, &document, args.out.as_deref())
}

/// Keeps the property model trait visible for downstream imports.
pub(crate) fn _trait_check(model: &dyn Predictor) -> &str {
    model.name()
}
