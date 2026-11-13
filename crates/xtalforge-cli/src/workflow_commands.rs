//! The `run` command: declarative workflows over many inputs (Task 16.1.2).
//!
//! The command validates the document, expands `--inputs`, runs the shared
//! workflow runner and prints the report. A failed item never stops the
//! others; the process exits with code 3 when at least one item failed, so a
//! shell loop can retry after fixing the inputs.

use std::path::{Path, PathBuf};

use xtalforge::agent::workflow::{
    RunOptions, Workflow, WorkflowError, WorkflowReport, expand_inputs, run_workflow, schema_text,
};

use clap::Args;

use crate::{Cli, CliError};

/// Arguments of `xtalforge run`.
#[derive(Debug, Args)]
pub(crate) struct RunArgs {
    /// Workflow file (JSON; a `.yaml` file must contain JSON).
    #[arg(required_unless_present = "schema")]
    pub(crate) workflow: Option<String>,
    /// Input file, directory or pattern; repeatable.
    #[arg(long = "inputs", value_name = "PATTERN")]
    pub(crate) inputs: Vec<String>,
    /// Items processed at once.
    #[arg(long, default_value_t = 1)]
    pub(crate) jobs: usize,
    /// Checkpoint file; an existing one is resumed.
    #[arg(long, value_name = "PATH")]
    pub(crate) checkpoint: Option<String>,
    /// Write the JSON report here as well as to standard output.
    #[arg(long, value_name = "PATH")]
    pub(crate) out: Option<String>,
    /// Ignore an existing checkpoint and start over.
    #[arg(long)]
    pub(crate) fresh: bool,
    /// Validate the workflow and print its steps without running it.
    #[arg(long)]
    pub(crate) validate: bool,
    /// Print the published workflow JSON Schema and exit.
    #[arg(long)]
    pub(crate) schema: bool,
}

pub(crate) fn run(cli: &Cli, args: &RunArgs) -> Result<(), CliError> {
    let RunArgs {
        workflow,
        inputs,
        jobs,
        checkpoint,
        out,
        fresh,
        validate,
        schema,
    } = args;
    let workflow = workflow.as_deref();
    let inputs = inputs.as_slice();
    let jobs = *jobs;
    let checkpoint = checkpoint.as_deref();
    let out = out.as_deref();
    let fresh = *fresh;
    let validate = *validate;
    let schema = *schema;
    if schema {
        let text = schema_text().map_err(json_error)?;
        println!("{text}");
        return Ok(());
    }
    let Some(path) = workflow else {
        return Err(CliError::input(
            "a workflow file is required; try `xtalforge run --schema`",
        ));
    };
    let workflow = Workflow::from_path(path).map_err(workflow_error)?;
    if validate {
        let summary = workflow.validate().map_err(workflow_error)?;
        if cli.json {
            println!("{}", serde_json::to_string(&summary).map_err(json_error)?);
        } else {
            println!(
                "{} (workflow_version {})",
                summary.name, summary.workflow_version
            );
            for (index, step) in summary.steps.iter().enumerate() {
                let target = step
                    .verb
                    .clone()
                    .or_else(|| step.path.clone())
                    .unwrap_or_default();
                println!("  {}. {:<10} {target}", index + 1, step.kind.as_str());
            }
        }
        return Ok(());
    }
    if inputs.is_empty() {
        return Err(CliError::input(
            "--inputs is required; give a file, a directory or a pattern",
        ));
    }
    if jobs == 0 {
        return Err(CliError::input("--jobs must be at least 1"));
    }
    let expanded = expand_inputs(inputs).map_err(workflow_error)?;
    if expanded.is_empty() {
        return Err(CliError::input(format!(
            "no input matched {}",
            inputs.join(", ")
        )));
    }
    let report = run_workflow(
        &workflow,
        &RunOptions {
            inputs: expanded,
            jobs,
            checkpoint: checkpoint.map(PathBuf::from),
            resume: !fresh,
        },
    )
    .map_err(workflow_error)?;
    if let Some(path) = out {
        write_report(path, &report)?;
    }
    if cli.json {
        println!("{}", serde_json::to_string(&report).map_err(json_error)?);
    } else {
        print_report(&report);
    }
    if report.failed > 0 {
        return Err(CliError {
            code: "workflow_failed",
            category: "workflow",
            message: format!(
                "{} of {} items failed; fix them and rerun with the same --checkpoint to resume",
                report.failed, report.inputs
            ),
        });
    }
    Ok(())
}

fn write_report(path: &str, report: &WorkflowReport) -> Result<(), CliError> {
    let text = serde_json::to_string_pretty(report).map_err(json_error)?;
    if let Some(parent) = Path::new(path)
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        std::fs::create_dir_all(parent)
            .map_err(|error| CliError::new("io", format!("{}: {error}", parent.display())))?;
    }
    std::fs::write(path, text.as_bytes())
        .map_err(|error| CliError::new("io", format!("{path}: {error}")))
}

fn print_report(report: &WorkflowReport) {
    println!(
        "{}: {} ok, {} failed, {} resumed of {} inputs ({} job{})",
        report.name,
        report.succeeded,
        report.failed,
        report.resumed,
        report.inputs,
        report.jobs,
        if report.jobs == 1 { "" } else { "s" }
    );
    for item in &report.items {
        match &item.error {
            Some(problem) => {
                let step = problem.step.as_deref().unwrap_or("-");
                println!(
                    "  failed {} at {step}: {} ({})",
                    item.input, problem.message, problem.code
                );
            }
            None if item.resumed => println!("  resumed {}", item.input),
            None => println!("  ok      {}", item.input),
        }
    }
}

fn workflow_error(error: WorkflowError) -> CliError {
    CliError {
        code: error.code(),
        category: error.category(),
        message: error.to_string(),
    }
}

fn json_error(error: serde_json::Error) -> CliError {
    CliError::new("json_error", error.to_string())
}
