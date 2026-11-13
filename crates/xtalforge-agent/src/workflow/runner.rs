//! The workflow runner: per-item sessions, isolation, order and checkpoints
//! (Task 16.1.2).
//!
//! Every input gets a fresh [`Session`], so a document name, an undo history
//! or an edit from one item can never leak into another. The report is built
//! in the input order regardless of the worker count, and an item failure
//! stops only that item.
//!
//! A checkpoint is a JSON file written after every completed item; on resume
//! the items already marked `ok` are skipped and failed items run again. The
//! file records the workflow fingerprint and the model version, so a
//! checkpoint from another workflow is refused instead of mixing results.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use serde_json::{Map, Value, json};

use super::WorkflowError;
use super::checkpoint::{CHECKPOINT_VERSION, WorkflowCheckpoint, fingerprint};
use super::exec::{RunContext, run_pending, write_checkpoint};
use super::model::{Problem, StepKind, StepSummary, Workflow};
use super::template::{TemplateContext, render_template};
use crate::{Answer, Session};

/// How a run processes its inputs.
#[derive(Debug, Clone, PartialEq)]
pub struct RunOptions {
    /// Input files, directories or expanded paths.
    pub inputs: Vec<PathBuf>,
    /// Items processed at once; `0` and `1` both mean sequential.
    pub jobs: usize,
    /// Checkpoint file; written after every item and resumed when it exists.
    pub checkpoint: Option<PathBuf>,
    /// When `false`, an existing checkpoint is ignored and overwritten.
    pub resume: bool,
}

impl Default for RunOptions {
    fn default() -> Self {
        Self {
            inputs: Vec::new(),
            jobs: 1,
            checkpoint: None,
            resume: true,
        }
    }
}

/// Whether one input succeeded.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ItemStatus {
    /// Every step ran.
    Ok,
    /// One step failed; the item stopped there.
    Failed,
}

/// The answer of one step.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StepReport {
    /// Effective step label.
    pub step: String,
    /// Step kind.
    pub kind: StepKind,
    /// The agent answer as JSON.
    pub answer: Option<Value>,
}

/// What happened to one input.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ItemReport {
    /// Input path as it was expanded.
    pub input: String,
    /// Overall outcome.
    pub status: ItemStatus,
    /// True when a checkpoint supplied this result instead of a run.
    pub resumed: bool,
    /// Completed steps, in workflow order.
    pub steps: Vec<StepReport>,
    /// The failure, when the item stopped early.
    pub error: Option<Problem>,
}

/// The result of one workflow invocation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WorkflowReport {
    /// Workflow model version.
    pub workflow_version: u32,
    /// Workflow name.
    pub name: String,
    /// Stable fingerprint of the document.
    pub fingerprint: String,
    /// Effective worker count.
    pub jobs: usize,
    /// Number of inputs after expansion, sorting and deduplication.
    pub inputs: usize,
    /// Items that ran to completion or were already complete.
    pub succeeded: usize,
    /// Items with a failed step.
    pub failed: usize,
    /// Items skipped because the checkpoint had them as `ok`.
    pub resumed: usize,
    /// One entry per input, in deterministic input order.
    pub items: Vec<ItemReport>,
}

/// Runs a workflow over every input and returns the report.
///
/// # Errors
///
/// Returns the validation errors of the document, checkpoint mismatches and
/// filesystem errors. An input that fails during processing is reported in
/// [`WorkflowReport::items`], not returned as an error.
pub fn run_workflow(
    workflow: &Workflow,
    options: &RunOptions,
) -> Result<WorkflowReport, WorkflowError> {
    let summary = workflow.validate()?;
    let fingerprint = fingerprint(workflow)?;
    let mut inputs = options.inputs.clone();
    inputs.sort();
    inputs.dedup();
    // The browser build has no threads; it runs sequentially, which keeps the
    // report and the checkpoint identical.
    let jobs = if cfg!(target_arch = "wasm32") {
        1
    } else {
        options.jobs.max(1).min(inputs.len().max(1))
    };

    let mut results: BTreeMap<usize, ItemReport> = BTreeMap::new();
    if let (Some(path), true) = (&options.checkpoint, options.resume)
        && path.exists()
    {
        let checkpoint = WorkflowCheckpoint::read(path)?;
        if checkpoint.checkpoint_version != CHECKPOINT_VERSION {
            return Err(WorkflowError::Checkpoint {
                path: path.display().to_string(),
                reason: format!(
                    "layout version {} is not supported; this build writes {CHECKPOINT_VERSION}",
                    checkpoint.checkpoint_version
                ),
            });
        }
        if checkpoint.workflow_version != workflow.workflow_version {
            return Err(WorkflowError::Checkpoint {
                path: path.display().to_string(),
                reason: format!(
                    "workflow_version {} does not match the document ({})",
                    checkpoint.workflow_version, workflow.workflow_version
                ),
            });
        }
        if checkpoint.name != workflow.name || checkpoint.fingerprint != fingerprint {
            return Err(WorkflowError::Checkpoint {
                path: path.display().to_string(),
                reason: format!(
                    "belongs to workflow {:?} ({}); use `--fresh` or another checkpoint file",
                    checkpoint.name, checkpoint.fingerprint
                ),
            });
        }
        for item in checkpoint.items {
            if item.status != ItemStatus::Ok {
                continue;
            }
            if let Some((index, _)) = inputs
                .iter()
                .enumerate()
                .find(|(_, path)| path.to_string_lossy() == item.input)
                && !results.contains_key(&index)
            {
                let mut item = item;
                item.resumed = true;
                results.insert(index, item);
            }
        }
    }

    let work: Vec<(usize, &PathBuf)> = inputs
        .iter()
        .enumerate()
        .filter(|(index, _)| !results.contains_key(index))
        .collect();
    let context = RunContext {
        steps: &summary.steps,
        raw_steps: &workflow.steps,
        jobs,
        checkpoint: &options.checkpoint,
        workflow,
        fingerprint: &fingerprint,
        inputs: &inputs,
    };
    if !work.is_empty() {
        run_pending(&context, &work, &mut results)?;
    }
    if options.checkpoint.is_some() {
        write_checkpoint(&context, &results)?;
    }

    let items = ordered_items(&inputs, &results);
    let succeeded = items
        .iter()
        .filter(|item| item.status == ItemStatus::Ok)
        .count();
    let failed = items.len() - succeeded;
    let resumed = items.iter().filter(|item| item.resumed).count();
    Ok(WorkflowReport {
        workflow_version: workflow.workflow_version,
        name: workflow.name.clone(),
        fingerprint,
        jobs,
        inputs: inputs.len(),
        succeeded,
        failed,
        resumed,
        items,
    })
}

pub(super) struct ResolvedStep<'a> {
    pub(super) summary: &'a StepSummary,
    pub(super) arguments: Option<&'a Map<String, Value>>,
}

pub(super) fn run_item(resolved: &[ResolvedStep<'_>], path: &Path, index: usize) -> ItemReport {
    let input = path.to_string_lossy().into_owned();
    let mut session = Session::new();
    let mut steps = Vec::with_capacity(resolved.len());
    for step in resolved {
        match run_step(&mut session, step, path, index) {
            Ok(answer) => match answer.to_json() {
                Ok(value) => steps.push(StepReport {
                    step: step.summary.name.clone(),
                    kind: step.summary.kind,
                    answer: Some(value),
                }),
                Err(error) => {
                    return failed(
                        input,
                        steps,
                        step.summary.name.clone(),
                        error.code(),
                        error.to_string(),
                    );
                }
            },
            Err(problem) => {
                return failed(
                    input,
                    steps,
                    step.summary.name.clone(),
                    &problem.code,
                    problem.message,
                );
            }
        }
    }
    ItemReport {
        input,
        status: ItemStatus::Ok,
        resumed: false,
        steps,
        error: None,
    }
}

fn run_step(
    session: &mut Session,
    step: &ResolvedStep<'_>,
    input: &Path,
    index: usize,
) -> Result<Answer, Problem> {
    let summary = step.summary;
    let context = TemplateContext::new(input, index + 1);
    match summary.kind {
        StepKind::Read => {
            let template = summary.path.as_deref().unwrap_or("{{input}}");
            let path = render_template(template, &context).map_err(|error| error.problem())?;
            let mut arguments = json!({"path": path, "name": summary.name.clone()});
            if let Some(format) = &summary.format
                && let Some(object) = arguments.as_object_mut()
            {
                object.insert("format".to_owned(), Value::String(format.clone()));
            }
            session.run("load", &arguments).map_err(agent_problem)
        }
        StepKind::Write => {
            let template = summary.path.as_deref().unwrap_or_default();
            let path = render_template(template, &context).map_err(|error| error.problem())?;
            ensure_parent(Path::new(&path))?;
            let mut arguments = json!({"path": path});
            if let Some(format) = &summary.format
                && let Some(object) = arguments.as_object_mut()
            {
                object.insert("format".to_owned(), Value::String(format.clone()));
            }
            session.run("save", &arguments).map_err(agent_problem)
        }
        _ => {
            let verb = summary.verb.as_deref().unwrap_or_default();
            let arguments = step
                .arguments
                .map(|map| Value::Object(map.clone()))
                .unwrap_or_else(|| json!({}));
            session.run(verb, &arguments).map_err(agent_problem)
        }
    }
}

fn ensure_parent(path: &Path) -> Result<(), Problem> {
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        std::fs::create_dir_all(parent).map_err(|error| Problem {
            code: "io".to_owned(),
            message: format!("{}: {error}", parent.display()),
            path: Some(parent.display().to_string()),
            step: None,
        })?;
    }
    Ok(())
}

fn agent_problem(error: crate::Error) -> Problem {
    Problem {
        code: error.code().to_owned(),
        message: error.to_string(),
        path: None,
        step: None,
    }
}

fn failed(
    input: String,
    steps: Vec<StepReport>,
    step: String,
    code: &str,
    message: String,
) -> ItemReport {
    ItemReport {
        input,
        status: ItemStatus::Failed,
        resumed: false,
        steps,
        error: Some(Problem {
            code: code.to_owned(),
            message,
            path: None,
            step: Some(step),
        }),
    }
}

fn ordered_items(inputs: &[PathBuf], results: &BTreeMap<usize, ItemReport>) -> Vec<ItemReport> {
    inputs
        .iter()
        .enumerate()
        .filter_map(|(index, _)| results.get(&index).cloned())
        .collect()
}

/// Stored items never carry the run-specific `resumed` marker; only the
/// report of the run that skipped them does.
pub(super) fn checkpoint_items(
    inputs: &[PathBuf],
    results: &BTreeMap<usize, ItemReport>,
) -> Vec<ItemReport> {
    ordered_items(inputs, results)
        .into_iter()
        .map(|mut item| {
            item.resumed = false;
            item
        })
        .collect()
}
