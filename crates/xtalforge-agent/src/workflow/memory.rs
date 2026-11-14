//! In-memory workflow execution for hosts without a filesystem (Task 17.2.5).
//!
//! The browser cannot read or write real paths. This runner takes the input
//! documents as a map of `path -> text`, executes the same workflow model and
//! verbs as [`run_workflow`](super::run_workflow), and returns the written
//! documents as a map. No checkpoint is produced because there is no disk to
//! resume from; the report otherwise keeps the same shape and ordering.

use std::collections::BTreeMap;
use std::path::Path;

use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use super::checkpoint::fingerprint as workflow_fingerprint;
use super::template::{TemplateContext, render_template};
use super::{
    ItemReport, ItemStatus, Problem, StepKind, StepReport, Workflow, WorkflowError, WorkflowReport,
};
use crate::Session;

/// Report and written documents from an in-memory run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MemoryWorkflowReport {
    /// The usual workflow report.
    pub report: WorkflowReport,
    /// Written documents keyed by their rendered path.
    pub outputs: BTreeMap<String, String>,
}

/// Runs a workflow over in-memory documents.
///
/// `files` maps the input paths named in the workflow to their text content.
/// A `read` step takes the text from that map (or from its rendered `path`,
/// when the step names a different key); a `write` step stores the document in
/// the returned `outputs` instead of touching the filesystem. Every other step
/// runs through the ordinary agent verb catalogue.
///
/// # Errors
///
/// Returns the validation errors of the workflow or a malformed step.
pub fn run_workflow_memory(
    workflow: &Workflow,
    files: &BTreeMap<String, String>,
) -> Result<MemoryWorkflowReport, WorkflowError> {
    let summary = workflow.validate()?;
    let fingerprint = workflow_fingerprint(workflow)?;
    let inputs: Vec<String> = files.keys().cloned().collect();
    let mut outputs = BTreeMap::new();
    let mut items = Vec::with_capacity(inputs.len());
    for (index, input) in inputs.iter().enumerate() {
        let mut session = Session::new();
        let mut steps = Vec::with_capacity(summary.steps.len());
        let mut failure = None;
        for (raw, step) in workflow.steps.iter().zip(summary.steps.iter()) {
            let context = TemplateContext::new(Path::new(input), index + 1);
            let outcome: Result<crate::Answer, Problem> = (|| match step.kind {
                StepKind::Read => {
                    let template = step.path.as_deref().unwrap_or("{{input}}");
                    let path =
                        render_template(template, &context).map_err(|error| error.problem())?;
                    let text = files.get(&path).ok_or_else(|| Problem {
                        code: "io".to_owned(),
                        message: format!("in-memory input {path:?} was not supplied"),
                        path: Some(path.clone()),
                        step: Some(step.name.clone()),
                    })?;
                    let mut arguments = json!({"text": text, "name": step.name.clone()});
                    if let Some(format) = &step.format
                        && let Some(object) = arguments.as_object_mut()
                    {
                        object.insert("format".to_owned(), Value::String(format.clone()));
                    }
                    session.run("load", &arguments).map_err(agent_problem)
                }
                StepKind::Write => {
                    let template = step.path.as_deref().unwrap_or("output");
                    let path =
                        render_template(template, &context).map_err(|error| error.problem())?;
                    let format = step
                        .format
                        .clone()
                        .or_else(|| {
                            Path::new(&path)
                                .extension()
                                .and_then(|extension| extension.to_str())
                                .map(str::to_owned)
                        })
                        .ok_or_else(|| Problem {
                            code: "argument".to_owned(),
                            message: format!(
                                "output path {path:?} has no extension; give a format"
                            ),
                            path: Some(path.clone()),
                            step: Some(step.name.clone()),
                        })?;
                    let structure = session
                        .document(None)
                        .map(|document| document.structure().clone())
                        .map_err(agent_problem)?;
                    let text = xtalforge_io::registry::write_with(&structure, &format)
                        .map_err(core_problem)?;
                    let bytes = text.len();
                    outputs.insert(path.clone(), text);
                    Ok(crate::Answer::Saved {
                        document: session.active_name().unwrap_or_default().to_owned(),
                        path,
                        format,
                        bytes,
                    })
                }
                _ => {
                    let verb = step.verb.as_deref().unwrap_or_default();
                    let arguments = raw
                        .arguments
                        .as_ref()
                        .map(|map| Value::Object(map.clone()))
                        .unwrap_or_else(|| json!({}));
                    session.run(verb, &arguments).map_err(agent_problem)
                }
            })();
            match outcome {
                Ok(answer) => match answer.to_json() {
                    Ok(value) => steps.push(StepReport {
                        step: step.name.clone(),
                        kind: step.kind,
                        answer: Some(value),
                    }),
                    Err(error) => {
                        failure = Some(Problem {
                            code: error.code().to_owned(),
                            message: error.to_string(),
                            path: None,
                            step: Some(step.name.clone()),
                        });
                        break;
                    }
                },
                Err(problem) => {
                    failure = Some(problem);
                    break;
                }
            }
        }
        let status = if failure.is_some() {
            ItemStatus::Failed
        } else {
            ItemStatus::Ok
        };
        items.push(ItemReport {
            input: input.clone(),
            status,
            resumed: false,
            steps,
            error: failure,
        });
    }
    let succeeded = items
        .iter()
        .filter(|item| item.status == ItemStatus::Ok)
        .count();
    let report = WorkflowReport {
        workflow_version: workflow.workflow_version,
        name: workflow.name.clone(),
        fingerprint,
        jobs: 1,
        inputs: items.len(),
        succeeded,
        failed: items.len() - succeeded,
        resumed: 0,
        items,
    };
    Ok(MemoryWorkflowReport { report, outputs })
}

fn agent_problem(error: crate::Error) -> Problem {
    Problem {
        code: error.code().to_owned(),
        message: error.to_string(),
        path: None,
        step: None,
    }
}

fn core_problem(error: xtalforge_core::Error) -> Problem {
    Problem {
        code: error.code().to_owned(),
        message: error.to_string(),
        path: None,
        step: None,
    }
}
