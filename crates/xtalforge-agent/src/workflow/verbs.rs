//! Workflow verbs exposed to agents and the MCP server (Task 16.1.3).
//!
//! `workflow_schema` returns the published JSON Schema, `workflow_validate`
//! checks a document and returns a structured verdict, and `workflow_run`
//! executes a document over a list of inputs with the same runner as the CLI.

use std::path::PathBuf;

use serde_json::Value;

use super::{RunOptions, Workflow, expand_inputs, run_workflow, schema_json};
use crate::answers::Answer;
use crate::run::optional_string;
use crate::{Error, Session};

/// `workflow_schema`: the JSON Schema a client uses to build a document.
pub(crate) fn schema_answer(_session: &mut Session, _arguments: &Value) -> Result<Answer, Error> {
    Ok(Answer::WorkflowSchema {
        workflow_version: super::WORKFLOW_VERSION,
        schema: schema_json(),
    })
}

/// `workflow_validate`: parses and validates without running.
pub(crate) fn validate_answer(_session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    match workflow_argument(arguments) {
        Ok(workflow) => {
            let summary = workflow.validate().map_err(Error::Workflow)?;
            Ok(Answer::WorkflowValidated {
                valid: true,
                summary: Some(summary),
                problem: None,
            })
        }
        Err(Error::Workflow(error)) => Ok(Answer::WorkflowValidated {
            valid: false,
            summary: None,
            problem: Some(error.problem()),
        }),
        Err(other) => Err(other),
    }
}

/// `workflow_run`: executes a document over the expanded inputs.
pub(crate) fn run_answer(_session: &mut Session, arguments: &Value) -> Result<Answer, Error> {
    let workflow = workflow_argument(arguments)?;
    let patterns = input_patterns(arguments)?;
    let inputs = expand_inputs(patterns).map_err(Error::Workflow)?;
    if inputs.is_empty() {
        return Err(Error::Argument {
            name: "inputs".to_owned(),
            reason: "no file matched the inputs".to_owned(),
        });
    }
    let jobs = crate::run::integer_or(arguments, "jobs", 1)?;
    if jobs == 0 {
        return Err(Error::Argument {
            name: "jobs".to_owned(),
            reason: "must be at least 1".to_owned(),
        });
    }
    let jobs = usize::try_from(jobs).map_err(|_| Error::Argument {
        name: "jobs".to_owned(),
        reason: "does not fit in a machine word".to_owned(),
    })?;
    let checkpoint = optional_string(arguments, "checkpoint")?
        .map(PathBuf::from)
        .filter(|path| !path.as_os_str().is_empty());
    let resume = optional_boolean(arguments, "resume")?.unwrap_or(true);
    let report = run_workflow(
        &workflow,
        &RunOptions {
            inputs,
            jobs,
            checkpoint,
            resume,
        },
    )
    .map_err(Error::Workflow)?;
    Ok(Answer::WorkflowRun { report })
}

fn workflow_argument(arguments: &Value) -> Result<Workflow, Error> {
    let inline = arguments.get("workflow");
    let text = arguments.get("text");
    let path = arguments.get("path");
    let given =
        usize::from(inline.is_some()) + usize::from(text.is_some()) + usize::from(path.is_some());
    if given == 0 {
        return Err(Error::Argument {
            name: "workflow".to_owned(),
            reason: "give one of `workflow` (object), `text` or `path`".to_owned(),
        });
    }
    if given > 1 {
        return Err(Error::Argument {
            name: "workflow".to_owned(),
            reason: "give only one of `workflow`, `text` or `path`".to_owned(),
        });
    }
    if let Some(value) = inline {
        if !value.is_object() {
            return Err(Error::Argument {
                name: "workflow".to_owned(),
                reason: "must be a JSON object".to_owned(),
            });
        }
        return Workflow::from_value(value.clone()).map_err(Error::Workflow);
    }
    if let Some(value) = text {
        let Some(text) = value.as_str() else {
            return Err(Error::Argument {
                name: "text".to_owned(),
                reason: "must be a string".to_owned(),
            });
        };
        return Workflow::from_text(text).map_err(Error::Workflow);
    }
    let Some(value) = path else {
        return Err(Error::Argument {
            name: "path".to_owned(),
            reason: "is required".to_owned(),
        });
    };
    let Some(path) = value.as_str() else {
        return Err(Error::Argument {
            name: "path".to_owned(),
            reason: "must be a string".to_owned(),
        });
    };
    Workflow::from_path(path).map_err(Error::Workflow)
}

fn input_patterns(arguments: &Value) -> Result<Vec<String>, Error> {
    let Some(value) = arguments.get("inputs") else {
        return Err(Error::Argument {
            name: "inputs".to_owned(),
            reason: "is required".to_owned(),
        });
    };
    let Some(items) = value.as_array() else {
        return Err(Error::Argument {
            name: "inputs".to_owned(),
            reason: "must be an array of paths or patterns".to_owned(),
        });
    };
    let mut patterns = Vec::with_capacity(items.len());
    for item in items {
        let Some(text) = item.as_str() else {
            return Err(Error::Argument {
                name: "inputs".to_owned(),
                reason: "every entry must be a string".to_owned(),
            });
        };
        patterns.push(text.to_owned());
    }
    if patterns.is_empty() {
        return Err(Error::Argument {
            name: "inputs".to_owned(),
            reason: "must contain at least one path or pattern".to_owned(),
        });
    }
    Ok(patterns)
}

fn optional_boolean(arguments: &Value, name: &str) -> Result<Option<bool>, Error> {
    match arguments.get(name) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(value)) => Ok(Some(*value)),
        Some(_) => Err(Error::Argument {
            name: name.to_owned(),
            reason: "must be a boolean".to_owned(),
        }),
    }
}
