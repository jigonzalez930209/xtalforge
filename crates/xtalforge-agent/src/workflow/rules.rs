//! Validation rules for a workflow document (Task 16.1.1).
//!
//! The rules mirror the published JSON Schema and check every verb against the
//! live catalogue, so a step can only name an operation this build has.

use std::collections::BTreeSet;

use super::WorkflowError;
use super::model::{Step, StepKind, StepSummary, WORKFLOW_VERSION, Workflow, WorkflowSummary};
use super::template::validate_template;

/// Verbs that describe session plumbing rather than a workflow step.
pub(super) const RESERVED_VERBS: &[&str] = &[
    "capabilities",
    "list_documents",
    "load",
    "save",
    "render",
    "close",
    "undo",
    "redo",
    "history",
    "changes",
    "workflow_schema",
    "workflow_validate",
    "workflow_run",
];

/// Validates a whole document and returns its effective summary.
///
/// # Errors
///
/// Returns [`WorkflowError::Version`] or [`WorkflowError::Invalid`] with the
/// JSON path of the offending value.
pub(super) fn validate(workflow: &Workflow) -> Result<WorkflowSummary, WorkflowError> {
    if workflow.workflow_version != WORKFLOW_VERSION {
        return Err(WorkflowError::Version {
            found: u64::from(workflow.workflow_version),
            supported: WORKFLOW_VERSION,
        });
    }
    if workflow.name.trim().is_empty() {
        return Err(invalid("name", "must not be empty"));
    }
    if workflow.steps.is_empty() {
        return Err(invalid("steps", "needs at least one step"));
    }
    let first_is_read = workflow
        .steps
        .first()
        .is_some_and(|step| step.kind == StepKind::Read);
    if !first_is_read {
        return Err(invalid(
            "steps[0].kind",
            "the first step must be `read` so every item starts from an input file",
        ));
    }
    let mut labels = BTreeSet::new();
    let mut steps = Vec::with_capacity(workflow.steps.len());
    for (index, step) in workflow.steps.iter().enumerate() {
        let label = step
            .name
            .clone()
            .unwrap_or_else(|| format!("step-{}", index + 1));
        let label_path = format!("steps[{index}].name");
        if label.trim().is_empty()
            || !label
                .chars()
                .all(|character| character.is_ascii_alphanumeric() || "_.-".contains(character))
        {
            return Err(invalid(
                label_path,
                "must contain only ASCII letters, digits, `_`, `.` and `-`",
            ));
        }
        if !labels.insert(label.clone()) {
            return Err(invalid(
                format!("steps[{index}].name"),
                format!("duplicate step name {label:?}"),
            ));
        }
        steps.push(summarize_step(step, index, &label)?);
    }
    Ok(WorkflowSummary {
        workflow_version: workflow.workflow_version,
        name: workflow.name.clone(),
        description: workflow.description.clone(),
        steps,
    })
}

fn summarize_step(step: &Step, index: usize, label: &str) -> Result<StepSummary, WorkflowError> {
    let base = format!("steps[{index}]");
    match step.kind {
        StepKind::Read => {
            forbid_verb(step, &base)?;
            let path = step.path.clone().unwrap_or_else(|| "{{input}}".to_owned());
            validate_template(&path, &format!("{base}.path"))?;
            if let Some(format) = &step.format {
                check_format(format, &format!("{base}.format"), false)?;
            }
            Ok(StepSummary {
                name: label.to_owned(),
                kind: step.kind,
                verb: None,
                path: Some(path),
                format: step.format.clone(),
            })
        }
        StepKind::Write => {
            forbid_verb(step, &base)?;
            let Some(path) = &step.path else {
                return Err(invalid(
                    format!("{base}.path"),
                    "a `write` step needs a path",
                ));
            };
            validate_template(path, &format!("{base}.path"))?;
            if let Some(format) = &step.format {
                check_format(format, &format!("{base}.format"), true)?;
            }
            Ok(StepSummary {
                name: label.to_owned(),
                kind: step.kind,
                verb: None,
                path: Some(path.clone()),
                format: step.format.clone(),
            })
        }
        kind => {
            if step.path.is_some() {
                return Err(invalid(
                    format!("{base}.path"),
                    "only `read` and `write` steps take a path",
                ));
            }
            if step.format.is_some() {
                return Err(invalid(
                    format!("{base}.format"),
                    "only `read` and `write` steps take a format",
                ));
            }
            let Some(name) = &step.verb else {
                return Err(invalid(
                    format!("{base}.verb"),
                    format!("a `{}` step needs a verb from the catalogue", kind.as_str()),
                ));
            };
            let verb = crate::verbs::lookup(name).ok_or_else(|| {
                invalid(
                    format!("{base}.verb"),
                    format!("unknown verb {name:?}; call `capabilities` for the catalogue"),
                )
            })?;
            if RESERVED_VERBS.contains(&verb.name) {
                return Err(invalid(
                    format!("{base}.verb"),
                    format!("verb {name:?} is session plumbing and cannot be a workflow step"),
                ));
            }
            match kind {
                StepKind::Transform if !verb.mutation => {
                    return Err(invalid(
                        format!("{base}.kind"),
                        format!("{name:?} does not mutate; use `analyse` or `compute`"),
                    ));
                }
                StepKind::Analyse | StepKind::Compute if verb.mutation => {
                    return Err(invalid(
                        format!("{base}.kind"),
                        format!("{name:?} mutates the document; use `transform`"),
                    ));
                }
                _ => {}
            }
            Ok(StepSummary {
                name: label.to_owned(),
                kind,
                verb: Some(verb.name.to_owned()),
                path: None,
                format: None,
            })
        }
    }
}

fn forbid_verb(step: &Step, base: &str) -> Result<(), WorkflowError> {
    if step.verb.is_some() {
        return Err(invalid(
            format!("{base}.verb"),
            "only `transform`, `analyse` and `compute` steps take a verb",
        ));
    }
    if step.arguments.is_some() {
        return Err(invalid(
            format!("{base}.arguments"),
            "only `transform`, `analyse` and `compute` steps take arguments",
        ));
    }
    Ok(())
}

fn check_format(format: &str, path: &str, write: bool) -> Result<(), WorkflowError> {
    let Some(spec) = xtalforge_io::registry::find(format) else {
        return Err(invalid(path, format!("unknown format {format:?}")));
    };
    let supported = if write {
        spec.write.is_some()
    } else {
        spec.read.is_some()
    };
    if !supported {
        let direction = if write { "written" } else { "read" };
        return Err(invalid(
            path,
            format!("format {:?} cannot be {direction}", spec.name),
        ));
    }
    Ok(())
}

fn invalid(path: impl Into<String>, reason: impl Into<String>) -> WorkflowError {
    WorkflowError::Invalid {
        path: path.into(),
        reason: reason.into(),
    }
}
