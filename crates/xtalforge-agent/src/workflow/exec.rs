//! Sequential and thread-based item execution (Task 16.1.2).
//!
//! One worker count drives both paths: `jobs` of one (and the browser build)
//! runs in order on the calling thread, larger counts share the work through
//! an atomic cursor. Results are merged by input index, so completion order
//! never changes the report order.

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc;

use super::WorkflowError;
use super::checkpoint::{CHECKPOINT_VERSION, WorkflowCheckpoint};
use super::model::{StepSummary, Workflow};
use super::runner::{ItemReport, ResolvedStep, checkpoint_items, run_item};

pub(super) struct RunContext<'a> {
    pub(super) steps: &'a [StepSummary],
    pub(super) raw_steps: &'a [super::Step],
    pub(super) jobs: usize,
    pub(super) checkpoint: &'a Option<PathBuf>,
    pub(super) workflow: &'a Workflow,
    pub(super) fingerprint: &'a str,
    pub(super) inputs: &'a [PathBuf],
}

impl RunContext<'_> {
    fn resolved(&self) -> Vec<ResolvedStep<'_>> {
        self.steps
            .iter()
            .zip(self.raw_steps.iter())
            .map(|(summary, raw)| ResolvedStep {
                summary,
                arguments: raw.arguments.as_ref(),
            })
            .collect()
    }
}

pub(super) fn write_checkpoint(
    context: &RunContext<'_>,
    results: &BTreeMap<usize, ItemReport>,
) -> Result<(), WorkflowError> {
    let Some(path) = context.checkpoint else {
        return Ok(());
    };
    WorkflowCheckpoint {
        checkpoint_version: CHECKPOINT_VERSION,
        workflow_version: context.workflow.workflow_version,
        name: context.workflow.name.clone(),
        fingerprint: context.fingerprint.to_owned(),
        items: checkpoint_items(context.inputs, results),
    }
    .write(path)
}

pub(super) fn run_pending(
    context: &RunContext<'_>,
    work: &[(usize, &PathBuf)],
    results: &mut BTreeMap<usize, ItemReport>,
) -> Result<(), WorkflowError> {
    if context.jobs <= 1 {
        return run_sequential(context, work, results);
    }
    run_parallel(context, work, results)
}

fn run_sequential(
    context: &RunContext<'_>,
    work: &[(usize, &PathBuf)],
    results: &mut BTreeMap<usize, ItemReport>,
) -> Result<(), WorkflowError> {
    let resolved = context.resolved();
    for (index, path) in work {
        results.insert(*index, run_item(&resolved, path, *index));
        write_checkpoint(context, results)?;
    }
    Ok(())
}

fn run_parallel(
    context: &RunContext<'_>,
    work: &[(usize, &PathBuf)],
    results: &mut BTreeMap<usize, ItemReport>,
) -> Result<(), WorkflowError> {
    let resolved = context.resolved();
    let next = AtomicUsize::new(0);
    let (sender, receiver) = mpsc::channel::<(usize, ItemReport)>();
    let mut checkpoint_error = None;
    std::thread::scope(|scope| {
        for _ in 0..context.jobs.min(work.len()) {
            let sender = sender.clone();
            let next = &next;
            let resolved = &resolved;
            scope.spawn(move || {
                loop {
                    let slot = next.fetch_add(1, Ordering::Relaxed);
                    let Some((index, path)) = work.get(slot) else {
                        break;
                    };
                    let report = run_item(resolved, path, *index);
                    if sender.send((*index, report)).is_err() {
                        break;
                    }
                }
            });
        }
        drop(sender);
        while let Ok((index, report)) = receiver.recv() {
            results.insert(index, report);
            if checkpoint_error.is_none()
                && let Err(error) = write_checkpoint(context, results)
            {
                checkpoint_error = Some(error);
            }
        }
    });
    match checkpoint_error {
        Some(error) => Err(error),
        None => Ok(()),
    }
}
