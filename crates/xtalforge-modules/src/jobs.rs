//! Module jobs and reports (Task 5.10).

use std::collections::VecDeque;

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

/// Job state.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum JobStatus {
    /// Waiting in the queue.
    Queued,
    /// Running.
    Running,
    /// Finished.
    Done,
    /// Failed.
    Failed,
    /// Cancelled.
    Cancelled,
}

/// A queued module invocation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Job {
    /// Stable id.
    pub id: u64,
    /// Module name.
    pub module: String,
    /// JSON parameters.
    pub parameters: serde_json::Value,
    /// State.
    pub status: JobStatus,
}

/// FIFO job queue.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct JobQueue {
    jobs: VecDeque<Job>,
    next: u64,
}

impl JobQueue {
    /// Empty queue.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a job and returns its id.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty module name.
    pub fn submit(&mut self, module: &str, parameters: serde_json::Value) -> Result<u64, Error> {
        if module.trim().is_empty() {
            return Err(Error::InvalidInput("a job needs a module name".to_owned()));
        }
        let id = self.next;
        self.next += 1;
        self.jobs.push_back(Job {
            id,
            module: module.to_owned(),
            parameters,
            status: JobStatus::Queued,
        });
        Ok(id)
    }

    /// Pops the next queued job and marks it running.
    pub fn start_next(&mut self) -> Option<Job> {
        let mut job = self.jobs.pop_front()?;
        job.status = JobStatus::Running;
        self.jobs.push_back(job.clone());
        Some(job)
    }

    /// Marks a running job finished.
    pub fn finish(&mut self, id: u64) -> bool {
        self.set_status(id, JobStatus::Done)
    }

    /// Marks a running job failed.
    pub fn fail(&mut self, id: u64) -> bool {
        self.set_status(id, JobStatus::Failed)
    }

    /// Cancels a queued job.
    pub fn cancel(&mut self, id: u64) -> bool {
        self.set_status(id, JobStatus::Cancelled)
    }

    /// All jobs in submission order.
    #[must_use]
    pub fn jobs(&self) -> Vec<&Job> {
        self.jobs.iter().collect()
    }

    fn set_status(&mut self, id: u64, status: JobStatus) -> bool {
        if let Some(job) = self.jobs.iter_mut().find(|job| job.id == id) {
            job.status = status;
            true
        } else {
            false
        }
    }
}
