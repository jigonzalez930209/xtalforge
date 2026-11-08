//! Integration tests for jobs (Task 5.10).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::json;
use xtalforge_modules::{JobQueue, JobStatus};

#[test]
fn jobs_move_through_the_queue() {
    let mut queue = JobQueue::new();
    let first = queue.submit("pxrd", json!({"step": 0.02})).unwrap();
    let second = queue.submit("porosity", json!({"spacing": 0.2})).unwrap();
    assert_eq!(queue.jobs().len(), 2);
    let job = queue.start_next().unwrap();
    assert_eq!(job.id, first);
    assert!(queue.finish(first));
    assert!(queue.fail(second));
    assert!(
        queue
            .jobs()
            .iter()
            .any(|job| job.id == first && job.status == JobStatus::Done)
    );
    assert!(
        queue
            .jobs()
            .iter()
            .any(|job| job.id == second && job.status == JobStatus::Failed)
    );
    assert!(!queue.finish(999));
    assert!(queue.submit("", json!({})).is_err());
    assert!(!queue.cancel(999));
}
