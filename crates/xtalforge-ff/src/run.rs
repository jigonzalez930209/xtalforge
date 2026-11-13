//! Cancellation and progress reporting for long force-field runs (Task 14.2.2).
//!
//! Every optimizer, scan and molecular-dynamics run accepts an optional
//! [`RunControl`]: a cancellation flag that another thread can set and a
//! callback that receives one [`ProgressEvent`] per recorded step.  The
//! defaults are "never cancel" and "no callback", so a run through
//! [`crate::fire`] or [`crate::relaxed_scan`] behaves exactly as before.
//!
//! A cancelled run is not an error: the optimizer, scan or trajectory returns
//! the work completed so far and sets its `cancelled` flag, which keeps the
//! last geometry usable.  The flag is polled once per step, so cancellation
//! latency is one force evaluation.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use serde::{Deserialize, Serialize};

/// A cancellation flag shared between a run and the thread that stops it.
///
/// [`CancelToken::flag`] returns the raw flag, which is what
/// [`RunControl::with_cancel`] and the external-process runner expect.
#[derive(Debug, Clone, Default)]
pub struct CancelToken {
    flag: Arc<AtomicBool>,
}

impl CancelToken {
    /// A new token that is not cancelled.
    #[must_use]
    pub fn new() -> Self {
        Self {
            flag: Arc::new(AtomicBool::new(false)),
        }
    }

    /// Requests cancellation; already-running steps finish before the run
    /// notices.
    pub fn cancel(&self) {
        self.flag.store(true, Ordering::SeqCst);
    }

    /// Clears the request so the token can be reused by another run.
    pub fn reset(&self) {
        self.flag.store(false, Ordering::SeqCst);
    }

    /// Whether cancellation has been requested.
    #[must_use]
    pub fn is_cancelled(&self) -> bool {
        self.flag.load(Ordering::SeqCst)
    }

    /// Borrows the raw flag.
    #[must_use]
    pub fn flag(&self) -> &AtomicBool {
        &self.flag
    }

    /// Shares the flag with a calculator or another thread.
    #[must_use]
    pub fn shared(&self) -> Arc<AtomicBool> {
        Arc::clone(&self.flag)
    }
}

/// One progress record emitted during a run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum ProgressEvent {
    /// One optimizer iteration.
    Optimizer {
        /// Iteration number, starting at the first evaluated step.
        step: usize,
        /// Potential energy in eV.
        energy: f64,
        /// Largest force on a free site in eV/Angstrom.
        max_force: f64,
    },
    /// One finished point of a relaxed scan.
    Scan {
        /// Zero-based point index.
        point: usize,
        /// Requested coordinate value.
        target: f64,
        /// Coordinate value reached after relaxation.
        achieved: f64,
        /// Potential energy in eV, penalty excluded.
        energy: f64,
    },
    /// One sampled step of a molecular-dynamics run.
    Dynamics {
        /// Integration step.
        step: usize,
        /// Potential energy in eV.
        potential: f64,
        /// Kinetic energy in eV.
        kinetic: f64,
        /// Instantaneous temperature in kelvin.
        temperature: f64,
    },
}

/// Optional cancellation and progress reporting for one run.
///
/// The control borrows its flag and callback for the duration of the run.
/// Both are optional; `RunControl::new()` means "run to completion, report
/// nothing".
#[derive(Default)]
pub struct RunControl<'a> {
    cancel: Option<&'a AtomicBool>,
    progress: Option<&'a mut dyn FnMut(&ProgressEvent)>,
}

impl std::fmt::Debug for RunControl<'_> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("RunControl")
            .field("cancel", &self.cancel.is_some())
            .field("progress", &self.progress.is_some())
            .finish()
    }
}

impl<'a> RunControl<'a> {
    /// A control that never cancels and reports nothing.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Watches `cancel`.
    #[must_use]
    pub fn with_cancel(mut self, cancel: &'a AtomicBool) -> Self {
        self.cancel = Some(cancel);
        self
    }

    /// Watches a [`CancelToken`].
    #[must_use]
    pub fn with_token(self, token: &'a CancelToken) -> Self {
        self.with_cancel(token.flag())
    }

    /// Calls `progress` once per recorded step.
    #[must_use]
    pub fn with_progress(mut self, progress: &'a mut dyn FnMut(&ProgressEvent)) -> Self {
        self.progress = Some(progress);
        self
    }

    /// Whether the run should stop.
    #[must_use]
    pub fn is_cancelled(&self) -> bool {
        self.cancel.is_some_and(|flag| flag.load(Ordering::SeqCst))
    }

    /// Emits one event when a callback is installed.
    pub fn report(&mut self, event: ProgressEvent) {
        if let Some(progress) = self.progress.as_mut() {
            progress(&event);
        }
    }
}
