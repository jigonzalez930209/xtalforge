//! Cooperative cancellation and progress reporting for long powder runs
//! (Task 14.1.5).
//!
//! Every entry point that iterates for more than a fraction of a second
//! accepts a [`RunControl`]. The control carries an optional [`Cancellation`]
//! flag, checked between work units, and an optional [`ProgressCallback`].
//!
//! A cancelled run returns [`Error::InvalidInput`] with the message
//! `the run was cancelled`; the core error enum has no dedicated
//! cancellation variant, so the stable `code()` of a cancelled run is
//! `invalid_input`. This is the documented mapping until Phase 17 adds the
//! variant.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

/// A cloneable cooperative-cancellation flag.
///
/// Clones share the same flag, so a caller keeps one handle while a worker
/// receives another. The flag is a plain `AtomicBool`; it is available on
/// `wasm32` without threads.
#[derive(Debug, Clone, Default)]
pub struct Cancellation {
    flag: Arc<AtomicBool>,
}

impl Cancellation {
    /// Creates a flag that is not cancelled.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets the flag; the running operation stops at its next check.
    pub fn cancel(&self) {
        self.flag.store(true, Ordering::SeqCst);
    }

    /// Whether [`cancel`](Self::cancel) has been called.
    #[must_use]
    pub fn is_cancelled(&self) -> bool {
        self.flag.load(Ordering::SeqCst)
    }

    /// Clears the flag so the handle can be reused.
    pub fn reset(&self) {
        self.flag.store(false, Ordering::SeqCst);
    }
}

/// One progress event.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Progress {
    /// Name of the stage that produced the event.
    pub stage: String,
    /// Completed work units (function evaluations or iterations).
    pub done: usize,
    /// Current objective value (Rwp, chi-squared or energy), when known.
    pub value: Option<f64>,
}

/// Callback invoked with every progress event.
pub type ProgressCallback = Arc<dyn Fn(&Progress) + Send + Sync>;

/// Cancellation and progress controls shared by the powder entry points.
///
/// The default control never cancels and never reports.
#[derive(Clone, Default)]
pub struct RunControl {
    cancellation: Option<Cancellation>,
    progress: Option<ProgressCallback>,
}

impl std::fmt::Debug for RunControl {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("RunControl")
            .field("cancelled", &self.is_cancelled())
            .field("has_progress", &self.progress.is_some())
            .finish()
    }
}

impl PartialEq for RunControl {
    /// Two controls are equal when they carry the same cancellation state and
    /// both either have or lack a progress sink; the callback identities are
    /// not comparable.
    fn eq(&self, other: &Self) -> bool {
        self.is_cancelled() == other.is_cancelled()
            && self.progress.is_some() == other.progress.is_some()
    }
}

impl RunControl {
    /// Creates an inactive control.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Copy with a cancellation flag.
    #[must_use]
    pub fn with_cancellation(mut self, cancellation: Cancellation) -> Self {
        self.cancellation = Some(cancellation);
        self
    }

    /// Copy with a progress callback.
    #[must_use]
    pub fn with_progress<F>(mut self, progress: F) -> Self
    where
        F: Fn(&Progress) + Send + Sync + 'static,
    {
        self.progress = Some(Arc::new(progress));
        self
    }

    /// Whether a cancellation has been requested.
    #[must_use]
    pub fn is_cancelled(&self) -> bool {
        self.cancellation
            .as_ref()
            .is_some_and(Cancellation::is_cancelled)
    }

    /// Returns an error when the run has been cancelled.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] with `the run was cancelled` when the
    /// attached flag is set.
    pub fn check(&self) -> Result<(), Error> {
        if self.is_cancelled() {
            return Err(Error::InvalidInput("the run was cancelled".to_owned()));
        }
        Ok(())
    }

    /// Emits one progress event when a callback is attached.
    pub fn report(&self, stage: &str, done: usize, value: Option<f64>) {
        if let Some(progress) = &self.progress {
            progress(&Progress {
                stage: stage.to_owned(),
                done,
                value,
            });
        }
    }
}
