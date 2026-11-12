//! Thread-pool options shared by the surfaces (Task 11.1.2).

/// Worker-thread policy for the algorithms that have a parallel path.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ParallelOptions {
    /// Number of workers; `None` uses the Rayon default.
    pub threads: Option<usize>,
}

impl ParallelOptions {
    /// Creates options with an optional worker count.
    #[must_use]
    pub const fn new(threads: Option<usize>) -> Self {
        Self { threads }
    }

    /// Runs `function` inside a scoped pool with the requested worker count.
    ///
    /// Results are independent of the thread count because the parallel
    /// reductions keep a fixed order (AGENTS.md rule 6).
    ///
    /// # Errors
    ///
    /// Returns a message when the pool cannot be built.
    #[cfg(feature = "parallel")]
    pub fn install<T, F>(&self, function: F) -> Result<T, String>
    where
        F: FnOnce() -> T + Send,
        T: Send,
    {
        match self.threads {
            None => Ok(function()),
            Some(threads) => {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .map_err(|error| format!("could not build the thread pool: {error}"))?;
                Ok(pool.install(function))
            }
        }
    }

    /// Runs `function` directly when the parallel feature is off.
    ///
    /// # Errors
    ///
    /// Never returns an error in this build.
    #[cfg(not(feature = "parallel"))]
    pub fn install<T, F>(&self, function: F) -> Result<T, String>
    where
        F: FnOnce() -> T + Send,
        T: Send,
    {
        Ok(function())
    }
}

/// Installs the process-wide Rayon pool with `threads` workers.
///
/// Returns `true` when the pool was installed and `false` when this build has
/// no parallel feature. A pool can be installed only once per process.
///
/// # Errors
///
/// Returns a message when the pool was already installed or cannot be built.
#[cfg(feature = "parallel")]
pub fn set_global_threads(threads: usize) -> Result<bool, String> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .map(|()| true)
        .map_err(|error| format!("could not install the global thread pool: {error}"))
}

/// No-op when the parallel feature is off.
///
/// # Errors
///
/// Never returns an error in this build.
#[cfg(not(feature = "parallel"))]
pub fn set_global_threads(_threads: usize) -> Result<bool, String> {
    Ok(false)
}
