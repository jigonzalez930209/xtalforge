//! In-process quantum methods (Task 4.10).
//!
//! The scientific backend is sci-form's GFN0, GFN1 and GFN2-xTB
//! implementation, behind the optional `sci-form` feature.  That dependency is
//! **not vendored in this workspace**, so this build cannot link it: the
//! capability table is real and stable today, while [`adapter`] reports a typed
//! [`Error::InvalidInput`] instead of pretending to run a method.  The same
//! refusal is returned whether or not the feature is enabled; the message
//! differs so a user knows whether to rebuild or to install the vendored
//! backend.
//!
//! The capability table records what the sci-form adapter supports once the
//! backend is linked: all three methods accept periodic cells and return
//! forces, and none returns an analytic stress.  A stress request must be
//! refused by the caller rather than replaced by a finite-difference value
//! silently.

use serde::{Deserialize, Serialize};

use crate::Error;

/// One in-process semi-empirical method.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum QuantumMethod {
    /// GFN0-xTB.
    Gfn0,
    /// GFN1-xTB.
    Gfn1,
    /// GFN2-xTB.
    Gfn2,
}

impl QuantumMethod {
    /// Stable method name used by the CLI and bindings.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Gfn0 => "gfn0",
            Self::Gfn1 => "gfn1",
            Self::Gfn2 => "gfn2",
        }
    }
}

/// What one method can answer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct QuantumCapabilities {
    /// Accepts a periodic cell rather than a molecule only.
    pub periodic: bool,
    /// Returns forces.
    pub forces: bool,
    /// Returns an analytic stress tensor.
    pub stress: bool,
}

/// Capabilities of one method.
///
/// All three GFN methods take periodic cells and provide forces; the adapter
/// exposes no analytic stress for any of them.
#[must_use]
pub fn capabilities(method: QuantumMethod) -> QuantumCapabilities {
    match method {
        QuantumMethod::Gfn0 | QuantumMethod::Gfn1 | QuantumMethod::Gfn2 => QuantumCapabilities {
            periodic: true,
            forces: true,
            stress: false,
        },
    }
}

/// Zero-size handle for the optional in-process calculator.
///
/// The handle exists even when the backend is not linked so callers can code
/// against the final API shape.  [`adapter`] is the only constructor and
/// currently always fails.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct QuantumAdapter;

/// Builds the adapter for one method.
///
/// # Errors
///
/// Without the `sci-form` feature returns
/// `Error::InvalidInput("sci-form support is not compiled; rebuild with
/// --features sci-form")`.  With the feature returns the clear
/// `"sci-form backend is not linked in this build"` error because the
/// dependency is not vendored here.
pub fn adapter(_method: QuantumMethod) -> Result<QuantumAdapter, Error> {
    #[cfg(not(feature = "sci-form"))]
    {
        Err(Error::InvalidInput(
            "sci-form support is not compiled; rebuild with --features sci-form".to_owned(),
        ))
    }
    #[cfg(feature = "sci-form")]
    {
        Err(Error::InvalidInput(
            "sci-form backend is not linked in this build; the dependency is not vendored"
                .to_owned(),
        ))
    }
}
