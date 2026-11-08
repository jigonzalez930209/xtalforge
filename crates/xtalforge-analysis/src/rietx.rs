//! Read-only inventory of the rietx powder bridge (Task 7.1).
//!
//! This module is deliberately **not a port** of rietx. It records what the
//! Crystal-Builder bridge uses and the acceptance limits a later port must
//! meet, so the scope decision is visible in code and can be asserted by a
//! test.
//!
//! # Source and license
//!
//! * Bridge under study: `Crystal-Builder/xtal/powder/bridge.py`.
//! * Reference implementation: [rietx](https://github.com/yue-here/rietx) and
//!   its kernels package. Both are **MIT** licensed, so a port would be
//!   permitted; this crate still reimplements the methods from the published
//!   powder-diffraction literature and uses rietx only as a behavioural
//!   oracle.
//! * Attribution is recorded in `THIRD_PARTY.md`; this file carries the note
//!   in code so the inventory is self-describing.
//!
//! # Approximations
//!
//! The inventory is a static description. It does not import, link or call
//! rietx, and the acceptance check compares caller-supplied measured metrics
//! with the documented limits rather than running a reference fit.

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

/// Maximum accepted weighted profile R factor for a Pawley/Rietveld fit.
///
/// Task 7.1.4 asks for a stated margin around the reference result. Crystal-
/// Builder's benchmark fits are in the low tens of percent; 30 % leaves room
/// for the documented profile and scattering-factor approximations of this
/// crate while still rejecting a fit that has not converged.
pub const ACCEPTANCE_MAX_RWP: f64 = 0.30;

/// Maximum accepted absolute peak-position error, in degrees two-theta.
///
/// This is half of the default 0.1 degree sampling step: a line must land
/// inside the observed peak's own step, not merely near it.
pub const ACCEPTANCE_MAX_POSITION_ERROR: f64 = 0.05;

/// License note for rietx and rietx-kernels.
pub const RIETX_LICENSE_NOTE: &str = "rietx and rietx-kernels are MIT licensed \
     (https://github.com/yue-here/rietx); this crate reimplements the methods \
     and records rietx as a behavioural reference, see THIRD_PARTY.md.";

/// A read-only description of the rietx bridge surface.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RietxBridge {
    /// Name of the reference implementation.
    pub name: &'static str,
    /// Capabilities the bridge exposes, in pipeline order.
    pub capabilities: Vec<&'static str>,
    /// License of the reference implementation.
    pub license: &'static str,
    /// Where the bridge and its entry points live.
    pub source: &'static str,
}

/// Measured acceptance metrics of one benchmark fit.
///
/// Both fields are maxima over the reflections of a run: the fit's Rwp is a
/// single number, while the position error is the largest absolute difference
/// between an observed and a calculated line position in degrees two-theta.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AcceptanceMetrics {
    /// Largest weighted profile R factor of the run.
    pub max_rwp: f64,
    /// Largest absolute peak-position error in degrees two-theta.
    pub max_position_error: f64,
}

/// Returns the static inventory of the rietx bridge.
#[must_use]
pub fn rietx_inventory() -> RietxBridge {
    RietxBridge {
        name: "rietx",
        capabilities: vec![
            "pattern_io",
            "background",
            "profile",
            "le_bail",
            "pawley",
            "rietveld",
            "indexing",
        ],
        license: "MIT",
        source: "Crystal-Builder xtal/powder/bridge.py -> rietx \
                 (https://github.com/yue-here/rietx)",
    }
}

/// Checks measured metrics against the acceptance limits of Task 7.1.4.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when either metric is not finite or
/// negative, or when it exceeds [`ACCEPTANCE_MAX_RWP`] respectively
/// [`ACCEPTANCE_MAX_POSITION_ERROR`]. The message names the metric, its value
/// and the limit so a caller can report the failure without parsing.
pub fn check_acceptance(metrics: &AcceptanceMetrics) -> Result<(), Error> {
    if !metrics.max_rwp.is_finite() || metrics.max_rwp < 0.0 {
        return Err(Error::InvalidInput(format!(
            "acceptance max_rwp must be finite and non-negative, got {}",
            metrics.max_rwp
        )));
    }
    if metrics.max_rwp > ACCEPTANCE_MAX_RWP {
        return Err(Error::InvalidInput(format!(
            "acceptance max_rwp {:.4} exceeds the limit {ACCEPTANCE_MAX_RWP:.4}",
            metrics.max_rwp
        )));
    }
    if !metrics.max_position_error.is_finite() || metrics.max_position_error < 0.0 {
        return Err(Error::InvalidInput(format!(
            "acceptance max_position_error must be finite and non-negative, got {}",
            metrics.max_position_error
        )));
    }
    if metrics.max_position_error > ACCEPTANCE_MAX_POSITION_ERROR {
        return Err(Error::InvalidInput(format!(
            "acceptance max_position_error {:.4} degrees exceeds the limit {ACCEPTANCE_MAX_POSITION_ERROR:.4}",
            metrics.max_position_error
        )));
    }
    Ok(())
}
