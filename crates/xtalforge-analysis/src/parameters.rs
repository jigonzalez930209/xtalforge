//! Persistent powder parameter sets (Task 14.1.5).
//!
//! A [`ParameterSet`] is a named, versioned map of floating-point parameters
//! serialized as JSON. It is the format a caller saves after a refinement and
//! reloads before the next one, so a Pawley, Rietveld or instrument state
//! survives a process boundary on every surface (the save/load helpers need a
//! filesystem and are therefore not built for WebAssembly; the JSON helpers
//! work everywhere).
//!
//! The values are plain `f64`, sorted by key, and the schema version is
//! [`PARAMETER_SET_VERSION`]. A set with a different version is refused rather
//! than silently reinterpreted.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

use crate::instrument::{Caglioti, FcjAsymmetry, InstrumentProfile, Monochromator};
use crate::pawley::{PawleyOptions, PawleyResult};
use crate::refine::RefinedPhase;

/// Schema version written by [`ParameterSet::new`] and accepted on load.
pub const PARAMETER_SET_VERSION: u32 = 1;

/// A named, versioned map of refinement parameters.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ParameterSet {
    /// Schema version; see [`PARAMETER_SET_VERSION`].
    pub version: u32,
    /// Free-form name, for example the phase label.
    pub name: String,
    /// Parameters sorted by key.
    pub values: BTreeMap<String, f64>,
}

impl ParameterSet {
    /// Creates an empty set with the current schema version.
    #[must_use]
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            version: PARAMETER_SET_VERSION,
            name: name.into(),
            values: BTreeMap::new(),
        }
    }

    /// Inserts one finite parameter.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-finite value or an empty key.
    pub fn set(&mut self, key: impl Into<String>, value: f64) -> Result<(), Error> {
        let key = key.into();
        if key.trim().is_empty() {
            return Err(Error::InvalidInput(
                "a parameter key must not be empty".to_owned(),
            ));
        }
        if !value.is_finite() {
            return Err(Error::InvalidInput(format!(
                "parameter {key} must be finite, got {value}"
            )));
        }
        self.values.insert(key, value);
        Ok(())
    }

    /// Reads one parameter.
    #[must_use]
    pub fn get(&self, key: &str) -> Option<f64> {
        self.values.get(key).copied()
    }

    /// Merges the values of another set of the same schema version.
    ///
    /// Existing keys are overwritten.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a different schema version.
    pub fn merge(&mut self, other: &Self) -> Result<(), Error> {
        if other.version != PARAMETER_SET_VERSION {
            return Err(Error::InvalidInput(format!(
                "parameter set schema {} is not supported (expected {PARAMETER_SET_VERSION})",
                other.version
            )));
        }
        for (key, value) in &other.values {
            self.values.insert(key.clone(), *value);
        }
        Ok(())
    }

    /// Serializes the set as pretty JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the serializer fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string_pretty(self)
            .map_err(|error| Error::InvalidInput(format!("parameter set encoding: {error}")))
    }

    /// Parses a set from JSON and checks its schema version and values.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for invalid JSON, an unsupported schema
    /// version or a non-finite value.
    pub fn from_json(text: &str) -> Result<Self, Error> {
        let set: Self = serde_json::from_str(text)
            .map_err(|error| Error::InvalidInput(format!("parameter set decoding: {error}")))?;
        if set.version != PARAMETER_SET_VERSION {
            return Err(Error::InvalidInput(format!(
                "parameter set schema {} is not supported (expected {PARAMETER_SET_VERSION})",
                set.version
            )));
        }
        for (key, value) in &set.values {
            if !value.is_finite() {
                return Err(Error::InvalidInput(format!(
                    "parameter {key} must be finite"
                )));
            }
        }
        Ok(set)
    }

    /// Writes the set as JSON to a file.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the file cannot be written.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn save(&self, path: impl AsRef<std::path::Path>) -> Result<(), Error> {
        let path = path.as_ref();
        std::fs::write(path, self.to_json()?).map_err(|error| {
            Error::InvalidInput(format!("could not write {}: {error}", path.display()))
        })
    }

    /// Reads a set from a JSON file.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the file cannot be read or parsed.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn load(path: impl AsRef<std::path::Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|error| {
            Error::InvalidInput(format!("could not read {}: {error}", path.display()))
        })?;
        Self::from_json(&text)
    }

    /// Collects the parameters of a finished structural refinement.
    ///
    /// Every named parameter is copied unchanged, with `rwp` added.
    #[must_use]
    pub fn from_refined_phase(name: impl Into<String>, phase: &RefinedPhase) -> Self {
        let mut set = Self::new(name);
        for (key, value) in &phase.parameters {
            let _ = set.set(key.clone(), *value);
        }
        let _ = set.set("rwp", phase.rwp);
        set
    }

    /// Collects the parameters of a finished Pawley fit.
    ///
    /// Keys: `rwp`, `eta`, `fwhm` (the fitted size term), `strain`, `zero` and
    /// `cell.a` … `cell.gamma`.
    #[must_use]
    pub fn from_pawley(name: impl Into<String>, fit: &PawleyResult, eta: f64) -> Self {
        let mut set = Self::new(name);
        let _ = set.set("rwp", fit.rwp);
        let _ = set.set("eta", eta);
        let _ = set.set("fwhm", fit.size);
        let _ = set.set("strain", fit.strain);
        let _ = set.set("zero", fit.zero_shift);
        for (index, value) in fit.cell.iter().enumerate() {
            let key = [
                "cell.a",
                "cell.b",
                "cell.c",
                "cell.alpha",
                "cell.beta",
                "cell.gamma",
            ]
            .get(index)
            .copied()
            .unwrap_or("cell");
            let _ = set.set(key, *value);
        }
        set
    }

    /// Copies the stored `eta` and `fwhm` onto a Pawley option struct.
    ///
    /// Values outside the valid range are ignored.
    #[must_use]
    pub fn pawley_options(&self, base: &PawleyOptions) -> PawleyOptions {
        let mut options = *base;
        if let Some(eta) = self.get("eta")
            && (0.0..=1.0).contains(&eta)
        {
            options.eta = eta;
        }
        if let Some(fwhm) = self.get("fwhm")
            && fwhm.is_finite()
            && fwhm > 0.0
        {
            options.fwhm = fwhm;
        }
        options
    }

    /// Collects the parameters of an instrument profile.
    ///
    /// Keys: `profile.fwhm`, `profile.eta`, `caglioti.u`, `caglioti.v`,
    /// `caglioti.w`, `asymmetry.fcj1`, `asymmetry.fcj2` and
    /// `monochromator.polarization`.
    #[must_use]
    pub fn from_instrument_profile(name: impl Into<String>, profile: &InstrumentProfile) -> Self {
        let mut set = Self::new(name);
        let _ = set.set("profile.fwhm", profile.fwhm);
        let _ = set.set("profile.eta", profile.eta);
        if let Some(caglioti) = profile.caglioti {
            let _ = set.set("caglioti.u", caglioti.u);
            let _ = set.set("caglioti.v", caglioti.v);
            let _ = set.set("caglioti.w", caglioti.w);
        }
        if let Some(asymmetry) = profile.asymmetry {
            let _ = set.set("asymmetry.fcj1", asymmetry.fcj1);
            let _ = set.set("asymmetry.fcj2", asymmetry.fcj2);
        }
        let _ = set.set(
            "monochromator.polarization",
            profile.monochromator.polarization,
        );
        set
    }

    /// Applies the stored instrument keys onto a profile; missing keys keep the
    /// value of `base`.
    ///
    /// A Caglioti law is installed only when all three coefficients are
    /// present, and an asymmetry only when both coefficients are present.
    #[must_use]
    pub fn instrument_profile(&self, base: InstrumentProfile) -> InstrumentProfile {
        let mut profile = base;
        if let Some(value) = self.get("profile.fwhm")
            && value.is_finite()
            && value > 0.0
        {
            profile.fwhm = value;
        }
        if let Some(value) = self.get("profile.eta")
            && (0.0..=1.0).contains(&value)
        {
            profile.eta = value;
        }
        if let (Some(u), Some(v), Some(w)) = (
            self.get("caglioti.u"),
            self.get("caglioti.v"),
            self.get("caglioti.w"),
        ) {
            profile.caglioti = Some(Caglioti { u, v, w });
        }
        if let (Some(fcj1), Some(fcj2)) = (self.get("asymmetry.fcj1"), self.get("asymmetry.fcj2")) {
            profile.asymmetry = Some(FcjAsymmetry { fcj1, fcj2 });
        }
        if let Some(value) = self.get("monochromator.polarization")
            && (0.0..=1.0).contains(&value)
        {
            profile.monochromator = Monochromator {
                polarization: value,
            };
        }
        profile
    }
}
