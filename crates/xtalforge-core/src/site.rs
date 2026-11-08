//! Asymmetric-unit sites: label, element, coordinates, occupancy and props.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::Error;
use crate::elements::element_by_symbol;

/// One site of the asymmetric unit.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Site {
    label: String,
    element: String,
    frac: [f64; 3],
    occupancy: f64,
    u_iso: Option<f64>,
    u_aniso: Option<[f64; 6]>,
    charge: Option<f64>,
    disorder_group: Option<u32>,
    props: BTreeMap<String, Value>,
}

impl Site {
    /// Creates a fully occupied site with no displacement parameters.
    ///
    /// The element is normalized through [`crate::elements::element_by_symbol`]
    /// (`Fe2+` becomes `Fe`, `D` stays deuterium).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty label, an unknown element,
    /// non-finite coordinates or an occupancy outside `[0, 1]`.
    pub fn new(label: impl Into<String>, element: &str, frac: [f64; 3]) -> Result<Self, Error> {
        let label = label.into();
        if label.trim().is_empty() {
            return Err(Error::InvalidInput(
                "site label must not be empty".to_owned(),
            ));
        }
        if !frac.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "site coordinates must be finite".to_owned(),
            ));
        }
        let element = element_by_symbol(element)?.symbol.clone();
        Ok(Self {
            label,
            element,
            frac,
            occupancy: 1.0,
            u_iso: None,
            u_aniso: None,
            charge: None,
            disorder_group: None,
            props: BTreeMap::new(),
        })
    }

    /// Site label.
    #[must_use]
    pub fn label(&self) -> &str {
        &self.label
    }

    /// Canonical element symbol.
    #[must_use]
    pub fn element(&self) -> &str {
        &self.element
    }

    /// Fractional coordinates.
    #[must_use]
    pub fn frac(&self) -> [f64; 3] {
        self.frac
    }

    /// Occupancy in `[0, 1]`.
    #[must_use]
    pub fn occupancy(&self) -> f64 {
        self.occupancy
    }

    /// Isotropic displacement parameter, when present.
    #[must_use]
    pub fn u_iso(&self) -> Option<f64> {
        self.u_iso
    }

    /// Anisotropic displacement parameters `U11 U22 U33 U12 U13 U23`.
    #[must_use]
    pub fn u_aniso(&self) -> Option<[f64; 6]> {
        self.u_aniso
    }

    /// Formal charge, when present.
    #[must_use]
    pub fn charge(&self) -> Option<f64> {
        self.charge
    }

    /// Disorder assembly or group number.
    #[must_use]
    pub fn disorder_group(&self) -> Option<u32> {
        self.disorder_group
    }

    /// Open extension map.
    #[must_use]
    pub fn props(&self) -> &BTreeMap<String, Value> {
        &self.props
    }

    /// Whether the site is partially occupied.
    #[must_use]
    pub fn is_partial(&self) -> bool {
        self.occupancy < 1.0
    }

    /// Copy with a different occupancy.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the occupancy is not finite or
    /// outside `[0, 1]`.
    pub fn with_occupancy(mut self, occupancy: f64) -> Result<Self, Error> {
        if !occupancy.is_finite() || !(0.0..=1.0).contains(&occupancy) {
            return Err(Error::InvalidInput(format!(
                "occupancy must be in [0, 1], got {occupancy}"
            )));
        }
        self.occupancy = occupancy;
        Ok(self)
    }

    /// Copy with different fractional coordinates.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for non-finite coordinates.
    pub fn with_frac(mut self, frac: [f64; 3]) -> Result<Self, Error> {
        if !frac.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "site coordinates must be finite".to_owned(),
            ));
        }
        self.frac = frac;
        Ok(self)
    }

    /// Copy with an isotropic displacement parameter.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the value is not finite or is
    /// negative.
    pub fn with_u_iso(mut self, u_iso: f64) -> Result<Self, Error> {
        if !u_iso.is_finite() || u_iso < 0.0 {
            return Err(Error::InvalidInput(format!(
                "Uiso must be non-negative and finite, got {u_iso}"
            )));
        }
        self.u_iso = Some(u_iso);
        Ok(self)
    }

    /// Copy with anisotropic displacement parameters.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for non-finite values.
    pub fn with_u_aniso(mut self, u_aniso: [f64; 6]) -> Result<Self, Error> {
        if !u_aniso.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "anisotropic displacement parameters must be finite".to_owned(),
            ));
        }
        self.u_aniso = Some(u_aniso);
        Ok(self)
    }

    /// Copy with a formal charge.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-finite charge.
    pub fn with_charge(mut self, charge: f64) -> Result<Self, Error> {
        if !charge.is_finite() {
            return Err(Error::InvalidInput("site charge must be finite".to_owned()));
        }
        self.charge = Some(charge);
        Ok(self)
    }

    /// Copy with a disorder group.
    #[must_use]
    pub fn with_disorder_group(mut self, group: u32) -> Self {
        self.disorder_group = Some(group);
        self
    }

    /// Copy with one extension property.
    #[must_use]
    pub fn with_prop(mut self, key: impl Into<String>, value: Value) -> Self {
        self.props.insert(key.into(), value);
        self
    }
}
