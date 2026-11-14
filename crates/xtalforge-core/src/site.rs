//! Asymmetric-unit sites: label, element, coordinates, occupancy and props.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::Error;
use crate::elements::element_by_symbol;

/// A crystal axis along which a collinear magnetic moment is aligned.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum MagneticAxis {
    /// The `a` axis.
    A,
    /// The `b` axis.
    B,
    /// The `c` axis.
    C,
}

/// A magnetic moment in Bohr magnetons.
///
/// Components follow the magCIF `crystalaxis` convention: they are the
/// projections of the moment onto the unit vectors along the direct lattice
/// vectors `a`, `b` and `c`. They are therefore independent of the cell
/// lengths. [`MagneticMoment::Collinear`] stores one signed magnitude along a
/// single axis; [`MagneticMoment::NonCollinear`] stores all three components.
///
/// The value is additive site data: it is absent for non-magnetic sites, and
/// it does not change the serialized structure schema version; files written
/// before it existed still load with `None`. Transforms that copy sites keep
/// the moment verbatim. Re-expressing it in a changed crystal frame is part of
/// the magnetic-symmetry work (Task 19.1.2).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MagneticMoment {
    /// A collinear moment along one crystal axis.
    Collinear {
        /// Crystal axis of the moment.
        axis: MagneticAxis,
        /// Signed moment in Bohr magnetons.
        magnitude: f64,
    },
    /// A moment with components along the three crystal axes.
    NonCollinear {
        /// Moment components in Bohr magnetons along `a`, `b` and `c`.
        components: [f64; 3],
    },
}

impl MagneticMoment {
    /// Creates a collinear moment along one crystal axis.
    ///
    /// A zero magnitude has no axis and is stored as a zero
    /// [`MagneticMoment::NonCollinear`] vector.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the magnitude is not finite.
    pub fn collinear(axis: MagneticAxis, magnitude: f64) -> Result<Self, Error> {
        if !magnitude.is_finite() {
            return Err(Error::InvalidInput(
                "magnetic moment must be finite".to_owned(),
            ));
        }
        if magnitude == 0.0 {
            return Ok(Self::NonCollinear {
                components: [0.0; 3],
            });
        }
        Ok(Self::Collinear { axis, magnitude })
    }

    /// Creates a moment from crystal-axis components.
    ///
    /// Exactly one non-zero component is normalized to
    /// [`MagneticMoment::Collinear`] so that a magCIF round trip keeps the
    /// same representation.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for non-finite components.
    pub fn vector([x, y, z]: [f64; 3]) -> Result<Self, Error> {
        if ![x, y, z].iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "magnetic moment components must be finite".to_owned(),
            ));
        }
        let nonzero = usize::from(x != 0.0) + usize::from(y != 0.0) + usize::from(z != 0.0);
        match nonzero {
            1 if x != 0.0 => Self::collinear(MagneticAxis::A, x),
            1 if y != 0.0 => Self::collinear(MagneticAxis::B, y),
            1 => Self::collinear(MagneticAxis::C, z),
            _ => Ok(Self::NonCollinear {
                components: [x, y, z],
            }),
        }
    }

    /// Crystal-axis components in Bohr magnetons.
    #[must_use]
    pub fn components(&self) -> [f64; 3] {
        match *self {
            Self::Collinear { axis, magnitude } => match axis {
                MagneticAxis::A => [magnitude, 0.0, 0.0],
                MagneticAxis::B => [0.0, magnitude, 0.0],
                MagneticAxis::C => [0.0, 0.0, magnitude],
            },
            Self::NonCollinear { components } => components,
        }
    }

    /// Length of the moment vector in Bohr magnetons.
    #[must_use]
    pub fn magnitude(&self) -> f64 {
        let [x, y, z] = self.components();
        (x * x + y * y + z * z).sqrt()
    }

    /// Whether this is a collinear moment.
    #[must_use]
    pub fn is_collinear(&self) -> bool {
        matches!(self, Self::Collinear { .. })
    }

    /// Whether every component is zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        let [x, y, z] = self.components();
        x == 0.0 && y == 0.0 && z == 0.0
    }
}

/// One site of the asymmetric unit.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Site {
    label: String,
    element: String,
    frac: [f64; 3],
    occupancy: f64,
    /// Additive since Phase 19; absent in schema version 1 files written
    /// before magnetic moments existed.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    magnetic_moment: Option<MagneticMoment>,
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
            magnetic_moment: None,
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

    /// Magnetic moment in Bohr magnetons, when the site is magnetic.
    #[must_use]
    pub fn magnetic_moment(&self) -> Option<MagneticMoment> {
        self.magnetic_moment
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

    /// Copy with a different label.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty label.
    pub fn with_label(mut self, label: impl Into<String>) -> Result<Self, Error> {
        let label = label.into();
        if label.trim().is_empty() {
            return Err(Error::InvalidInput(
                "site label must not be empty".to_owned(),
            ));
        }
        self.label = label;
        Ok(self)
    }

    /// Copy with a different fractional position.
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

    /// Copy with a magnetic moment.
    #[must_use]
    pub fn with_magnetic_moment(mut self, moment: MagneticMoment) -> Self {
        self.magnetic_moment = Some(moment);
        self
    }

    /// Copy with a collinear magnetic moment along one crystal axis.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the magnitude is not finite.
    pub fn with_collinear_moment(
        mut self,
        axis: MagneticAxis,
        magnitude: f64,
    ) -> Result<Self, Error> {
        self.magnetic_moment = Some(MagneticMoment::collinear(axis, magnitude)?);
        Ok(self)
    }

    /// Copy with a non-collinear magnetic moment from crystal-axis components.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for non-finite components.
    pub fn with_moment_vector(mut self, components: [f64; 3]) -> Result<Self, Error> {
        self.magnetic_moment = Some(MagneticMoment::vector(components)?);
        Ok(self)
    }

    /// Copy without a magnetic moment.
    #[must_use]
    pub fn without_magnetic_moment(mut self) -> Self {
        self.magnetic_moment = None;
        self
    }

    /// Copy with a disorder group.
    #[must_use]
    pub fn with_disorder_group(mut self, group: u32) -> Self {
        self.disorder_group = Some(group);
        self
    }

    /// Copy without a disorder group.
    #[must_use]
    pub fn without_disorder_group(mut self) -> Self {
        self.disorder_group = None;
        self
    }

    /// Copy with one extension property.
    #[must_use]
    pub fn with_prop(mut self, key: impl Into<String>, value: Value) -> Self {
        self.props.insert(key.into(), value);
        self
    }

    /// Copy without one extension property.
    #[must_use]
    pub fn without_prop(mut self, key: &str) -> Self {
        self.props.remove(key);
        self
    }
}
