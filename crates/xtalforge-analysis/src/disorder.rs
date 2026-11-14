//! Partial occupancies and disorder groups as first-class data (Task 19.2.1).
//!
//! A disorder group is a set of sites that share a fractional position (or an
//! explicit [`Site::disorder_group`](xtalforge_core::Site::disorder_group)
//! number) and whose occupancies add up to at most one. The helpers in this
//! module read that information from the P1 expansion, so they see every
//! symmetry image exactly once.
//!
//! The composition, formula and density helpers are occupancy weighted. They
//! agree with [`Structure::formula`] and [`Structure::density`], which already
//! treat occupancies as first-class data; the extra value here is the grouped
//! view and the scattering-weight hook for custom powder-diffraction code.
//!
//! # PXRD weighting hook
//!
//! [`disorder_scattering_weights`] turns every disordered position into one
//! coherent amplitude `sum_i occ_i f_i(s)` plus an incoherent variance
//! `sum_i occ_i f_i(s)^2 - (sum_i occ_i f_i(s))^2`. X-ray, neutron and electron
//! amplitudes come from [`crate::scattering`]. The built-in
//! [`calculate_pxrd`](crate::calculate_pxrd) already applies the same
//! occupancy weighting internally; this entry point exists for callers that
//! build structure factors themselves.
//!
//! The module is deterministic: groups, positions and occupants are sorted by
//! key, element symbol and fractional coordinate.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use xtalforge_core::{Error, Structure, elements::element_by_symbol};

use crate::instrument::RadiationKind;
use crate::scattering::{electron_scattering_factor, neutron_scattering_length, xray_form_factor};

/// Relative tolerance used to decide whether a position is fully occupied.
const OCCUPANCY_TOLERANCE: f64 = 1e-9;

/// One species and its occupancy at a disordered position.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DisorderOccupant {
    /// Canonical element symbol.
    pub element: String,
    /// Occupancy summed over every P1 image with this element at the position.
    pub occupancy: f64,
}

/// One fractional position of a disorder group.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DisorderPosition {
    /// Fractional coordinates in the parent cell, wrapped into `[0, 1)`.
    pub frac: [f64; 3],
    /// Occupants sorted by element symbol.
    pub occupants: Vec<DisorderOccupant>,
}

impl DisorderPosition {
    /// Sum of the occupant occupancies.
    #[must_use]
    pub fn total_occupancy(&self) -> f64 {
        self.occupants
            .iter()
            .map(|occupant| occupant.occupancy)
            .sum()
    }

    /// Whether the position has more than one species or a vacancy.
    #[must_use]
    pub fn is_disordered(&self) -> bool {
        self.occupants.len() > 1 || self.total_occupancy() < 1.0 - OCCUPANCY_TOLERANCE
    }
}

/// A disorder group: one or more positions tied by a group number or by
/// sharing a fractional position.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DisorderGroup {
    /// Explicit disorder-group number, absent for position-based groups.
    pub disorder_group: Option<u32>,
    /// Positions of the group, sorted by fractional coordinate.
    pub positions: Vec<DisorderPosition>,
}

impl DisorderGroup {
    /// Sum of the occupant occupancies over every position.
    #[must_use]
    pub fn total_occupancy(&self) -> f64 {
        self.positions
            .iter()
            .map(DisorderPosition::total_occupancy)
            .sum()
    }

    /// Occupant occupancies aggregated over the positions of the group.
    ///
    /// This is the composition a concentration constraint uses when a group
    /// spans more than one split position.
    #[must_use]
    pub fn aggregate_occupants(&self) -> Vec<DisorderOccupant> {
        let mut totals: BTreeMap<String, f64> = BTreeMap::new();
        for position in &self.positions {
            for occupant in &position.occupants {
                *totals.entry(occupant.element.clone()).or_default() += occupant.occupancy;
            }
        }
        totals
            .into_iter()
            .map(|(element, occupancy)| DisorderOccupant { element, occupancy })
            .collect()
    }
}

/// Occupancy-weighted amount of one element in the cell.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ElementOccupancy {
    /// Canonical element symbol.
    pub element: String,
    /// Occupancy sum over every P1 image.
    pub occupancy: f64,
    /// Number of P1 atoms that contributed to the occupancy.
    pub atoms: usize,
}

/// Composition of a possibly disordered cell.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DisorderComposition {
    /// Elements in first-appearance order of the P1 expansion.
    pub elements: Vec<ElementOccupancy>,
    /// Disordered groups found in the same P1 expansion.
    pub disorder_groups: Vec<DisorderGroup>,
}

impl DisorderComposition {
    /// Occupancy-aware formula, for example `Zn0.5 O0.5`.
    #[must_use]
    pub fn formula(&self) -> String {
        self.elements
            .iter()
            .map(|element| {
                let count = element.occupancy;
                if (count - count.round()).abs() < 1e-9 {
                    let integer = count.round() as u64;
                    if integer == 1 {
                        element.element.clone()
                    } else {
                        format!("{}{integer}", element.element)
                    }
                } else {
                    format!("{}{count:.3}", element.element)
                }
            })
            .collect::<Vec<_>>()
            .join(" ")
    }

    /// Occupancy-weighted cell mass in unified atomic mass units.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when an element is missing from the
    /// element table.
    pub fn cell_mass(&self) -> Result<f64, Error> {
        let mut mass = 0.0;
        for element in &self.elements {
            mass += element_by_symbol(&element.element)?.mass * element.occupancy;
        }
        Ok(mass)
    }

    /// Density in g/cm3 from the cell mass and an external volume.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive or non-finite volume,
    /// or the errors of [`DisorderComposition::cell_mass`].
    pub fn density(&self, volume: f64) -> Result<f64, Error> {
        if !volume.is_finite() || volume <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "cell volume must be positive and finite, got {volume}"
            )));
        }
        Ok(self.cell_mass()? * 1.660_539_066_60 / volume)
    }
}

/// Collects the disordered groups of a structure from its P1 expansion.
///
/// Explicit `disorder_group` numbers are kept together even when the group
/// spans several split positions. A position without a group number becomes a
/// group when it carries more than one species or a vacancy. Fully occupied
/// single-species positions are omitted.
///
/// # Errors
///
/// Returns the errors of the P1 expansion.
pub fn disorder_groups(structure: &Structure) -> Result<Vec<DisorderGroup>, Error> {
    let p1 = structure.p1()?;
    let mut groups: BTreeMap<GroupKey, BTreeMap<String, PositionAccumulator>> = BTreeMap::new();
    for atom in &p1.atoms {
        let order = position_key(atom.frac);
        let key = atom
            .disorder_group
            .map_or_else(|| GroupKey::Implicit(order.clone()), GroupKey::Explicit);
        groups
            .entry(key)
            .or_default()
            .entry(order)
            .or_insert_with(|| PositionAccumulator::new(atom.frac))
            .add(&atom.element, atom.occupancy);
    }
    let mut output = Vec::new();
    for (key, positions) in groups {
        let positions: Vec<DisorderPosition> = positions
            .into_values()
            .map(PositionAccumulator::finish)
            .collect();
        let explicit = matches!(&key, GroupKey::Explicit(_));
        let disordered = explicit || positions.iter().any(|position| position.is_disordered());
        if !disordered {
            continue;
        }
        output.push(DisorderGroup {
            disorder_group: match key {
                GroupKey::Explicit(id) => Some(id),
                GroupKey::Implicit(_) => None,
            },
            positions,
        });
    }
    Ok(output)
}

/// Composition of the P1 cell with occupancies kept as fractional amounts.
///
/// # Errors
///
/// Returns the errors of the P1 expansion or of the element table.
pub fn disorder_composition(structure: &Structure) -> Result<DisorderComposition, Error> {
    let p1 = structure.p1()?;
    let mut elements: Vec<ElementOccupancy> = Vec::new();
    for atom in &p1.atoms {
        if let Some(entry) = elements
            .iter_mut()
            .find(|entry| entry.element == atom.element)
        {
            entry.occupancy += atom.occupancy;
            entry.atoms += 1;
        } else {
            elements.push(ElementOccupancy {
                element: atom.element.clone(),
                occupancy: atom.occupancy,
                atoms: 1,
            });
        }
    }
    Ok(DisorderComposition {
        elements,
        disorder_groups: disorder_groups(structure)?,
    })
}

/// Occupancy-aware formula of the whole cell.
///
/// # Errors
///
/// Returns the errors of [`disorder_composition`].
pub fn disorder_formula(structure: &Structure) -> Result<String, Error> {
    Ok(disorder_composition(structure)?.formula())
}

/// Occupancy-aware density of the whole cell in g/cm3.
///
/// # Errors
///
/// Returns the errors of [`disorder_composition`] or its density helper.
pub fn disorder_density(structure: &Structure) -> Result<f64, Error> {
    let composition = disorder_composition(structure)?;
    composition.density(structure.lattice().volume())
}

/// One occupant of a scattering position, with its unweighted amplitude.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScatteringOccupant {
    /// Canonical element symbol.
    pub element: String,
    /// Occupancy of this occupant at the position.
    pub occupancy: f64,
    /// Amplitude at the requested `s`, in electrons (X-ray) or Angstrom.
    pub amplitude: f64,
}

/// Coherent and incoherent weights of one disordered position.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DisorderScatterer {
    /// Fractional coordinates of the position.
    pub frac: [f64; 3],
    /// Occupants sorted by element symbol.
    pub occupants: Vec<ScatteringOccupant>,
    /// `sum_i occ_i f_i(s)`.
    pub coherent_amplitude: f64,
    /// `sum_i occ_i f_i(s)^2 - coherent_amplitude^2`.
    pub incoherent_intensity: f64,
}

/// Scattering weights of every position at one `s = sin(theta) / lambda`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DisorderScatteringWeights {
    /// Probe radiation kind.
    pub kind: RadiationKind,
    /// `s` in reciprocal Angstrom.
    pub s: f64,
    /// One entry per distinct fractional position, sorted by coordinate.
    pub scatterers: Vec<DisorderScatterer>,
    /// Sum of the position occupancies.
    pub total_occupancy: f64,
}

/// Builds the disorder-aware PXRD weighting hook.
///
/// The amplitude is the neutral-atom X-ray form factor, the coherent neutron
/// scattering length or the Mott-Bethe electron factor. Occupancies multiply
/// the amplitude; the incoherent term is the occupancy-weighted variance, a
/// measure of substitution contrast. The built-in powder engine already folds
/// the same coherent term into `F(hkl)`; this hook is for custom engines and
/// for reporting the neutron incoherent part.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a negative or non-finite `s`, an unknown
/// element or a missing neutron scattering length, or the errors of the P1
/// expansion.
pub fn disorder_scattering_weights(
    structure: &Structure,
    kind: RadiationKind,
    s: f64,
) -> Result<DisorderScatteringWeights, Error> {
    if !s.is_finite() || s < 0.0 {
        return Err(Error::InvalidInput(format!(
            "s must be non-negative and finite, got {s}"
        )));
    }
    let p1 = structure.p1()?;
    let mut positions: BTreeMap<String, PositionAccumulator> = BTreeMap::new();
    for atom in &p1.atoms {
        positions
            .entry(position_key(atom.frac))
            .or_insert_with(|| PositionAccumulator::new(atom.frac))
            .add(&atom.element, atom.occupancy);
    }
    let mut scatterers = Vec::with_capacity(positions.len());
    let mut total_occupancy = 0.0;
    for accumulator in positions.into_values() {
        let position = accumulator.finish();
        let mut occupants = Vec::with_capacity(position.occupants.len());
        let mut coherent = 0.0;
        let mut squared = 0.0;
        for occupant in &position.occupants {
            let amplitude = match kind {
                RadiationKind::XRay => xray_form_factor(&occupant.element, s)?,
                RadiationKind::Neutron => neutron_scattering_length(&occupant.element)?,
                RadiationKind::Electron => electron_scattering_factor(&occupant.element, s)?,
            };
            coherent += occupant.occupancy * amplitude;
            squared += occupant.occupancy * amplitude * amplitude;
            occupants.push(ScatteringOccupant {
                element: occupant.element.clone(),
                occupancy: occupant.occupancy,
                amplitude,
            });
        }
        let variance = squared - coherent * coherent;
        total_occupancy += position.total_occupancy();
        scatterers.push(DisorderScatterer {
            frac: position.frac,
            occupants,
            coherent_amplitude: coherent,
            incoherent_intensity: variance.max(0.0),
        });
    }
    Ok(DisorderScatteringWeights {
        kind,
        s,
        scatterers,
        total_occupancy,
    })
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum GroupKey {
    Explicit(u32),
    Implicit(String),
}

#[derive(Debug, Clone)]
struct PositionAccumulator {
    frac: [f64; 3],
    occupants: BTreeMap<String, f64>,
}

impl PositionAccumulator {
    fn new(frac: [f64; 3]) -> Self {
        Self {
            frac,
            occupants: BTreeMap::new(),
        }
    }

    fn add(&mut self, element: &str, occupancy: f64) {
        *self.occupants.entry(element.to_owned()).or_default() += occupancy;
    }

    fn finish(self) -> DisorderPosition {
        DisorderPosition {
            frac: self.frac,
            occupants: self
                .occupants
                .into_iter()
                .map(|(element, occupancy)| DisorderOccupant { element, occupancy })
                .collect(),
        }
    }
}

fn position_key(frac: [f64; 3]) -> String {
    let quantized = frac.map(|value| {
        let wrapped = value - value.floor();
        let wrapped = if wrapped >= 1.0 { 0.0 } else { wrapped };
        (wrapped * 100_000.0).round() as i64
    });
    let [first, second, third] = quantized;
    format!("{first},{second},{third}")
}
