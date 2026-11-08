//! Commands whose behaviour belongs to later phases (Task 3.10.7).
//!
//! The catalogue is complete here so the document can record every edit even
//! when the heavy algorithm lives in another crate: FF charges and hydrogen
//! placement come from Phase 4, interpenetration from Phase 5 and connection
//! points from Phase 6. Each wrapper stores plain data and keeps the command
//! undoable through the shared [`Command`] trait.

use xtalforge_core::{Error, Site, Structure};

use crate::command::Command;

/// Writes per-site charges onto the asymmetric unit (Task 4.3.4).
#[derive(Debug, Clone, PartialEq)]
pub struct SetCharges {
    charges: Vec<Option<f64>>,
}

impl SetCharges {
    /// Creates the command from one entry per asymmetric site.
    #[must_use]
    pub fn new(charges: Vec<Option<f64>>) -> Self {
        Self { charges }
    }
}

impl Command for SetCharges {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        if self.charges.len() != structure.sites().len() {
            return Err(Error::InvalidInput(format!(
                "expected {} charges for {} sites",
                structure.sites().len(),
                self.charges.len()
            )));
        }
        let sites: Vec<Site> = structure
            .sites()
            .iter()
            .zip(self.charges.iter())
            .map(|(site, charge)| {
                let mut site = site.clone();
                if let Some(value) = charge {
                    site = site.with_charge(*value)?;
                }
                Ok(site)
            })
            .collect::<Result<_, Error>>()?;
        rebuild(structure, sites)
    }

    fn description(&self) -> String {
        "Set atomic charges".to_owned()
    }
}

/// Duplicates every site at a fractional shift (Task 5.8.2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interpenetrate {
    shift: [f64; 3],
}

impl Interpenetrate {
    /// Creates the command from the fractional lattice shift.
    #[must_use]
    pub fn new(shift: [f64; 3]) -> Self {
        Self { shift }
    }
}

impl Command for Interpenetrate {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        if !self.shift.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "interpenetration shift must be finite".to_owned(),
            ));
        }
        let original = structure.sites().to_vec();
        let mut sites = original.clone();
        for site in &original {
            let frac = site.frac();
            let label = format!("{}'", site.label());
            sites.push(
                Site::new(
                    label,
                    site.element(),
                    [
                        (frac[0] + self.shift[0]).rem_euclid(1.0),
                        (frac[1] + self.shift[1]).rem_euclid(1.0),
                        (frac[2] + self.shift[2]).rem_euclid(1.0),
                    ],
                )?
                .with_occupancy(site.occupancy())?,
            );
        }
        rebuild(structure, sites)
    }

    fn description(&self) -> String {
        "Interpenetrate the structure".to_owned()
    }
}

/// Adds dummy `X` connection points (Task 6.13.4).
#[derive(Debug, Clone, PartialEq)]
pub struct MarkConnectionPoints {
    points: Vec<[f64; 3]>,
}

impl MarkConnectionPoints {
    /// Creates the command from fractional connection points.
    #[must_use]
    pub fn new(points: Vec<[f64; 3]>) -> Self {
        Self { points }
    }
}

impl Command for MarkConnectionPoints {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        if self.points.is_empty() {
            return Err(Error::InvalidInput(
                "at least one connection point is required".to_owned(),
            ));
        }
        let mut sites = structure.sites().to_vec();
        for (index, point) in self.points.iter().enumerate() {
            if !point.iter().all(|value| value.is_finite()) {
                return Err(Error::InvalidInput(
                    "connection points must be finite".to_owned(),
                ));
            }
            sites.push(Site::new(format!("X{}", index + 1), "X", *point)?);
        }
        rebuild(structure, sites)
    }

    fn description(&self) -> String {
        format!("Mark {} connection points", self.points.len())
    }
}

fn rebuild(structure: &Structure, sites: Vec<Site>) -> Result<Structure, Error> {
    let mut builder = xtalforge_core::StructureBuilder::new(structure.name(), *structure.lattice())
        .hall(structure.hall().to_owned());
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for site in sites {
        builder = builder.site(site);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}
