//! Crystal structures: a lattice, a space-group identity and asymmetric sites.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::Error;
use crate::Lattice;
use crate::elements::element_by_symbol;
use crate::site::Site;
use crate::symop::SymOp;

/// A crystal structure.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Structure {
    name: String,
    lattice: Lattice,
    hall: String,
    sites: Vec<Site>,
    #[serde(default)]
    operations: Vec<SymOp>,
    metadata: BTreeMap<String, Value>,
}

impl Structure {
    /// Structure name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Unit cell.
    #[must_use]
    pub fn lattice(&self) -> &Lattice {
        &self.lattice
    }

    /// Hall symbol of the space group; `P 1` for a P1 structure.
    #[must_use]
    pub fn hall(&self) -> &str {
        &self.hall
    }

    /// Asymmetric-unit sites.
    #[must_use]
    pub fn sites(&self) -> &[Site] {
        &self.sites
    }

    /// Explicit symmetry operations; empty means identity only.
    #[must_use]
    pub fn operations(&self) -> &[SymOp] {
        &self.operations
    }

    /// Open structure metadata.
    #[must_use]
    pub fn metadata(&self) -> &BTreeMap<String, Value> {
        &self.metadata
    }

    /// Whether the structure is P1.
    #[must_use]
    pub fn is_p1(&self) -> bool {
        let hall = self.hall.trim();
        hall.is_empty() || hall == "P 1" || hall == "P1"
    }

    /// Cartesian coordinates of the asymmetric sites.
    #[must_use]
    pub fn cart_sites(&self) -> Vec<[f64; 3]> {
        self.sites
            .iter()
            .map(|site| self.lattice.to_cart(site.frac()))
            .collect()
    }

    /// Cell contents as `(element, occupancy sum)` in first-appearance order.
    #[must_use]
    pub fn cell_contents(&self) -> Vec<(String, f64)> {
        let mut contents: Vec<(String, f64)> = Vec::new();
        for site in &self.sites {
            if let Some(entry) = contents
                .iter_mut()
                .find(|(element, _)| element == site.element())
            {
                entry.1 += site.occupancy();
            } else {
                contents.push((site.element().to_owned(), site.occupancy()));
            }
        }
        contents
    }

    /// Formula string, for example `Fe2 O3` or `C0.5`.
    #[must_use]
    pub fn formula(&self) -> String {
        self.cell_contents()
            .iter()
            .map(|(element, count)| {
                if (*count - count.round()).abs() < 1e-9 {
                    let count = count.round() as u64;
                    if count == 1 {
                        element.clone()
                    } else {
                        format!("{element}{count}")
                    }
                } else {
                    format!("{element}{count:.3}")
                }
            })
            .collect::<Vec<_>>()
            .join(" ")
    }

    /// Mass of the cell contents in unified atomic mass units.
    ///
    /// # Errors
    ///
    /// Returns an error if a site element is missing from the element table.
    pub fn cell_mass(&self) -> Result<f64, Error> {
        let mut mass = 0.0;
        for site in &self.sites {
            mass += site.occupancy() * element_by_symbol(site.element())?.mass;
        }
        Ok(mass)
    }

    /// Molar mass of the cell contents in g/mol.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`Structure::cell_mass`].
    pub fn molar_mass(&self) -> Result<f64, Error> {
        self.cell_mass()
    }

    /// Density in g/cm3, assuming the sites are the whole cell (P1).
    ///
    /// # Errors
    ///
    /// Returns the errors of [`Structure::cell_mass`].
    pub fn density(&self) -> Result<f64, Error> {
        Ok(self.cell_mass()? * 1.660_539_066_60 / self.lattice.volume())
    }

    /// Number of formula units in the cell.
    ///
    /// The current model stores the cell contents in the asymmetric unit, so a
    /// P1 structure has one formula unit. The P1 expansion in Task 1.6 refines
    /// this when the asymmetric unit is not the cell.
    #[must_use]
    pub fn z(&self) -> u32 {
        1
    }

    /// Sum of `occupancy * charge` over the sites.
    #[must_use]
    pub fn total_charge(&self) -> f64 {
        self.sites
            .iter()
            .map(|site| site.occupancy() * site.charge().unwrap_or(0.0))
            .sum()
    }

    /// Charge balance, the same as [`Structure::total_charge`].
    #[must_use]
    pub fn charge_balance(&self) -> f64 {
        self.total_charge()
    }

    /// Summary used by the CLI, Python and WASM `info` surfaces.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`Structure::density`].
    pub fn summary(&self) -> Result<StructureSummary, Error> {
        Ok(StructureSummary {
            name: self.name.clone(),
            formula: self.formula(),
            hall: self.hall.clone(),
            n_sites: self.sites.len(),
            volume: self.lattice.volume(),
            density: self.density()?,
            total_charge: self.total_charge(),
            is_p1: self.is_p1(),
        })
    }

    /// Copy with a new name.
    #[must_use]
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = name.into();
        self
    }

    /// Copy with a new lattice.
    #[must_use]
    pub fn with_lattice(mut self, lattice: Lattice) -> Self {
        self.lattice = lattice;
        self
    }

    /// Copy with one more site.
    #[must_use]
    pub fn with_site(mut self, site: Site) -> Self {
        self.sites.push(site);
        self
    }

    /// A non-periodic molecule in a padded P1 box.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for no points, a non-positive padding,
    /// or a point whose element is unknown.
    pub fn molecule(
        name: impl Into<String>,
        points: &[(String, [f64; 3])],
        padding: f64,
    ) -> Result<Self, Error> {
        if points.is_empty() {
            return Err(Error::InvalidInput(
                "a molecule needs at least one atom".to_owned(),
            ));
        }
        if !padding.is_finite() || padding <= 0.0 {
            return Err(Error::InvalidInput(
                "molecule padding must be positive and finite".to_owned(),
            ));
        }
        let [mut min_x, mut min_y, mut min_z] = [f64::INFINITY; 3];
        let [mut max_x, mut max_y, mut max_z] = [f64::NEG_INFINITY; 3];
        for (_, point) in points {
            let [x, y, z] = *point;
            min_x = min_x.min(x);
            min_y = min_y.min(y);
            min_z = min_z.min(z);
            max_x = max_x.max(x);
            max_y = max_y.max(y);
            max_z = max_z.max(z);
        }
        let [length_x, length_y, length_z] = [
            max_x - min_x + 2.0 * padding,
            max_y - min_y + 2.0 * padding,
            max_z - min_z + 2.0 * padding,
        ];
        let lattice = Lattice::orthorhombic(length_x, length_y, length_z)?;
        let [origin_x, origin_y, origin_z] = [min_x - padding, min_y - padding, min_z - padding];
        let mut builder = StructureBuilder::new(name, lattice);
        for (index, (element, point)) in points.iter().enumerate() {
            let [x, y, z] = *point;
            let frac = [
                (x - origin_x) / length_x,
                (y - origin_y) / length_y,
                (z - origin_z) / length_z,
            ];
            builder = builder.site(Site::new(
                format!("{}{}", element, index + 1),
                element,
                frac,
            )?);
        }
        builder.build()
    }
}

/// Version of the serialized structure schema.
pub const SCHEMA_VERSION: u32 = 1;

/// Versioned JSON envelope for a structure.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StructureDocument {
    /// Schema version, independent of the crate version.
    pub schema_version: u32,
    /// The structure payload.
    pub structure: Structure,
}

impl StructureDocument {
    /// Wraps a structure in the current schema version.
    #[must_use]
    pub fn new(structure: Structure) -> Self {
        Self {
            schema_version: SCHEMA_VERSION,
            structure,
        }
    }

    /// Serializes the document to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] if serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self)
            .map_err(|error| Error::InvalidInput(format!("could not serialize structure: {error}")))
    }

    /// Parses a document and checks its schema version.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for invalid JSON, unknown fields or an
    /// unsupported schema version.
    pub fn from_json(text: &str) -> Result<Self, Error> {
        let document: Self = serde_json::from_str(text)
            .map_err(|error| Error::InvalidInput(format!("invalid structure JSON: {error}")))?;
        if document.schema_version != SCHEMA_VERSION {
            return Err(Error::InvalidInput(format!(
                "unsupported structure schema_version {}, expected {SCHEMA_VERSION}",
                document.schema_version
            )));
        }
        Ok(document)
    }
}

impl Structure {
    /// Wraps this structure in the versioned document envelope.
    #[must_use]
    pub fn to_document(&self) -> StructureDocument {
        StructureDocument::new(self.clone())
    }
}

/// Parses a structure JSON document and returns its summary.
///
/// This is the Phase 1 bridge until the CIF reader lands in Task 1.15.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid JSON or a structure that cannot
/// produce a summary.
pub fn structure_info(text: &str) -> Result<StructureSummary, Error> {
    let structure: Structure = serde_json::from_str(text)
        .map_err(|error| Error::InvalidInput(format!("invalid structure JSON: {error}")))?;
    structure.summary()
}

/// Summary of a structure for the surfaces.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct StructureSummary {
    /// Structure name.
    pub name: String,
    /// Formula string.
    pub formula: String,
    /// Hall symbol.
    pub hall: String,
    /// Number of asymmetric sites.
    pub n_sites: usize,
    /// Cell volume in cubic Angstrom.
    pub volume: f64,
    /// Density in g/cm3.
    pub density: f64,
    /// Total charge.
    pub total_charge: f64,
    /// Whether the structure is P1.
    pub is_p1: bool,
}

/// Builder for a [`Structure`].
#[derive(Debug, Clone)]
pub struct StructureBuilder {
    name: String,
    lattice: Lattice,
    hall: String,
    sites: Vec<Site>,
    operations: Vec<SymOp>,
    metadata: BTreeMap<String, Value>,
}

impl StructureBuilder {
    /// Starts a P1 structure.
    #[must_use]
    pub fn new(name: impl Into<String>, lattice: Lattice) -> Self {
        Self {
            name: name.into(),
            lattice,
            hall: "P 1".to_owned(),
            sites: Vec::new(),
            operations: Vec::new(),
            metadata: BTreeMap::new(),
        }
    }

    /// Sets the Hall symbol.
    #[must_use]
    pub fn hall(mut self, hall: impl Into<String>) -> Self {
        self.hall = hall.into();
        self
    }

    /// Adds a site.
    #[must_use]
    pub fn site(mut self, site: Site) -> Self {
        self.sites.push(site);
        self
    }

    /// Adds an explicit symmetry operation.
    #[must_use]
    pub fn operation(mut self, operation: SymOp) -> Self {
        self.operations.push(operation);
        self
    }

    /// Adds one metadata entry.
    #[must_use]
    pub fn metadata(mut self, key: impl Into<String>, value: Value) -> Self {
        self.metadata.insert(key.into(), value);
        self
    }

    /// Builds the structure.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the Hall symbol is empty after
    /// trimming.
    pub fn build(self) -> Result<Structure, Error> {
        if self.hall.trim().is_empty() {
            return Err(Error::InvalidInput(
                "space-group Hall symbol must not be empty; use P 1".to_owned(),
            ));
        }
        Ok(Structure {
            name: self.name,
            lattice: self.lattice,
            hall: self.hall,
            sites: self.sites,
            operations: self.operations,
            metadata: self.metadata,
        })
    }
}
