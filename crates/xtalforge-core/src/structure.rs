//! Crystal structures: a lattice, a space-group identity and asymmetric sites.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::Error;
use crate::Lattice;
use crate::bonding::{ExplicitBond, TopologyBond};
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
    #[serde(default)]
    explicit_bonds: Vec<ExplicitBond>,
    #[serde(default)]
    topology_bonds: Vec<TopologyBond>,
    metadata: BTreeMap<String, Value>,
}

mod properties;

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

    /// User bonds and removals applied after perception.
    #[must_use]
    pub fn explicit_bonds(&self) -> &[ExplicitBond] {
        &self.explicit_bonds
    }

    /// Topology edges, kept separate from chemical bonds.
    #[must_use]
    pub fn topology_bonds(&self) -> &[TopologyBond] {
        &self.topology_bonds
    }

    /// Open structure metadata.
    #[must_use]
    pub fn metadata(&self) -> &BTreeMap<String, Value> {
        &self.metadata
    }

    /// Whether the structure is P1: the Hall symbol is `P 1` (or empty) and
    /// every operation is the identity. A CIF may carry an explicit operation
    /// loop without a space-group symbol; such a file is not P1 even though
    /// its Hall symbol defaults to `P 1`.
    #[must_use]
    pub fn is_p1(&self) -> bool {
        let hall = self.hall.trim();
        (hall.is_empty() || hall == "P 1" || hall == "P1")
            && self.operations.iter().all(SymOp::is_identity)
    }

    /// Cartesian coordinates of the asymmetric sites.
    #[must_use]
    pub fn cart_sites(&self) -> Vec<[f64; 3]> {
        self.sites
            .iter()
            .map(|site| self.lattice.to_cart(site.frac()))
            .collect()
    }

    /// Structured warnings collected from the structure metadata.
    #[must_use]
    pub fn warnings(&self) -> Vec<StructureWarning> {
        let mut warnings = Vec::new();
        for (key, value) in &self.metadata {
            let Some(message) = value.as_str() else {
                continue;
            };
            let code = key
                .strip_suffix("_warning")
                .unwrap_or(key.as_str())
                .to_owned();
            let hint = if code.contains("symmetry") {
                Some("use --spacegroup SYMBOL to override the setting".to_owned())
            } else if code.contains("coincident") {
                Some("remove one of the overlapping atoms or lower --symprec".to_owned())
            } else {
                None
            };
            warnings.push(StructureWarning {
                code,
                message: message.to_owned(),
                atoms: Vec::new(),
                fields: vec![key.clone()],
                hint,
            });
        }
        warnings
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
            warnings: self.warnings(),
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
/// Accepts a bare structure object or the versioned `StructureDocument`
/// envelope that the CLI writes with `--json`. CIF text is read with the `io`
/// crate, not here.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid JSON or a structure that cannot
/// produce a summary.
pub fn structure_info(text: &str) -> Result<StructureSummary, Error> {
    let structure: Structure = if text.contains("\"schema_version\"") {
        StructureDocument::from_json(text)?.structure
    } else {
        serde_json::from_str(text)
            .map_err(|error| Error::InvalidInput(format!("invalid structure JSON: {error}")))?
    };
    structure.summary()
}

/// One structured warning attached to a structure.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct StructureWarning {
    /// Stable snake_case code.
    pub code: String,
    /// Human-readable message.
    pub message: String,
    /// P1 atom indices involved, when known.
    pub atoms: Vec<usize>,
    /// Fields involved, when known.
    pub fields: Vec<String>,
    /// Suggested fix, when one exists.
    pub hint: Option<String>,
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
    /// Structured reader and geometry warnings.
    pub warnings: Vec<StructureWarning>,
}

/// Builder for a [`Structure`].
#[derive(Debug, Clone)]
pub struct StructureBuilder {
    name: String,
    lattice: Lattice,
    hall: String,
    sites: Vec<Site>,
    operations: Vec<SymOp>,
    explicit_bonds: Vec<ExplicitBond>,
    topology_bonds: Vec<TopologyBond>,
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
            explicit_bonds: Vec::new(),
            topology_bonds: Vec::new(),
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

    /// Adds an explicit bond or removal.
    #[must_use]
    pub fn explicit_bond(mut self, bond: ExplicitBond) -> Self {
        self.explicit_bonds.push(bond);
        self
    }

    /// Adds a topology edge.
    #[must_use]
    pub fn topology_bond(mut self, bond: TopologyBond) -> Self {
        self.topology_bonds.push(bond);
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
            explicit_bonds: self.explicit_bonds,
            topology_bonds: self.topology_bonds,
            metadata: self.metadata,
        })
    }
}
