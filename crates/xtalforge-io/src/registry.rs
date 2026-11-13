//! Format registry, capabilities and loss warnings (Task 3.1).
//!
//! One list drives the CLI, the bindings, conversion and any future file
//! dialog. A format declares its name, extensions, optional reader/writer and
//! what it preserves, so a conversion can report the data it loses.

use xtalforge_core::{Error, Structure};

mod formats;

pub use formats::FORMATS;

/// What a format preserves.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Capabilities {
    /// Cell parameters.
    pub cell: bool,
    /// Space-group identity and operations.
    pub symmetry: bool,
    /// Partial occupancies.
    pub occupancy: bool,
    /// Anisotropic displacement parameters.
    pub anisotropic_adp: bool,
    /// Formal charges.
    pub charges: bool,
    /// Explicit bonds.
    pub bonds: bool,
    /// Site labels.
    pub labels: bool,
    /// Multiple frames.
    pub frames: bool,
}

/// A structure reader.
pub type Reader = fn(&str) -> Result<Structure, Error>;
/// A structure writer.
pub type Writer = fn(&Structure) -> Result<String, Error>;

/// One registered format.
#[derive(Clone, Copy, Debug)]
pub struct FormatSpec {
    /// Stable format name.
    pub name: &'static str,
    /// Lower-case extensions, with the leading dot.
    pub extensions: &'static [&'static str],
    /// Reader, when the format can be read.
    pub read: Option<Reader>,
    /// Writer, when the format can be written.
    pub write: Option<Writer>,
    /// Preserved features.
    pub capabilities: Capabilities,
}

/// Finds a format by name or extension (case-insensitive).
#[must_use]
pub fn find(name_or_extension: &str) -> Option<&'static FormatSpec> {
    let wanted = name_or_extension.trim().to_ascii_lowercase();
    let wanted = if wanted.starts_with('.') {
        wanted
    } else {
        format!(".{wanted}")
    };
    FORMATS.iter().find(|format| {
        format.name == wanted.trim_start_matches('.')
            || format
                .extensions
                .iter()
                .any(|extension| *extension == wanted)
    })
}

/// Reads text with an explicit format name or extension.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown format or a reader failure.
pub fn read_with(text: &str, hint: &str) -> Result<Structure, Error> {
    let format = find(hint)
        .ok_or_else(|| Error::InvalidInput(format!("unknown structure format: {hint}")))?;
    let reader = format
        .read
        .ok_or_else(|| Error::InvalidInput(format!("format {} is write-only", format.name)))?;
    reader(text)
}

/// Writes a structure with an explicit format name or extension.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown format or a writer failure.
pub fn write_with(structure: &Structure, hint: &str) -> Result<String, Error> {
    let format = find(hint)
        .ok_or_else(|| Error::InvalidInput(format!("unknown structure format: {hint}")))?;
    let writer = format
        .write
        .ok_or_else(|| Error::InvalidInput(format!("format {} is read-only", format.name)))?;
    writer(structure)
}

/// Detects the format from the content.
#[must_use]
pub fn sniff(text: &str) -> Option<&'static FormatSpec> {
    let trimmed = text.trim_start_matches('\u{feff}').trim_start();
    if text.contains("_atom_site.") || text.contains("_cell.") {
        return find("mmcif");
    }
    if trimmed.starts_with("data_")
        || text.contains("_atom_site_")
        || text.contains("_cell_length_")
    {
        return find("cif");
    }
    if trimmed.starts_with('{') && text.contains("@module") {
        return find("pmg-json");
    }
    if trimmed.starts_with("@<TRIPOS>") || trimmed.starts_with("@<tripos>") {
        return find("mol2");
    }
    if text.contains("<modeling") || text.contains("<vasprun") {
        return find("vasprun");
    }
    if (text.contains("TITEL") && text.contains("ions per type"))
        || (text.contains("POSITION") && text.contains("TOTAL-FORCE"))
    {
        return find("outcar");
    }
    if trimmed.starts_with("%BLOCK LATTICE") || trimmed.starts_with("%block lattice") {
        return find("castep");
    }
    if trimmed.starts_with("lattice_vector")
        || trimmed.starts_with("atom ")
        || trimmed.starts_with("atom_frac ")
    {
        return find("fhi-aims");
    }
    if text.contains("&FORCE_EVAL") || text.contains("&SUBSYS") {
        return find("cp2k");
    }
    if text.contains("&CONTROL") || text.contains("ATOMIC_POSITIONS") {
        return find("pwscf");
    }
    if text.contains("M  END") && (text.contains("V2000") || text.contains("V3000")) {
        return find("sdf");
    }
    if text.contains("natom")
        && (text.contains("xred") || text.contains("xcart"))
        && text.contains("acell")
    {
        return find("abinit");
    }
    let record = trimmed.get(0..6).unwrap_or("").trim();
    if matches!(
        record,
        "ATOM" | "HETATM" | "CRYST1" | "HEADER" | "TITLE" | "MODEL"
    ) {
        return find("pdb");
    }
    let mut lines = trimmed.lines().filter(|line| !line.trim().is_empty());
    let first = lines.next().unwrap_or("").trim();
    let second = lines.next().unwrap_or("").trim();
    if first.parse::<usize>().is_ok() && !second.is_empty() {
        if second.contains("Properties=") {
            return find("extxyz");
        }
        return find("xyz");
    }
    None
}

/// Names of the features a conversion loses.
#[must_use]
pub fn losses(from: Capabilities, to: Capabilities) -> Vec<&'static str> {
    let mut lost = Vec::new();
    for (kept, name) in [
        (from.symmetry && !to.symmetry, "space-group symmetry"),
        (from.occupancy && !to.occupancy, "partial occupancy"),
        (
            from.anisotropic_adp && !to.anisotropic_adp,
            "anisotropic ADPs",
        ),
        (from.charges && !to.charges, "formal charges"),
        (from.bonds && !to.bonds, "explicit bonds"),
        (from.labels && !to.labels, "site labels"),
        (from.frames && !to.frames, "multiple frames"),
    ] {
        if kept {
            lost.push(name);
        }
    }
    lost
}

/// What a file leaving the library should contain (Task 3.1.4).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExportPolicy {
    /// Write the full P1 cell instead of the asymmetric unit.
    pub expand_to_p1: bool,
    /// Keep hydrogen sites.
    pub keep_hydrogens: bool,
    /// Keep the dummy `X` connection points.
    pub keep_dummy: bool,
    /// Keep explicit bonds and removals.
    pub keep_explicit_bonds: bool,
}

impl Default for ExportPolicy {
    fn default() -> Self {
        Self {
            expand_to_p1: false,
            keep_hydrogens: true,
            keep_dummy: true,
            keep_explicit_bonds: true,
        }
    }
}

impl ExportPolicy {
    /// Applies the policy to a structure.
    ///
    /// # Errors
    ///
    /// Returns the errors of the P1 expansion.
    pub fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        let mut source = if self.expand_to_p1 {
            let p1 = structure.p1()?;
            let mut builder = xtalforge_core::StructureBuilder::new(structure.name(), p1.lattice);
            for atom in &p1.atoms {
                builder = builder.site(
                    xtalforge_core::Site::new(atom.label.clone(), &atom.element, atom.frac)?
                        .with_occupancy(atom.occupancy)?,
                );
            }
            builder.build()?
        } else {
            structure.clone()
        };
        if !self.keep_hydrogens || !self.keep_dummy || !self.keep_explicit_bonds {
            let sites: Vec<xtalforge_core::Site> = source
                .sites()
                .iter()
                .filter(|site| {
                    (self.keep_hydrogens || !matches!(site.element(), "H" | "D"))
                        && (self.keep_dummy || site.element() != "X")
                })
                .cloned()
                .collect();
            let mut builder =
                xtalforge_core::StructureBuilder::new(source.name(), *source.lattice())
                    .hall(source.hall().to_owned());
            for operation in source.operations() {
                builder = builder.operation(*operation);
            }
            for site in sites {
                builder = builder.site(site);
            }
            if self.keep_explicit_bonds {
                for bond in source.explicit_bonds() {
                    builder = builder.explicit_bond(*bond);
                }
            }
            for (key, value) in source.metadata() {
                builder = builder.metadata(key.clone(), value.clone());
            }
            source = builder.build()?;
        }
        Ok(source)
    }
}

/// A runtime format registry for plugins (Task 8.3.1).
///
/// [`FORMATS`] is the built-in list compiled into the crate. A front end or a
/// plugin can copy it and add formats without editing the library:
///
/// ```
/// use xtalforge_io::registry::{FormatRegistry, FormatSpec, Capabilities};
///
/// let mut registry = FormatRegistry::new();
/// registry.register(FormatSpec {
///     name: "twoline",
///     extensions: &[".2l"],
///     read: Some(|text| xtalforge_io::read_xyz(text)),
///     write: None,
///     capabilities: Capabilities::default(),
/// })
/// .unwrap();
/// assert!(registry.find(".2l").is_some());
/// ```
#[derive(Debug, Clone)]
pub struct FormatRegistry {
    formats: Vec<FormatSpec>,
}

impl Default for FormatRegistry {
    fn default() -> Self {
        Self::new()
    }
}

impl FormatRegistry {
    /// A registry with the built-in formats.
    #[must_use]
    pub fn new() -> Self {
        Self {
            formats: FORMATS.to_vec(),
        }
    }

    /// An empty registry.
    #[must_use]
    pub fn empty() -> Self {
        Self {
            formats: Vec::new(),
        }
    }

    /// Registers a format.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty name, a duplicate name or
    /// an extension already claimed by another format.
    pub fn register(&mut self, spec: FormatSpec) -> Result<(), Error> {
        if spec.name.trim().is_empty() {
            return Err(Error::InvalidInput(
                "a format needs a non-empty name".to_owned(),
            ));
        }
        if spec.extensions.is_empty() {
            return Err(Error::InvalidInput(format!(
                "format {} needs at least one extension",
                spec.name
            )));
        }
        if self.formats.iter().any(|format| format.name == spec.name) {
            return Err(Error::InvalidInput(format!(
                "format {} is already registered",
                spec.name
            )));
        }
        for extension in spec.extensions {
            if self
                .formats
                .iter()
                .any(|format| format.extensions.contains(extension))
            {
                return Err(Error::InvalidInput(format!(
                    "extension {extension} is already registered"
                )));
            }
        }
        self.formats.push(spec);
        Ok(())
    }

    /// Finds a format by name or extension (case-insensitive).
    #[must_use]
    pub fn find(&self, name_or_extension: &str) -> Option<&FormatSpec> {
        let wanted = name_or_extension.trim().to_ascii_lowercase();
        let wanted = if wanted.starts_with('.') {
            wanted
        } else {
            format!(".{wanted}")
        };
        self.formats.iter().find(|format| {
            format.name == wanted.trim_start_matches('.')
                || format
                    .extensions
                    .iter()
                    .any(|extension| *extension == wanted)
        })
    }

    /// Every registered format, built-ins first.
    #[must_use]
    pub fn formats(&self) -> &[FormatSpec] {
        &self.formats
    }

    /// Reads text with an explicit format name or extension.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an unknown or write-only format,
    /// or the reader error.
    pub fn read_with(&self, text: &str, hint: &str) -> Result<Structure, Error> {
        let format = self
            .find(hint)
            .ok_or_else(|| Error::InvalidInput(format!("unknown structure format: {hint}")))?;
        let reader = format
            .read
            .ok_or_else(|| Error::InvalidInput(format!("format {} is write-only", format.name)))?;
        reader(text)
    }

    /// Writes a structure with an explicit format name or extension.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an unknown or read-only format, or
    /// the writer error.
    pub fn write_with(&self, structure: &Structure, hint: &str) -> Result<String, Error> {
        let format = self
            .find(hint)
            .ok_or_else(|| Error::InvalidInput(format!("unknown structure format: {hint}")))?;
        let writer = format
            .write
            .ok_or_else(|| Error::InvalidInput(format!("format {} is read-only", format.name)))?;
        writer(structure)
    }
}
