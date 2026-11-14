//! Named parser callbacks for compute steps (Task 18.2/18.4).

use std::collections::BTreeMap;

use crate::Error;
use crate::result::CalculationResult;

/// A parser callback: engine text to a typed result.
pub type ParserFn = fn(&str) -> Result<CalculationResult, Error>;

/// Named parser callbacks used by compute steps.
///
/// The registry is the serializable half of the integration point: a
/// [`super::ComputeStep`] stores the parser *name*, and the caller supplies the
/// registry (built-ins or custom callbacks) to [`super::run_job`]. This keeps
/// the step serializable for workflows on every surface.
#[derive(Debug, Clone, Default)]
pub struct ParserRegistry {
    parsers: BTreeMap<String, ParserFn>,
}

impl ParserRegistry {
    /// An empty registry.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// A registry with every parser of Task 18.2 registered.
    #[must_use]
    pub fn with_builtins() -> Self {
        let mut registry = Self::new();
        registry.register("vasp_outcar", crate::parse::parse_vasp_outcar);
        registry.register("vasp_vasprun", crate::parse::parse_vasp_vasprun);
        registry.register("vasp_oszicar", crate::parse::parse_vasp_oszicar);
        registry.register("vasp_eigenval", crate::parse::parse_vasp_eigenval);
        registry.register("vasp_doscar", crate::parse::parse_vasp_doscar);
        registry.register("vasp_procar", crate::parse::parse_vasp_procar);
        registry.register("qe", crate::parse::parse_qe);
        registry.register("qe_xml", crate::parse::parse_qe_xml);
        registry.register("qe_bands", crate::parse::parse_qe_bands);
        registry.register("qe_projwfc", crate::parse::parse_qe_projwfc);
        registry.register("cp2k", crate::parse::parse_cp2k);
        registry.register("lammps", crate::parse::parse_lammps);
        registry.register("lammps_log", crate::parse::parse_lammps_log);
        registry.register("lammps_dump", crate::parse::parse_lammps_dump);
        registry.register("gulp", crate::parse::parse_gulp);
        registry
    }

    /// Registers a parser, returning the callback it replaced.
    pub fn register(&mut self, name: impl Into<String>, parser: ParserFn) -> Option<ParserFn> {
        self.parsers.insert(name.into(), parser)
    }

    /// Looks up a parser by name.
    #[must_use]
    pub fn parser(&self, name: &str) -> Option<ParserFn> {
        self.parsers.get(name).copied()
    }

    /// Registered parser names in sorted order.
    #[must_use]
    pub fn names(&self) -> Vec<&str> {
        self.parsers.keys().map(String::as_str).collect()
    }

    /// Runs a registered parser.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Missing`] for an unregistered name and the parser's
    /// own error for malformed output.
    pub fn parse(&self, name: &str, text: &str) -> Result<CalculationResult, Error> {
        let parser = self
            .parser(name)
            .ok_or_else(|| Error::Missing(format!("parser {name:?} is not registered")))?;
        parser(text)
    }
}
