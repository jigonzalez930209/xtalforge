//! DFTB+ HSD input and `detailed.out` parsing (Task 4.9.2).

use serde::{Deserialize, Serialize};
use xtalforge_core::{Structure, element_info};

use super::text::{find_line, find_number, finite, numbers, push_line, validate_token};
use super::{
    EngineResult, HARTREE_PER_BOHR_TO_EV_PER_ANGSTROM, HARTREE_PER_BOHR3_TO_GPA, HARTREE_TO_EV,
};
use crate::Error;

/// DFTB+ input options.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DftbOptions {
    /// HSD Hamiltonian block name, normally `DFTB`.
    pub method: String,
    /// Maximum angular momentum written for every element (`s`, `p`, `d`...).
    pub max_angular_momentum: String,
    /// k-point mesh along the three reciprocal axes.
    pub kpoints: [u32; 3],
    /// Whether self-consistent charge (SCC) is enabled.
    pub scc: bool,
}

impl Default for DftbOptions {
    fn default() -> Self {
        Self {
            method: "DFTB".to_owned(),
            max_angular_momentum: "p".to_owned(),
            kpoints: [1, 1, 1],
            scc: true,
        }
    }
}

impl DftbOptions {
    /// Checks the options that end up as HSD tokens.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty or unsafe method name,
    /// empty angular momentum, or a zero k-point count.
    pub fn validate(&self) -> Result<(), Error> {
        validate_token("DFTB method", &self.method)?;
        validate_token("DFTB maximum angular momentum", &self.max_angular_momentum)?;
        if self.kpoints.contains(&0) {
            return Err(Error::InvalidInput(
                "DFTB k-point counts must all be positive".to_owned(),
            ));
        }
        Ok(())
    }
}

/// Writes a DFTB+ HSD input for `structure`.
///
/// The geometry is written as `GenFormat` type `S F` in fractional
/// coordinates; every element gets the configured maximum angular momentum,
/// and the Slater-Koster files use `Prefix = "slakos/"`,
/// `Separator = "-"`, `Suffix = ".skf"`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, an empty structure or
/// a dummy atom, and a core error when the P1 expansion or an element lookup
/// fails.
pub fn write_dftb_hsd(structure: &Structure, options: &DftbOptions) -> Result<String, Error> {
    options.validate()?;
    let cell = structure.p1()?;
    if cell.n_atoms() == 0 {
        return Err(Error::InvalidInput(
            "cannot write a DFTB+ input for an empty structure".to_owned(),
        ));
    }
    let mut types = std::collections::BTreeMap::new();
    for atom in &cell.atoms {
        if element_info(&atom.element)?.is_dummy() {
            return Err(Error::InvalidInput(
                "external engines do not accept the dummy X species".to_owned(),
            ));
        }
        let next = types.len() + 1;
        types.entry(atom.element.clone()).or_insert(next);
    }
    let mut out = String::new();
    push_line(
        &mut out,
        format!("# xtalforge DFTB+ input for {}", structure.name()),
    );
    push_line(&mut out, "Geometry = GenFormat {");
    push_line(&mut out, format!("  {}  S F", types.len()));
    let atom_types = cell
        .atoms
        .iter()
        .map(|atom| types.get(&atom.element).copied().unwrap_or(0).to_string())
        .collect::<Vec<_>>()
        .join(" ");
    push_line(&mut out, format!("  {}  {atom_types}", cell.n_atoms()));
    for (index, atom) in cell.atoms.iter().enumerate() {
        let atom_type = types.get(&atom.element).copied().unwrap_or(0);
        push_line(
            &mut out,
            format!(
                "  {}  {atom_type}  {:.10}  {:.10}  {:.10}",
                index + 1,
                atom.frac[0],
                atom.frac[1],
                atom.frac[2]
            ),
        );
    }
    push_line(&mut out, "}");
    push_line(&mut out, format!("Hamiltonian = {} {{", options.method));
    push_line(
        &mut out,
        format!("  SCC = {}", if options.scc { "Yes" } else { "No" }),
    );
    push_line(&mut out, "  MaxAngularMomentum = {");
    for element in types.keys() {
        push_line(
            &mut out,
            format!("    {element} = \"{}\"", options.max_angular_momentum),
        );
    }
    push_line(&mut out, "  }");
    push_line(&mut out, "  SlaterKosterFiles = Type2FileNames {");
    push_line(&mut out, "    Prefix = \"slakos/\"");
    push_line(&mut out, "    Separator = \"-\"");
    push_line(&mut out, "    Suffix = \".skf\"");
    push_line(&mut out, "  }");
    push_line(&mut out, "  KPointsAndWeights = SupercellFolding {");
    push_line(
        &mut out,
        format!(
            "    {} 0 0\n    0 {} 0\n    0 0 {}\n    0.5 0.5 0.5",
            options.kpoints[0], options.kpoints[1], options.kpoints[2]
        ),
    );
    push_line(&mut out, "  }");
    push_line(&mut out, "}");
    Ok(out)
}

/// Parses the parts of a DFTB+ `detailed.out` this crate consumes.
///
/// `Total energy:` is in Hartree, `Total Forces` in Hartree/Bohr and
/// `Total stress tensor` in GPa or Hartree/Bohr^3 as marked by its header.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the energy line is missing, when a
/// present stress block is incomplete, or when a parsed value is not finite.
pub fn parse_dftb_output(text: &str) -> Result<EngineResult, Error> {
    let lines: Vec<&str> = text.lines().collect();
    let energy = find_number(&lines, "Total energy:")
        .ok_or_else(|| Error::InvalidInput("DFTB output has no 'Total energy:' line".to_owned()))?
        * HARTREE_TO_EV;
    finite("DFTB energy", energy)?;
    let forces = if let Some(index) = find_line(&lines, "Total Forces") {
        super::text::parse_vectors(&lines, index + 1, HARTREE_PER_BOHR_TO_EV_PER_ANGSTROM)
    } else {
        Vec::new()
    };
    let stress = if let Some(index) = find_line(&lines, "Total stress tensor") {
        Some(parse_stress(&lines, index)?)
    } else {
        None
    };
    Ok(EngineResult {
        energy,
        forces,
        stress,
    })
}

fn parse_stress(lines: &[&str], index: usize) -> Result<[[f64; 3]; 3], Error> {
    let header = lines.get(index).copied().unwrap_or_default();
    let scale = if header.to_ascii_lowercase().contains("gpa") {
        1.0
    } else {
        HARTREE_PER_BOHR3_TO_GPA
    };
    let mut values = Vec::new();
    for line in lines.iter().skip(index + 1) {
        if line.trim().is_empty() && !values.is_empty() {
            break;
        }
        for value in numbers(line) {
            if values.len() < 9 {
                values.push(value);
            }
        }
        if values.len() >= 9 {
            break;
        }
    }
    if values.len() != 9 {
        return Err(Error::InvalidInput(format!(
            "DFTB stress tensor has {} numbers, expected 9",
            values.len()
        )));
    }
    let mut tensor = [[0.0; 3]; 3];
    for (row, values) in values.chunks_exact(3).enumerate() {
        for (column, value) in values.iter().enumerate() {
            if let Some(slot) = tensor.get_mut(row).and_then(|row| row.get_mut(column)) {
                *slot = value * scale;
            }
        }
    }
    finite("DFTB stress", tensor.iter().flatten().copied().sum())?;
    Ok(tensor)
}
