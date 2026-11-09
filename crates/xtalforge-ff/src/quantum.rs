//! In-process quantum methods (Task 4.10).
//!
//! The scientific backend is sci-form's GFN0, GFN1 and GFN2-xTB
//! implementation, behind the optional `sci-form` feature. sci-form is a
//! published MIT crate (2.0.0), so the adapter is real: [`adapter`] builds a
//! [`Calculator`] that evaluates the requested method.
//!
//! # Capability limits
//!
//! The sci-form GFN methods are **molecular**: they take element numbers and
//! Cartesian positions, with no lattice. The adapter therefore refuses a
//! structure that carries space-group symmetry, because the atoms of a
//! periodic crystal are not one isolated molecule; build a molecule with
//! [`xtalforge_core::Structure::molecule`] and the P1 cell is evaluated as an
//! isolated cluster. The cell itself is never used by the energy.
//!
//! GFN2 returns forces, computed by sci-form's central-difference gradient.
//! GFN0 and GFN1 return an energy only and refuse `energy_forces` rather than
//! report zero forces. None of the three returns an analytic stress.
//!
//! The methods are neutral-molecule solvers: the adapter passes charge `0`
//! to GFN0 and the GFN1/GFN2 entry points take no charge. Occupancies other
//! than one, dummy `X` atoms and periodic space groups are refused.

use serde::{Deserialize, Serialize};

use crate::Error;

use crate::calculator::Calculator;
use xtalforge_core::Structure;

#[cfg(feature = "sci-form")]
use xtalforge_core::elements::element_by_symbol;

/// One in-process semi-empirical method.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum QuantumMethod {
    /// GFN0-xTB.
    Gfn0,
    /// GFN1-xTB.
    Gfn1,
    /// GFN2-xTB.
    Gfn2,
}

impl QuantumMethod {
    /// Stable method name used by the CLI and bindings.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Gfn0 => "gfn0",
            Self::Gfn1 => "gfn1",
            Self::Gfn2 => "gfn2",
        }
    }
}

/// What one method can answer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct QuantumCapabilities {
    /// Accepts a periodic cell rather than a molecule only.
    pub periodic: bool,
    /// Returns forces.
    pub forces: bool,
    /// Returns an analytic stress tensor.
    pub stress: bool,
}

/// Capabilities of one method.
///
/// sci-form's GFN methods are molecular only. GFN2 provides forces through a
/// central-difference gradient; GFN0 and GFN1 provide an energy only. None
/// exposes an analytic stress.
#[must_use]
pub fn capabilities(method: QuantumMethod) -> QuantumCapabilities {
    QuantumCapabilities {
        periodic: false,
        forces: method == QuantumMethod::Gfn2,
        stress: false,
    }
}

/// Handle for the in-process calculator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QuantumAdapter {
    method: QuantumMethod,
}

impl Default for QuantumAdapter {
    fn default() -> Self {
        Self {
            method: QuantumMethod::Gfn2,
        }
    }
}

impl QuantumAdapter {
    /// Method this adapter evaluates.
    #[must_use]
    pub const fn method(&self) -> QuantumMethod {
        self.method
    }
}

/// Builds the adapter for one method.
///
/// # Errors
///
/// Without the `sci-form` feature returns
/// `Error::InvalidInput("sci-form support is not compiled; rebuild with
/// --features sci-form")`.
pub fn adapter(method: QuantumMethod) -> Result<QuantumAdapter, Error> {
    #[cfg(not(feature = "sci-form"))]
    {
        let _ = method;
        Err(Error::InvalidInput(
            "sci-form support is not compiled; rebuild with --features sci-form".to_owned(),
        ))
    }
    #[cfg(feature = "sci-form")]
    {
        Ok(QuantumAdapter { method })
    }
}

impl Calculator for QuantumAdapter {
    fn name(&self) -> &'static str {
        self.method.name()
    }

    fn energy(&self, structure: &Structure) -> Result<f64, Error> {
        #[cfg(not(feature = "sci-form"))]
        let result = {
            let _ = structure;
            Err(Error::InvalidInput(
                "sci-form support is not compiled; rebuild with --features sci-form".to_owned(),
            ))
        };
        #[cfg(feature = "sci-form")]
        let result = {
            let (elements, positions) = molecular_input(structure)?;
            match self.method {
                QuantumMethod::Gfn0 => sci_form::xtb::solve_gfn0(&elements, &positions, 0.0)
                    .map(|result| result.total_energy)
                    .map_err(|message| Error::Quantum(format!("gfn0: {message}"))),
                QuantumMethod::Gfn1 => sci_form::compute_gfn1(&elements, &positions)
                    .map(|result| result.total_energy)
                    .map_err(|error| Error::Quantum(format!("gfn1: {error}"))),
                QuantumMethod::Gfn2 => sci_form::compute_gfn2(&elements, &positions)
                    .map(|result| result.total_energy)
                    .map_err(|error| Error::Quantum(format!("gfn2: {error}"))),
            }
        };
        result
    }

    fn energy_forces(&self, structure: &Structure) -> Result<(f64, Vec<[f64; 3]>), Error> {
        #[cfg(not(feature = "sci-form"))]
        let result = {
            let _ = structure;
            Err(Error::InvalidInput(
                "sci-form support is not compiled; rebuild with --features sci-form".to_owned(),
            ))
        };
        #[cfg(feature = "sci-form")]
        let result = {
            if self.method != QuantumMethod::Gfn2 {
                return Err(Error::InvalidInput(format!(
                    "{} returns an energy only; GFN2 is the method with forces",
                    self.method.name()
                )));
            }
            let (elements, positions) = molecular_input(structure)?;
            let result = sci_form::compute_xtb_gradient(&elements, &positions)
                .map_err(|error| Error::Quantum(format!("gfn2 gradient: {error}")))?;
            let forces = result
                .gradients
                .iter()
                .map(|gradient| gradient.map(|value| -value))
                .collect();
            Ok((result.energy, forces))
        };
        result
    }

    fn provides_forces(&self) -> bool {
        self.method == QuantumMethod::Gfn2
    }

    fn provides_stress(&self) -> bool {
        false
    }
}

/// Converts a P1 molecular structure to sci-form's element and position arrays.
#[cfg(feature = "sci-form")]
fn molecular_input(structure: &Structure) -> Result<(Vec<u8>, Vec<[f64; 3]>), Error> {
    if !structure.is_p1() {
        return Err(Error::InvalidInput(
            "sci-form GFN methods are molecular only; build the molecule with \
             Structure::molecule (P1) instead of a space group"
                .to_owned(),
        ));
    }
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "sci-form GFN methods need at least one atom".to_owned(),
        ));
    }
    let mut elements = Vec::with_capacity(p1.atoms.len());
    let mut positions = Vec::with_capacity(p1.atoms.len());
    for atom in &p1.atoms {
        if (atom.occupancy - 1.0).abs() > 1e-8 {
            return Err(Error::InvalidInput(format!(
                "sci-form GFN methods need full occupancy; atom {} has {}",
                atom.label, atom.occupancy
            )));
        }
        let element = element_by_symbol(&atom.element)?;
        if element.z == 0 {
            return Err(Error::InvalidInput(format!(
                "sci-form GFN methods do not support the dummy species at atom {}",
                atom.label
            )));
        }
        elements.push(element.z);
        positions.push(p1.lattice.to_cart(atom.frac));
    }
    Ok((elements, positions))
}
