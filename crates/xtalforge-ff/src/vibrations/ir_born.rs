//! Born-effective-charge IR intensities (Task 19.4.1).
//!
//! [`ir_spectrum`](super::ir_spectrum) uses fixed scalar QEq charges; this
//! module adds the tensor path:
//!
//! ```text
//! p_a = sum_i sum_b Z*_{i,a b} u_{i,b},
//! I_m = C |p|^2 / sum_i m_i |u_i|^2.
//! ```
//!
//! `Z*` is the Born effective-charge tensor in units of the elementary charge,
//! `u` a Cartesian normal-mode displacement in Angstrom, `m` an atomic mass in
//! amu and `C` the same `(e / sqrt(amu))^2` conversion the fixed-charge path
//! uses, so both functions report km/mol on one scale.  The tensors can come
//! from a periodic engine through [`BornCharges::from_text`], or be
//! approximated from diagonal QEq or EQeq charges through
//! [`BornCharges::qeq`] and [`BornCharges::eqeq`].  The diagonal shortcut is
//! exactly the fixed-charge model; a full engine tensor also carries the charge
//! flow that the fixed-charge path misses.

use serde::{Deserialize, Serialize};
use xtalforge_core::{Structure, element_info};

use crate::Error;
use crate::charges::{ChargeOptions, Charges, eqeq_charges, qeq_charges};

use super::constants::E2_PER_AMU_TO_KM_PER_MOL;
use super::ir::{IrOptions, IrPeak, IrSpectrum, MAX_IR_POINTS};
use super::modes::VibrationalModes;
use super::spectrum::broaden;
use super::tensors::parse_tensor_file;

/// Born effective charges of one P1 cell.
///
/// `per_atom[i]` is the `3 x 3` tensor of atom `i`, indexed as
/// `tensor[polarisation axis][displacement axis]`, the Cartesian convention of
/// VASP and Quantum ESPRESSO:
///
/// ```text
/// Z*_{i,a b} = d p_a / d u_{i,b}
/// ```
///
/// in units of the elementary charge.  The tensors are supplied by an engine
/// or approximated from diagonal QEq or EQeq charges.  Unlike the fixed-charge
/// model, a full tensor carries the charge that flows when a bond stretches.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BornCharges {
    /// One `3 x 3` tensor per P1 atom.
    pub per_atom: Vec<[[f64; 3]; 3]>,
}

impl BornCharges {
    /// Builds Born charges from explicit per-atom tensors.
    #[must_use]
    pub fn from_tensors(per_atom: Vec<[[f64; 3]; 3]>) -> Self {
        Self { per_atom }
    }

    /// Builds isotropic tensors `q_i * I` from scalar charges.
    ///
    /// This is exactly the fixed-charge dipole model expressed as a Born
    /// tensor, so it is the natural bridge from QEq and EQeq.
    #[must_use]
    pub fn from_charges(charges: &Charges) -> Self {
        let mut per_atom = Vec::with_capacity(charges.per_atom.len());
        for charge in &charges.per_atom {
            let mut tensor = [[0.0; 3]; 3];
            for (axis, row) in tensor.iter_mut().enumerate() {
                if let Some(value) = row.get_mut(axis) {
                    *value = *charge;
                }
            }
            per_atom.push(tensor);
        }
        Self { per_atom }
    }

    /// Reads tensors from a simple text or JSON tensor file.
    ///
    /// The file formats are documented on [`parse_tensor_file`](super::parse_tensor_file).  A full
    /// `OUTCAR` or Quantum ESPRESSO output is not parsed: export the tensor
    /// block, one `3 x 3` tensor per atom in atom order.
    ///
    /// # Errors
    ///
    /// Returns the parser errors of [`parse_tensor_file`](super::parse_tensor_file).
    pub fn from_text(text: &str) -> Result<Self, Error> {
        Ok(Self {
            per_atom: parse_tensor_file(text)?.tensors,
        })
    }

    /// QEq charges as isotropic Born tensors, with default QEq options.
    ///
    /// The structure must be the P1 cell the modes describe; a non-P1
    /// structure expands to more atoms than the modes.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`qeq_charges`].
    pub fn qeq(structure: &Structure) -> Result<Self, Error> {
        let charges = qeq_charges(structure, &ChargeOptions::default())?;
        Ok(Self::from_charges(&charges))
    }

    /// EQeq charges as isotropic Born tensors, with default EQeq options.
    ///
    /// The structure must be the P1 cell the modes describe.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`eqeq_charges`].
    pub fn eqeq(structure: &Structure) -> Result<Self, Error> {
        let charges = eqeq_charges(structure, &ChargeOptions::eqeq())?;
        Ok(Self::from_charges(&charges))
    }

    /// Number of atoms.
    #[must_use]
    pub fn len(&self) -> usize {
        self.per_atom.len()
    }

    /// Whether no tensors are stored.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.per_atom.is_empty()
    }

    /// Checks the atom count and rejects non-finite components.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the tensor count differs from
    /// `n_atoms` or a component is not finite.
    pub fn validate(&self, n_atoms: usize) -> Result<(), Error> {
        if self.per_atom.len() != n_atoms {
            return Err(Error::InvalidInput(format!(
                "got {} Born tensors for {n_atoms} atoms",
                self.per_atom.len()
            )));
        }
        for (index, tensor) in self.per_atom.iter().enumerate() {
            if tensor.iter().flatten().any(|value| !value.is_finite()) {
                return Err(Error::InvalidInput(format!(
                    "Born tensor {index} has a non-finite component"
                )));
            }
        }
        Ok(())
    }

    /// Removes the mean tensor component (the Born-charge acoustic sum rule).
    ///
    /// A uniform translation then produces no dipole derivative.  The shift is
    /// not applied automatically because it changes a non-neutral cell; call
    /// it explicitly when the tensors do not already sum to zero.
    #[must_use]
    pub fn neutralised(&self) -> Self {
        let count = self.per_atom.len();
        if count == 0 {
            return self.clone();
        }
        let mut mean = [[0.0; 3]; 3];
        for tensor in &self.per_atom {
            for (row, values) in tensor.iter().enumerate() {
                for (column, value) in values.iter().enumerate() {
                    if let Some(slot) = mean.get_mut(row).and_then(|row| row.get_mut(column)) {
                        *slot += value / count as f64;
                    }
                }
            }
        }
        let per_atom = self
            .per_atom
            .iter()
            .map(|tensor| {
                let mut shifted = [[0.0; 3]; 3];
                for (row, values) in tensor.iter().enumerate() {
                    for (column, value) in values.iter().enumerate() {
                        let average = mean
                            .get(row)
                            .and_then(|row| row.get(column))
                            .copied()
                            .unwrap_or(0.0);
                        if let Some(slot) = shifted.get_mut(row).and_then(|row| row.get_mut(column))
                        {
                            *slot = value - average;
                        }
                    }
                }
                shifted
            })
            .collect();
        Self { per_atom }
    }
}

/// Builds an IR spectrum from Gamma-point modes and Born effective charges.
///
/// This is the engine-grade counterpart of [`ir_spectrum`](super::ir_spectrum):
/// the dipole derivative of mode `m` is
/// `p_a = sum_i sum_b Z*_{i,a b} u_{i,b}`, and the stick intensity is
///
/// ```text
/// I_m = C |p|^2 / sum_i m_i |u_i|^2
/// ```
///
/// in km/mol, exactly the normalisation of the fixed-charge path.  The
/// structure and the modes must describe the same P1 atoms.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, a non-P1 structure, a
/// mode/atom mismatch, a tensor count mismatch, a non-finite tensor, or a
/// non-positive mass.
pub fn ir_spectrum_with_born(
    structure: &Structure,
    modes: &VibrationalModes,
    born: &BornCharges,
    options: &IrOptions,
) -> Result<IrSpectrum, Error> {
    options.validate()?;
    if !structure.is_p1() {
        return Err(Error::InvalidInput(
            "IR intensities need a P1 structure; expand the asymmetric unit first".to_owned(),
        ));
    }
    let count = structure.sites().len();
    if modes.n_atoms() != count {
        return Err(Error::InvalidInput(format!(
            "modes describe {} atoms but the structure has {count}",
            modes.n_atoms()
        )));
    }
    born.validate(count)?;

    let mut masses = Vec::with_capacity(count);
    for site in structure.sites() {
        let element = element_info(site.element())?;
        if element.mass <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "element {} has no positive mass",
                site.element()
            )));
        }
        masses.push(element.mass);
    }

    let mut peaks = Vec::new();
    for (mode, frequency) in modes.real_modes() {
        let Some(displacements) = modes.modes.get(mode) else {
            continue;
        };
        if displacements.len() != count {
            return Err(Error::InvalidInput(format!(
                "mode {mode} has {} atoms, expected {count}",
                displacements.len()
            )));
        }
        let mut dipole = [0.0; 3];
        let mut mass_norm = 0.0;
        for (index, (point, mass)) in displacements.iter().zip(&masses).enumerate() {
            let length_squared = point[0] * point[0] + point[1] * point[1] + point[2] * point[2];
            mass_norm += mass * length_squared;
            let Some(tensor) = born.per_atom.get(index) else {
                continue;
            };
            for (axis, row) in tensor.iter().enumerate() {
                let mut derivative = 0.0;
                for (beta, value) in point.iter().enumerate() {
                    derivative += row.get(beta).copied().unwrap_or(0.0) * value;
                }
                if let Some(slot) = dipole.get_mut(axis) {
                    *slot += derivative;
                }
            }
        }
        if mass_norm <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "mode {mode} has zero mass-weighted norm"
            )));
        }
        let intensity = E2_PER_AMU_TO_KM_PER_MOL
            * (dipole[0] * dipole[0] + dipole[1] * dipole[1] + dipole[2] * dipole[2])
            / mass_norm;
        peaks.push(IrPeak {
            mode,
            frequency,
            intensity,
        });
    }

    peaks.sort_by(|first, second| {
        first
            .frequency
            .partial_cmp(&second.frequency)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(first.mode.cmp(&second.mode))
    });
    let sticks = peaks
        .iter()
        .map(|peak| (peak.frequency, peak.intensity))
        .collect::<Vec<_>>();
    let points = broaden(
        &sticks,
        options.fwhm,
        options.step,
        options.broadening,
        MAX_IR_POINTS,
        "IR",
    )?;
    Ok(IrSpectrum { points, peaks })
}
