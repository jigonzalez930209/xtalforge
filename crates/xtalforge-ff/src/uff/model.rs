//! The prepared UFF model: typing, topology and pair list built once.
//!
//! Ported from Crystal-Builder `xtal/ff/uff/calculator.py`. Construction is
//! the expensive half -- atom typing, angle and torsion enumeration over the
//! periodic bond graph, the van der Waals pair list -- and then the geometry
//! moves under a fixed topology. [`crate::Calculator`] builds one model per
//! call; optimisation code can hold a model and call [`UffModel::energy_at`]
//! repeatedly.
//!
//! Dummy `X` atoms are absent from the model. [`UffModel::n_atoms`] is the
//! number of real atoms; forces returned by the model have one entry per P1
//! atom including a zero for each dummy.

use xtalforge_core::{Structure, perceive_bonds};

use crate::Error;
use crate::calculator::CalculatorOptions;
use crate::uff::geometry::{AtomView, Geometry, local_bonds, local_rings};
use crate::uff::math::virial_to_stress;
use crate::uff::params::KCAL_PER_MOL_TO_EV;
use crate::uff::topology::Topology;
use crate::uff::typer::{Typing, assign};

/// A UFF calculation prepared for one structure and one set of options.
#[derive(Debug, Clone)]
pub struct UffModel {
    view: AtomView,
    typing: Typing,
    topology: Topology,
    options: CalculatorOptions,
    volume: f64,
}

impl UffModel {
    /// Builds the model: P1 expansion, bond graph, typing and topology.
    ///
    /// # Errors
    ///
    /// Returns a typing error for an element without UFF parameters, a core
    /// error when the P1 cell or bond graph cannot be built, and a topology
    /// error for inconsistent geometry or parameters.
    pub fn build(structure: &Structure, options: CalculatorOptions) -> Result<Self, Error> {
        options.validate()?;
        let p1 = structure.p1()?;
        let view = AtomView::from_p1(structure, &p1)?;
        let graph = perceive_bonds(structure, options.bond_rules)?;
        let bonds = local_bonds(&graph, &view);
        let geometry = Geometry::new(view.n_atoms(), &bonds);
        let rings = local_rings(structure, &view)?;
        let mut typing = assign(&view, &bonds, &geometry, &rings, options.parameter_set)?;
        if bonds.is_empty() && view.n_atoms() > 1 {
            typing.warnings.push(
                "no bonds were perceived, so this is a gas of atoms held together by nothing but van der Waals"
                    .to_owned(),
            );
        }
        let topology = Topology::build(
            structure,
            &view,
            &bonds,
            &geometry,
            &typing,
            options.vdw,
            options.vdw_cutoff,
        )?;
        let volume = p1.lattice.volume();
        Ok(Self {
            view,
            typing,
            topology,
            options,
            volume,
        })
    }

    /// Builds a model with default options.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`UffModel::build`].
    pub fn from_structure(structure: &Structure) -> Result<Self, Error> {
        Self::build(structure, CalculatorOptions::default())
    }

    /// Options this model was built with.
    #[must_use]
    pub fn options(&self) -> &CalculatorOptions {
        &self.options
    }

    /// Number of real (non-dummy) atoms.
    #[must_use]
    pub fn n_atoms(&self) -> usize {
        self.view.n_atoms()
    }

    /// Number of P1 atoms including dummies.
    #[must_use]
    pub fn n_p1_atoms(&self) -> usize {
        self.view.n_full
    }

    /// Indices in the P1 cell of the atoms this model evaluates.
    #[must_use]
    pub fn active_atoms(&self) -> &[usize] {
        &self.view.to_full
    }

    /// Cartesian positions at which the model was built, keyed like
    /// [`UffModel::active_atoms`].
    #[must_use]
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.view.positions
    }

    /// UFF type names in active-atom order.
    #[must_use]
    pub fn types(&self) -> Vec<&str> {
        self.typing.names()
    }

    /// Full typing, with reasons and confidences.
    #[must_use]
    pub fn typing(&self) -> &Typing {
        &self.typing
    }

    /// Warnings worth showing before the energy is trusted.
    #[must_use]
    pub fn warnings(&self) -> &[String] {
        &self.typing.warnings
    }

    /// Number of bond terms.
    #[must_use]
    pub fn n_bonds(&self) -> usize {
        self.topology.bonds.len()
    }

    /// Number of angle terms.
    #[must_use]
    pub fn n_angles(&self) -> usize {
        self.topology.angles.len()
    }

    /// Number of torsion terms.
    #[must_use]
    pub fn n_torsions(&self) -> usize {
        self.topology.torsions.len()
    }

    /// Number of inversion terms.
    #[must_use]
    pub fn n_inversions(&self) -> usize {
        self.topology.inversions.len()
    }

    /// Number of van der Waals pairs.
    #[must_use]
    pub fn n_vdw_pairs(&self) -> usize {
        self.topology.vdw.len()
    }

    /// Energy in eV at the built geometry.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`UffModel::energy_at`].
    pub fn energy(&self) -> Result<f64, Error> {
        self.energy_at(&self.view.positions)
    }

    /// Energy in eV at these active-atom positions.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the position count does not match
    /// the model, or a topology error from a term.
    pub fn energy_at(&self, positions: &[[f64; 3]]) -> Result<f64, Error> {
        let (energy, _) = self.raw(positions, None)?;
        Ok(energy)
    }

    /// Energy in eV and forces in eV/Angstrom at the built geometry.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`UffModel::energy_forces_at`].
    pub fn energy_forces(&self) -> Result<(f64, Vec<[f64; 3]>), Error> {
        self.energy_forces_at(&self.view.positions)
    }

    /// Energy in eV and forces in eV/Angstrom at these positions.
    ///
    /// The force vector has one entry per P1 atom; dummy `X` atoms carry a
    /// zero.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the position count does not match
    /// the model, or a topology error from a term.
    pub fn energy_forces_at(&self, positions: &[[f64; 3]]) -> Result<(f64, Vec<[f64; 3]>), Error> {
        let (energy, gradient) = self.raw(positions, None)?;
        let forces: Vec<[f64; 3]> = gradient
            .iter()
            .map(|value| {
                let [x, y, z] = *value;
                [
                    -x * KCAL_PER_MOL_TO_EV,
                    -y * KCAL_PER_MOL_TO_EV,
                    -z * KCAL_PER_MOL_TO_EV,
                ]
            })
            .collect();
        Ok((energy, self.view.full_forces(&forces)))
    }

    /// Stress tensor in GPa at the built geometry.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`UffModel::stress_at`].
    pub fn stress(&self) -> Result<[[f64; 3]; 3], Error> {
        self.stress_at(&self.view.positions)
    }

    /// Stress tensor in GPa at these positions.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a zero-volume cell or a position
    /// count that does not match the model.
    pub fn stress_at(&self, positions: &[[f64; 3]]) -> Result<[[f64; 3]; 3], Error> {
        if self.volume <= 0.0 {
            return Err(Error::InvalidInput(
                "stress needs a positive cell volume".to_owned(),
            ));
        }
        let mut virial = [[0.0; 3]; 3];
        let _ = self.raw(positions, Some(&mut virial))?;
        Ok(virial_to_stress(virial, self.volume))
    }

    fn raw(
        &self,
        positions: &[[f64; 3]],
        virial: Option<&mut [[f64; 3]; 3]>,
    ) -> Result<(f64, Vec<[f64; 3]>), Error> {
        if positions.len() != self.view.n_atoms() {
            return Err(Error::InvalidInput(format!(
                "expected {} positions, got {}",
                self.view.n_atoms(),
                positions.len()
            )));
        }
        let (energy, gradient) =
            self.topology
                .energy_and_gradient(positions, self.view.matrix, virial)?;
        Ok((energy * KCAL_PER_MOL_TO_EV, gradient))
    }
}
