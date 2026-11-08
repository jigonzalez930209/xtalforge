//! The optimizer's variables: asymmetric-unit sites and a symmetric strain.
//!
//! Ported from the core idea of Crystal-Builder `xtal/ff/optimize.py`.  The
//! calculator works on the P1 cell, because that is what has an energy; an
//! optimizer that moved those atoms independently would break the space group
//! on its first step and leave no structure to write back.  The variables here
//! are therefore the sites of the asymmetric unit, the P1 cell is regenerated
//! from them at every evaluation, and the gradient is carried back through the
//! operations that generated it.
//!
//! **Documented subset.**  `RelaxMode::Asymmetric` keeps the asymmetric unit
//! and projects each site's displacement onto what its site symmetry leaves
//! invariant.  The strain is averaged over the point group, so a cubic cell
//! keeps one free strain and a triclinic one keeps six; held cell quantities
//! and external pressure are Task 4.4.2 work not yet in this subset.
//! `RelaxMode::P1` requires a P1 structure and gives every site three free
//! coordinates.

use xtalforge_core::{Lattice, Structure, SymOp};

use super::tensor::{symmetrize, tensor_of, w_of};
use crate::Error;
use crate::math::{add, identity3, mat3_mul, mat3_transpose, norm, scale, sub, vec_mat};
use crate::optimize::RelaxMode;
use crate::structure_util::with_frac;

use xtalforge_core::P1Cell;

/// Tolerance, in Angstrom, for deciding that an operation maps a site onto
/// itself.
pub(crate) const SPECIAL_POSITION_TOL: f64 = 0.05;

/// Degrees of freedom prepared for one structure.
#[derive(Debug, Clone)]
pub(crate) struct Dof {
    base: Structure,
    mode: RelaxMode,
    cell: bool,
    base_lattice: Lattice,
    base_inverse: [[f64; 3]; 3],
    n_sites: usize,
    cell_scale: f64,
    parents: Vec<usize>,
    rotations: Vec<[[f64; 3]; 3]>,
    multiplicities: Vec<usize>,
    projectors: Vec<[[f64; 3]; 3]>,
    fixed: Vec<bool>,
    point_group: Vec<[[f64; 3]; 3]>,
    start_cart: Vec<[f64; 3]>,
}

impl Dof {
    /// Prepares variables for one structure and set of options.
    pub(crate) fn new(
        structure: &Structure,
        mode: RelaxMode,
        cell: bool,
        fixed_atoms: &[usize],
    ) -> Result<Self, Error> {
        if mode == RelaxMode::P1 && !structure.is_p1() {
            return Err(Error::InvalidInput(
                "RelaxMode::P1 needs a P1 structure; use RelaxMode::Asymmetric to keep the space group"
                    .to_owned(),
            ));
        }
        let n_sites = structure.sites().len();
        if n_sites == 0 {
            return Err(Error::InvalidInput(
                "cannot optimize a structure with no sites".to_owned(),
            ));
        }
        let mut fixed = vec![false; n_sites];
        for index in fixed_atoms {
            let Some(slot) = fixed.get_mut(*index) else {
                return Err(Error::InvalidInput(format!(
                    "fixed atom index {index} is outside 0..{n_sites}"
                )));
            };
            *slot = true;
        }
        if fixed.iter().all(|value| *value) && !cell {
            return Err(Error::InvalidInput(
                "every site is fixed, so there is nothing to relax".to_owned(),
            ));
        }
        let base_lattice = *structure.lattice();
        let volume = base_lattice.volume().abs();
        if volume <= 0.0 {
            return Err(Error::InvalidInput(
                "optimization needs a cell with positive volume".to_owned(),
            ));
        }
        let base_inverse = base_lattice.inverse();
        let cell_scale = volume.cbrt();
        let start_cart: Vec<[f64; 3]> = structure
            .sites()
            .iter()
            .map(|site| base_lattice.to_cart(site.frac()))
            .collect();

        let (parents, rotations, multiplicities) = mapping(structure, mode, &base_inverse)?;
        let projectors = projectors(structure, mode, &base_inverse, &base_lattice);
        let point_group = point_group(structure, &base_inverse, &base_lattice);

        Ok(Self {
            base: structure.clone(),
            mode,
            cell,
            base_lattice,
            base_inverse,
            n_sites,
            cell_scale,
            parents,
            rotations,
            multiplicities,
            projectors,
            fixed,
            point_group,
            start_cart,
        })
    }

    /// Whether the cell is a variable.
    pub(crate) fn cell(&self) -> bool {
        self.cell
    }

    /// Number of variables.
    pub(crate) fn n_variables(&self) -> usize {
        self.n_sites * 3 + if self.cell { 6 } else { 0 }
    }

    /// Initial variable vector.
    pub(crate) fn start(&self) -> Vec<f64> {
        let mut values = Vec::with_capacity(self.n_variables());
        for position in &self.start_cart {
            values.extend_from_slice(position);
        }
        values.resize(self.n_variables(), 0.0);
        values
    }

    /// Cartesian site positions at these variables.
    pub(crate) fn sites_of(&self, variables: &[f64]) -> Result<Vec<[f64; 3]>, Error> {
        if variables.len() < self.n_sites * 3 {
            return Err(Error::InvalidInput(
                "optimizer variable vector is shorter than the site block".to_owned(),
            ));
        }
        let mut sites = Vec::with_capacity(self.n_sites);
        for index in 0..self.n_sites {
            let offset = index * 3;
            let x = *variables.get(offset).unwrap_or(&0.0);
            let y = *variables.get(offset + 1).unwrap_or(&0.0);
            let z = *variables.get(offset + 2).unwrap_or(&0.0);
            sites.push([x, y, z]);
        }
        Ok(sites)
    }

    /// Six strain variables at these variables.
    fn strain_variables(&self, variables: &[f64]) -> [f64; 6] {
        if !self.cell {
            return [0.0; 6];
        }
        let base = self.n_sites * 3;
        let mut values = [0.0; 6];
        for (index, slot) in values.iter_mut().enumerate() {
            *slot = variables.get(base + index).copied().unwrap_or(0.0) / self.cell_scale;
        }
        values
    }

    /// Symmetric strain tensor, projected onto the point group when the
    /// asymmetric unit is the variable.
    pub(crate) fn strain_of(&self, variables: &[f64]) -> [[f64; 3]; 3] {
        if !self.cell {
            return [[0.0; 3]; 3];
        }
        let tensor = tensor_of(self.strain_variables(variables));
        if self.mode == RelaxMode::Asymmetric {
            self.project_strain(tensor)
        } else {
            tensor
        }
    }

    /// Lattice at these variables.
    pub(crate) fn lattice_of(&self, variables: &[f64]) -> Result<Lattice, Error> {
        if !self.cell {
            return Ok(self.base_lattice);
        }
        self.base_lattice
            .strained(self.strain_of(variables))
            .map_err(Error::from)
    }

    /// The structure the calculator sees at these variables.
    pub(crate) fn structure_of(&self, variables: &[f64]) -> Result<Structure, Error> {
        let sites = self.sites_of(variables)?;
        let mut frac = Vec::with_capacity(self.n_sites);
        for site in sites {
            frac.push(vec_mat(site, self.base_inverse));
        }
        with_frac(&self.base, self.lattice_of(variables)?, &frac)
    }

    /// Carries P1 forces back onto the site variables.
    pub(crate) fn site_gradient(&self, forces: &[[f64; 3]]) -> Result<Vec<[f64; 3]>, Error> {
        if forces.len() != self.parents.len() {
            return Err(Error::InvalidInput(format!(
                "expected {} P1 forces, got {}",
                self.parents.len(),
                forces.len()
            )));
        }
        let mut gradient = vec![[0.0; 3]; self.n_sites];
        for (index, force) in forces.iter().enumerate() {
            let parent = self.parents.get(index).copied().unwrap_or(0);
            let rotation = self.rotations.get(index).copied().unwrap_or(identity3());
            let contribution = vec_mat(scale(*force, -1.0), mat3_transpose(rotation));
            if let Some(slot) = gradient.get_mut(parent) {
                *slot = add(*slot, contribution);
            }
        }
        for (index, slot) in gradient.iter_mut().enumerate() {
            let projector = self.projectors.get(index).copied().unwrap_or(identity3());
            *slot = vec_mat(*slot, projector);
            if self.fixed.get(index).copied().unwrap_or(false) {
                *slot = [0.0; 3];
            }
        }
        Ok(gradient)
    }

    /// Strain gradient from the stress, in eV per variable.
    ///
    /// `stress` is in eV/Angstrom^3, the same convention as
    /// [`crate::Calculator::stress`] after converting its GPa.  The derivative
    /// of the energy with respect to a strain is `V * stress`, and the
    /// projected strain means the gradient is projected the same way.
    pub(crate) fn strain_gradient(&self, variables: &[f64], stress: [[f64; 3]; 3]) -> [f64; 6] {
        if !self.cell {
            return [0.0; 6];
        }
        let volume = self
            .lattice_of(variables)
            .map(|lattice| lattice.volume().abs())
            .unwrap_or(0.0);
        let tensor = crate::math::mat3_scale(stress, volume);
        let tensor = if self.mode == RelaxMode::Asymmetric {
            self.project_strain(tensor)
        } else {
            tensor
        };
        let mut values = w_of(tensor);
        for value in &mut values {
            *value /= self.cell_scale;
        }
        values
    }

    /// Flattens site gradients and, when present, a strain gradient.
    pub(crate) fn flatten(&self, site: &[[f64; 3]], strain: Option<[f64; 6]>) -> Vec<f64> {
        let mut values = Vec::with_capacity(self.n_variables());
        for position in site {
            values.extend_from_slice(position);
        }
        if self.cell {
            values.extend_from_slice(&strain.unwrap_or([0.0; 6]));
        }
        values
    }

    /// Projects a step direction onto the allowed subspace.
    pub(crate) fn project_direction(&self, direction: &mut [f64]) {
        for index in 0..self.n_sites {
            let offset = index * 3;
            let vector = [
                direction.get(offset).copied().unwrap_or(0.0),
                direction.get(offset + 1).copied().unwrap_or(0.0),
                direction.get(offset + 2).copied().unwrap_or(0.0),
            ];
            let projector = self.projectors.get(index).copied().unwrap_or(identity3());
            let projected = vec_mat(vector, projector);
            for (component, value) in projected.iter().enumerate() {
                if let Some(slot) = direction.get_mut(offset + component) {
                    *slot = *value;
                }
            }
            if self.fixed.get(index).copied().unwrap_or(false) {
                for component in 0..3 {
                    if let Some(slot) = direction.get_mut(offset + component) {
                        *slot = 0.0;
                    }
                }
            }
        }
        if self.cell {
            let base = self.n_sites * 3;
            let mut strain = [0.0; 6];
            for (index, slot) in strain.iter_mut().enumerate() {
                *slot = direction.get(base + index).copied().unwrap_or(0.0);
            }
            if self.mode == RelaxMode::Asymmetric {
                strain = w_of(self.project_strain(tensor_of(strain)));
            }
            for (index, value) in strain.iter().enumerate() {
                if let Some(slot) = direction.get_mut(base + index) {
                    *slot = *value;
                }
            }
        }
    }

    /// Largest force on any free site, in eV/Angstrom.
    ///
    /// A site of multiplicity `m` collects `m` atoms' worth of gradient, so
    /// the number is divided by `m` to mean the same thing on every Wyckoff
    /// position.  Cell variables are not part of this number; a cell is
    /// converged on the energy change instead.
    pub(crate) fn max_force(&self, gradient: &[[f64; 3]]) -> f64 {
        let mut maximum: f64 = 0.0;
        for (index, value) in gradient.iter().enumerate() {
            if self.fixed.get(index).copied().unwrap_or(false) {
                continue;
            }
            let multiplicity = self.multiplicities.get(index).copied().unwrap_or(1).max(1);
            let force = norm(*value) / multiplicity as f64;
            maximum = maximum.max(force);
        }
        maximum
    }

    /// Projects a tensor onto the point group's invariant strains.
    fn project_strain(&self, tensor: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
        if self.point_group.is_empty() {
            return symmetrize(tensor);
        }
        let mut total = [[0.0; 3]; 3];
        for rotation in &self.point_group {
            let transpose = mat3_transpose(*rotation);
            let value = mat3_mul(mat3_mul(transpose, tensor), *rotation);
            total = crate::math::mat3_add(total, value);
        }
        let count = self.point_group.len() as f64;
        symmetrize(crate::math::mat3_scale(total, 1.0 / count))
    }
}

/// Parent site, Cartesian rotation and multiplicity of each P1 atom.
type Mapping = (Vec<usize>, Vec<[[f64; 3]; 3]>, Vec<usize>);

/// Builds the P1-atom-to-site mapping for the two modes.
fn mapping(
    structure: &Structure,
    mode: RelaxMode,
    inverse: &[[f64; 3]; 3],
) -> Result<Mapping, Error> {
    if mode == RelaxMode::P1 {
        return Ok((
            (0..structure.sites().len()).collect(),
            vec![identity3(); structure.sites().len()],
            vec![1_usize; structure.sites().len()],
        ));
    }
    let cell = P1Cell::expand(structure, SPECIAL_POSITION_TOL)?;
    let operations = operations_of(structure);
    let matrix = structure.lattice().matrix();
    let mut parents = Vec::with_capacity(cell.n_atoms());
    let mut rotations = Vec::with_capacity(cell.n_atoms());
    let mut multiplicities = vec![0_usize; structure.sites().len()];
    for atom in &cell.atoms {
        let operation = operations
            .get(atom.op_idx)
            .copied()
            .unwrap_or_else(SymOp::identity);
        parents.push(atom.site_idx);
        rotations.push(cartesian_rotation(operation, inverse, &matrix));
        if let Some(slot) = multiplicities.get_mut(atom.site_idx) {
            *slot += 1;
        }
    }
    for value in &mut multiplicities {
        *value = (*value).max(1);
    }
    Ok((parents, rotations, multiplicities))
}

/// Symmetry operations of a structure, with identity when none are stored.
fn operations_of(structure: &Structure) -> Vec<SymOp> {
    if structure.operations().is_empty() {
        vec![SymOp::identity()]
    } else {
        structure.operations().to_vec()
    }
}

/// One (3,3) projector per site: what its site symmetry leaves free.
fn projectors(
    structure: &Structure,
    mode: RelaxMode,
    inverse: &[[f64; 3]; 3],
    lattice: &Lattice,
) -> Vec<[[f64; 3]; 3]> {
    let count = structure.sites().len();
    if mode == RelaxMode::P1 {
        return vec![identity3(); count];
    }
    let operations = operations_of(structure);
    let matrix = lattice.matrix();
    let mut result = Vec::with_capacity(count);
    for site in structure.sites() {
        let mut total = [[0.0; 3]; 3];
        let mut found = 0_usize;
        for operation in &operations {
            let image = operation.apply(site.frac());
            let shift = sub(image, site.frac());
            let wrapped = [
                shift[0] - shift[0].round(),
                shift[1] - shift[1].round(),
                shift[2] - shift[2].round(),
            ];
            if norm(vec_mat(wrapped, matrix)) < SPECIAL_POSITION_TOL {
                let rotation = cartesian_rotation(*operation, inverse, &matrix);
                total = crate::math::mat3_add(total, rotation);
                found += 1;
            }
        }
        result.push(if found == 0 {
            identity3()
        } else {
            crate::math::mat3_scale(total, 1.0 / found as f64)
        });
    }
    result
}

/// Cartesian rotations of every operation, used for strain projection.
fn point_group(
    structure: &Structure,
    inverse: &[[f64; 3]; 3],
    lattice: &Lattice,
) -> Vec<[[f64; 3]; 3]> {
    let operations = operations_of(structure);
    let matrix = lattice.matrix();
    operations
        .iter()
        .map(|operation| cartesian_rotation(*operation, inverse, &matrix))
        .collect()
}

/// `M^-1 R^T M`: the Cartesian rotation an operation applies to a row vector.
fn cartesian_rotation(
    operation: SymOp,
    inverse: &[[f64; 3]; 3],
    matrix: &[[f64; 3]; 3],
) -> [[f64; 3]; 3] {
    let rotation = operation.rotation();
    let cartesian = [
        [
            f64::from(rotation[0][0]),
            f64::from(rotation[1][0]),
            f64::from(rotation[2][0]),
        ],
        [
            f64::from(rotation[0][1]),
            f64::from(rotation[1][1]),
            f64::from(rotation[2][1]),
        ],
        [
            f64::from(rotation[0][2]),
            f64::from(rotation[1][2]),
            f64::from(rotation[2][2]),
        ],
    ];
    mat3_mul(mat3_mul(*inverse, cartesian), *matrix)
}
