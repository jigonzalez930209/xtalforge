//! Symmetry mapping and projectors used by [`super::dof::Dof`].
//!
//! These helpers turn the stored operations of a structure into the data the
//! optimizer needs: which asymmetric site generated each P1 atom, the
//! Cartesian rotation that carries a force back, the multiplicity of every
//! site, the projector onto the displacement a site symmetry allows, and the
//! point group used to average the strain.

use xtalforge_core::{Lattice, P1Cell, Structure, SymOp};

use crate::Error;
use crate::math::{identity3, mat3_mul, norm, sub, vec_mat};

use super::RelaxMode;
use super::dof::SPECIAL_POSITION_TOL;

/// Parent site, Cartesian rotation and multiplicity of each P1 atom.
pub(crate) type Mapping = (Vec<usize>, Vec<[[f64; 3]; 3]>, Vec<usize>);

/// Builds the P1-atom-to-site mapping for the two modes.
pub(crate) fn mapping(
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
pub(crate) fn operations_of(structure: &Structure) -> Vec<SymOp> {
    if structure.operations().is_empty() {
        vec![SymOp::identity()]
    } else {
        structure.operations().to_vec()
    }
}

/// One (3,3) projector per site: what its site symmetry leaves free.
pub(crate) fn projectors(
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
pub(crate) fn point_group(
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
