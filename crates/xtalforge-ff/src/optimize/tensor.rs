//! Symmetric strain tensors in six metric components.
//!
//! A symmetric 3x3 strain has six independent numbers.  Stored as
//! `[e_xx, e_yy, e_zz, sqrt(2) e_yz, sqrt(2) e_xz, sqrt(2) e_xy]`, their
//! Euclidean norm is the tensor's Frobenius norm, which makes the point-group
//! average an orthogonal projector.  The strain variables of
//! [`super::dof::Dof`] use exactly this basis.

/// Symmetric tensor from six metric components.
pub(crate) fn tensor_of(values: [f64; 6]) -> [[f64; 3]; 3] {
    let root = std::f64::consts::SQRT_2;
    [
        [values[0], values[5] / root, values[4] / root],
        [values[5] / root, values[1], values[3] / root],
        [values[4] / root, values[3] / root, values[2]],
    ]
}

/// Six metric components of a symmetric tensor.
pub(crate) fn w_of(tensor: [[f64; 3]; 3]) -> [f64; 6] {
    let root = std::f64::consts::SQRT_2;
    [
        tensor[0][0],
        tensor[1][1],
        tensor[2][2],
        root * tensor[1][2],
        root * tensor[0][2],
        root * tensor[0][1],
    ]
}

/// Symmetric part of a tensor.
pub(crate) fn symmetrize(tensor: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [
            tensor[0][0],
            0.5 * (tensor[0][1] + tensor[1][0]),
            0.5 * (tensor[0][2] + tensor[2][0]),
        ],
        [
            0.5 * (tensor[1][0] + tensor[0][1]),
            tensor[1][1],
            0.5 * (tensor[1][2] + tensor[2][1]),
        ],
        [
            0.5 * (tensor[2][0] + tensor[0][2]),
            0.5 * (tensor[2][1] + tensor[1][2]),
            tensor[2][2],
        ],
    ]
}
