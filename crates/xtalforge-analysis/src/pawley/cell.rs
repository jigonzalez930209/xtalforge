//! Metric family and free-parameter packing for the Pawley fit (Task 7.6.2).

/// Metric family of the refined cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Family {
    Cubic,
    Tetragonal,
    Hexagonal,
    Orthorhombic,
    Monoclinic,
    Triclinic,
}

impl Family {
    /// Detects the family from the cell numbers, using a small tolerance.
    pub(crate) fn from_parameters(parameters: [f64; 6]) -> Self {
        let [a, b, c, alpha, beta, gamma] = parameters;
        let close = |first: f64, second: f64| (first - second).abs() <= 0.01;
        let right = |angle: f64| (angle - 90.0).abs() <= 0.01;
        if close(a, b) && close(b, c) && right(alpha) && right(beta) && right(gamma) {
            Self::Cubic
        } else if close(a, b) && right(alpha) && right(beta) && right(gamma) {
            Self::Tetragonal
        } else if close(a, b) && right(alpha) && right(beta) && close(gamma, 120.0) {
            Self::Hexagonal
        } else if right(alpha) && right(beta) && right(gamma) {
            Self::Orthorhombic
        } else if right(alpha) && right(gamma) {
            Self::Monoclinic
        } else {
            Self::Triclinic
        }
    }

    /// Indices of the cell numbers the family leaves free.
    pub(crate) fn free_indices(self) -> &'static [usize] {
        match self {
            Self::Cubic => &[0],
            Self::Tetragonal | Self::Hexagonal => &[0, 2],
            Self::Orthorhombic => &[0, 1, 2],
            Self::Monoclinic => &[0, 1, 2, 4],
            Self::Triclinic => &[0, 1, 2, 3, 4, 5],
        }
    }

    /// Packs the free numbers of a full cell.
    pub(crate) fn pack(self, parameters: [f64; 6]) -> Vec<f64> {
        self.free_indices()
            .iter()
            .filter_map(|index| parameters.get(*index).copied())
            .collect()
    }

    /// Restores a full cell from the free values, enforcing the family ties.
    pub(crate) fn unpack(self, _base: [f64; 6], free: &[f64]) -> [f64; 6] {
        let value = |index: usize| free.get(index).copied().unwrap_or(0.0);
        match self {
            Self::Cubic => {
                let a = value(0);
                [a, a, a, 90.0, 90.0, 90.0]
            }
            Self::Tetragonal => {
                let a = value(0);
                [a, a, value(1), 90.0, 90.0, 90.0]
            }
            Self::Hexagonal => {
                let a = value(0);
                [a, a, value(1), 90.0, 90.0, 120.0]
            }
            Self::Orthorhombic => [value(0), value(1), value(2), 90.0, 90.0, 90.0],
            Self::Monoclinic => [value(0), value(1), value(2), 90.0, value(3), 90.0],
            Self::Triclinic => [value(0), value(1), value(2), value(3), value(4), value(5)],
        }
    }
}
