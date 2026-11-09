//! Supercells and origin shifts (Task 1.7).

use crate::limits::{Limits, check_supercell};
use crate::spatial::{FractionalGrid, perpendicular_widths};
use crate::{Error, Lattice, Site, Structure, StructureBuilder};

impl Structure {
    /// Builds an `na x nb x nc` supercell as a P1 structure.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for zero repeats, or
    /// [`Error::LimitExceeded`] when the estimate exceeds the hard limit.
    pub fn supercell(&self, na: u32, nb: u32, nc: u32) -> Result<Self, Error> {
        self.supercell_with_limits(na, nb, nc, &Limits::default())
    }

    /// Builds an `na x nb x nc` supercell with an explicit limits profile.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for zero repeats, or
    /// [`Error::LimitExceeded`] when the estimate exceeds the profile's hard
    /// limit.
    pub fn supercell_with_limits(
        &self,
        na: u32,
        nb: u32,
        nc: u32,
        limits: &Limits,
    ) -> Result<Self, Error> {
        if na == 0 || nb == 0 || nc == 0 {
            return Err(Error::InvalidInput(
                "supercell repeats must be positive".to_owned(),
            ));
        }
        self.supercell_matrix_with_limits(
            [
                [f64::from(na), 0.0, 0.0],
                [0.0, f64::from(nb), 0.0],
                [0.0, 0.0, f64::from(nc)],
            ],
            limits,
        )
    }

    /// Builds a supercell from a positive integer transformation matrix.
    ///
    /// Rows of `p` are the new vectors in units of the old ones.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for non-integer, singular or negative
    /// determinant matrices, or [`Error::LimitExceeded`] for a refused size.
    pub fn supercell_matrix(&self, p: [[f64; 3]; 3]) -> Result<Self, Error> {
        self.supercell_matrix_with_limits(p, &Limits::default())
    }

    /// Builds a supercell with an explicit limits profile.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for non-integer, singular or negative
    /// determinant matrices, or [`Error::LimitExceeded`] for a refused size.
    pub fn supercell_matrix_with_limits(
        &self,
        p: [[f64; 3]; 3],
        limits: &Limits,
    ) -> Result<Self, Error> {
        for value in p.iter().flatten() {
            if !value.is_finite() || (value - value.round()).abs() > 1e-9 {
                return Err(Error::InvalidInput(
                    "supercell matrix must contain integers".to_owned(),
                ));
            }
        }
        let verdict = check_supercell(limits, self.sites().len() as u64, p)?;
        if verdict.refused() {
            return Err(Error::LimitExceeded(verdict));
        }
        let lattice = self.lattice().transform(p)?;
        let inverse = Lattice::new(p)?.inverse();
        let max_inverse = inverse
            .iter()
            .flatten()
            .map(|value| value.abs())
            .fold(0.0, f64::max);
        let range = (max_inverse * 4.0).ceil() as i64 + 2;
        let widths = perpendicular_widths(&lattice)?;
        let mut grid = FractionalGrid::new(widths, 1e-6);
        let mut sites: Vec<Site> = Vec::new();
        for site in self.sites() {
            for i in -range..=range {
                for j in -range..=range {
                    for k in -range..=range {
                        let candidate = [
                            site.frac()[0] + i as f64,
                            site.frac()[1] + j as f64,
                            site.frac()[2] + k as f64,
                        ];
                        let frac = multiply_row(candidate, inverse);
                        if !frac.iter().all(|value| (-1e-9..1.0).contains(value)) {
                            continue;
                        }
                        let mut duplicate = false;
                        grid.for_each_candidate(site.element(), frac, |existing_index| {
                            let Some(existing) = sites.get(existing_index) else {
                                return true;
                            };
                            if lattice.distance(existing.frac(), frac) <= 1e-6 {
                                duplicate = true;
                                return false;
                            }
                            true
                        });
                        if duplicate {
                            continue;
                        }
                        let index = sites.len();
                        sites.push(site.clone().with_frac(frac)?);
                        grid.insert(site.element(), frac, index);
                    }
                }
            }
        }
        let mut builder = StructureBuilder::new(self.name(), lattice);
        for site in sites {
            builder = builder.site(site);
        }
        for (key, value) in self.metadata() {
            builder = builder.metadata(key.clone(), value.clone());
        }
        builder.build()
    }

    /// Shifts the origin by a fractional vector and rewraps the sites.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for non-finite shifts.
    pub fn origin_shift(&self, delta: [f64; 3]) -> Result<Self, Error> {
        if !delta.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "origin shift must be finite".to_owned(),
            ));
        }
        let mut builder = StructureBuilder::new(self.name(), *self.lattice());
        for site in self.sites() {
            let frac = [
                site.frac()[0] + delta[0],
                site.frac()[1] + delta[1],
                site.frac()[2] + delta[2],
            ];
            builder = builder.site(site.clone().with_frac(Lattice::wrap(frac))?);
        }
        for (key, value) in self.metadata() {
            builder = builder.metadata(key.clone(), value.clone());
        }
        builder.build()
    }
}

fn multiply_row(vector: [f64; 3], matrix: [[f64; 3]; 3]) -> [f64; 3] {
    let [x, y, z] = vector;
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    [
        x * a + y * d + z * g,
        x * b + y * e + z * h,
        x * c + y * f + z * i,
    ]
}
