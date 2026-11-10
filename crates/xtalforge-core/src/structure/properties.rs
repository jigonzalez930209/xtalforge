//! Properties of the whole cell: contents, formula, mass, density, Z and charge.

use super::{Error, Site, Structure, element_by_symbol};

impl Structure {
    /// Atoms of the cell as `(element, occupancy, charge)`: the P1 expansion for
    /// a structure with symmetry operations, the sites themselves for a P1
    /// structure or when the expansion fails.
    fn cell_atoms(&self) -> Vec<(String, f64, f64)> {
        if !self.is_p1()
            && let Ok(p1) = self.p1()
        {
            return p1
                .atoms
                .iter()
                .map(|atom| {
                    let charge = self
                        .sites
                        .get(atom.site_idx)
                        .and_then(Site::charge)
                        .unwrap_or(0.0);
                    (atom.element.clone(), atom.occupancy, charge)
                })
                .collect();
        }
        self.sites
            .iter()
            .map(|site| {
                (
                    site.element().to_owned(),
                    site.occupancy(),
                    site.charge().unwrap_or(0.0),
                )
            })
            .collect()
    }

    /// Cell contents as `(element, occupancy sum)` in first-appearance order,
    /// counted on the whole cell (the P1 expansion when the structure has
    /// symmetry operations).
    #[must_use]
    pub fn cell_contents(&self) -> Vec<(String, f64)> {
        let mut contents: Vec<(String, f64)> = Vec::new();
        for (element, occupancy, _) in self.cell_atoms() {
            if let Some(entry) = contents.iter_mut().find(|(name, _)| *name == element) {
                entry.1 += occupancy;
            } else {
                contents.push((element, occupancy));
            }
        }
        contents
    }

    /// Formula string, for example `Fe2 O3` or `C0.5`.
    #[must_use]
    pub fn formula(&self) -> String {
        self.cell_contents()
            .iter()
            .map(|(element, count)| {
                if (*count - count.round()).abs() < 1e-9 {
                    let count = count.round() as u64;
                    if count == 1 {
                        element.clone()
                    } else {
                        format!("{element}{count}")
                    }
                } else {
                    format!("{element}{count:.3}")
                }
            })
            .collect::<Vec<_>>()
            .join(" ")
    }

    /// Mass of the whole cell contents in unified atomic mass units.
    ///
    /// # Errors
    ///
    /// Returns an error if a site element is missing from the element table.
    pub fn cell_mass(&self) -> Result<f64, Error> {
        let mut mass = 0.0;
        for (element, occupancy, _) in self.cell_atoms() {
            mass += occupancy * element_by_symbol(&element)?.mass;
        }
        Ok(mass)
    }

    /// Molar mass of the cell contents in g/mol.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`Structure::cell_mass`].
    pub fn molar_mass(&self) -> Result<f64, Error> {
        self.cell_mass()
    }

    /// Density in g/cm3 of the whole cell (the P1 expansion when the structure
    /// has symmetry operations).
    ///
    /// # Errors
    ///
    /// Returns the errors of [`Structure::cell_mass`].
    pub fn density(&self) -> Result<f64, Error> {
        Ok(self.cell_mass()? * 1.660_539_066_60 / self.lattice.volume())
    }

    /// Number of formula units in the cell: the greatest common divisor of the
    /// element counts of the whole cell, or 1 when a count is not an integer.
    #[must_use]
    pub fn z(&self) -> u32 {
        let mut divisor: u64 = 0;
        for (_, count) in self.cell_contents() {
            if count < 0.5 || (count - count.round()).abs() > 1e-6 {
                return 1;
            }
            let mut a = divisor;
            let mut b = count.round() as u64;
            while b != 0 {
                (a, b) = (b, a % b);
            }
            divisor = a;
        }
        u32::try_from(divisor.max(1)).unwrap_or(1)
    }

    /// Sum of `occupancy * charge` over the atoms of the whole cell.
    #[must_use]
    pub fn total_charge(&self) -> f64 {
        self.cell_atoms()
            .iter()
            .map(|(_, occupancy, charge)| occupancy * charge)
            .sum()
    }

    /// Charge balance, the same as [`Structure::total_charge`].
    #[must_use]
    pub fn charge_balance(&self) -> f64 {
        self.total_charge()
    }
}
