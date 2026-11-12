//! Open metal sites (Task 10.2.4).

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure, element_info};

/// Options for [`open_metal_sites`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OpenMetalOptions {
    /// Neighbor cutoff in Angstrom.
    pub cutoff: f64,
    /// Minimum Voronoi solid angle for a coordinated atom.
    pub solid_angle: f64,
    /// Probe radius used for the accessibility test.
    pub probe_radius: f64,
}

impl Default for OpenMetalOptions {
    fn default() -> Self {
        Self {
            cutoff: 3.5,
            solid_angle: 0.05,
            probe_radius: 1.2,
        }
    }
}

/// One coordinatively unsaturated metal.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OpenMetalSite {
    /// P1 atom index.
    pub atom: usize,
    /// Metal element.
    pub element: String,
    /// Coordination number found.
    pub coordination: usize,
    /// Typical coordination number for the element family.
    pub expected_coordination: usize,
    /// Missing neighbors with respect to the typical number.
    pub missing: usize,
    /// Free direction in Cartesian coordinates.
    pub direction: [f64; 3],
    /// Distance from the metal to the nearest framework surface along the free
    /// direction, in Angstrom.
    pub clearance: f64,
    /// Whether a probe of `probe_radius` reaches the free direction.
    pub accessible: bool,
}

/// Finds coordinatively unsaturated metals.
///
/// A metal is open when its Voronoi coordination number is below the typical
/// value for its element family and its coordination sphere leaves a free
/// direction. Accessibility is the distance-grid clearance along that
/// direction.
///
/// # Errors
///
/// Returns the errors of the P1 expansion, coordination analysis or distance
/// grid.
pub fn open_metal_sites(
    structure: &Structure,
    options: &OpenMetalOptions,
) -> Result<Vec<OpenMetalSite>, Error> {
    let p1 = structure.p1()?;
    let p1_structure = crate::platform::util::p1_structure(&p1)?;
    let grid = xtalforge_analysis::distance_grid(
        &p1_structure,
        &xtalforge_analysis::default_radii(&p1_structure)?,
        0.5,
    )?;
    let mut results = Vec::new();
    for (index, atom) in p1.atoms.iter().enumerate() {
        let info = match element_info(&atom.element) {
            Ok(info) => info,
            Err(_) => continue,
        };
        if !info.is_metal {
            continue;
        }
        let shell = xtalforge_core::coordination(
            &p1_structure,
            index,
            xtalforge_core::CoordinationStrategy::Voronoi {
                cutoff: options.cutoff,
                solid_angle_min: options.solid_angle,
            },
        )?;
        let expected = expected_coordination(&atom.element);
        if shell.number >= expected {
            continue;
        }
        let centre = p1.lattice.to_cart(atom.frac);
        let mut sum = [0.0; 3];
        for neighbour in &shell.neighbors {
            let Some(other) = p1.atoms.get(neighbour.index) else {
                continue;
            };
            let cart = p1.lattice.to_cart(other.frac);
            let delta = [
                cart[0]
                    + f64::from(neighbour.image[0]) * p1.lattice.matrix()[0][0]
                    + f64::from(neighbour.image[1]) * p1.lattice.matrix()[1][0]
                    + f64::from(neighbour.image[2]) * p1.lattice.matrix()[2][0]
                    - centre[0],
                cart[1]
                    + f64::from(neighbour.image[0]) * p1.lattice.matrix()[0][1]
                    + f64::from(neighbour.image[1]) * p1.lattice.matrix()[1][1]
                    + f64::from(neighbour.image[2]) * p1.lattice.matrix()[2][1]
                    - centre[1],
                cart[2]
                    + f64::from(neighbour.image[0]) * p1.lattice.matrix()[0][2]
                    + f64::from(neighbour.image[1]) * p1.lattice.matrix()[1][2]
                    + f64::from(neighbour.image[2]) * p1.lattice.matrix()[2][2]
                    - centre[2],
            ];
            sum[0] += delta[0];
            sum[1] += delta[1];
            sum[2] += delta[2];
        }
        let norm = (sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]).sqrt();
        if norm < 1e-9 {
            continue;
        }
        let direction = [-sum[0] / norm, -sum[1] / norm, -sum[2] / norm];
        let probe = [
            centre[0] + direction[0] * 2.0,
            centre[1] + direction[1] * 2.0,
            centre[2] + direction[2] * 2.0,
        ];
        let clearance = sample_grid(&grid, &p1_structure, probe);
        results.push(OpenMetalSite {
            atom: index,
            element: atom.element.clone(),
            coordination: shell.number,
            expected_coordination: expected,
            missing: expected.saturating_sub(shell.number),
            direction,
            clearance,
            accessible: clearance >= options.probe_radius,
        });
    }
    Ok(results)
}

fn expected_coordination(element: &str) -> usize {
    match element {
        "Li" | "Na" | "K" | "Rb" | "Cs" => 6,
        "Ca" | "Sr" | "Ba" | "Mg" => 7,
        "La" | "Ce" | "Pr" | "Nd" | "Sm" | "Eu" | "Gd" | "Tb" | "Dy" | "Ho" | "Er" | "Tm"
        | "Yb" | "Lu" | "Y" | "Sc" => 8,
        _ => 6,
    }
}

fn sample_grid(
    grid: &xtalforge_analysis::DistanceGrid,
    structure: &Structure,
    point: [f64; 3],
) -> f64 {
    let frac = structure.lattice().to_frac(point);
    let frac = [
        frac[0].rem_euclid(1.0),
        frac[1].rem_euclid(1.0),
        frac[2].rem_euclid(1.0),
    ];
    let [n_x, n_y, n_z] = grid.grid.shape;
    let x = ((frac[0] * n_x as f64).floor() as usize).min(n_x - 1);
    let y = ((frac[1] * n_y as f64).floor() as usize).min(n_y - 1);
    let z = ((frac[2] * n_z as f64).floor() as usize).min(n_z - 1);
    grid.grid.get(x, y, z).unwrap_or(f64::INFINITY)
}
