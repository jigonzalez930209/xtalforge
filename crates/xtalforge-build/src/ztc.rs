//! Schwarzite level-set surfaces: gyroid and Schwarz P (Task 6.8).
//!
//! Both fields are the periodic trigonometric approximations of triply
//! periodic minimal surfaces. [`ztc_surface_mesh`] samples the gyroid, the
//! surface used for zeolite-templated carbon models, and
//! [`schwarzite_surface_mesh`] samples Schwarz P. The samples are handed to
//! the `xtalforge_analysis::isosurface` extractor, so the mesh is periodic and
//! deterministic for a given grid.

use std::f64::consts::TAU;

use xtalforge_analysis::{TriangleMesh, isosurface};
use xtalforge_core::{Error, ExplicitBond, Lattice, Site, Structure, StructureBuilder};

use crate::carbon_recipe::editor::remesh;
use crate::carbon_recipe::field::LevelSet;
use crate::carbon_recipe::graph::dualise;
use crate::carbon_recipe::mesh::weld;
use crate::carbon_recipe::settle::settle;
use xtalforge_io::grid::Grid3;

/// Default number of samples along each cube axis.
pub const DEFAULT_SAMPLES: usize = 16;
/// Default cubic cell side in Angstrom.
pub const DEFAULT_CELL: f64 = 12.0;
/// Minimum spacing between carbon sites of the former vertex-thinning model.
///
/// Kept for API stability; the builders now place one carbon per triangle of
/// the remeshed surface and relax the bonds to 1.42 Angstrom.
pub const MIN_CARBON_DISTANCE: f64 = 1.0;

/// Grid and level used to mesh a schwarzite surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ZtcOptions {
    /// Number of samples along each cube axis; at least two.
    pub samples: usize,
    /// Cubic cell side in Angstrom; positive and finite.
    pub cell: f64,
    /// Field value that defines the surface.
    pub isovalue: f64,
}

impl Default for ZtcOptions {
    fn default() -> Self {
        Self {
            samples: DEFAULT_SAMPLES,
            cell: DEFAULT_CELL,
            isovalue: 0.0,
        }
    }
}

/// Meshes the gyroid level set of [`ZtcOptions`].
///
/// The gyroid field is `sin x cos y + sin y cos z + sin z cos x`, sampled on a
/// periodic cubic grid that spans one period.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for fewer than two samples, a non-positive
/// or non-finite cell, a non-finite isovalue, or a grid that overflows.
pub fn ztc_surface_mesh(options: &ZtcOptions) -> Result<TriangleMesh, Error> {
    mesh(options, gyroid)
}

/// Meshes the Schwarz P level set of [`ZtcOptions`].
///
/// The Schwarz P field is `cos x + cos y + cos z`, sampled on the same
/// periodic cubic grid as [`ztc_surface_mesh`].
///
/// # Errors
///
/// Returns the errors of [`ztc_surface_mesh`].
pub fn schwarzite_surface_mesh(options: &ZtcOptions) -> Result<TriangleMesh, Error> {
    mesh(options, schwarz_p)
}

/// Builds a carbon structure from the gyroid mesh of [`ZtcOptions`].
///
/// Mesh vertices are retained as carbon sites when they are at least
/// [`MIN_CARBON_DISTANCE`] Angstrom from every site already kept, which thins
/// the densely duplicated marching-tetrahedra vertices into a usable starting
/// model.
///
/// # Errors
///
/// Returns the errors of [`ztc_surface_mesh`] or [`Error::InvalidInput`] when
/// the mesh has no vertex.
pub fn build_ztc(options: &ZtcOptions) -> Result<Structure, Error> {
    minimum_cell(options, MIN_GYROID_CELL, "gyroid")?;
    let options = &carbon_grid(options)?;
    let mesh = ztc_surface_mesh(options)?;
    structure_from_mesh("ztc", options, &mesh, Surface::Gyroid)
}

/// Builds a carbon structure from the Schwarz P mesh of [`ZtcOptions`].
///
/// The thinning rule is the one of [`build_ztc`].
///
/// # Errors
///
/// Returns the errors of [`schwarzite_surface_mesh`] or
/// [`Error::InvalidInput`] when the mesh has no vertex.
pub fn build_schwarzite(options: &ZtcOptions) -> Result<Structure, Error> {
    minimum_cell(options, MIN_SCHWARZ_P_CELL, "Schwarz P")?;
    let options = &carbon_grid(options)?;
    let mesh = schwarzite_surface_mesh(options)?;
    structure_from_mesh("schwarzite", options, &mesh, Surface::SchwarzP)
}

fn mesh(options: &ZtcOptions, field: fn(f64, f64, f64) -> f64) -> Result<TriangleMesh, Error> {
    let grid = sample(options, field)?;
    isosurface(&grid, options.isovalue)
}

fn sample(options: &ZtcOptions, field: fn(f64, f64, f64) -> f64) -> Result<Grid3, Error> {
    if options.samples < 2 {
        return Err(Error::InvalidInput(
            "a schwarzite grid needs at least two samples per axis".to_owned(),
        ));
    }
    if !options.cell.is_finite() || options.cell <= 0.0 {
        return Err(Error::InvalidInput(
            "the schwarzite cell side must be positive and finite".to_owned(),
        ));
    }
    if !options.isovalue.is_finite() {
        return Err(Error::InvalidInput(
            "the schwarzite isovalue must be finite".to_owned(),
        ));
    }
    // One extra sample closes the period: the marching tetrahedra stop one
    // sample before the end of each axis, so without the plane at `cell` the
    // surface would have an open seam on the cell faces.
    let points = options
        .samples
        .checked_add(1)
        .ok_or_else(|| Error::InvalidInput("the schwarzite grid dimensions overflow".to_owned()))?;
    let total = points
        .checked_mul(points)
        .and_then(|value| value.checked_mul(points))
        .ok_or_else(|| Error::InvalidInput("the schwarzite grid dimensions overflow".to_owned()))?;
    let samples = options.samples;
    let scale = TAU / samples as f64;
    let step = options.cell / samples as f64;
    let mut data = Vec::with_capacity(total);
    for z in 0..=samples {
        for y in 0..=samples {
            for x in 0..=samples {
                // A sample exactly on the isovalue (the gyroid vanishes on
                // many grid points of an even grid) produces no triangle and
                // leaves a hole; nudge it to one side.
                let value = field(x as f64 * scale, y as f64 * scale, z as f64 * scale);
                let value = if (value - options.isovalue).abs() < 1e-9 {
                    options.isovalue + 1e-9
                } else {
                    value
                };
                data.push(value);
            }
        }
    }
    Grid3::new(
        [0.0; 3],
        [[step, 0.0, 0.0], [0.0, step, 0.0], [0.0, 0.0, step]],
        [samples + 1; 3],
        data,
        true,
    )
}

fn structure_from_mesh(
    name: &str,
    options: &ZtcOptions,
    mesh: &TriangleMesh,
    surface: Surface,
) -> Result<Structure, Error> {
    let lattice = Lattice::cubic(options.cell)?;
    let welded = weld(
        &mesh.vertices,
        &mesh.triangles,
        &lattice,
        [options.samples + 1; 3],
    )?;
    let field = TrigField {
        surface,
        cell: options.cell,
    };
    let even = remesh(&welded, &field, options.isovalue, 12, DUAL_EDGE)?;
    let graph = settle(&dualise(&even)?);
    let mut builder = StructureBuilder::new(name, lattice);
    for (index, center) in graph.centers.iter().enumerate() {
        builder = builder.site(Site::new(format!("C{}", index + 1), "C", *center)?);
    }
    for bond in &graph.bonds {
        builder = builder.explicit_bond(ExplicitBond {
            i: bond.a,
            j: bond.b,
            image: bond.image,
            order: 1.0,
            remove: false,
        });
    }
    builder.build()
}

/// Grid spacing, in Angstrom, from which the carbon remesh converges to an
/// even honeycomb.
const CARBON_GRID_SPACING: f64 = 1.2;

/// The options with the sampling the carbon builders start from.
///
/// The caller's `samples` is still validated, so the builders reject the same
/// options as the mesh functions.
fn carbon_grid(options: &ZtcOptions) -> Result<ZtcOptions, Error> {
    if options.samples < 2 {
        return Err(Error::InvalidInput(
            "a schwarzite grid needs at least two samples per axis".to_owned(),
        ));
    }
    let samples = if options.cell.is_finite() && options.cell > 0.0 {
        ((options.cell / CARBON_GRID_SPACING).round() as usize).clamp(4, 200)
    } else {
        options.samples
    };
    Ok(ZtcOptions {
        samples,
        ..*options
    })
}

/// Smallest gyroid cell, in Angstrom, whose channels hold an sp2 sheet: at
/// 12 Angstrom 98 % of the carbons are three-coordinate, at 10 only 36 %.
pub const MIN_GYROID_CELL: f64 = 12.0;
/// Smallest Schwarz P cell, in Angstrom (96 % three-coordinate at 8).
pub const MIN_SCHWARZ_P_CELL: f64 = 8.0;

fn minimum_cell(options: &ZtcOptions, minimum: f64, surface: &str) -> Result<(), Error> {
    if options.cell.is_finite() && options.cell < minimum {
        return Err(Error::InvalidInput(format!(
            "a {surface} carbon needs a cell of at least {minimum} Angstrom; at {} Angstrom the channels are too narrow for three-coordinate carbon",
            options.cell
        )));
    }
    Ok(())
}

/// Mesh edge whose dual is a 1.42 Angstrom honeycomb (`1.42 sqrt(3)`).
const DUAL_EDGE: f64 = 2.459;

/// Which triply periodic surface a carbon is built on.
#[derive(Debug, Clone, Copy)]
enum Surface {
    Gyroid,
    SchwarzP,
}

/// The trigonometric field of [`Surface`] over a cubic cell.
struct TrigField {
    surface: Surface,
    cell: f64,
}

impl LevelSet for TrigField {
    fn value_gradient(&self, frac: [f64; 3]) -> (f64, [f64; 3]) {
        let [x, y, z] = frac.map(|value| value * TAU);
        let (value, gradient) = match self.surface {
            Surface::Gyroid => (
                gyroid(x, y, z),
                [
                    x.cos() * y.cos() - z.sin() * x.sin(),
                    y.cos() * z.cos() - x.sin() * y.sin(),
                    z.cos() * x.cos() - y.sin() * z.sin(),
                ],
            ),
            Surface::SchwarzP => (schwarz_p(x, y, z), [-x.sin(), -y.sin(), -z.sin()]),
        };
        (value, gradient.map(|component| component * TAU / self.cell))
    }
}

fn gyroid(x: f64, y: f64, z: f64) -> f64 {
    x.sin() * y.cos() + y.sin() * z.cos() + z.sin() * x.cos()
}

fn schwarz_p(x: f64, y: f64, z: f64) -> f64 {
    x.cos() + y.cos() + z.cos()
}
