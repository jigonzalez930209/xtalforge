//! Schwarzite level-set surfaces: gyroid and Schwarz P (Task 6.8).
//!
//! Both fields are the periodic trigonometric approximations of triply
//! periodic minimal surfaces. [`ztc_surface_mesh`] samples the gyroid, the
//! surface used for zeolite-templated carbon models, and
//! [`schwarzite_surface_mesh`] samples Schwarz P. The samples are handed to
//! [`xtalforge_analysis::isosurface`], so the mesh is periodic and
//! deterministic for a given grid.

use std::f64::consts::TAU;

use xtalforge_analysis::{TriangleMesh, isosurface};
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};
use xtalforge_io::grid::Grid3;

/// Default number of samples along each cube axis.
pub const DEFAULT_SAMPLES: usize = 16;
/// Default cubic cell side in Angstrom.
pub const DEFAULT_CELL: f64 = 12.0;
/// Minimum spacing between carbon sites when a mesh becomes a structure.
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
    let mesh = ztc_surface_mesh(options)?;
    structure_from_mesh("ztc", options, &mesh)
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
    let mesh = schwarzite_surface_mesh(options)?;
    structure_from_mesh("schwarzite", options, &mesh)
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
    let total = options
        .samples
        .checked_mul(options.samples)
        .and_then(|value| value.checked_mul(options.samples))
        .ok_or_else(|| Error::InvalidInput("the schwarzite grid dimensions overflow".to_owned()))?;
    let samples = options.samples;
    let scale = TAU / samples as f64;
    let step = options.cell / samples as f64;
    let mut data = Vec::with_capacity(total);
    for z in 0..samples {
        for y in 0..samples {
            for x in 0..samples {
                data.push(field(x as f64 * scale, y as f64 * scale, z as f64 * scale));
            }
        }
    }
    Grid3::new(
        [0.0; 3],
        [[step, 0.0, 0.0], [0.0, step, 0.0], [0.0, 0.0, step]],
        [samples; 3],
        data,
        true,
    )
}

fn structure_from_mesh(
    name: &str,
    options: &ZtcOptions,
    mesh: &TriangleMesh,
) -> Result<Structure, Error> {
    let lattice = Lattice::cubic(options.cell)?;
    let mut builder = StructureBuilder::new(name, lattice);
    let mut kept: Vec<[f64; 3]> = Vec::new();
    for vertex in &mesh.vertices {
        let frac = [
            vertex[0] / options.cell,
            vertex[1] / options.cell,
            vertex[2] / options.cell,
        ];
        let close = kept
            .iter()
            .any(|other| lattice.distance(frac, *other) < MIN_CARBON_DISTANCE);
        if close {
            continue;
        }
        kept.push(frac);
        builder = builder.site(Site::new(format!("C{}", kept.len()), "C", frac)?);
    }
    if kept.is_empty() {
        return Err(Error::InvalidInput(
            "the schwarzite mesh has no vertex to place".to_owned(),
        ));
    }
    builder.build()
}

fn gyroid(x: f64, y: f64, z: f64) -> f64 {
    x.sin() * y.cos() + y.sin() * z.cos() + z.sin() * x.cos()
}

fn schwarz_p(x: f64, y: f64, z: f64) -> f64 {
    x.cos() + y.cos() + z.cos()
}
