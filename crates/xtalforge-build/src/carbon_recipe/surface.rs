//! Cell solve for a carbon recipe (Task 6.8.1).
//!
//! A sheet's area grows as the square of the cell and the volume as the
//! cube, so for a fixed net, radius ratio and coverage there is one scale
//! at which the carbon the sheet can hold is the density asked for. The
//! distortions are in Angstrom and do not scale, so the correction is a
//! fixed point that settles in a few rounds.

#![allow(clippy::indexing_slicing, clippy::too_many_arguments)] // bounded solve context

use xtalforge_core::Error;

use super::field::{GaussianField, WIDTH_FRACTION, level_for};
use super::mesh::Mesh;
use super::net::{BEAD_SPACING, Distortion, net_of, repeated, skeleton};
use super::{ATOMS_PER_GCC, Rng};

/// A solved sheet and the carbon count it must hold.
pub(crate) struct SolvedSurface {
    pub(crate) mesh: Mesh,
    pub(crate) edge_length: f64,
    pub(crate) carbons: usize,
    pub(crate) segments: Vec<(usize, usize, [f64; 3], [f64; 3])>,
    pub(crate) vertices: Vec<[f64; 3]>,
}

/// Solves the cell for the requested carbon density.
pub(crate) fn solve_surface(
    recipe_net: &str,
    repeat: [usize; 3],
    density: f64,
    radius_ratio: f64,
    coverage: f64,
    sigma_vertex: f64,
    sigma_edge: f64,
    rng: &mut Rng,
) -> Result<SolvedSurface, Error> {
    if !coverage.is_finite() || !(0.0..=1.0).contains(&coverage) || coverage == 0.0 {
        return Err(Error::InvalidInput(
            "carbon coverage must be above zero and at most one".to_owned(),
        ));
    }
    if !density.is_finite() || density <= 0.0 {
        return Err(Error::InvalidInput(
            "carbon density must be positive and finite".to_owned(),
        ));
    }
    let net = repeated(&net_of(recipe_net)?, repeat)?;
    let distortion = Distortion {
        vertex: (0..net.frac.len()).map(|_| rng.normal3()).collect(),
        edge: (0..net.edges.len()).map(|_| rng.normal3()).collect(),
    };
    let mut scale = 13.0;
    let mut solved = None;
    for _round in 0..8 {
        let frame = skeleton(&net, scale, &distortion, sigma_vertex, sigma_edge)?;
        let radius = radius_ratio * frame.edge_length;
        if !radius.is_finite() || radius <= 0.0 {
            return Err(Error::InvalidInput(
                "carbon radius must be positive and finite".to_owned(),
            ));
        }
        let width = WIDTH_FRACTION * radius;
        let level = level_for(radius, width, BEAD_SPACING);
        let field = GaussianField::new(&frame, width)?;
        let mesh = field.march_natural(level)?;
        let mesh = super::editor::remesh(&mesh, &field, level, 6, super::editor::TARGET_EDGE)?;
        let held = mesh.triangles() as f64 * coverage;
        let wanted = density * ATOMS_PER_GCC * mesh.lattice.volume();
        let carbons = wanted.round().max(1.0) as usize;
        let matrix = frame.lattice.matrix();
        let mut segments = Vec::with_capacity(frame.edges.len());
        for edge in &frame.edges {
            let start = frame.lattice.to_cart(frame.frac[edge.i]);
            let end = super::net::add(
                frame.lattice.to_cart(frame.frac[edge.j]),
                super::net::transform(edge.image.map(f64::from), matrix),
            );
            segments.push((edge.i, edge.j, start, end));
        }
        solved = Some(SolvedSurface {
            mesh,
            edge_length: frame.edge_length,
            carbons,
            segments,
            vertices: frame
                .frac
                .iter()
                .map(|frac| frame.lattice.to_cart(*frac))
                .collect(),
        });
        let ratio = held / wanted;
        if !ratio.is_finite() || ratio <= 0.0 {
            break;
        }
        if (ratio - 1.0).abs() < 0.01 {
            break;
        }
        scale *= ratio;
    }
    solved.ok_or_else(|| Error::InvalidInput("the carbon scale did not solve".to_owned()))
}
