//! Disordered carbon frameworks from a recipe (Task 6.8).
//!
//! The route is a convolution surface around a periodic net, a periodic
//! triangle mesh of one of its level sets, a cut into a connected ribbon
//! per layer, and a dualised carbon graph with the requested edge
//! termination. The mesh modules are private; the recipe builder is the
//! public entry point.

mod build;
mod clash;
pub(crate) mod editor;
mod editor_ops;
pub(crate) mod field;
pub(crate) mod graph;
pub(crate) mod mesh;
mod net;
mod relax;
pub(crate) mod settle;
mod surface;
mod terminate;

pub use build::{CarbonBuild, CarbonRecipe, build_carbon};

use self::field::GaussianField;
use xtalforge_core::Error;

/// Carbon atoms per cubic Angstrom in one g/cm3.
pub(crate) const ATOMS_PER_GCC: f64 = 0.602_214_08 / 12.011;

/// Summary of the surface stage, used by the recipe builder and tests.
#[derive(Debug, Clone, PartialEq)]
pub struct CarbonSurface {
    /// Vertices after welding and remeshing.
    pub vertices: usize,
    /// Triangles after remeshing.
    pub triangles: usize,
    /// Euler characteristic of the welded mesh.
    pub euler: i64,
    /// Surface area in square Angstrom.
    pub area: f64,
    /// Mean net edge in Angstrom.
    pub edge_length: f64,
    /// Net vertices in the supercell.
    pub net_vertices: usize,
    /// Net edges in the supercell.
    pub net_edges: usize,
    /// Beads in the convolution field.
    pub beads: usize,
    /// Connected triangle pieces.
    pub pieces: usize,
    /// Carbons on the closed dual.
    pub carbons: usize,
    /// Directions in which the closed dual percolates.
    pub periodicity: usize,
}

/// Builds, remeshes and dualises one level set of a net.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unsupported net, a bad scale,
/// a bad repeat or a level that produces no sheet.
pub fn carbon_surface(
    net_name: &str,
    repeat: [usize; 3],
    scale: f64,
    radius_ratio: f64,
    seed: u64,
) -> Result<CarbonSurface, Error> {
    let net = net::net_of(net_name)?;
    let net = net::repeated(&net, repeat)?;
    let mut rng = Rng::new(seed);
    let distortion = net::Distortion {
        vertex: (0..net.frac.len()).map(|_| rng.normal3()).collect(),
        edge: (0..net.edges.len()).map(|_| rng.normal3()).collect(),
    };
    let skeleton = net::skeleton(&net, scale, &distortion, 0.0, 0.0)?;
    let radius = radius_ratio * skeleton.edge_length;
    let width = field::WIDTH_FRACTION * radius;
    let level = field::level_for(radius, width, net::BEAD_SPACING);
    let field = GaussianField::new(&skeleton, width)?;
    let mesh = field.march_natural(level)?;
    let mesh = editor::remesh(&mesh, &field, level, 6, editor::TARGET_EDGE)?;
    let graph = graph::dualise(&mesh)?;
    Ok(CarbonSurface {
        vertices: mesh.frac.len(),
        triangles: mesh.triangles(),
        euler: mesh.euler(),
        area: mesh.area(),
        edge_length: skeleton.edge_length,
        net_vertices: skeleton.frac.len(),
        net_edges: skeleton.edges.len(),
        beads: skeleton.bead_edge.len(),
        pieces: mesh
            .components()
            .iter()
            .collect::<std::collections::BTreeSet<_>>()
            .len(),
        carbons: graph.len(),
        periodicity: graph.periodicity(),
    })
}

/// Small deterministic xoshiro256++ generator for the carbon recipe.
#[derive(Debug, Clone)]
pub(crate) struct Rng {
    state: [u64; 4],
}

impl Rng {
    pub(crate) fn new(seed: u64) -> Self {
        let mut state = [0u64; 4];
        let mut value = seed;
        for slot in &mut state {
            value = value.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = value;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            *slot = z ^ (z >> 31);
        }
        Self { state }
    }

    pub(crate) fn next(&mut self) -> u64 {
        let result = self.state[0]
            .wrapping_add(self.state[3])
            .rotate_left(23)
            .wrapping_add(self.state[0]);
        let t = self.state[1] << 17;
        self.state[2] ^= self.state[0];
        self.state[3] ^= self.state[1];
        self.state[1] ^= self.state[2];
        self.state[0] ^= self.state[3];
        self.state[2] ^= t;
        self.state[3] = self.state[3].rotate_left(45);
        result
    }

    pub(crate) fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }

    pub(crate) fn normal(&mut self) -> f64 {
        let first = self.unit().max(1e-12);
        let second = self.unit();
        (-2.0 * first.ln()).sqrt() * (std::f64::consts::TAU * second).cos()
    }

    pub(crate) fn normal3(&mut self) -> [f64; 3] {
        [self.normal(), self.normal(), self.normal()]
    }
}
