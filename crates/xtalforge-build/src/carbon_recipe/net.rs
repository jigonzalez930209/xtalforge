//! Periodic nets used by the carbon recipe builder (Task 6.8).
//!
//! Only the nets the recipes need are placed here. `dia` is the diamond
//! net in its primitive rhombohedral cell (two vertices, four edges) and
//! `pcu` is the primitive cubic net with three self-edges. Both are
//! scaled so the mean edge is one Angstrom; [`repeated`] grows the
//! quotient into a supercell and [`skeleton`] places the beads the
//! convolution surface is built from.

#![allow(clippy::indexing_slicing)] // bounded local arrays

use xtalforge_core::{Error, Lattice};

/// One oriented net edge: `i -> j` plus the lattice image of `j`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct NetEdge {
    /// First vertex.
    pub(crate) i: usize,
    /// Second vertex.
    pub(crate) j: usize,
    /// Lattice image of `j`.
    pub(crate) image: [i32; 3],
}

/// A net in a cell whose mean edge is one Angstrom.
#[derive(Debug, Clone)]
pub(crate) struct Net {
    /// Cell.
    pub(crate) lattice: Lattice,
    /// Vertices in fractional coordinates.
    pub(crate) frac: Vec<[f64; 3]>,
    /// Oriented edges.
    pub(crate) edges: Vec<NetEdge>,
}

/// The net with the given name, or a typed error for an unsupported one.
pub(crate) fn net_of(name: &str) -> Result<Net, Error> {
    match name {
        "dia" => {
            // Conventional cubic diamond cell: eight vertices and sixteen
            // nearest-neighbour edges, with the cell scaled so the mean edge
            // is one Angstrom.
            let a = 4.0 / 3.0_f64.sqrt();
            let lattice = Lattice::cubic(a)?;
            let mut frac = Vec::with_capacity(8);
            for vertex in [
                [0.0, 0.0, 0.0],
                [0.0, 0.5, 0.5],
                [0.5, 0.0, 0.5],
                [0.5, 0.5, 0.0],
                [0.25, 0.25, 0.25],
                [0.25, 0.75, 0.75],
                [0.75, 0.25, 0.75],
                [0.75, 0.75, 0.25],
            ] {
                frac.push(vertex);
            }
            let edges = [
                (0, 4, [0, 0, 0]),
                (0, 5, [0, -1, -1]),
                (0, 6, [-1, 0, -1]),
                (0, 7, [-1, -1, 0]),
                (1, 4, [0, 0, 0]),
                (1, 5, [0, 0, 0]),
                (1, 6, [-1, 0, 0]),
                (1, 7, [-1, 0, 0]),
                (2, 4, [0, 0, 0]),
                (2, 5, [0, -1, 0]),
                (2, 6, [0, 0, 0]),
                (2, 7, [0, -1, 0]),
                (3, 4, [0, 0, 0]),
                (3, 5, [0, 0, -1]),
                (3, 6, [0, 0, -1]),
                (3, 7, [0, 0, 0]),
            ]
            .into_iter()
            .map(|(i, j, image)| NetEdge { i, j, image })
            .collect();
            Ok(Net {
                lattice,
                frac,
                edges,
            })
        }
        "pcu" => Ok(Net {
            lattice: Lattice::cubic(1.0)?,
            frac: vec![[0.0, 0.0, 0.0]],
            edges: vec![
                NetEdge {
                    i: 0,
                    j: 0,
                    image: [1, 0, 0],
                },
                NetEdge {
                    i: 0,
                    j: 0,
                    image: [0, 1, 0],
                },
                NetEdge {
                    i: 0,
                    j: 0,
                    image: [0, 0, 1],
                },
            ],
        }),
        other => Err(Error::InvalidInput(format!(
            "carbon recipes support the dia and pcu nets, not {other:?}"
        ))),
    }
}

/// The net repeated over a supercell.
pub(crate) fn repeated(net: &Net, repeat: [usize; 3]) -> Result<Net, Error> {
    if repeat.contains(&0) || repeat.iter().any(|value| *value > 8) {
        return Err(Error::InvalidInput(
            "a carbon repeat must be between 1 and 8 along every axis".to_owned(),
        ));
    }
    let cells = cell_indices(repeat);
    let index = |cell: [usize; 3]| -> usize {
        cell[0] * repeat[1] * repeat[2] + cell[1] * repeat[2] + cell[2]
    };
    let mut frac = Vec::with_capacity(net.frac.len() * cells.len());
    for cell in &cells {
        for vertex in &net.frac {
            frac.push([
                (vertex[0] + cell[0] as f64) / repeat[0] as f64,
                (vertex[1] + cell[1] as f64) / repeat[1] as f64,
                (vertex[2] + cell[2] as f64) / repeat[2] as f64,
            ]);
        }
    }
    let n = net.frac.len();
    let mut edges = Vec::with_capacity(net.edges.len() * cells.len());
    for (cell_index, cell) in cells.iter().enumerate() {
        for edge in &net.edges {
            let far = [
                cell[0] as i32 + edge.image[0],
                cell[1] as i32 + edge.image[1],
                cell[2] as i32 + edge.image[2],
            ];
            let home = [
                far[0].rem_euclid(repeat[0] as i32) as usize,
                far[1].rem_euclid(repeat[1] as i32) as usize,
                far[2].rem_euclid(repeat[2] as i32) as usize,
            ];
            edges.push(NetEdge {
                i: edge.i + n * cell_index,
                j: edge.j + n * index(home),
                image: [
                    (far[0] - home[0] as i32) / repeat[0] as i32,
                    (far[1] - home[1] as i32) / repeat[1] as i32,
                    (far[2] - home[2] as i32) / repeat[2] as i32,
                ],
            });
        }
    }
    let matrix = net.lattice.matrix();
    let super_matrix = [
        [
            matrix[0][0] * repeat[0] as f64,
            matrix[0][1] * repeat[0] as f64,
            matrix[0][2] * repeat[0] as f64,
        ],
        [
            matrix[1][0] * repeat[1] as f64,
            matrix[1][1] * repeat[1] as f64,
            matrix[1][2] * repeat[1] as f64,
        ],
        [
            matrix[2][0] * repeat[2] as f64,
            matrix[2][1] * repeat[2] as f64,
            matrix[2][2] * repeat[2] as f64,
        ],
    ];
    Ok(Net {
        lattice: Lattice::new(super_matrix)?,
        frac,
        edges,
    })
}

fn cell_indices(repeat: [usize; 3]) -> Vec<[usize; 3]> {
    let mut cells = Vec::with_capacity(repeat[0] * repeat[1] * repeat[2]);
    for a in 0..repeat[0] {
        for b in 0..repeat[1] {
            for c in 0..repeat[2] {
                cells.push([a, b, c]);
            }
        }
    }
    cells
}

/// The skeleton of one build: the distorted net in its cell and the
/// beads of the convolution field.
#[derive(Debug, Clone)]
pub(crate) struct Skeleton {
    /// Scaled cell.
    pub(crate) lattice: Lattice,
    /// Distorted vertices in fractional coordinates.
    pub(crate) frac: Vec<[f64; 3]>,
    /// Edges, unchanged by the distortion.
    pub(crate) edges: Vec<NetEdge>,
    /// Beads in Cartesian coordinates, wrapped into the cell.
    pub(crate) beads: Vec<[f64; 3]>,
    /// Mean edge length in Angstrom.
    pub(crate) edge_length: f64,
    /// Which edge each bead belongs to.
    pub(crate) bead_edge: Vec<usize>,
}

/// Bead spacing along an edge.
pub(crate) const BEAD_SPACING: f64 = 0.3;

/// Draws the distortion once so the scale can change without moving it.
#[derive(Debug, Clone)]
pub(crate) struct Distortion {
    /// One standard-normal displacement per vertex.
    pub(crate) vertex: Vec<[f64; 3]>,
    /// One standard-normal bow per edge.
    pub(crate) edge: Vec<[f64; 3]>,
}

/// Places the net at `scale` Angstrom per unit edge.
pub(crate) fn skeleton(
    net: &Net,
    scale: f64,
    distortion: &Distortion,
    sigma_vertex: f64,
    sigma_edge: f64,
) -> Result<Skeleton, Error> {
    if !scale.is_finite() || scale <= 0.0 {
        return Err(Error::InvalidInput(
            "the carbon scale must be positive and finite".to_owned(),
        ));
    }
    let matrix = net.lattice.matrix();
    let scaled = Lattice::new(matrix.map(|row| row.map(|value| value * scale)))?;
    let mut frac = net.frac.clone();
    if sigma_vertex > 0.0 {
        for (index, vertex) in frac.iter_mut().enumerate() {
            let draw = distortion.vertex.get(index).copied().unwrap_or([0.0; 3]);
            let offset = scaled.to_frac(draw.map(|value| value * sigma_vertex));
            for axis in 0..3 {
                vertex[axis] += offset[axis];
            }
        }
    }
    let cart = scaled.to_cart_many(&frac);
    let scaled_matrix = scaled.matrix();
    let mut beads = Vec::new();
    let mut bead_edge = Vec::new();
    let mut lengths = Vec::with_capacity(net.edges.len());
    for (index, edge) in net.edges.iter().enumerate() {
        let start = cart[edge.i];
        let end = add(
            cart[edge.j],
            transform(edge.image.map(f64::from), scaled_matrix),
        );
        let axis = sub(end, start);
        let length = norm(axis);
        lengths.push(length);
        let bow = if sigma_edge > 0.0 && length > 1e-9 {
            let draw = distortion.edge.get(index).copied().unwrap_or([0.0; 3]);
            let unit = axis.map(|value| value / length);
            let along = dot(draw, unit);
            sub(draw, unit.map(|value| value * along)).map(|value| value * sigma_edge)
        } else {
            [0.0; 3]
        };
        let middle = add(
            [
                0.5 * (start[0] + end[0]),
                0.5 * (start[1] + end[1]),
                0.5 * (start[2] + end[2]),
            ],
            bow.map(|value| 2.0 * value),
        );
        let count = (length / BEAD_SPACING).ceil().max(1.0) as usize;
        for step in 0..count {
            let t = step as f64 / count as f64;
            let one = 1.0 - t;
            let point = [
                one * one * start[0] + 2.0 * one * t * middle[0] + t * t * end[0],
                one * one * start[1] + 2.0 * one * t * middle[1] + t * t * end[1],
                one * one * start[2] + 2.0 * one * t * middle[2] + t * t * end[2],
            ];
            let wrapped = scaled.to_frac(point).map(|value| value.rem_euclid(1.0));
            beads.push(scaled.to_cart(wrapped));
            bead_edge.push(index);
        }
    }
    let mean = if lengths.is_empty() {
        0.0
    } else {
        lengths.iter().sum::<f64>() / lengths.len() as f64
    };
    Ok(Skeleton {
        lattice: scaled,
        frac,
        edges: net.edges.clone(),
        beads,
        edge_length: mean,
        bead_edge,
    })
}

pub(crate) fn add(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[0] + second[0],
        first[1] + second[1],
        first[2] + second[2],
    ]
}

pub(crate) fn sub(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[0] - second[0],
        first[1] - second[1],
        first[2] - second[2],
    ]
}

pub(crate) fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    first[0] * second[0] + first[1] * second[1] + first[2] * second[2]
}

pub(crate) fn norm(value: [f64; 3]) -> f64 {
    dot(value, value).sqrt()
}

pub(crate) fn transform(vector: [f64; 3], matrix: [[f64; 3]; 3]) -> [f64; 3] {
    [
        vector[0] * matrix[0][0] + vector[1] * matrix[1][0] + vector[2] * matrix[2][0],
        vector[0] * matrix[0][1] + vector[1] * matrix[1][1] + vector[2] * matrix[2][1],
        vector[0] * matrix[0][2] + vector[1] * matrix[1][2] + vector[2] * matrix[2][2],
    ]
}
