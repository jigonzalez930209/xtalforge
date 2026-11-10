//! Honeycomb sheets and nanotubes (Tasks 6.5 and 6.6).

const SQRT_3: f64 = 1.732_050_807_568_877_2;

use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

/// Nearest-neighbour C-C distance in graphene, Angstrom.
pub const C_C: f64 = 1.42;
/// Graphene lattice constant, Angstrom.
pub const A_GRAPHENE: f64 = C_C * SQRT_3;

/// Builds an `nx x ny` graphene sheet with `layers` AA-stacked layers.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for zero repeats or layers.
pub fn graphene(
    nx: u32,
    ny: u32,
    layers: u32,
    interlayer: f64,
    vacuum: f64,
) -> Result<Structure, Error> {
    if nx == 0 || ny == 0 || layers == 0 {
        return Err(Error::InvalidInput(
            "graphene repeats and layers must be positive".to_owned(),
        ));
    }
    if !interlayer.is_finite() || interlayer <= 0.0 || !vacuum.is_finite() || vacuum < 0.0 {
        return Err(Error::InvalidInput(
            "graphene spacing and vacuum must be positive".to_owned(),
        ));
    }
    let a = A_GRAPHENE;
    let height = f64::from(layers) * interlayer + vacuum;
    let lattice = Lattice::new([
        [a * f64::from(nx), 0.0, 0.0],
        [
            a * f64::from(ny) / 2.0,
            a * f64::from(ny) * SQRT_3 / 2.0,
            0.0,
        ],
        [0.0, 0.0, height],
    ])?;
    let mut builder = StructureBuilder::new("graphene", lattice);
    for layer in 0..layers {
        let z = (f64::from(layer) * interlayer + 1.0) / height;
        for i in 0..nx {
            for j in 0..ny {
                // Primitive cell (i, j) of the nx x ny supercell: the two atoms
                // sit at (0, 0) and (1/3, 1/3) of the primitive vectors
                // a1 = (a, 0) and a2 = (a/2, a sqrt(3)/2), so the supercell
                // fractional coordinates are (i + u) / nx and (j + v) / ny.
                let positions = [
                    [f64::from(i) / f64::from(nx), f64::from(j) / f64::from(ny)],
                    [
                        (f64::from(i) + 1.0 / 3.0) / f64::from(nx),
                        (f64::from(j) + 1.0 / 3.0) / f64::from(ny),
                    ],
                ];
                for (index, [x, y]) in positions.iter().enumerate() {
                    builder = builder.site(Site::new(
                        format!("C{}_{}_{}", layer, i * ny + j, index),
                        "C",
                        [*x, *y, z],
                    )?);
                }
            }
        }
    }
    builder.build()
}

/// Rolls a graphene sheet into a `(n, m)` nanotube of `length` cells.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for `(0, 0)`, zero length or bad vacuum.
pub fn nanotube(n: i32, m: i32, length: u32, vacuum: f64) -> Result<Structure, Error> {
    if n == 0 && m == 0 {
        return Err(Error::InvalidInput(
            "a nanotube needs a nonzero chiral vector".to_owned(),
        ));
    }
    if length == 0 {
        return Err(Error::InvalidInput(
            "nanotube length must be positive".to_owned(),
        ));
    }
    if !vacuum.is_finite() || vacuum < 0.0 {
        return Err(Error::InvalidInput(
            "nanotube vacuum must be non-negative".to_owned(),
        ));
    }
    let [ni, mi] = [f64::from(n), f64::from(m)];
    let a1 = [A_GRAPHENE, 0.0];
    let a2 = [A_GRAPHENE / 2.0, A_GRAPHENE * SQRT_3 / 2.0];
    let chiral = [ni * a1[0] + mi * a2[0], ni * a1[1] + mi * a2[1]];
    let circumference = (chiral[0] * chiral[0] + chiral[1] * chiral[1]).sqrt();
    let radius = circumference / (2.0 * std::f64::consts::PI);
    let translation = find_translation(n, m, a1, a2)?;
    let length_vector = [
        f64::from(length) * translation[0],
        f64::from(length) * translation[1],
    ];
    let t_one = (translation[0] * translation[0] + translation[1] * translation[1]).sqrt();
    let t_norm = (length_vector[0] * length_vector[0] + length_vector[1] * length_vector[1]).sqrt();
    if t_norm < 1e-9 {
        return Err(Error::InvalidInput(
            "nanotube translation vector is degenerate".to_owned(),
        ));
    }
    // The tube axis is y; the cross-section box is 2R + vacuum on x and z.
    let width = 2.0 * radius + vacuum;
    let lattice = Lattice::new([[width, 0.0, 0.0], [0.0, t_norm, 0.0], [0.0, 0.0, width]])?;
    let mut builder = StructureBuilder::new("nanotube", lattice);
    let count = (ni.abs() + mi.abs()).ceil() as i32 * 4 + 16;
    let eps = 1e-6;
    let mut index = 0;
    for u in -count..=count {
        for v in -count..=count {
            for basis in [[0.0, 0.0], [A_GRAPHENE / 2.0, A_GRAPHENE / (2.0 * SQRT_3)]] {
                let x = f64::from(u) * a1[0] + f64::from(v) * a2[0] + basis[0];
                let y = f64::from(u) * a1[1] + f64::from(v) * a2[1] + basis[1];
                // Coordinates of the sheet atom along the chiral vector and
                // along the tube translation, as fractions of the tube cell.
                let around = (x * chiral[0] + y * chiral[1]) / (circumference * circumference);
                let along =
                    (x * translation[0] + y * translation[1]) / (t_one * t_one * f64::from(length));
                let inside = |value: f64| value > -eps && value < 1.0 - eps;
                if !inside(around) || !inside(along) {
                    continue;
                }
                let angle = 2.0 * std::f64::consts::PI * around;
                let frac_x = (radius * angle.sin() + width / 2.0) / width;
                let frac_y = along.max(0.0);
                let frac_z = (radius * angle.cos() + width / 2.0) / width;
                builder = builder.site(Site::new(
                    format!("C{}", index + 1),
                    "C",
                    [frac_x, frac_y, frac_z],
                )?);
                index += 1;
            }
        }
    }
    if index == 0 {
        return Err(Error::InvalidInput(
            "the nanotube contains no atoms".to_owned(),
        ));
    }
    builder.build()
}

fn find_translation(n: i32, m: i32, a1: [f64; 2], a2: [f64; 2]) -> Result<[f64; 2], Error> {
    // d_R = gcd(2n + m, 2m + n); T = ((2m + n) a1 - (2n + m) a2) / d_R.
    let d = gcd(2 * n + m, 2 * m + n).max(1);
    let t1 = (2 * m + n) / d;
    let t2 = -(2 * n + m) / d;
    let vector = [
        f64::from(t1) * a1[0] + f64::from(t2) * a2[0],
        f64::from(t1) * a1[1] + f64::from(t2) * a2[1],
    ];
    if vector[0].is_finite() && vector[1].is_finite() {
        Ok(vector)
    } else {
        Err(Error::InvalidInput(
            "could not build the nanotube translation".to_owned(),
        ))
    }
}

fn gcd(mut first: i32, mut second: i32) -> i32 {
    while second != 0 {
        let remainder = first % second;
        first = second;
        second = remainder;
    }
    first.abs()
}
