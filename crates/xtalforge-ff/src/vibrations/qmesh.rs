//! Supercell force constants, Monkhorst-Pack q-mesh and phonon DOS.

use std::collections::BTreeMap;
use std::f64::consts::TAU;

use serde::{Deserialize, Serialize};
use xtalforge_core::{P1Cell, Structure, element_info};

use crate::Error;
use crate::calculator::Calculator;

use super::constants::CM_PER_SQRT_EV_PER_A2_AMU;
use super::hessian::{HessianOptions, hessian};

type Block = [[f64; 3]; 3];
type BlockMap = BTreeMap<[i32; 3], Block>;
type BlockGrid = Vec<Vec<BlockMap>>;

/// Force constants folded onto the primitive cell.
#[derive(Debug, Clone)]
pub(crate) struct ForceConstants {
    pub(crate) n_prim: usize,
    pub(crate) blocks: BlockGrid,
    masses: Vec<f64>,
}

/// Phonon density of states on a linear frequency grid.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PhononDos {
    /// Bin centres in cm^-1, ascending.
    pub frequencies: Vec<f64>,
    /// Normalised density, `sum(density * spacing) = 1`.
    pub density: Vec<f64>,
    /// Gaussian smearing in cm^-1.
    pub sigma: f64,
}

/// Options for the phonon DOS histogram.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DosOptions {
    /// Number of bins.
    pub bins: usize,
    /// Gaussian smearing in cm^-1.
    pub sigma: f64,
    /// Lower edge; the data minimum when absent.
    pub min: Option<f64>,
    /// Upper edge; the data maximum when absent.
    pub max: Option<f64>,
}

impl Default for DosOptions {
    fn default() -> Self {
        Self {
            bins: 200,
            sigma: 5.0,
            min: None,
            max: None,
        }
    }
}

impl DosOptions {
    /// Checks the histogram options.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for zero bins, a non-positive sigma or
    /// a reversed range.
    pub fn validate(&self) -> Result<(), Error> {
        if self.bins == 0 {
            return Err(Error::InvalidInput(
                "phonon DOS needs at least one bin".to_owned(),
            ));
        }
        if !self.sigma.is_finite() || self.sigma <= 0.0 {
            return Err(Error::InvalidInput(
                "phonon DOS sigma must be positive and finite".to_owned(),
            ));
        }
        if let (Some(min), Some(max)) = (self.min, self.max)
            && min >= max
        {
            return Err(Error::InvalidInput(
                "phonon DOS min must be below max".to_owned(),
            ));
        }
        Ok(())
    }
}

pub(crate) fn force_constants(
    calculator: &dyn Calculator,
    structure: &Structure,
    repeats: [u32; 3],
    options: &HessianOptions,
) -> Result<ForceConstants, Error> {
    let [na, nb, nc] = repeats;
    if na == 0 || nb == 0 || nc == 0 {
        return Err(Error::InvalidInput(
            "phonon supercell repeats must be positive".to_owned(),
        ));
    }
    let primitive = structure.p1()?;
    let supercell = structure.supercell(na, nb, nc)?;
    let hessian = hessian(calculator, &supercell, options)?;
    let mapping = map_supercell(&primitive, &supercell, repeats)?;
    let cells = na as usize * nb as usize * nc as usize;
    let n = primitive.atoms.len();
    let mut blocks: BlockGrid = vec![vec![BTreeMap::new(); n]; n];
    for (i, (p, first)) in mapping.iter().enumerate() {
        for (j, (q, second)) in mapping.iter().enumerate() {
            let shift = wrap_shift(
                [
                    second[0] - first[0],
                    second[1] - first[1],
                    second[2] - first[2],
                ],
                repeats,
            );
            let mut block = [[0.0; 3]; 3];
            for (alpha, row) in block.iter_mut().enumerate() {
                for (beta, cell) in row.iter_mut().enumerate() {
                    *cell = hessian.value(3 * i + alpha, 3 * j + beta).unwrap_or(0.0);
                }
            }
            if let Some(map) = blocks.get_mut(*p).and_then(|row| row.get_mut(*q)) {
                let entry = map.entry(shift).or_insert([[0.0; 3]; 3]);
                for (alpha, row) in entry.iter_mut().enumerate() {
                    for (beta, cell) in row.iter_mut().enumerate() {
                        let value = block
                            .get(alpha)
                            .and_then(|row| row.get(beta))
                            .copied()
                            .unwrap_or(0.0);
                        *cell += value / cells as f64;
                    }
                }
            }
        }
    }
    symmetrize(&mut blocks);
    acoustic_sum_rule(&mut blocks);
    let mut masses = Vec::with_capacity(n);
    for atom in &primitive.atoms {
        masses.push(element_info(&atom.element)?.mass);
    }
    Ok(ForceConstants {
        n_prim: n,
        blocks,
        masses,
    })
}

/// Monkhorst-Pack q-points in fractional reciprocal coordinates.
pub(crate) fn q_points(grid: [u32; 3]) -> Vec<[f64; 3]> {
    let [na, nb, nc] = grid;
    let mut points = Vec::with_capacity((na as usize) * (nb as usize) * (nc as usize));
    for a in 0..na {
        for b in 0..nb {
            for c in 0..nc {
                points.push([monkhorst(a, na), monkhorst(b, nb), monkhorst(c, nc)]);
            }
        }
    }
    points
}

fn monkhorst(index: u32, count: u32) -> f64 {
    if count <= 1 {
        0.0
    } else {
        (f64::from(index) + 0.5) / f64::from(count) - 0.5
    }
}

pub(crate) fn frequencies(constants: &ForceConstants, q: [f64; 3]) -> Result<Vec<f64>, Error> {
    let n = constants.n_prim;
    let size = 3 * n;
    let mut real = vec![vec![0.0; size]; size];
    let mut imaginary = vec![vec![0.0; size]; size];
    for p in 0..n {
        for q_atom in 0..n {
            let Some(map) = constants.blocks.get(p).and_then(|row| row.get(q_atom)) else {
                continue;
            };
            for (shift, block) in map {
                let phase = TAU
                    * (q[0] * f64::from(shift[0])
                        + q[1] * f64::from(shift[1])
                        + q[2] * f64::from(shift[2]));
                let cosine = phase.cos();
                let sine = phase.sin();
                for (alpha, row) in block.iter().enumerate() {
                    for (beta, value) in row.iter().enumerate() {
                        let row_index = 3 * p + alpha;
                        let column = 3 * q_atom + beta;
                        add(&mut real, row_index, column, value * cosine);
                        add(&mut imaginary, row_index, column, value * sine);
                    }
                }
            }
        }
    }
    for (row_index, row) in real.iter_mut().enumerate() {
        for (column, value) in row.iter_mut().enumerate() {
            let first = constants.masses.get(row_index / 3).copied().unwrap_or(1.0);
            let second = constants.masses.get(column / 3).copied().unwrap_or(1.0);
            *value /= (first * second).sqrt();
        }
    }
    for (row_index, row) in imaginary.iter_mut().enumerate() {
        for (column, value) in row.iter_mut().enumerate() {
            let first = constants.masses.get(row_index / 3).copied().unwrap_or(1.0);
            let second = constants.masses.get(column / 3).copied().unwrap_or(1.0);
            *value /= (first * second).sqrt();
        }
    }
    let mut embedded = vec![vec![0.0; 2 * size]; 2 * size];
    for (row, values) in real.iter().enumerate() {
        for (column, value) in values.iter().enumerate() {
            set(&mut embedded, row, column, *value);
            set(
                &mut embedded,
                row,
                size + column,
                -imaginary
                    .get(row)
                    .and_then(|values| values.get(column))
                    .copied()
                    .unwrap_or(0.0),
            );
            set(
                &mut embedded,
                size + row,
                column,
                imaginary
                    .get(row)
                    .and_then(|values| values.get(column))
                    .copied()
                    .unwrap_or(0.0),
            );
            set(&mut embedded, size + row, size + column, *value);
        }
    }
    let (eigenvalues, _) = crate::math::symmetric_eigen(&embedded).ok_or_else(|| {
        Error::InvalidInput("the phonon dynamical matrix could not be diagonalized".to_owned())
    })?;
    // The real embedding duplicates every Hermitian eigenvalue once.
    let unique = eigenvalues
        .into_iter()
        .step_by(2)
        .take(size)
        .collect::<Vec<_>>();
    if unique.len() != size {
        return Err(Error::InvalidInput(
            "the phonon dynamical matrix lost an eigenvalue".to_owned(),
        ));
    }
    Ok(unique
        .into_iter()
        .map(|value| value.signum() * (value.abs().sqrt() * CM_PER_SQRT_EV_PER_A2_AMU))
        .collect())
}

pub(crate) fn dos(frequencies: &[Vec<f64>], options: &DosOptions) -> Result<PhononDos, Error> {
    options.validate()?;
    let mut all = frequencies.iter().flatten().copied().collect::<Vec<_>>();
    if all.is_empty() {
        return Err(Error::InvalidInput(
            "phonon DOS needs at least one frequency".to_owned(),
        ));
    }
    all.sort_by(f64::total_cmp);
    let first = all.first().copied().unwrap_or(0.0);
    let last = all.last().copied().unwrap_or(first + 1.0);
    let min = options.min.unwrap_or(first);
    let max = options.max.unwrap_or(last.max(min + 1.0));
    if min >= max {
        return Err(Error::InvalidInput("phonon DOS range is empty".to_owned()));
    }
    let bins = options.bins;
    let step = (max - min) / bins as f64;
    let mut density = vec![0.0; bins];
    let sigma2 = 2.0 * options.sigma * options.sigma;
    for value in &all {
        for (index, cell) in density.iter_mut().enumerate() {
            let center = min + (index as f64 + 0.5) * step;
            let offset = center - value;
            *cell += (-offset * offset / sigma2).exp();
        }
    }
    let integral = density.iter().sum::<f64>() * step;
    if integral > 0.0 {
        for value in &mut density {
            *value /= integral;
        }
    }
    Ok(PhononDos {
        frequencies: (0..bins)
            .map(|index| min + (index as f64 + 0.5) * step)
            .collect(),
        density,
        sigma: options.sigma,
    })
}

fn map_supercell(
    primitive: &P1Cell,
    supercell: &Structure,
    repeats: [u32; 3],
) -> Result<Vec<(usize, [i32; 3])>, Error> {
    let [na, nb, nc] = repeats;
    let matrix = [
        [f64::from(na), 0.0, 0.0],
        [0.0, f64::from(nb), 0.0],
        [0.0, 0.0, f64::from(nc)],
    ];
    let mut mapping = Vec::with_capacity(supercell.sites().len());
    for site in supercell.sites() {
        let frac = multiply_row(site.frac(), matrix);
        let mut best: Option<(f64, usize)> = None;
        for (index, atom) in primitive.atoms.iter().enumerate() {
            let distance = primitive.lattice.distance(frac, atom.frac);
            if best.is_none_or(|(value, _)| distance < value) {
                best = Some((distance, index));
            }
        }
        let Some((distance, index)) = best else {
            return Err(Error::InvalidInput(
                "phonon supercell has no primitive atom".to_owned(),
            ));
        };
        if distance > 1e-3 {
            return Err(Error::InvalidInput(
                "phonon supercell atom does not map to the primitive cell".to_owned(),
            ));
        }
        let atom = primitive
            .atoms
            .get(index)
            .ok_or_else(|| Error::InvalidInput("primitive index out of range".to_owned()))?;
        let shift = [
            (frac[0] - atom.frac[0]).round() as i32,
            (frac[1] - atom.frac[1]).round() as i32,
            (frac[2] - atom.frac[2]).round() as i32,
        ];
        mapping.push((index, shift));
    }
    Ok(mapping)
}

fn symmetrize(blocks: &mut BlockGrid) {
    let n = blocks.len();
    let mut updates = Vec::new();
    for p in 0..n {
        for q in 0..n {
            let Some(map) = blocks.get(p).and_then(|row| row.get(q)) else {
                continue;
            };
            for (shift, block) in map {
                let opposite = [-shift[0], -shift[1], -shift[2]];
                let other = blocks
                    .get(q)
                    .and_then(|row| row.get(p))
                    .and_then(|map| map.get(&opposite))
                    .copied()
                    .unwrap_or([[0.0; 3]; 3]);
                let other = transpose(other);
                let mut average = [[0.0; 3]; 3];
                for (alpha, row) in average.iter_mut().enumerate() {
                    for (beta, cell) in row.iter_mut().enumerate() {
                        let first = block
                            .get(alpha)
                            .and_then(|row| row.get(beta))
                            .copied()
                            .unwrap_or(0.0);
                        let second = other
                            .get(alpha)
                            .and_then(|row| row.get(beta))
                            .copied()
                            .unwrap_or(0.0);
                        *cell = 0.5 * (first + second);
                    }
                }
                updates.push((p, q, *shift, opposite, average));
            }
        }
    }
    for (p, q, shift, opposite, average) in updates {
        if let Some(map) = blocks.get_mut(p).and_then(|row| row.get_mut(q)) {
            map.insert(shift, average);
            map.insert(opposite, transpose(average));
        }
    }
}

fn acoustic_sum_rule(blocks: &mut BlockGrid) {
    let n = blocks.len();
    for p in 0..n {
        let mut total = [[0.0; 3]; 3];
        if let Some(row) = blocks.get(p) {
            for map in row {
                for block in map.values() {
                    for (alpha, row) in block.iter().enumerate() {
                        for (beta, value) in row.iter().enumerate() {
                            if let Some(cell) =
                                total.get_mut(alpha).and_then(|row| row.get_mut(beta))
                            {
                                *cell += value;
                            }
                        }
                    }
                }
            }
        }
        if let Some(map) = blocks.get_mut(p).and_then(|row| row.get_mut(p)) {
            let onsite = map.entry([0, 0, 0]).or_insert([[0.0; 3]; 3]);
            for (alpha, row) in onsite.iter_mut().enumerate() {
                for (beta, cell) in row.iter_mut().enumerate() {
                    *cell -= total
                        .get(alpha)
                        .and_then(|row| row.get(beta))
                        .copied()
                        .unwrap_or(0.0);
                }
            }
        }
    }
}

fn transpose(block: Block) -> Block {
    let mut out = [[0.0; 3]; 3];
    for (alpha, row) in out.iter_mut().enumerate() {
        for (beta, cell) in row.iter_mut().enumerate() {
            *cell = block
                .get(beta)
                .and_then(|row| row.get(alpha))
                .copied()
                .unwrap_or(0.0);
        }
    }
    out
}

fn add(matrix: &mut [Vec<f64>], row: usize, column: usize, value: f64) {
    if let Some(cell) = matrix.get_mut(row).and_then(|row| row.get_mut(column)) {
        *cell += value;
    }
}

fn set(matrix: &mut [Vec<f64>], row: usize, column: usize, value: f64) {
    if let Some(cell) = matrix.get_mut(row).and_then(|row| row.get_mut(column)) {
        *cell = value;
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

fn wrap_shift(mut shift: [i32; 3], repeats: [u32; 3]) -> [i32; 3] {
    for (value, count) in shift.iter_mut().zip(repeats) {
        let count = count as i32;
        if count > 1 {
            *value = value.rem_euclid(count);
            if *value > count / 2 {
                *value -= count;
            }
        }
    }
    shift
}
