//! Deterministic block evaluation shared by the UFF terms (Task 15.2.3).
//!
//! Every term computes one [`Part`] per term index: an energy, the virial
//! increments it makes, and the gradient scatters it makes, all in exactly the
//! order the sequential loop used to perform them. Blocks of term indices are
//! computed on worker threads when the `parallel` feature is on, then applied
//! to the accumulator **in block order, term by term**. Floating-point
//! addition therefore follows the same association as the sequential path:
//! the energy is added in term order and every gradient entry receives its
//! contributions in term order. Results are bit-identical for every thread
//! count, and the sequential reference is the same code with parallelism
//! disabled.

use crate::Error;
use crate::uff::math::add_virial;

/// Term indices evaluated per block.
pub(crate) const BLOCK: usize = 1_024;
#[cfg(feature = "parallel")]
const MAX_WINDOW: usize = 8;
/// Maximum scatters and virial increments one term makes.
const CAPACITY: usize = 4;

/// One gradient contribution.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Scatter {
    /// Atom index.
    pub atom: usize,
    /// Gradient value.
    pub value: [f64; 3],
}

/// One term's contribution, with the accumulation order preserved.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Part {
    /// Energy of this term in kcal/mol.
    pub energy: f64,
    n_virial: usize,
    virial: [([f64; 3], [f64; 3]); CAPACITY],
    n_scatter: usize,
    scatter: [Scatter; CAPACITY],
}

impl Part {
    /// A term that contributes nothing.
    pub(crate) const fn zero() -> Self {
        Self {
            energy: 0.0,
            n_virial: 0,
            virial: [([0.0; 3], [0.0; 3]); CAPACITY],
            n_scatter: 0,
            scatter: [Scatter {
                atom: 0,
                value: [0.0; 3],
            }; CAPACITY],
        }
    }

    /// Records one `sum d (x) g` virial increment.
    pub(crate) fn push_virial(&mut self, displacement: [f64; 3], gradient: [f64; 3]) {
        if let Some(slot) = self.virial.get_mut(self.n_virial) {
            *slot = (displacement, gradient);
            self.n_virial += 1;
        }
    }

    /// Records one `gradient[atom] += value` scatter.
    pub(crate) fn push_scatter(&mut self, atom: usize, value: [f64; 3]) {
        if let Some(slot) = self.scatter.get_mut(self.n_scatter) {
            *slot = Scatter { atom, value };
            self.n_scatter += 1;
        }
    }
}

/// Evaluates `count` terms and accumulates them in index order.
///
/// `compute` must be pure apart from reading its inputs. `parallel` is `true`
/// for the default path and `false` for the sequential reference.
///
/// # Errors
///
/// Returns the first error from `compute` or a gradient index outside
/// `gradient`.
pub(crate) fn evaluate<F>(
    count: usize,
    parallel: bool,
    compute: F,
    energy: &mut f64,
    gradient: &mut [[f64; 3]],
    virial: &mut [[f64; 3]; 3],
) -> Result<(), Error>
where
    F: Fn(usize) -> Result<Part, Error> + Sync,
{
    let window = window(parallel);
    let mut start = 0;
    while start < count {
        let ranges: Vec<(usize, usize)> = (0..window)
            .filter_map(|block| {
                let block_start = start.saturating_add(block.saturating_mul(BLOCK));
                if block_start >= count {
                    return None;
                }
                let end = block_start.saturating_add(BLOCK).min(count);
                Some((block_start, end))
            })
            .collect();
        let blocks = compute_blocks(&ranges, &compute, parallel)?;
        for block in blocks {
            for part in block {
                apply(&part, energy, gradient, virial)?;
            }
        }
        start = start.saturating_add(window.saturating_mul(BLOCK));
    }
    Ok(())
}

#[cfg(feature = "parallel")]
fn window(parallel: bool) -> usize {
    if parallel {
        rayon::current_num_threads().clamp(2, MAX_WINDOW)
    } else {
        1
    }
}

#[cfg(not(feature = "parallel"))]
fn window(_parallel: bool) -> usize {
    1
}

#[cfg(feature = "parallel")]
fn compute_blocks<F>(
    ranges: &[(usize, usize)],
    compute: &F,
    parallel: bool,
) -> Result<Vec<Vec<Part>>, Error>
where
    F: Fn(usize) -> Result<Part, Error> + Sync,
{
    if parallel && ranges.len() > 1 {
        use rayon::prelude::*;
        ranges
            .par_iter()
            .map(|(start, end)| (*start..*end).map(compute).collect::<Result<Vec<_>, _>>())
            .collect()
    } else {
        ranges
            .iter()
            .map(|(start, end)| (*start..*end).map(compute).collect::<Result<Vec<_>, _>>())
            .collect()
    }
}

#[cfg(not(feature = "parallel"))]
fn compute_blocks<F>(
    ranges: &[(usize, usize)],
    compute: &F,
    _parallel: bool,
) -> Result<Vec<Vec<Part>>, Error>
where
    F: Fn(usize) -> Result<Part, Error> + Sync,
{
    ranges
        .iter()
        .map(|(start, end)| (*start..*end).map(compute).collect::<Result<Vec<_>, _>>())
        .collect()
}

fn apply(
    part: &Part,
    energy: &mut f64,
    gradient: &mut [[f64; 3]],
    virial: &mut [[f64; 3]; 3],
) -> Result<(), Error> {
    *energy += part.energy;
    for (displacement, value) in part.virial.iter().take(part.n_virial) {
        add_virial(virial, *displacement, *value);
    }
    for scatter in part.scatter.iter().take(part.n_scatter) {
        let Some(slot) = gradient.get_mut(scatter.atom) else {
            return Err(Error::Topology(format!(
                "atom index {} is out of range",
                scatter.atom
            )));
        };
        *slot = crate::uff::math::add(*slot, scatter.value);
    }
    Ok(())
}
