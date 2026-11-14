//! Small deterministic helpers of the multicomponent SQS search.

use xtalforge_core::{Error, Structure, StructureBuilder};

use super::MulticomponentSqsOptions;

pub(super) fn host_without_position(p1: &Structure, frac: [f64; 3]) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new(p1.name(), *p1.lattice());
    for site in p1.sites() {
        if !same_position(site.frac(), frac) {
            builder = builder.site(site.clone());
        }
    }
    builder.build()
}

fn same_position(first: [f64; 3], second: [f64; 3]) -> bool {
    let delta = [
        first[0] - second[0],
        first[1] - second[1],
        first[2] - second[2],
    ];
    delta.iter().all(|value| {
        let wrapped = value - value.round();
        wrapped.abs() < 1e-5
    })
}

pub(super) fn allocate_counts(
    fractions: &[f64],
    sites: usize,
    options: &MulticomponentSqsOptions,
) -> Option<Vec<usize>> {
    let mut counts = vec![0usize; fractions.len()];
    let mut remainders: Vec<(f64, usize)> = Vec::new();
    let mut assigned = 0usize;
    for (index, fraction) in fractions.iter().enumerate() {
        if *fraction <= options.concentration_tolerance {
            continue;
        }
        let exact = fraction * sites as f64;
        let floor = exact.floor() as usize;
        if let Some(slot) = counts.get_mut(index) {
            *slot = floor;
        }
        assigned += floor;
        remainders.push((exact - floor as f64, index));
    }
    remainders.sort_by(|first, second| {
        second
            .0
            .partial_cmp(&first.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| first.1.cmp(&second.1))
    });
    let mut remaining = sites.saturating_sub(assigned);
    for (_, index) in &remainders {
        if remaining == 0 {
            break;
        }
        if let Some(slot) = counts.get_mut(*index) {
            *slot += 1;
            remaining -= 1;
        }
    }
    for (fraction, count) in fractions.iter().zip(&counts) {
        if *fraction > options.concentration_tolerance && *count == 0 {
            return None;
        }
    }
    Some(counts)
}

pub(super) struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub(super) fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    pub(super) fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    pub(super) fn shuffle<T>(&mut self, values: &mut [T]) {
        for index in (1..values.len()).rev() {
            let other = (self.next_u64() % (index as u64 + 1)) as usize;
            values.swap(index, other);
        }
    }
}
