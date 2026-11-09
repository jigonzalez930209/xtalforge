//! Process-level benchmark for the linked-cell periodic neighbor search.
//!
//! Run with `cargo bench -p xtalforge-core --bench neighbors`. The output is
//! the median of five runs for 10^3, 10^4 and 10^5 atoms at a fixed density.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::time::{Duration, Instant};

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder, neighbors};

fn main() {
    for atoms in [1_000_usize, 10_000, 100_000] {
        let structure = random_structure(atoms);
        let median = median_time(|| {
            let found = neighbors(&structure, 3.0).unwrap();
            std::hint::black_box(found);
        });
        println!(
            "{atoms:>7} atoms  cutoff 3.0 A  {:.1} ms  ({} neighbors)",
            median.as_secs_f64() * 1e3,
            neighbors(&structure, 3.0).unwrap().len()
        );
    }
}

fn random_structure(atoms: usize) -> Structure {
    // About 0.05 atoms per cubic Angstrom, a typical framework density.
    let length = (atoms as f64 * 20.0).cbrt();
    let mut builder = StructureBuilder::new("bench", Lattice::cubic(length).unwrap());
    let mut seed = 0x1234_5678_9abc_def0_u64;
    for index in 0..atoms {
        let frac = [random(&mut seed), random(&mut seed), random(&mut seed)];
        builder = builder.site(Site::new(format!("C{index}"), "C", frac).unwrap());
    }
    builder.build().unwrap()
}

fn random(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    ((*seed >> 11) as f64) / ((1_u64 << 53) as f64)
}

fn median_time(mut run: impl FnMut()) -> Duration {
    let mut samples = Vec::with_capacity(5);
    for _ in 0..5 {
        let start = Instant::now();
        run();
        samples.push(start.elapsed());
    }
    samples.sort_unstable();
    samples[2]
}
