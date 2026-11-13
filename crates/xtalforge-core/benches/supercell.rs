//! Process-level benchmark for a 10^6-atom supercell (Task 15.1.1).
//!
//! Run with `cargo bench -p xtalforge-core --bench supercell`. The output is
//! the median of three in-memory builds; the source is a synthetic eight-atom
//! cubic cell, so the measurement does not depend on an input parser or on
//! serialising one million atoms.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::time::Instant;

use xtalforge_core::limits::Limits;
use xtalforge_core::{Lattice, Profile, Site, Structure, StructureBuilder};

const REPEAT: u32 = 50;
const ATOM_COUNT: usize = 1_000_000;

fn main() {
    let structure = silicon_cell();
    let repetitions = std::env::var("XTALFORGE_BENCH_REPS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(5)
        .max(1);
    let mut samples = Vec::with_capacity(repetitions);
    for _ in 0..repetitions {
        let start = Instant::now();
        let built = structure
            .supercell_with_limits(REPEAT, REPEAT, REPEAT, &Limits::new(Profile::WarnOnly))
            .unwrap();
        std::hint::black_box(&built);
        assert_eq!(built.sites().len(), ATOM_COUNT);
        samples.push(start.elapsed());
    }
    samples.sort_unstable();
    println!(
        "supercell atoms={ATOM_COUNT} median_ms={:.3}",
        samples[repetitions / 2].as_secs_f64() * 1e3
    );
    println!(
        "samples_ms={}",
        samples
            .iter()
            .map(|sample| format!("{:.3}", sample.as_secs_f64() * 1e3))
            .collect::<Vec<_>>()
            .join(",")
    );
}

fn silicon_cell() -> Structure {
    let lattice = Lattice::cubic(5.430_941).unwrap();
    let mut builder = StructureBuilder::new("Si", lattice);
    for index in 0_u32..8 {
        let frac = [
            f64::from(index & 1) * 0.5,
            f64::from((index >> 1) & 1) * 0.5,
            f64::from((index >> 2) & 1) * 0.5,
        ];
        builder = builder.site(Site::new(format!("Si{index}"), "Si", frac).unwrap());
    }
    builder.build().unwrap()
}
