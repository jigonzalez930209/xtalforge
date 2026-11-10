//! Process-level benchmark for PXRD and porosity (Task 8.6.1).
//!
//! Run with `cargo bench -p xtalforge-analysis --bench pxrd`.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::time::{Duration, Instant};

use xtalforge_analysis::{PxrdOptions, calculate_pxrd, default_radii, distance_grid};
use xtalforge_io::read_cif;

fn median(mut samples: Vec<Duration>) -> Duration {
    samples.sort_unstable();
    samples[samples.len() / 2]
}

fn main() {
    let text = include_str!("../../../tests/corpus/cod/MOF-5.cif");
    let structure = read_cif(text).unwrap();

    let mut samples = Vec::new();
    for _ in 0..10 {
        let start = Instant::now();
        let pattern = calculate_pxrd(&structure, &PxrdOptions::default()).unwrap();
        std::hint::black_box(pattern);
        samples.push(start.elapsed());
    }
    println!("MOF-5 PXRD {:.1} ms", median(samples).as_secs_f64() * 1e3);

    let radii = default_radii(&structure).unwrap();
    let mut samples = Vec::new();
    for _ in 0..3 {
        let start = Instant::now();
        let grid = distance_grid(&structure, &radii, 0.8).unwrap();
        std::hint::black_box(grid.accessible_volume(1.2).unwrap());
        samples.push(start.elapsed());
    }
    println!(
        "MOF-5 porosity (0.8 A) {:.1} ms",
        median(samples).as_secs_f64() * 1e3
    );
}
