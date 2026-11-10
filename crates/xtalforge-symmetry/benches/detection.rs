//! Process-level benchmark for symmetry detection (Task 8.6.1).
//!
//! Run with `cargo bench -p xtalforge-symmetry --bench detection`.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::time::{Duration, Instant};

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_symmetry::{SymmetryOptions, detect_symmetry};

fn grid(repeats: usize) -> Structure {
    let length = 4.0 * repeats as f64;
    let mut builder = StructureBuilder::new("grid", Lattice::cubic(length).unwrap());
    let mut index = 0;
    for x in 0..repeats {
        for y in 0..repeats {
            for z in 0..repeats {
                let frac = [
                    x as f64 / repeats as f64,
                    y as f64 / repeats as f64,
                    z as f64 / repeats as f64,
                ];
                builder = builder.site(Site::new(format!("C{index}"), "C", frac).unwrap());
                index += 1;
            }
        }
    }
    builder.build().unwrap()
}

fn median(mut samples: Vec<Duration>) -> Duration {
    samples.sort_unstable();
    samples[samples.len() / 2]
}

fn main() {
    let options = SymmetryOptions {
        symprec: 1e-3,
        ..SymmetryOptions::default()
    };
    for repeats in [2_usize, 3, 4] {
        let structure = grid(repeats);
        let mut samples = Vec::new();
        for _ in 0..7 {
            let start = Instant::now();
            let dataset = detect_symmetry(&structure, options).unwrap();
            std::hint::black_box(dataset);
            samples.push(start.elapsed());
        }
        println!(
            "{:>4} atoms  symmetry {:.2} ms",
            structure.sites().len(),
            median(samples).as_secs_f64() * 1e3
        );
    }
}
