//! Process-level benchmark for UFF typing and energy (Task 8.6.1).
//!
//! Run with `cargo bench -p xtalforge-ff --bench uff`.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::time::{Duration, Instant};

use xtalforge_ff::{CalculatorOptions, UffModel};
use xtalforge_io::read_cif;

fn median(mut samples: Vec<Duration>) -> Duration {
    samples.sort_unstable();
    samples[samples.len() / 2]
}

fn main() {
    let text = include_str!("../../../tests/corpus/cod/MOF-5.cif");
    let structure = read_cif(text).unwrap();
    let mut samples = Vec::new();
    let mut energy = 0.0;
    for _ in 0..5 {
        let start = Instant::now();
        let model = UffModel::build(&structure, CalculatorOptions::default()).unwrap();
        energy = model.energy().unwrap();
        std::hint::black_box(&model);
        samples.push(start.elapsed());
    }
    println!(
        "{:>5} sites  UFF build+energy {:.1} ms  ({energy:.3} eV)",
        structure.sites().len(),
        median(samples).as_secs_f64() * 1e3
    );

    let methane = "5\nmethane\nC 0 0 0\nH 0.63 0.63 0.63\nH 0.63 -0.63 -0.63\nH -0.63 0.63 -0.63\nH -0.63 -0.63 0.63\n";
    let methane = xtalforge_io::read_xyz(methane).unwrap();
    let mut samples = Vec::new();
    for _ in 0..100 {
        let start = Instant::now();
        let model = UffModel::build(&methane, CalculatorOptions::default()).unwrap();
        std::hint::black_box(model.energy().unwrap());
        samples.push(start.elapsed());
    }
    println!(
        "    5 atoms  methane build+energy {:.3} ms",
        median(samples).as_secs_f64() * 1e3
    );
}
