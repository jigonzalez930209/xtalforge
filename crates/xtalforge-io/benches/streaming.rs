//! Process-level benchmark for the streaming readers (Task 15.2.2).
//!
//! Run with `cargo bench -p xtalforge-io --bench streaming`. The numbers are
//! the median of five runs for an in-memory multi-frame XYZ and a large CIF.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::io::Cursor;
use std::time::{Duration, Instant};

use xtalforge_io::{XyzStream, read_cif, read_cif_stream};

fn median(mut samples: Vec<Duration>) -> Duration {
    samples.sort_unstable();
    samples[samples.len() / 2]
}

fn median_duration(mut run: impl FnMut()) -> Duration {
    let mut samples = Vec::with_capacity(5);
    for _ in 0..5 {
        let start = Instant::now();
        run();
        samples.push(start.elapsed());
    }
    median(samples)
}

fn xyz_text(frames: usize, atoms: usize) -> String {
    let mut text = String::new();
    for frame in 0..frames {
        text.push_str(&format!("{atoms}\nframe {frame}\n"));
        for atom in 0..atoms {
            text.push_str(&format!(
                "C {:.6} {:.6} {:.6}\n",
                (atom % 100) as f64 * 0.1,
                (atom % 97) as f64 * 0.1,
                (atom % 89) as f64 * 0.1
            ));
        }
    }
    text
}

fn cif_text(atoms: usize) -> String {
    let mut text = String::from(
        "data_bench\n_cell_length_a 100\n_cell_length_b 100\n_cell_length_c 100\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n",
    );
    for atom in 0..atoms {
        text.push_str(&format!(
            "C{atom} C {:.6} {:.6} {:.6}\n",
            (atom % 1000) as f64 / 1000.0,
            (atom % 991) as f64 / 991.0,
            (atom % 983) as f64 / 983.0
        ));
    }
    text
}

fn main() {
    for (frames, atoms) in [(10_usize, 10_000_usize), (100, 1_000)] {
        let text = xyz_text(frames, atoms);
        let elapsed = median_duration(|| {
            let mut stream = XyzStream::new(Cursor::new(text.as_bytes()));
            let mut count = 0;
            while let Some(frame) = stream.next_frame().unwrap() {
                count += frame.len();
            }
            std::hint::black_box(count);
        });
        println!(
            "{frames:>4} frames x {atoms:>6} atoms  XYZ stream  {:.1} MB/s",
            text.len() as f64 / elapsed.as_secs_f64() / 1e6
        );
    }

    let text = cif_text(20_000);
    let whole = median_duration(|| {
        std::hint::black_box(read_cif(&text).unwrap());
    });
    let stream = median_duration(|| {
        std::hint::black_box(read_cif_stream(Cursor::new(text.as_bytes())).unwrap());
    });
    let megabytes = text.len() as f64 / 1e6;
    println!(
        "20,000 atoms  CIF whole {:.1} ms  CIF stream {:.1} ms  ({megabytes:.1} MB)",
        whole.as_secs_f64() * 1e3,
        stream.as_secs_f64() * 1e3
    );
}
