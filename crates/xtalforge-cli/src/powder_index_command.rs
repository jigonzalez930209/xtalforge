//! `xtalforge powder index` (Task 7.5 and 17.1.7).

use serde_json::json;
use xtalforge::analysis::{
    IndexingOptions, PeakOptions, find_peaks, index_cubic, index_pattern, indexing_agreement,
};

use crate::{Cli, CliError, input};

pub(crate) fn pattern(path: &str) -> Result<Vec<xtalforge::io::PatternPoint>, CliError> {
    let text = input::read_source(path).map_err(CliError::core)?;
    xtalforge::io::read_xy(&text).map_err(CliError::core)
}

#[allow(clippy::too_many_arguments)] // one flag per documented search switch
pub(crate) fn index(
    cli: &Cli,
    path: &str,
    wavelength: f64,
    max_index: i32,
    tolerance: f64,
    general: bool,
    background: bool,
    refine_peaks: bool,
) -> Result<(), CliError> {
    let points = pattern(path)?;
    // `--background` also applies the noise threshold, so a noisy pattern on a
    // raised baseline is not indexed on its noise maxima.
    let peak_options = PeakOptions {
        refine: refine_peaks,
        ..if background {
            PeakOptions::noise_limited(&points)
        } else {
            PeakOptions::default()
        }
    };
    let peaks = find_peaks(&points, &peak_options).map_err(CliError::analysis)?;
    let cubic =
        index_cubic(&peaks, wavelength, max_index, tolerance).map_err(CliError::analysis)?;
    let candidates = if general {
        Some(
            index_pattern(
                &peaks,
                wavelength,
                &IndexingOptions {
                    max_index,
                    tolerance,
                    ..IndexingOptions::default()
                },
            )
            .map_err(CliError::analysis)?,
        )
    } else {
        None
    };
    let agreement = if general {
        Some(
            indexing_agreement(
                &peaks,
                wavelength,
                &IndexingOptions {
                    max_index,
                    tolerance,
                    ..IndexingOptions::default()
                },
                0.01,
            )
            .map_err(CliError::analysis)?,
        )
    } else {
        None
    };
    if cli.json {
        println!(
            "{}",
            json!({
                "peaks": peaks,
                "index": cubic,
                "candidates": candidates,
                "agreement": agreement,
            })
        );
        return Ok(());
    }
    println!("{} peaks", peaks.len());
    match cubic {
        Some(index) => {
            println!(
                "a = {:.6} A, completeness {:.2}, mean error {:.4} deg",
                index.a, index.completeness, index.mean_error
            );
        }
        None => println!("no cubic candidate"),
    }
    if let Some(candidates) = candidates {
        for (rank, candidate) in candidates.candidates.iter().take(5).enumerate() {
            println!(
                "{}. {} a={:.4} b={:.4} c={:.4} beta={:.2} M20={:.1} F_N={:.1} coverage={:.2}",
                rank + 1,
                candidate.system.name(),
                candidate.cell.a,
                candidate.cell.b,
                candidate.cell.c,
                candidate.cell.beta,
                candidate.m20,
                candidate.f_n,
                candidate.completeness
            );
        }
    }
    if let Some(agreement) = agreement {
        println!(
            "agreement: {} ({})",
            if agreement.agree { "yes" } else { "no" },
            agreement.message
        );
    }
    Ok(())
}
