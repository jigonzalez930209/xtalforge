//! Automatic peak search, indexing, Le Bail and Rietveld run (Task 7.10).
//!
//! The pipeline copies the order of `powder/auto.py`: find the lines, index the
//! metric families, build the leading-cell table with a Pawley fit of every
//! leading candidate, extract Le Bail intensities, then fit the phase. Every
//! completed step is recorded in [`AutomaticResult::steps`] in execution
//! order, so a CLI or a run folder can replay the decision without guessing
//! what happened.
//!
//! # Approximations
//!
//! * The radiation is fixed to Cu K alpha (`1.54184 Angstrom`) because this
//!   entry point has no radiation argument; callers with another source should
//!   use [`crate::indexing::index_cubic`] and [`crate::lebail::le_bail`]
//!   directly.
//! * The indexed cells are the ranked candidates of the five metric families.
//!   The leading-cell table fits a Pawley model of each of the first
//!   [`AutomaticOptions::leading_cells`] candidates and ranks the table by
//!   the Pawley Rwp, which is the reference `b_indy -> c_paw` order. Space
//!   groups are not assigned from systematic absences; each row keeps its
//!   metric family.
//! * The Rietveld step uses [`crate::rietveld::fit_phase`], which refines the
//!   scale and a constant background, not the structural parameters; the
//!   structural refinement lives in [`crate::refine`] and is a separate step.

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Lattice, Structure};
use xtalforge_io::xy::PatternPoint;

use crate::indexing::{IndexedCell, IndexingOptions, index_pattern};
use crate::lebail::{LebailOptions, LebailResult, le_bail};
use crate::pawley::{PawleyOptions, PawleyResult, pawley_refine};
use crate::peaks::{Peak, PeakOptions, find_peaks};
use crate::pxrd::Radiation;
use crate::rietveld::{RietveldOptions, RietveldPhase, fit_phase};

/// Options for [`run_automatic`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AutomaticOptions {
    /// Maximum number of detected peaks handed to indexing.
    pub max_peaks: usize,
    /// Position tolerance in degrees two-theta for cubic indexing.
    pub index_tolerance: f64,
    /// Minimum fraction of peaks that must index for a cell to be kept.
    pub min_completeness: f64,
    /// Leading cells that receive a Pawley fit in the ranking table.
    pub leading_cells: usize,
}

impl Default for AutomaticOptions {
    fn default() -> Self {
        Self {
            max_peaks: 30,
            index_tolerance: 0.1,
            min_completeness: 0.3,
            leading_cells: 3,
        }
    }
}

/// One row of the automatic leading-cell table.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LeadingCell {
    /// Indexing candidate: metric family, cell, peaks and figures of merit.
    pub index: IndexedCell,
    /// Pawley Rwp of this cell, when the fit ran.
    pub pawley_rwp: Option<f64>,
    /// Pawley cell esds, when the fit ran.
    pub pawley_cell_esds: Option<[f64; 6]>,
}

/// Result of one automatic run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AutomaticResult {
    /// Detected peaks, ordered by increasing two-theta.
    pub peaks: Vec<Peak>,
    /// Best indexed cell, if any candidate passed `min_completeness`.
    pub index: Option<IndexedCell>,
    /// Le Bail extraction.
    pub lebail: Option<LebailResult>,
    /// Pawley fit of the best indexed cell.
    pub pawley: Option<PawleyResult>,
    /// Leading cells ranked by Pawley fit quality.
    pub leading_cells: Vec<LeadingCell>,
    /// Rietveld scale/background phase fit.
    pub phase: Option<RietveldPhase>,
    /// Completed steps, in execution order.
    pub steps: Vec<String>,
}

/// Runs peaks, general indexing, the leading-cell table, Le Bail and a phase
/// fit on one pattern.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for fewer than three points, non-finite or
/// non-monotone observed data, invalid options, fewer than three detected
/// peaks, or an error from one of the underlying steps.
pub fn run_automatic(
    structure: &Structure,
    observed: &[PatternPoint],
    options: &AutomaticOptions,
) -> Result<AutomaticResult, Error> {
    validate_observed(observed)?;
    if options.max_peaks < 2 {
        return Err(Error::InvalidInput(
            "automatic run needs at least two peaks to index".to_owned(),
        ));
    }
    if !options.index_tolerance.is_finite() || options.index_tolerance <= 0.0 {
        return Err(Error::InvalidInput(
            "automatic index tolerance must be positive and finite".to_owned(),
        ));
    }
    if !options.min_completeness.is_finite() || !(0.0..=1.0).contains(&options.min_completeness) {
        return Err(Error::InvalidInput(format!(
            "automatic minimum completeness must be in [0, 1], got {}",
            options.min_completeness
        )));
    }
    if options.leading_cells == 0 {
        return Err(Error::InvalidInput(
            "automatic run needs at least one leading cell".to_owned(),
        ));
    }
    let maximum = observed
        .iter()
        .map(|point| point.intensity)
        .fold(0.0_f64, f64::max);
    let step = median_step(observed);
    let mut peaks = find_peaks(
        observed,
        &PeakOptions {
            min_height: 0.01 * maximum,
            min_distance: (2.0 * step).max(0.01),
            smoothing: 1,
            ..PeakOptions::default()
        },
    )?;
    if peaks.len() < 3 {
        return Err(Error::InvalidInput(format!(
            "automatic indexing needs at least three peaks, found {}",
            peaks.len()
        )));
    }
    peaks.sort_by(|left, right| {
        right
            .intensity
            .total_cmp(&left.intensity)
            .then_with(|| left.two_theta.total_cmp(&right.two_theta))
    });
    peaks.truncate(options.max_peaks);
    peaks.sort_by(|left, right| left.two_theta.total_cmp(&right.two_theta));

    let mut steps = vec!["peaks".to_owned()];
    let wavelength = Radiation::CuKa.wavelength()?;
    let indexing = index_pattern(
        &peaks,
        wavelength,
        &IndexingOptions {
            max_index: 10,
            tolerance: options.index_tolerance,
            min_completeness: options.min_completeness,
            ..IndexingOptions::default()
        },
    )?;
    let index = indexing.best().cloned();
    steps.push("index".to_owned());
    let mut leading_cells = Vec::new();
    let mut pawley = None;
    let mut pawley_error = None;
    for (position, candidate) in indexing
        .candidates
        .iter()
        .take(options.leading_cells)
        .enumerate()
    {
        let cell = Lattice::from_parameters(
            candidate.cell.a,
            candidate.cell.b,
            candidate.cell.c,
            candidate.cell.alpha,
            candidate.cell.beta,
            candidate.cell.gamma,
        )?;
        let pawley_structure = structure.clone().with_lattice(cell);
        match pawley_refine(
            &pawley_structure,
            observed,
            &PawleyOptions {
                background_order: 2,
                ..PawleyOptions::default()
            },
        ) {
            Ok(fit) => {
                if position == 0 {
                    pawley = Some(fit.clone());
                }
                leading_cells.push(LeadingCell {
                    index: candidate.clone(),
                    pawley_rwp: Some(fit.rwp),
                    pawley_cell_esds: Some(fit.cell_esds),
                });
            }
            Err(error) => {
                if pawley_error.is_none() {
                    pawley_error = Some(error.message());
                }
                leading_cells.push(LeadingCell {
                    index: candidate.clone(),
                    pawley_rwp: None,
                    pawley_cell_esds: None,
                });
            }
        }
    }
    leading_cells.sort_by(
        |first, second| match (first.pawley_rwp, second.pawley_rwp) {
            (Some(left), Some(right)) => left.total_cmp(&right),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => second.index.m20.total_cmp(&first.index.m20),
        },
    );
    steps.push("leading_cells".to_owned());
    let lebail = le_bail(structure, observed, &LebailOptions::default())?;
    steps.push("le_bail".to_owned());
    if pawley.is_some() {
        steps.push("pawley".to_owned());
    } else if let Some(message) = pawley_error {
        // A cell that indexes the peaks but predicts too many reflections
        // must not sink the rest of the run.
        steps.push(format!("pawley_skipped: {message}"));
    }
    let first = observed.first().map_or(5.0, |point| point.two_theta);
    let last = observed.last().map_or(80.0, |point| point.two_theta);
    let phase = fit_phase(
        structure,
        observed,
        &RietveldOptions {
            two_theta_range: [first.max(0.0), last.min(180.0)],
            step: 0.02,
            ..RietveldOptions::default()
        },
    )?;
    steps.push("rietveld".to_owned());
    Ok(AutomaticResult {
        peaks,
        index,
        lebail: Some(lebail),
        pawley,
        leading_cells,
        phase: Some(phase),
        steps,
    })
}

fn validate_observed(observed: &[PatternPoint]) -> Result<(), Error> {
    if observed.len() < 3 {
        return Err(Error::InvalidInput(
            "automatic run needs at least three observed points".to_owned(),
        ));
    }
    for (index, point) in observed.iter().enumerate() {
        if !point.two_theta.is_finite() || !point.intensity.is_finite() {
            return Err(Error::InvalidInput(format!(
                "observed point {index} must be finite"
            )));
        }
        if let Some(previous) = index.checked_sub(1).and_then(|value| observed.get(value))
            && point.two_theta <= previous.two_theta
        {
            return Err(Error::InvalidInput(format!(
                "observed 2 theta must be strictly increasing at point {index}"
            )));
        }
    }
    Ok(())
}

fn median_step(observed: &[PatternPoint]) -> f64 {
    let mut steps = Vec::new();
    for pair in observed.windows(2) {
        if let [first, second] = pair
            && second.two_theta > first.two_theta
        {
            steps.push(second.two_theta - first.two_theta);
        }
    }
    steps.sort_by(f64::total_cmp);
    steps
        .get(steps.len() / 2)
        .copied()
        .unwrap_or(0.02)
        .max(0.001)
}
