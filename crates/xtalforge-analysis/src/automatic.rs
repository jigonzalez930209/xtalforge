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
//!   For each of the first [`AutomaticOptions::leading_cells`] candidates the
//!   extinction classes of its crystal system are fitted with Le Bail and
//!   ranked by BIC against the absence-free lattice class; the leading
//!   [`AutomaticOptions::classes`] classes are then Pawley fitted. The table
//!   is ranked by the Pawley Rwp, which is the reference `b_indy -> c_paw`
//!   order, and every row names its class and the groups that share it.
//! * The Rietveld step uses [`crate::rietveld::fit_phase`], which refines the
//!   scale and a constant background, not the structural parameters; the
//!   structural refinement lives in [`crate::refine`] and is a separate step.

#![allow(clippy::indexing_slicing)] // bounded local row slices

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};
use xtalforge_io::xy::PatternPoint;

use crate::extinction::{ExtinctionClass, SpaceGroupMatch, extinction_classes};
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
    /// Leading cells that receive an extinction-class ranking.
    pub leading_cells: usize,
    /// Leading surviving classes Pawley fitted per cell.
    pub classes: usize,
}

impl Default for AutomaticOptions {
    fn default() -> Self {
        Self {
            max_peaks: 30,
            index_tolerance: 0.1,
            min_completeness: 0.6,
            leading_cells: 5,
            classes: 3,
        }
    }
}

/// One row of the automatic leading-cell table.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LeadingCell {
    /// Indexing candidate: metric family, cell, peaks and figures of merit.
    pub index: IndexedCell,
    /// Rank of the indexing cell (1 is best).
    pub cell_rank: usize,
    /// Rank of the extinction class within its cell (1 is best).
    pub class_rank: usize,
    /// Hermann-Mauguin symbol of the class representative.
    pub class_symbol: String,
    /// BIC of the class minus the BIC of the absence-free lattice class.
    pub delta_bic: f64,
    /// Pawley Rwp of this (cell, class), when the fit ran.
    pub pawley_rwp: Option<f64>,
    /// Pawley cell esds, when the fit ran.
    pub pawley_cell_esds: Option<[f64; 6]>,
    /// Groups that share this extinction class.
    pub space_groups: Vec<SpaceGroupMatch>,
}

/// Two Pawley fits whose Rwp are this close are a tie, and the indexing and
/// class order breaks it.
const RWP_TIE: f64 = 0.01;

/// Ranks the table: Rwp first, ties in indexing and class order, failed last.
fn ranked_rows(mut rows: Vec<LeadingCell>) -> Vec<LeadingCell> {
    rows.sort_by(|first, second| {
        first
            .pawley_rwp
            .unwrap_or(f64::INFINITY)
            .total_cmp(&second.pawley_rwp.unwrap_or(f64::INFINITY))
    });
    let mut out = Vec::with_capacity(rows.len());
    let mut cursor = 0usize;
    while cursor < rows.len() {
        let Some(best) = rows[cursor].pawley_rwp else {
            break;
        };
        let limit = best * (1.0 + RWP_TIE);
        let start = cursor;
        while cursor < rows.len() && rows[cursor].pawley_rwp.is_some_and(|value| value <= limit) {
            cursor += 1;
        }
        let mut group = rows[start..cursor].to_vec();
        group.sort_by_key(|row| (row.cell_rank, row.class_rank));
        out.extend(group);
    }
    let mut failed = rows[cursor..].to_vec();
    failed.sort_by_key(|row| (row.cell_rank, row.class_rank));
    out.extend(failed);
    out
}

/// Median FWHM of the detected lines, used by the class fits.
fn median_fwhm(peaks: &[Peak]) -> f64 {
    let mut widths: Vec<f64> = peaks
        .iter()
        .map(|peak| peak.fwhm)
        .filter(|value| value.is_finite() && *value > 0.0)
        .collect();
    widths.sort_by(|first, second| first.total_cmp(second));
    widths
        .get(widths.len() / 2)
        .copied()
        .unwrap_or(0.1)
        .clamp(0.01, 0.5)
}

/// A structure with one general-position atom in the class's group.
fn class_structure(cell: Lattice, class: &ExtinctionClass, name: &str) -> Result<Structure, Error> {
    let mut builder =
        StructureBuilder::new(name, cell).hall(class.representative.hall_symbol.clone());
    for operation in &class.representative.operations {
        builder = builder.operation(*operation);
    }
    builder
        .site(Site::new("C1", "C", [0.137, 0.281, 0.419])?)
        .build()
}

/// BIC of one Le Bail fit, using the weighted residual sum.
fn profile_bic(result: &LebailResult, observed: usize) -> f64 {
    let n = observed.max(1) as f64;
    let rss = result.weighted_sum.max(1e-12);
    let parameters = (result.intensities.len() + 2) as f64;
    n * (rss / n).ln() + parameters * n.ln()
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
    if options.classes == 0 {
        return Err(Error::InvalidInput(
            "automatic run needs at least one extinction class".to_owned(),
        ));
    }
    let mut peaks = find_peaks(observed, &PeakOptions::noise_limited(observed))?;
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
    let fwhm = median_fwhm(&peaks);
    for (cell_rank, candidate) in indexing
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
        let hkl: Vec<[i32; 3]> = candidate.peaks.iter().map(|peak| peak.hkl).collect();
        let classes = extinction_classes(candidate.system, &hkl)?;
        let reference = classes
            .iter()
            .find(|class| class.absences == 0)
            .and_then(|class| class_structure(cell, class, "reference").ok())
            .and_then(|structure| {
                le_bail(
                    &structure,
                    observed,
                    &LebailOptions {
                        fwhm,
                        ..LebailOptions::default()
                    },
                )
                .ok()
            })
            .map(|result| profile_bic(&result, observed.len()))
            .unwrap_or(0.0);
        let mut alive = Vec::new();
        for class in classes.into_iter().filter(|class| !class.refuted) {
            let Ok(structure) = class_structure(cell, &class, "class") else {
                continue;
            };
            let Ok(fit) = le_bail(
                &structure,
                observed,
                &LebailOptions {
                    fwhm,
                    ..LebailOptions::default()
                },
            ) else {
                continue;
            };
            let delta_bic = profile_bic(&fit, observed.len()) - reference;
            alive.push((delta_bic, class));
        }
        alive.sort_by(|first, second| first.0.total_cmp(&second.0));
        for (class_rank, (delta_bic, class)) in alive.into_iter().take(options.classes).enumerate()
        {
            let Ok(structure) = class_structure(cell, &class, "pawley") else {
                continue;
            };
            match pawley_refine(
                &structure,
                observed,
                &PawleyOptions {
                    fwhm,
                    background_order: 2,
                    ..PawleyOptions::default()
                },
            ) {
                Ok(fit) => {
                    if pawley
                        .as_ref()
                        .is_none_or(|best: &PawleyResult| fit.rwp < best.rwp)
                    {
                        pawley = Some(fit.clone());
                    }
                    leading_cells.push(LeadingCell {
                        index: candidate.clone(),
                        cell_rank: cell_rank + 1,
                        class_rank: class_rank + 1,
                        class_symbol: class.symbol,
                        delta_bic,
                        pawley_rwp: Some(fit.rwp),
                        pawley_cell_esds: Some(fit.cell_esds),
                        space_groups: class.groups,
                    });
                }
                Err(error) => {
                    if pawley_error.is_none() {
                        pawley_error = Some(error.message());
                    }
                    leading_cells.push(LeadingCell {
                        index: candidate.clone(),
                        cell_rank: cell_rank + 1,
                        class_rank: class_rank + 1,
                        class_symbol: class.symbol,
                        delta_bic,
                        pawley_rwp: None,
                        pawley_cell_esds: None,
                        space_groups: class.groups,
                    });
                }
            }
        }
    }
    leading_cells = ranked_rows(leading_cells);
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
