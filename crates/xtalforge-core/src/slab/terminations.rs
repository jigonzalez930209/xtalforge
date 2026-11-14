//! Termination cuts and the slabs they produce (Task 17.3.2).
//!
//! The automatic cuts are the midpoints between consecutive atomic planes in
//! the surface basis, which is the behaviour of Task 9.4. [`TerminationOptions`]
//! adds explicit cut shifts and cuts inside every interlayer gap, and records
//! how many layers below the top the adsorption-site search reaches.

use crate::{Error, Lattice, Structure, make_slab, surface_basis};

use super::{Slab, gcd, multiply_row};

/// Height tolerance, in units of one interplanar spacing, that merges two
/// atomic planes into one layer by default.
const DEFAULT_LAYER_TOLERANCE: f64 = 1e-4;

/// Distance below which two cut shifts are the same cut.
const CUT_TOLERANCE: f64 = 1e-9;

/// Options for slab termination cuts and adsorption-site depth.
///
/// The defaults reproduce the Task 9.4 behaviour: cuts at the midpoints
/// between consecutive atomic planes and adsorption sites on the top layer
/// only.
#[derive(Debug, Clone, PartialEq)]
pub struct TerminationOptions {
    /// Extra cut shifts in units of one interplanar spacing.
    ///
    /// Each finite value is wrapped into `[0, 1)` and added to the automatic
    /// cuts. Empty by default.
    pub cuts: Vec<f64>,
    /// Add cuts at the quarter and three-quarter points of every gap between
    /// consecutive atomic planes, in addition to its midpoint.
    pub include_intra_layer: bool,
    /// Height tolerance, in units of one interplanar spacing, used to merge
    /// atomic planes into one layer when the cuts are generated.
    pub layer_tolerance: f64,
    /// Atomic layers below the top that produce adsorption sites.
    ///
    /// `1` keeps the top layer only, which is the default; `2` adds the layer
    /// immediately below it, and so on.
    pub site_layers: u32,
}

impl Default for TerminationOptions {
    fn default() -> Self {
        Self {
            cuts: Vec::new(),
            include_intra_layer: false,
            layer_tolerance: DEFAULT_LAYER_TOLERANCE,
            site_layers: 1,
        }
    }
}

/// One termination: the cut shift and the slab it produced.
#[derive(Debug, Clone, PartialEq)]
pub struct SlabTermination {
    /// Cut shift in units of one interplanar spacing, wrapped into `[0, 1)`.
    pub shift: f64,
    /// The slab built at that shift.
    pub slab: Slab,
}

/// All distinct termination cuts of one `(h k l)` surface.
#[derive(Debug, Clone, PartialEq)]
pub struct SlabTerminationReport {
    /// Reduced Miller index.
    pub hkl: [i32; 3],
    /// Terminations sorted by cut shift.
    pub terminations: Vec<SlabTermination>,
}

/// Generates the deterministic cut shifts of one `(h k l)` surface.
///
/// The automatic cuts sit at the midpoint of every gap between consecutive
/// atomic planes, including the gap that wraps across the cell boundary.
/// [`TerminationOptions::include_intra_layer`] adds the quarter and
/// three-quarter points of every gap; [`TerminationOptions::cuts`] adds the
/// caller's shifts. The result is sorted and free of duplicates.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty structure, invalid options or
/// a plane with no layers, or the errors of [`surface_basis`].
pub fn slab_cut_shifts(
    structure: &Structure,
    hkl: [i32; 3],
    options: &TerminationOptions,
) -> Result<Vec<f64>, Error> {
    validate(options)?;
    let basis = surface_basis(hkl)?;
    let heights = plane_heights(structure, &basis, options.layer_tolerance)?;
    let Some(&first_height) = heights.first() else {
        return Err(Error::InvalidInput(
            "the structure has no layers in this direction".to_owned(),
        ));
    };
    let mut cuts = Vec::with_capacity(heights.len() * 3 + options.cuts.len());
    for (index, &lower) in heights.iter().enumerate() {
        let upper = heights
            .get(index + 1)
            .copied()
            .unwrap_or(first_height + 1.0);
        let gap = upper - lower;
        if !gap.is_finite() || gap <= 0.0 {
            continue;
        }
        cuts.push(((lower + upper) / 2.0).rem_euclid(1.0));
        if options.include_intra_layer {
            cuts.push((lower + 0.25 * gap).rem_euclid(1.0));
            cuts.push((lower + 0.75 * gap).rem_euclid(1.0));
        }
    }
    for cut in &options.cuts {
        cuts.push(cut.rem_euclid(1.0));
    }
    cuts.sort_by(|first, second| {
        first
            .partial_cmp(second)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    cuts.dedup_by(|first, second| (*first - *second).abs() <= CUT_TOLERANCE);
    Ok(cuts)
}

/// Builds every termination of one `(h k l)` surface with default options.
///
/// # Errors
///
/// Returns the errors of [`slab_terminations_with_options`].
pub fn slab_terminations(
    structure: &Structure,
    hkl: [i32; 3],
    layers: u32,
    vacuum: f64,
) -> Result<SlabTerminationReport, Error> {
    slab_terminations_with_options(
        structure,
        hkl,
        layers,
        vacuum,
        &TerminationOptions::default(),
    )
}

/// Builds every termination of one `(h k l)` surface with explicit options.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty structure, a zero layer count,
/// a negative vacuum or invalid options, or the errors of [`make_slab`].
pub fn slab_terminations_with_options(
    structure: &Structure,
    hkl: [i32; 3],
    layers: u32,
    vacuum: f64,
    options: &TerminationOptions,
) -> Result<SlabTerminationReport, Error> {
    if layers == 0 {
        return Err(Error::InvalidInput(
            "a termination study needs at least one layer".to_owned(),
        ));
    }
    if !vacuum.is_finite() || vacuum < 0.0 {
        return Err(Error::InvalidInput(
            "the vacuum cannot be negative".to_owned(),
        ));
    }
    let shifts = slab_cut_shifts(structure, hkl, options)?;
    let mut terminations = Vec::with_capacity(shifts.len());
    for shift in shifts {
        let slab = make_slab(structure, hkl, layers, vacuum, shift)?;
        terminations.push(SlabTermination { shift, slab });
    }
    let gcd_value = gcd(gcd(hkl[0].abs(), hkl[1].abs()), hkl[2].abs()).max(1);
    Ok(SlabTerminationReport {
        hkl: [hkl[0] / gcd_value, hkl[1] / gcd_value, hkl[2] / gcd_value],
        terminations,
    })
}

/// Fractional heights of the atomic planes in the surface basis.
fn plane_heights(
    structure: &Structure,
    basis: &[[i32; 3]; 3],
    tolerance: f64,
) -> Result<Vec<f64>, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput("the structure has no atoms".to_owned()));
    }
    let basis_float = basis.map(|row| row.map(f64::from));
    let inverse = Lattice::new(basis_float)?.inverse();
    let mut heights: Vec<f64> = p1
        .atoms
        .iter()
        .map(|atom| multiply_row(atom.frac, inverse)[2].rem_euclid(1.0))
        .collect();
    heights.sort_by(|first, second| {
        first
            .partial_cmp(second)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    heights.dedup_by(|first, second| (*first - *second).abs() <= tolerance);
    Ok(heights)
}

pub(crate) fn validate(options: &TerminationOptions) -> Result<(), Error> {
    if options.cuts.iter().any(|cut| !cut.is_finite()) {
        return Err(Error::InvalidInput(
            "every cut shift must be finite".to_owned(),
        ));
    }
    if !options.layer_tolerance.is_finite() || options.layer_tolerance <= 0.0 {
        return Err(Error::InvalidInput(
            "the layer tolerance must be positive and finite".to_owned(),
        ));
    }
    if options.site_layers == 0 {
        return Err(Error::InvalidInput(
            "the adsorption search needs at least one layer".to_owned(),
        ));
    }
    Ok(())
}
