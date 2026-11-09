# Powder refinement

The powder layer lives in `xtalforge-analysis`: measured patterns, peaks,
indexing, Le Bail/Pawley extraction and a least-squares structural refinement.

## rietx inventory

Task 7.1 is a study, not a port. `rietx_inventory()` returns a `RietxBridge`
with the documented capabilities (pattern IO, background, profile, Le Bail,
Pawley, Rietveld, indexing), the MIT license and the source note. The acceptance
metrics are `AcceptanceMetrics { max_rwp, max_position_error }` and
`check_acceptance` rejects values outside them. `RIETX_LICENSE_NOTE`,
`ACCEPTANCE_MAX_RWP` and `ACCEPTANCE_MAX_POSITION_ERROR` are the public
constants. The exact rietx sources are not vendored; the native implementation
is a simplified model documented in the module comments and the roadmap's
known-gaps section.

## Instrument

`Instrument { geometry, zero_shift, displacement, ka2_ratio, polarization,
sample_radius }` applies a zero shift and the sample-displacement correction,
and `doublet` splits a reflection into the Kα1/Kα2 pair. `Geometry` is
`BraggBrentano` or `DebyeScherrer`. `doublet_checked` validates the ratio.

## Le Bail

`le_bail(structure, observed, &LebailOptions)` iterates the observed-to-
calculated ratio inside every reflection window and returns `LebailResult {
intensities, intensity_esds, positions, pattern, r_factor, rwp, chi_squared }`.
The intensity esds are the diagonal of `(A^T W A)^-1` scaled by the reduced
chi-squared; above 128 reflections the covariance diagonal is used instead of
the exact inverse. See [Le Bail extraction](../algorithms/lebail.md).

## Pawley

`pawley_refine(structure, observed, &PawleyOptions)` fits the whole pattern
with one free intensity per allowed reflection, a Chebyshev background and
optional free cell, zero-shift, size and strain parameters. It returns
`PawleyResult { reflections, cell, cell_esds, zero_shift, zero_shift_esd, size,
strain, rp, rwp, rexp, chi_squared, pattern, converged, iterations }`; every
`PawleyReflection` carries `hkl`, `d`, `two_theta`, `intensity` and `esd`. The
linear intensities and their covariance are solved exactly for each non-linear
point and the cell/zero/width model is refined by variable projection. See
[Pawley fit](../algorithms/pawley.md) for the conventions and limits.

## Whole-range and structural refinement

`refine_phase(structure, observed, &RietveldOptions, &RefineOptions)` runs the
McCusker stages: the scale and background analytically, then the zero shift and
free cell numbers, the profile width, and the structure. `RefineOptions`
selects them with `refine_positions`, `refine_adps`, `refine_occupancies`,
`refine_zero`, `refine_cell`, `refine_width`, `refine_texture`, `texture_axis`
and `doublet`. Positions move only along the directions
`free_directions(operations, frac)` allows for each Wyckoff site, the cell ties
of the metric family are restored after every move, and a
`March-Dollase` coefficient is refined about `texture_axis`.

It refuses to refine more parameters than the observed peaks support and
returns `RefinedPhase { structure, rwp, parameters, parameter_esds, stages,
converged, iterations }`. `parameter_esds` holds the covariance-based esds of
the free parameters, and `stages` is the order that actually ran.
`combined_objective(energy, rwp, weight)` is the energy-plus-fit objective of
Task 7.9.1. See [Rietveld refinement](../algorithms/rietveld.md) for the
stage details, the Wyckoff rules and the documented deviations.

## Cells and extinction

`free_cell_parameters(number)` returns the free cell parameters for the
crystal system of a space-group number (triclinic 6, monoclinic 4, orthorhombic
3, tetragonal/hexagonal 2, cubic 1); `TrigonalSetting` and
`free_cell_parameters_with_setting` choose the trigonal setting, and
`CELL_PARAMETER_NAMES` holds the names. `extinction_symbol(number)` returns
the Bravais extinction letter for the 14 Bravais types.

## Automatic run

`run_automatic(structure, observed, &AutomaticOptions)` finds the peaks, runs
the general metric-family indexing, builds a leading-cell table, extracts Le
Bail intensities, fits a Pawley model of the best indexed cell and fits the
phase, recording every step in order in `AutomaticResult { peaks, index,
leading_cells, lebail, pawley, phase, steps }`. `leading_cells` holds the first
`AutomaticOptions::leading_cells` indexing candidates (default 3), each with its
`IndexedCell`, Pawley `pawley_rwp` and `pawley_cell_esds`; the table is sorted
by Pawley Rwp, so the best fit is first even when the indexing ranking differs.
The `steps` list contains `peaks`, `index`, `leading_cells`, `le_bail`,
`pawley` and `rietveld` (a Pawley fit the cell makes impossible is recorded as
`pawley_skipped`). It refuses fewer than three points or peaks and a zero
`leading_cells`. `IndexedPeak` is one indexed reflection.

## Indexing

`index_pattern(peaks, wavelength, &IndexingOptions)` solves the cubic,
tetragonal, hexagonal, orthorhombic and monoclinic families, refines each
candidate and ranks them by coverage, symmetry and M20, with F_N and the
`coverage` fraction reported. `figures_of_merit`, `predict_reflections` and
`match_peaks` expose the individual pieces. The definitions and the ranking
rules are on the [indexing algorithm page](../algorithms/indexing.md).

## Surfaces

`xtalforge powder {index,lebail,pawley,rietveld,automatic,cells,rietx}`;
`powder index --general` prints the ranked metric families,
`--background`/`--refine-peaks` control peak detection, `powder pawley`
prints the full fit with Rwp and esds, and `powder rietveld --refine` runs the
McCusker stages (`--occupancies` and `--texture` add those parameters). Python `powder_index`, `powder_pawley`,
`powder_automatic`, `powder_refine`, `powder_cells`, `rietx`; WASM
`powder_index_info`, `powder_pawley_info`, `powder_automatic_info`,
`powder_refine_info`, `powder_cells_info`, `rietx_inventory_info`.
