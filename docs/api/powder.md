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

## Whole-range and structural refinement

`refine_phase(structure, observed, &RietveldOptions, &RefineOptions)` fits the
scale and background analytically, then iterates a finite-difference
least-squares step over the requested fractional coordinates, a global
isotropic ADP multiplier and occupancies. It refuses to refine more parameters
than the observed peaks support, and returns
`RefinedPhase { structure, rwp, parameters, converged, iterations }`.
`combined_objective(energy, rwp, weight)` is the energy-plus-fit objective of
Task 7.9.1.

## Cells and extinction

`free_cell_parameters(number)` returns the free cell parameters for the
crystal system of a space-group number (triclinic 6, monoclinic 4, orthorhombic
3, tetragonal/hexagonal 2, cubic 1); `TrigonalSetting` and
`free_cell_parameters_with_setting` choose the trigonal setting, and
`CELL_PARAMETER_NAMES` holds the names. `extinction_symbol(number)` returns
the Bravais extinction letter for the 14 Bravais types.

## Automatic run

`run_automatic(structure, observed, &AutomaticOptions)` finds the peaks, runs
the general metric-family indexing, extracts Le Bail intensities and fits the
phase, recording every step in order in `AutomaticResult { peaks, index,
lebail, phase, steps }`. It refuses fewer than three points or peaks.
`IndexedPeak` is one indexed reflection.

## Indexing

`index_pattern(peaks, wavelength, &IndexingOptions)` solves the cubic,
tetragonal, hexagonal, orthorhombic and monoclinic families, refines each
candidate and ranks them by coverage, symmetry and M20, with F_N and the
`coverage` fraction reported. `figures_of_merit`, `predict_reflections` and
`match_peaks` expose the individual pieces. The definitions and the ranking
rules are on the [indexing algorithm page](../algorithms/indexing.md).

## Surfaces

`xtalforge powder {index,lebail,pawley,rietveld,automatic,cells,rietx}`;
`powder index --general` prints the ranked metric families and
`--background`/`--refine-peaks` control peak detection. Python
`powder_index`, `powder_automatic`, `powder_refine`, `powder_cells`, `rietx`;
WASM `powder_index_info`, `powder_automatic_info`, `powder_refine_info`,
`powder_cells_info`, `rietx_inventory_info`.
