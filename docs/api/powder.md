# Powder refinement

The powder layer lives in `xtalforge-analysis`: measured patterns, peaks,
indexing, Le Bail/Pawley extraction and a least-squares structural refinement.

## rietx inventory

Task 7.1 is a study, not a port. `rietx_inventory()` returns a `RietxBridge`
with the documented capabilities (pattern IO, background, profile, Le Bail,
Pawley, Rietveld, indexing), the MIT license and the source note. The acceptance
metrics are `AcceptanceMetrics { max_rwp, max_position_error }` and
`check_acceptance` rejects values outside them. The exact rietx sources are not
vendored; `xtalforge-powder` will port what is needed when that crate starts.

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
3, tetragonal/hexagonal 2, cubic 1; trigonal has both settings).
`extinction_symbol(number)` returns the Bravais extinction letter for the 14
Bravais types.

## Automatic run

`run_automatic(structure, observed, &AutomaticOptions)` finds the peaks, indexes
a cubic cell, extracts Le Bail intensities and fits the phase, recording every
step in order. It refuses fewer than three points or peaks and reports the
cubic-only limitation in the module documentation.

## Surfaces

`xtalforge powder {index,lebail,pawley,rietveld,automatic,cells,rietx}`;
Python `powder_automatic`, `powder_refine`, `powder_cells`, `rietx`; WASM
`powder_automatic_info`, `powder_refine_info`, `powder_cells_info`,
`rietx_inventory_info`.
