# Powder patterns and porosity

`xtalforge-analysis` starts with calculated PXRD (Task 5.1) and the periodic
surface-distance grid with porosity numbers (Tasks 5.2/5.3).

## PXRD

The Rust entry point is `calculate_pxrd`; Python exposes `pxrd`, WASM
`pxrd_info` and the CLI `xtalforge pxrd`.

```rust
use xtalforge_analysis::{PxrdOptions, Radiation, calculate_pxrd};

let pattern = calculate_pxrd(&structure, &PxrdOptions {
    radiation: Radiation::CuKa,
    ..PxrdOptions::default()
})?;
```

`Radiation` carries the wavelengths (CuKa, CuKa1, MoKa, CoKa, FeKa, CrKa,
AgKa). `calculate_pxrd` enumerates hkl reflections in the requested 2-theta
range, computes `d = 1/|hkl @ reciprocal|`, structure factors from the P1
atoms, Lorentz-polarization and multiplicity, and samples a Gaussian profile.
`Pattern { points, reflections }` carries the sampled curve and the peak list;
`Reflection { hkl, d, two_theta, multiplicity, intensity }`.

## Pareto front

`pareto_front(points)` returns the indices of the non-dominated `(energy,
residual)` points, minimising both objectives.

## Rietveld phase model

`fit_phase(structure, observed, &RietveldOptions)` simulates the phase from the
structure and PXRD engine, then fits its scale and background, returning
`RietveldPhase { scale, background, rwp, pattern }`. The structural
refinement is `refine_phase(structure, observed, &RietveldOptions,
&RefineOptions)` in the [powder refinement](./powder) page; it runs the
McCusker stages (zero/cell, width, then Wyckoff coordinates, a global isotropic
ADP multiplier, occupancies and texture) and returns
`RefinedPhase { structure, rwp, parameters, parameter_esds, stages, converged,
iterations }`. `free_directions(operations, frac)` is the Wyckoff-allowed
subspace helper.

## Le Bail

`le_bail(structure, observed, &LebailOptions)` partitions the observed
intensity among overlapping reflections and iterates the profile sum, returning
`LebailResult { intensities, positions, pattern, r_factor, rwp, chi_squared,
intensity_esds }`. The esds come from the inverse of the profile normal matrix
scaled by the reduced chi-squared; above 128 reflections only the covariance
diagonal is used. See [Le Bail extraction](../algorithms/lebail.md).

## Indexing

`index_pattern(peaks, wavelength, &IndexingOptions)` searches the cubic,
tetragonal, hexagonal, orthorhombic and monoclinic metric families and returns
`IndexingResult { candidates }`, best first. Each `IndexedCell` carries the
family, the refined `CellMetric { a, b, c, alpha, beta, gamma }`, the indexed
peaks, `completeness`, `mean_error`, the de Wolff `m20`, the Smith-Snyder
`f_n` and the prediction `coverage`. The search is bounded by `max_index`,
`max_volume`, `min_completeness` and `max_evaluations`.

`figures_of_merit(&cell, peaks, wavelength, tolerance)` evaluates one cell and
`predict_reflections(&cell, max_q)` / `match_peaks(...)` expose the prediction
and the matching. `index_cubic(peaks, wavelength, max_index, tolerance)` is
the original cubic-only search and keeps returning `CubicIndex`.

Run `xtalforge powder index PATTERN --general` for the ranked table. See
[Powder indexing](../algorithms/indexing.md) for the solve, the figures of
merit and the ranking rules.

## Profiles and background

`Profile { kind, center, fwhm, eta, height }` evaluates Gaussian, Lorentzian
and pseudo-Voigt peaks and samples a range. `Background::{Constant, Linear,
Quadratic}` are the background models.

## Least squares

`least_squares(initial, &FitOptions, residuals)` is a Levenberg-Marquardt fit
with finite-difference Jacobians. `FitResult { parameters, residuals,
chi_squared, iterations, converged }` is the result.

## Measured patterns and peaks

`xtalforge_io::read_xy` reads a measured pattern; `find_peaks(&points,
&PeakOptions)` finds local maxima with optional smoothing, a height threshold
and a minimum separation, returning `Peak { two_theta, intensity, fwhm }`.
`PeakOptions` also turns on a SNIP background subtraction
(`background`, `background_window`) and a per-line Gaussian plus linear
background fit (`refine`, `fit_window`). `background_envelope(&points, window)`
returns the SNIP envelope and `fit_peaks(&points, &positions, &options)` fits
named positions, which is what an indexing run wants to refine.

## Nets

`quotient_graph(&CgdEntry, tolerance)` reduces a CGD net to nodes, edges and
lattice shifts; `NetInvariants` reports node/edge counts, cycle rank and
coordination, and `shift_map` groups shifts by node pair. `NetIndex` indexes a CGD file and searches by invariants; `components` finds
connected components (interpenetration candidates) and `interpenetrate` duplicates
a net at a shift. `RcsrIndex::from_subset()` identifies the curated 17-net
subset by counts, cycle rank and coordination; the full RCSR catalogue and its
canonical keys are not bundled.

## Isosurfaces

`isosurface(&Grid3, value)` extracts a `TriangleMesh { vertices, triangles }`
with marching tetrahedra (six tetrahedra per cell, no ambiguous marching-cube
cases). It rejects grids with fewer than two samples along an axis.

## Distance grid

`distance_grid(structure, radii, spacing)` returns a `DistanceGrid` wrapping
`xtalforge_io::grid::Grid3`. The value is the distance to the nearest atom
surface (negative inside an atom), using minimum images. `contact_surface`,
`accessible_volume(probe)`, `pore_volume`, `surface_area`,
`largest_included_diameter` and `pore_diameter(bin_width)` are the porosity
numbers; `PoreSizeBin` is one histogram bin. `default_radii` uses the element
van der Waals radii. `PatternPoint` is the shared measured-pattern point type,
and `ProfileKind` selects Gaussian, Lorentzian or pseudo-Voigt.

All operations return typed errors for non-positive spacing, missing radii or
invalid probe values.

## Zeo++, RCSR and overlays

`run_zeo(structure, &ZeoOptions)` probes `zeo++` and parses its `Channels`
output (`parse_zeo_channels`, `parse_zeo_volume`); without the binary it
returns a typed error. `write_zeo_inputs` writes the CSSR and the network
arguments.

`RcsrIndex::from_subset()` holds a curated 17-net subset (the full RCSR
database is an external bundle); `identify` matches by node/edge counts, cycle
rank and coordination and `validate` checks the subset. `QuotientGraph`,
`NetEdge`, `NetInvariants` and `RcsrEntry` are the net types. `merge_nets`
joins two quotient graphs and `overlay_net` turns one into an `Overlay` with an
`OverlayStyle`; `OverlaySet` collects overlays for a viewer.

`xtalforge analysis {porosity,nets,zeo,overlay}`, Python `porosity`,
`net_invariants`, `zeo_parse` and the WASM equivalents expose the same
operations.
