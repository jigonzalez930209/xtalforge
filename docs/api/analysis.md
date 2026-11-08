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
`RietveldPhase { scale, background, rwp, pattern }`. Atomic positions and ADPs
are the caller's model; a full refinement loop is later Task 7.8.

## Le Bail

`le_bail(structure, observed, &LebailOptions)` partitions the observed
intensity among overlapping reflections and iterates the profile sum, returning
`LebailResult { intensities, positions, pattern, r_factor }`.

## Cubic indexing

`index_cubic(peaks, wavelength, max_index, tolerance)` searches cubic `a`,
indexes each peak and returns `CubicIndex { a, peaks, completeness,
mean_error }`. It is the first indexing engine; non-cubic cells and figures of
merit are later Task 7.5 work.

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

## Nets

`quotient_graph(&CgdEntry, tolerance)` reduces a CGD net to nodes, edges and
lattice shifts; `NetInvariants` reports node/edge counts, cycle rank and
coordination, and `shift_map` groups shifts by node pair. `NetIndex` indexes a CGD file and searches by invariants; `components` finds
connected components (interpenetration candidates) and `interpenetrate` duplicates
a net at a shift. RCSR identification and
canonical keys are later Task 5.6 work.

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
numbers. `default_radii` uses the element van der Waals radii.

All operations return typed errors for non-positive spacing, missing radii or
invalid probe values.
