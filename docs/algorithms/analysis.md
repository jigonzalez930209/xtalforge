# Analysis: patterns, peaks, grids, surfaces and nets

This page covers the analysis operations that do not have a dedicated page:
the calculated pattern, peak finding, the distance grid, isosurfaces, net
merge/overlay, the Pareto front and the Zeo++ bridge. Powder refinement has its
own pages ([indexing](./indexing), [Pawley](./pawley), [Le Bail](./lebail),
[Rietveld](./rietveld)).

## Calculated PXRD

![Bragg's law and the calculated pattern](/diagrams/bragg-pxrd.svg)


For every allowed reflection `(hkl)` with d-spacing `d`:

```
2 theta = 2 asin(lambda / (2 d))
I(hkl)  = |F(hkl)|^2 * multiplicity * LP(2 theta)
I(2 theta) = sum_hkl I(hkl) * G(2 theta - 2 theta_hkl; FWHM)
```

`F` is the structure factor from the site positions, occupancies and atomic
scattering factors; `LP` is the Lorentz-polarisation factor; `G` is a
normalised Gaussian. The pattern is sampled on a regular `2 theta` grid. The
answer is normalised to 100 in the API; the unnormalised reflection intensities
are kept in `Pattern::reflections`.

Validation: the rutile pattern from Crystal-Builder/rietx and the NaCl
calculated pattern (`tests/golden/phase-7`, corpus tests). The profile is a
model, not an instrument: the FWHM is an input (one Gaussian of unit area, no
Caglioti width) and no axial divergence or sample broadening is applied. An
isotropic `U_iso` enters as `exp(-8 pi^2 U s^2)`; anisotropic parameters and
anomalous scattering are ignored, and the scattering-factor fit has no core
tail, so absolute intensities and the high-angle falloff are approximate while
positions, multiplicities and absences are exact.

## Peaks

`find_peaks` reports local maxima above `min_height` and `min_distance`, with
optional smoothing, SNIP background subtraction and a per-line fit; the defaults
(`min_height = 0`, no smoothing, no background) report every maximum, so noisy
data should use `PeakOptions::noise_limited` (background subtracted, lines above
1 % of the strongest and five times the noise); `fit_peaks` refines position,
height, width and a background term. The peak list is the input to indexing.
The algorithms are deterministic and do not use a random initialisation.

## Distance grid

The cell is sampled with `max(1, round(L_i / spacing))` points along each axis,
so anisotropic cells get a similar resolution. For each sample the signed
distance to the surface of the nearest atom is stored
(`centre distance - van der Waals radius`, minimum image over the periodic
images): negative inside an atom, zero on its surface. From the grid:

- accessible volume: count of samples whose stored distance is at least the
  probe radius, times the volume of one sample (converges from below, with a
  boundary-layer error of about one spacing);
- surface area: grid faces whose two ends have opposite sign (outside versus
  inside a van der Waals atom) times the face area. It is the voxelised area of
  the van der Waals surface, not of the probe-accessible boundary, and
  overestimates a smooth surface at coarse spacing;
- largest included diameter: twice the largest stored distance;
- pore-size bins: histogram of `2 * distance` at every free sample with a
  caller-supplied bin width (a geometric distribution, not a Zeo++ result).

Cost is `prod(round(L_i / spacing)) * N_atoms`: every sample is compared with
every atom, which is exact but meant for small and medium cells. The grid
refuses a size above the budget before allocation; the generic
`check_grid` estimate uses the same rounding with a minimum of 8 points per
axis. Validation: the CoRE MOF subset (HKUST-1 and ZIF-8) checks a bounded
accessible fraction and a positive surface area; the reference values are in
`tests/external/core_mof/README.md`.

## Isosurfaces

`isosurface` runs **marching tetrahedra** (six tetrahedra per grid cell) over a
periodic scalar grid; there are no ambiguous marching-cube cases. The four
crossings of a 2-2 sign pattern are ordered around their quadrilateral, so
neighbouring tetrahedra share the same diagonal and the surface is watertight;
vertices are welded by fractional position and the periodic boundaries close by
wraparound. `TriangleMesh { vertices, triangles }` holds Cartesian vertices and
index triples; it stores no normals. Used by the gyroid/Schwarzite carbon
builders.

## Nets

`quotient_graph` contracts a structure to its net: metal clusters and
connection points become vertices, linkers become edges, and the periodic
images become edge offsets. `fingerprint` computes the coordination sequence
and point symbol; `RcsrIndex` matches them against the embedded RCSR 2019-06-01
snapshot (2,931 entries). `merge_nets` unions two quotient graphs;
`overlay_net` writes a display overlay; `pareto_front` keeps the non-dominated
points of a set of candidates.

The RCSR match returns a unique name, an ambiguous set or the nearest
candidates: a fingerprint can be shared, and the API says so instead of
guessing.

## Zeo++ bridge

`write_zeo_inputs` writes the CSSR structure and a radii file;
`run_zeo` launches the external binary and `parse_zeo_channels`/
`parse_zeo_volume` parse its output. Without the binary the writer and parser
still work and the tests use recorded output. This is an optional engine, not a
native implementation.

## Validation

- `crates/xtalforge-analysis/tests/` covers PXRD, peaks, grids, isosurfaces,
  nets, merge/overlay, Pareto, Zeo parsing and the automatic run.
- The CoRE MOF external test covers real frameworks.
- The RCSR index build validates 2,909 keyed entries against their canonical
  keys; the remaining 22 are unkeyed by design.
