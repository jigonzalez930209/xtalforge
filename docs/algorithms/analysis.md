# Analysis: patterns, peaks, grids, surfaces and nets

This page covers the analysis operations that do not have a dedicated page:
the calculated pattern, peak finding, the distance grid, isosurfaces, net
merge/overlay, the Pareto front and the Zeo++ bridge. Powder refinement has its
own pages ([indexing](./indexing), [Pawley](./pawley), [Le Bail](./lebail),
[Rietveld](./rietveld)).

## Calculated PXRD

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
model, not an instrument: the FWHM is an input and no axial divergence or
sample broadening is applied.

## Peaks

`find_peaks` smooths the pattern, estimates a background envelope and reports
local maxima above a prominence threshold; `fit_peaks` refines position,
height, width and a background term. The peak list is the input to indexing.
The algorithms are deterministic and do not use a random initialisation.

## Distance grid

The cell is sampled at approximately the requested spacing. For each sample the
distance to the nearest atom surface is stored (van der Waals radius subtracted,
minimum image over the periodic images). From the grid:

- accessible volume: count of samples whose distance is at least the probe
  radius, times the sample volume;
- surface area: area of the boundary between accessible and inaccessible
  samples;
- largest included diameter: twice the largest stored distance;
- pore-size bins: histogram of the distance values with a bin width.

Cost is `prod(ceil(L_i / spacing)) * N_atoms`; the limits check refuses a grid
above the budget before allocation. Validation: the CoRE MOF subset (HKUST-1
and ZIF-8) checks a bounded accessible fraction and a positive surface area;
the reference values are in `tests/external/core_mof/README.md`.

## Isosurfaces

`isosurface` runs marching cubes over a periodic scalar grid. The mesh is a
triangle list in Cartesian coordinates with normals; the periodic boundaries
are closed by sampling with wraparound. Used by the gyroid/Schwarzite carbon
builders and exposed as `TriangleMesh`.

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
