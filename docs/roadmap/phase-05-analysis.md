# Phase 5 — Analysis

**Goal:** calculated powder diffraction, the distance grid and everything read
off it (accessible volume, surface area, pore sizes, accessibility of voids),
periodic isosurfaces, the Zeo++ bridge, nets as periodic graphs with canonical
keys and RCSR identification, net search, merged-net decomposition,
interpenetration generation and detection, overlays, and the module system
that turns analyses into runs with reports.

**Reference:** `Crystal-Builder/xtal/analysis/` (pxrd 655 lines, grid,
voids 481, isosurface, porosity 1022 (Zeo++ result parsing and records),
components, topology 1484, rcsr 1249, netsearch, interpenetrate 612, overlays,
data/ with the RCSR snapshots), `xtal/modules/` (registry, job, record, report
987, stub, pxrd, zeopp 995, poregrid, net), `xtal/references.py`; sci-form
`src/materials/{cell.rs,scattering_factors.rs}`, `src/surface/sasa.rs`.

Crates: `xtalforge-analysis`, `xtalforge-modules`.

---

## Task 5.1 — Calculated powder diffraction

- [ ] 5.1.1 Reflections
  Enumerate (hkl) up to a maximum 2θ for the chosen wavelength; d-spacings;
  systematic absences from the space group operations (not from hard-coded
  rules); multiplicities by grouping symmetry-equivalent reflections (Friedel
  pairs merged).
- [ ] 5.1.2 Structure factors
  Atomic form factors from Cromer-Mann/Waasmaier-Kirfel coefficients (port
  sci-form `scattering_factors.rs`, coefficients to `data/`), occupancies,
  isotropic and anisotropic Debye-Waller factors, |F|².
- [ ] 5.1.3 Intensity corrections and profiles
  Port `analysis/pxrd.py`, whose science came from the same author's
  `DataPlotter` application: Lorentz-polarization factor (with the
  monochromator angle option), multiplicity, Debye-Waller attenuation,
  pseudo-Voigt profile with the Caglioti width relation (U, V, W) and mixing
  parameter; continuous pattern on a 2θ grid; tick marks list with (hkl), d,
  2θ, intensity.
- [ ] 5.1.4 Radiation
  Cu Kα1, Cu Kα1+Kα2 doublet with ratio, Mo, Co, Cr, Fe, Ag and arbitrary
  wavelength; output in 2θ, d or Q.
- [ ] 5.1.5 Module and outputs
  Port `modules/pxrd.py`: the calculation as a module run that leaves the
  pattern as `.xy` (Task 3.5.3) and a reflection table in the run folder, plus
  a report with a `Curve` block (Task 5.10.4).
- [ ] 5.1.6 Validation
  Peak positions and normalized intensities against Crystal-Builder golden
  patterns and published patterns (NaCl, Si, α-quartz, HKUST-1).

## Task 5.2 — Distance grid

- [ ] 5.2.1 Surface-distance field
  Port `analysis/grid.py`: on a regular grid in fractional coordinates, the
  distance from each point to the surface of the nearest atom (centre
  distance minus that atom's radius). Its zero set is the van der Waals
  surface; its `r` level is where the centre of a probe of radius `r` can sit.
  Crystal-Builder computes this itself because Zeo++ 0.3's grid outputs fail
  on real MOFs (documented in the source). Use the neighbor machinery of
  Task 1.10 instead of a k-d tree; parallel over grid slabs with a
  deterministic result.
- [ ] 5.2.2 Radii sets and limits
  Radii sets (UFF van der Waals, CCDC, Zeo++ default radii from
  `tests/data/zeo_radii.cc` once its license is checked); grid spacing with
  the memory estimate of Task 1.3 (a 0.2 Angstrom grid of MFU-4l is 3.8 M
  points, about 2.2 GB in Python).

## Task 5.3 — Accessibility of voids

- [ ] 5.3.1 Accessible versus inaccessible space
  Port `analysis/voids.py` (481 lines): which of the space a probe fits into
  it can reach from outside, by periodic connected-component labelling of the
  grid region where the probe fits, and testing which components percolate
  through the cell (channel dimensionality 0D to 3D).
- [ ] 5.3.2 Pore volume, surface area, pore sizes from the grid
  Port `modules/poregrid.py`: the faster grid-based twins of Zeo++ `-sa` and
  `-vol`/`-volpo` (accessible surface area per volume and per mass, accessible
  and probe-occupiable pore volume, void fraction), and the largest included
  sphere and largest free sphere from the distance grid.
- [ ] 5.3.3 Validation
  Compare with Zeo++ reference values in `Crystal-Builder/tests/data/zeopp_reference.json`
  with tolerances documented per quantity.

## Task 5.4 — Isosurfaces

- [ ] 5.4.1 Periodic marching tetrahedra
  Port `analysis/isosurface.py`: grid plus level plus lattice to points and
  triangles in fractional coordinates. **Marching tetrahedra, not marching
  cubes**: 14 of the 256 cube patterns are ambiguous and wrong choices leave
  holes; tetrahedra have no ambiguous case. Periodic closure across cell
  faces, vertex welding, normals. sci-form's marching cubes is not used for
  this reason (keep it in mind only for non-periodic orbital surfaces).
- [ ] 5.4.2 Mesh output
  Mesh type (vertices, indices, normals) with JSON, OBJ and binary buffer
  output, typed arrays in WASM.

## Task 5.5 — Zeo++ bridge

- [ ] 5.5.1 Result records and parsers
  Port `analysis/porosity.py` (1022 lines; split): parsers and result records
  for Zeo++'s outputs (the three diameters in `.res` — largest included
  sphere, largest free sphere, largest included sphere along the free path —
  surface area, volumes, pore size distribution), with the sentences that say
  what each number means, because these names are often reported wrongly.
  Parsers are tested on captured outputs, so no binary is needed in CI.
- [ ] 5.5.2 Running Zeo++
  Port `modules/zeopp.py` (995 lines; split) on top of the process runner of
  Task 4.9.1: the `network` binary, CSSR input (Task 3.3.2), probe and gas
  defaults, channel and accessibility options. Excluded from WASM.

## Task 5.6 — Nets as periodic graphs

- [ ] 5.6.1 Data and licensing
  The RCSR snapshots shipped by Crystal-Builder (`analysis/data/RCSRnets-2019-06-01.cgd.gz`,
  `rcsr-2019-06-01.json.gz`, `rcsr-transitivity-2026-09-24.json.gz`): verify
  the RCSR terms of use and record them in `THIRD_PARTY.md` before copying
  them into `data/`. If redistribution is not allowed, download on first use
  with a checksum.
- [ ] 5.6.2 Quotient graph
  Port `analysis/topology.py` (1484 lines; split): a net is a finite quotient
  graph whose edges carry lattice translations, up to relabelling of vertices
  and change of basis in GL(d, Z). Built from a structure by contracting
  building units (metal clusters and organic linkers) into vertices, or from a
  CGD file. Only connectivity matters, so the net is unchanged by setting
  changes, supercells and relaxations.
- [ ] 5.6.3 Invariants and canonical key
  Coordination sequences, vertex symbols, point symbols, periodicity
  dimension, and a canonical key (canonical labelling of the quotient graph
  with a normalized translation basis) so equal nets give equal keys.
- [ ] 5.6.4 RCSR identification
  Port `analysis/rcsr.py` (1249 lines): name a net against the RCSR (2665 of
  the 2679 3-periodic nets are distinguished uniquely by the Crystal-Builder
  method; the remaining ones and nets not in the RCSR are reported as such,
  with the nearest candidates). Build step that precomputes the keys of the
  database into a data file (`python -m xtal.analysis.rcsr build` in the
  original).
- [ ] 5.6.5 Validation
  Known frameworks from the corpus identified correctly (HKUST-1 `tbo`, MOF-5
  `pcu`, UiO-66 `fcu`, ZIF-8 `sod`, MOF-74 `msf`; expected names taken from
  the Crystal-Builder golden data, not from memory);
  Systre test data in `Crystal-Builder/tests/data/systre` and `rcsr`.

## Task 5.7 — Net search and merged nets

- [ ] 5.7.1 Net search
  Port `analysis/netsearch.py`: find nets by what a chemist knows:
  coordination ("a 3,6 net"), space group number, transitivity (vertex, edge,
  face, tile counts, from the RCSR data, unknown where the RCSR has none) and
  name. Shared by the MOF builder and the net builder so both lists answer the
  same query identically.
- [ ] 5.7.2 Merged-net components
  Port `analysis/components.py`: break a merged net (Eddaoudi and O'Keeffe;
  two or three edge-transitive nets sharing vertices) into the nets that each
  kind of edge spans.

## Task 5.8 — Interpenetration

- [ ] 5.8.1 Detection
  Count independent subnets of a framework, their translation and symmetry
  relations, and classify by Blatov's scheme: Class Ia (copies related by a
  translation; one net's lattice is an index-n sublattice of the array's),
  Class Ib, Class II (related by a non-translational symmetry operation).
- [ ] 5.8.2 Generation
  Port `analysis/interpenetrate.py` (612 lines): generate an n-fold
  interpenetrated framework by placing copies at coset representatives
  (`pcu`, `dia`, `srs` are the labyrinth nets of the P, D and G minimal
  surfaces, which is why they interpenetrate two-fold readily), check the
  result with the detector, report clashes. The `Interpenetrate` command
  (Task 3.10.7) applies it.

## Task 5.9 — Overlays

- [ ] 5.9.1 Overlay data
  Port `analysis/overlays.py`: what a run can draw over a crystal that is not
  part of it: atomic charges (as per-atom scalars) and an orbital's two lobes
  (two isosurfaces of a cube grid), plus the pore network surface. Output is
  view data (Task 3.14.4), not rendering.

## Task 5.10 — Modules, jobs, run records, reports

- [ ] 5.10.1 Module registry
  Port `modules/registry.py`: what can be run, as data (module name, label,
  actions, parameters from Task 3.8, whether it needs a structure, an
  availability check). The CLI `modules` and `run` commands, the Python API
  and the agent tools all read it.
- [ ] 5.10.2 Jobs
  Port `modules/job.py`: a module's `run` gets one `Job` argument (structure,
  parameters, run folder, cancellation token, progress sink, log) and returns
  a `JobResult` (sentence, optional new structure, report). A `Job` can grow
  fields without breaking existing modules.
- [ ] 5.10.3 Run records
  Port `modules/record.py`: the module that runs writes its run folder
  (inputs, outputs, `run.log`), shared with the workspace (Task 3.12.2).
- [ ] 5.10.4 Reports
  Port `modules/report.py` (987 lines; split): a report is presentation, not
  results; every number in it was already written to the run folder. A few
  frozen block types (key values, table, histogram, curve) and renderers to
  Markdown, JSON and HTML.
- [ ] 5.10.5 Test stub module
  Port `modules/stub.py`: a module that does nothing, so the registry, the
  parameter forms of front ends, cancellation and run folders can be tested.
- [ ] 5.10.6 References
  Port `references.py`: one small type for citations of methods, programs and
  databases, attached to calculators, modules and builders and printed in
  reports and docs.

## Task 5.11 — Surfaces and docs

- [ ] 5.11.1 CLI
  `pxrd`, `grid`, `porosity` (grid-based), `zeopp`, `isosurface`, `net
  {extract,identify,search,components}`, `interpenetration {detect,generate}`,
  `modules`, `run <module> <action>`, `report`.
- [ ] 5.11.2 Python and WASM
  Same functions; NumPy and typed arrays for patterns, grids and meshes; WASM
  excludes Zeo++; long computations report progress and can be cancelled.
- [ ] 5.11.3 Docs and demos
  Algorithm pages: PXRD (with all corrections), distance grid, accessibility
  labelling, marching tetrahedra, quotient graphs and canonical keys, RCSR
  identification, Blatov classes. Demos: PXRD of a pasted CIF; accessible
  surface of a small MOF; net identification.
- [ ] 5.11.4 Parity and release
  Parity matrix; changelog; MINOR release at phase completion.
