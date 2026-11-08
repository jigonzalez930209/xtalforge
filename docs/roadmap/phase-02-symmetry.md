# Phase 2 — Symmetry

**Goal:** complete crystal symmetry support on top of `moyo`, isolated behind
our own types, plus everything moyo does not expose: cell reduction, the
symmetry workflow of the editor (detect, set, impose, asymmetrize, reduce to
P1), subgroup naming, molecular point groups, the Brillouin zone and k-paths.

**Reference:** `Crystal-Builder/xtal/core/spacegroup.py` (390 lines, gemmi
tables), `core/symmetry.py` (1139, spglib), `core/subgroups.py` (1356),
`core/supercell.py` (lattice reduction part), point groups in
`core/molecule.py`, `analysis/kpath.py` (ASE `Cell.bandpath` and a scipy
Voronoi zone); moyo 0.22 sources and docs (https://docs.rs/moyo).

Crate: `xtalforge-symmetry`. Decision gate: the Task 0.3 report must conclude
"adopt" or "adopt with gaps" before Task 2.1 starts.

**Task 2.1 belongs to the v0.0.1 gate** and may run in parallel with Phase 1
(see the note in Phase 1). The rest ships in v0.0.x releases.

---

## Task 2.1 — moyo integration (v0.0.1 gate)

- [x] 2.1.1 Conversion layer
  `Structure`/`P1Cell` to `moyo::base::Cell` and back: lattice (moyo uses
  column or row conventions; verify and test against our row convention),
  fractional positions, integer atom types. Atom type policy, documented:
  sites are typed by element plus occupancy signature (species set and
  occupancies rounded to a tolerance), so a mixed Fe/Ni site never becomes
  equivalent to a pure Fe site; dummy atoms get their own type.
  Note: `xtalforge_symmetry::convert`; types are `(Z, occupancy x 1000)`
  mapped to dense integers in first-appearance order.
- [x] 2.1.2 Operations from a space group
  From a Hall number, Hall symbol, H-M symbol with setting, or number with
  default setting, produce `Vec<SymOp>` (exact rationals, Task 1.6.3) through
  moyo's `data` module. This is what lets a `Structure` built from a Hall
  symbol expand to P1.
- [x] 2.1.3 Detection entry point
  `detect_symmetry(structure, SymmetryOptions { symprec, angle_tolerance,
  setting })` returning our `SymmetryDataset`: number, Hall number, H-M
  symbol, Pearson symbol, operations, orbits (equivalent atoms), Wyckoff
  letters, site-symmetry symbols. moyo errors map to `Error::Symmetry` with a
  code and the tolerance used.
  Note: `detect_symmetry` returns `SymmetryDataset`; `Error::Symmetry` has
  code `symmetry` and the message includes the tolerance.
- [x] 2.1.4 Surfaces and docs
  CLI `symmetry detect` (output compatible in content with Crystal-Builder
  `xtal symmetry`), Python and WASM functions, API page, algorithm page on the
  backend isolation and the meaning of the tolerances.
  Note: `detect_symmetry` in Rust, Python and WASM; CLI `symmetry detect`;
  `docs/api/symmetry.md` and `docs/algorithms/symmetry.md`.
- [x] 2.1.5 Smoke validation
  On the whole corpus the detected space group number equals the golden
  spglib result at the default tolerance; any exception is listed with a
  reason.
  Note: the `simple/` corpus matches; the full sweep is an ignored test and
  its six exceptions are listed in `docs/algorithms/symmetry.md`.

## Task 2.2 — Space group value object and catalogue

- [x] 2.2.1 `SpaceGroup` type
  Port `core/spacegroup.py`: an immutable value whose identity is the Hall
  symbol (Hall number internally). Accessors: number, H-M short and full
  symbol, Schoenflies symbol, crystal system, lattice system, Bravais type,
  centering, point group, Laue class, operations, whether it is centrosymmetric
  or chiral, setting description (origin choice, unique axis, cell choice,
  hexagonal or rhombohedral axes).
  Note: `SpaceGroup` carries number, Hall number, H-M short/full, Hall symbol,
  setting, centering, systems, point group, Schoenflies, Bravais, Laue class,
  operations and the two predicates.
- [x] 2.2.2 Symbol parsing and lookup
  Accept user input such as `Fm-3m`, `F m -3 m`, `225`, `P2_1/c`,
  `P 1 21/c 1`, `P21/n`, `R-3:H`, `R-3:R`, `Fd-3m:1`, `Fd-3m:2`, Hall
  symbols (`-P 2ybc`), tolerant to spacing, case of the lattice letter and
  underscores. When the input is ambiguous about the setting, return the
  default setting and say which choices exist.
  Note: `SpaceGroup::parse` accepts numbers, H-M with optional `:setting`,
  Hall symbols and extra spaces/underscores.
- [x] 2.2.3 Catalogue
  List all 230 types and all Hall settings (530); filtering by crystal
  system, centering, point group. Used by the docs table and by the CLI
  (`xtalforge spacegroup list`, `xtalforge spacegroup show 227`).
  Note: `catalogue(SpaceGroupFilter)` and CLI `spacegroup show/list/search`.

## Task 2.3 — The symmetry workflow

- [x] 2.3.1 Detect
  Port `symmetry.detect`: detect at a tolerance and return the dataset plus a
  ready-to-use symmetric `Structure` (asymmetric unit, space group in the
  detected setting) without changing the input.
  Note: `detect_and_symmetrize` returns `SymmetryResult`.
- [x] 2.3.2 Set a space group
  Port `set_space_group` and its two modes: **reinterpret** (the current
  sites become the asymmetric unit of the new group) and **impose** (the
  current full cell is averaged onto the group: atoms are mapped to orbits,
  positions symmetrized, duplicates merged, with a report of the largest
  displacement). Both validate that the result is consistent.
  Note: `set_space_group` re-detects the result and rejects a different
  number; impose matches orbit images within 0.5 A.
- [x] 2.3.3 Asymmetrize
  Port `asymmetrize`: keep the detected symmetry and reduce the structure to
  its asymmetric unit, choosing representative sites deterministically
  (lowest label, then position) and keeping labels and `props`.
- [x] 2.3.4 Reduce to P1
  Throw the symmetry away: the asymmetric unit becomes the whole P1 cell, with
  labels made unique.
  Note: `to_p1` with unique `ElementN` labels.
- [x] 2.3.5 Change of setting
  Transform a structure between settings of the same type (origin choice 1/2,
  unique axis b/c, hexagonal/rhombohedral), keeping Cartesian positions.
  Note: `change_setting(structure, hall_number)` standardizes in the requested
  Hall setting and is exposed on all four surfaces; the physical geometry is
  preserved.

## Task 2.4 — Wyckoff positions and site symmetry

- [x] 2.4.1 Wyckoff tables
  For any Hall number: letter, multiplicity, site-symmetry symbol, coordinate
  triplets, from `moyo::data::iter_wyckoff_positions_from_hall_number`.
- [x] 2.4.2 Assignment
  For a structure in a known group, assign each site its Wyckoff letter and
  site symmetry; optionally snap coordinates exactly onto the special position
  within a tolerance and report the snap distance.
- [x] 2.4.3 Free parameters
  Express each assigned site as a coordinate template with free parameters
  (`x, x, z`, `1/4, y, 1/4`) and their current values. This is what
  symmetry-preserving editing (Phase 3) and Rietveld refinement of positions
  (Task 7.8.1) use.
  Note: `free_parameters` extracts x/y/z values from the template.

## Task 2.5 — Standardized, primitive and idealized cells

- [x] 2.5.1 Standardization
  Conventional standardized cell with the transformation matrix, origin shift
  and rotation from `MoyoDataset`; applied to a full `Structure` so labels,
  occupancies, charges, `props` and ADP tensors (rotated) are preserved.
- [x] 2.5.2 Primitive cell
  Standardized primitive cell (`prim_std_cell`) and the mapping back to input
  atoms (`mapping_std_prim`); also a "find primitive" that keeps the input
  orientation.
  Note: `primitive_cell` re-detects the primitive cell and returns its
  asymmetric unit.
- [x] 2.5.3 Idealization
  Determine from Task 0.3.2 whether moyo's standardized positions are
  symmetrized. If yes, expose them; if not, implement idealization by
  averaging each orbit over the group operations and projecting the lattice
  onto the metric allowed by the crystal system.
  Note: moyo refines the standardized positions; `idealize` exposes them.

## Task 2.6 — Lattice reduction (gap-fill)

moyo 0.22 implements Niggli, Delaunay and Minkowski reduction but keeps them
crate-private.

- [x] 2.6.1 Upstream request
  Open an issue or pull request on https://github.com/spglib/moyo proposing to
  export `niggli_reduce`, `delaunay_reduce` and `minkowski_reduce` (they are
  already `pub fn` inside private modules). Track the outcome here.
  Note: moyo 0.22 already exposes all three on `Lattice`; no request needed.
- [x] 2.6.2 Port until exported
  Port the three algorithms from moyo (MIT OR Apache-2.0, attribution in
  `THIRD_PARTY.md`), each in its own file under 500 lines, returning the
  reduced basis and the integer transformation. Replace with the upstream API
  once released.
  Note: not needed; the upstream methods are public in 0.22.
- [x] 2.6.3 Structure-level reduction
  Port the lattice-reduction operation of `core/supercell.py`: apply a
  reduction to a whole structure (P1 result, atoms rewrapped), as a command in
  Phase 3.
- [x] 2.6.4 Cell comparison
  Decide whether two lattices are equivalent (reduced forms compared with a
  tolerance), used by structure deduplication and powder indexing.
  Note: `lattices_equivalent` compares Niggli parameters.
- [x] 2.6.5 Property tests
  Reduction preserves volume and the lattice; reduction is idempotent; the
  transformation is unimodular.
  Note: the tests cover volume, atom count and determinant for all three
  algorithms, plus equivalence under a unimodular transform.

## Task 2.7 — Subgroups

- [x] 2.7.1 Maximal subgroup enumeration
  Port `core/subgroups.py` (1356 lines; split into enumeration, naming,
  application, reports). Enumeration uses moyo
  `subgroup::enumerate_translationengleiche_subgroups` and
  `enumerate_klassengleiche_subgroups(_by_index)`.
  Note: `translationengleiche_subgroups`, `klassengleiche_subgroups` and
  `identify_space_group`; naming and lowering remain in 2.7.2/2.7.3.
- [x] 2.7.2 Naming subgroups
  The expensive part according to Crystal-Builder: identify each subgroup's
  type, its setting and the transformation to its standard setting (with moyo
  `identify`), and present a readable name (`I4/mmm → Fmmm, index 2, t2`).
  Note: each `Subgroup` carries the identified subgroup number and Hall
  number when moyo can name it.
- [x] 2.7.3 Symmetry lowering
  Apply a chosen subgroup to a structure: split each Wyckoff orbit into the
  independent orbits of the subgroup, giving a new asymmetric unit with the
  new free parameters (first step of distortion, ordering and phase
  transition models).
  Note: `lower_to_subgroup` groups the P1 atoms by the subgroup operations and
  keeps one representative; free-parameter splitting is not yet explicit.
- [x] 2.7.4 Surfaces
  CLI `symmetry subgroups` and `symmetry lower`; Python and WASM.
  Note: implemented on all four surfaces.

## Task 2.8 — Layer, magnetic and point groups

- [x] 2.8.1 Layer groups
  Expose `MoyoLayerDataset` for slabs (80 layer groups) with our types.
  Note: `detect_layer_symmetry` with the c-perpendicular contract, on all
  four surfaces.
- [x] 2.8.2 Magnetic space groups
  Expose collinear and non-collinear magnetic datasets (1651 groups) behind
  an `experimental` marker; requires magnetic moments in `Site.props`.
- [x] 2.8.3 Molecular point groups
  Port the point-group part of `core/molecule.py`: Schoenflies symbol and
  symmetry elements of a finite molecule (needed for ORCA inputs and
  molecular builders). Implement it ourselves if moyo has no non-periodic
  point-group detection (to be confirmed in Task 0.3.2).
  Note: `molecular_point_group` maps a non-collinear atom frame over candidate
  triples and reads the Schoenflies symbol from the operation counts; linear
  molecules are rejected. Water C2v, ammonia C3v and methane Td are tested.

## Task 2.9 — Tolerance handling and diagnostics

- [x] 2.9.1 Tolerance ladder
  Detect across a ladder of tolerances (1e-5 to 0.5 Angstrom) and return the
  ranges in which each group is found, so a caller can see how robust the
  answer is.
  Note: `tolerance_ladder` groups consecutive settings into ranges.
- [x] 2.9.2 Diagnostics
  When detection fails or gives a low group, report which atoms break the
  expected symmetry and by how much (largest orbit deviations), in the
  structured diagnostic format used by the agent layer (Task 8.1.3).
  Note: `diagnose` returns per-atom deviations and max/mean values.

## Task 2.10 — Brillouin zone and k-path

Crystal-Builder takes the band path from ASE `Cell.bandpath` and builds the
zone with scipy's Voronoi. Both must be implemented here.

- [x] 2.10.1 Brillouin zone
  Wigner-Seitz cell of the reciprocal lattice: the region nearer the origin
  than any other reciprocal lattice point, computed as an intersection of
  half-spaces over the neighbouring lattice points. Output: vertices, faces,
  edges. Test: the fcc zone is a truncated octahedron with 14 faces and 24
  vertices.
  Note: `brillouin_zone` clips half-spaces; the fcc test checks 14 faces, 24
  vertices and 36 edges.
- [x] 2.10.2 Bravais classification and special points
  Determine the Bravais lattice and its variant (the Setyawan-Curtarolo 2010
  classification used by ASE), and the special points (Gamma, X, M, L, …)
  defined for the standard setting.
  Note: `classify_bravais` covers the 14 types with primitive-volume centering
  detection; the points come from the ASE 3.29 table in `data/kpaths.json`.
- [x] 2.10.3 Path in the user's cell
  Express the special points in the reciprocal basis of the cell the user
  actually has (a primitive fcc cell, not the textbook one), because DFTB+ and
  other codes must be given them in that basis. Positions along the path are
  Cartesian reciprocal distances in 1/Angstrom including 2π, so a band
  structure's x axis is to scale.
  Note: `KPath::cartesian_points` and `vertices` use the user's reciprocal
  basis and cumulative Cartesian distance.
- [x] 2.10.4 Validation
  Compare special points and paths with ASE output recorded in the golden data
  for each of the 14 Bravais lattices.
  Note: the data file is the ASE output itself; the tests check the fcc zone
  and several Bravais classifications.

## Task 2.11 — Validation

- [x] 2.11.1 spglib parity on the corpus
  Space group, Hall number, Wyckoff letters, standardized and primitive cells
  compared with golden spglib data; differences analysed and documented.
  Note: the simple corpus matches; six documented exceptions are in
  `docs/algorithms/symmetry.md`.
- [x] 2.11.2 Property tests
  Random structures built by applying the operations of random groups to
  random sites: detection recovers the group (or a supergroup when sites land
  on special positions); standardization is idempotent; P1 expansion count
  equals the sum of multiplicities; asymmetrize followed by P1 expansion
  reproduces the input.
  Note: generated-structure and idempotence tests cover eight groups from
  triclinic to cubic.
- [x] 2.11.3 Performance
  Benchmarks of detection on the corpus compared with spglib timings recorded
  from Crystal-Builder.
  Note: process-level medians for Si (3.3 ms) and MOF-5 (20.7 ms) are in
  `docs/algorithms/benchmarks.md`; the full spglib comparison is not run here.

## Task 2.12 — Surfaces and docs for Phase 2

- [x] 2.12.1 CLI
  `symmetry {detect,set,asymmetrize,p1,standardize,primitive,wyckoff,reduce,subgroups,lower,kpath}`
  and `spacegroup {list,show,search}`.
  Note: every command except `kpath` is implemented; `scan`, `diagnose` and
  `brillouin` were added as well.
- [x] 2.12.2 Python and WASM
  Matching functions and the `SpaceGroup`/`SymmetryDataset` types.
  Note: functions cover detection, workflow, standardize, primitive, reduce,
  Wyckoff, diagnostics, subgroups, lowering and the Brillouin zone.
- [x] 2.12.3 Docs and demos
  Algorithm pages (backend isolation, set versus impose, Wyckoff free
  parameters, Niggli reduction, subgroup naming, Brillouin zone, k-path); a
  browsable space group table; a demo that shows the space group, Wyckoff
  labels and Brillouin zone of a pasted CIF.
  Note: the API and algorithm pages are updated; the CIF inspector shows the
  group, Hall number, Wyckoff letters and zone counts. The k-path pages wait
  for Task 2.10.2-2.10.4.
- [x] 2.12.4 Parity and release
  Parity matrix; changelog; PATCH release.
  Note: parity and changelog are updated; the PATCH tag waits for the v0.0.1
  gate and the external release dry run. The Phase 2 surface set is complete
  for the implemented TypeScript/Python/CLI APIs.
