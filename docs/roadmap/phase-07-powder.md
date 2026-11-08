# Phase 7 — Powder pattern refinement

**Goal:** the refinement workbench of Crystal-Builder without its GUI: a
measured pattern and its radiation, peak finding and whole-range peak
refinement, indexing with extinction-symbol determination, Pawley, Rietveld,
Rietveld with force-field energies (as Materials Studio Reflex does), a Pareto
sweep over the energy weight, persistent parameter sets, and the automatic run
(peaks to ranked Pawley fits). The workflow follows TOPAS terminology
(`xo_Is`, `b_indy`, `c_paw`, `d_riet`), as Crystal-Builder does.

**Reference:** `Crystal-Builder/xtal/powder/` (data 233, peaks 604, index
398, cell 97, pawley 316, rietveld 247, energy 572, pareto 210, parameters
517, auto 342, bridge 1268), `xtal/modules/powder.py` (1893). Crystal-Builder
delegates the numerical core to **rietx** (https://github.com/yue-here/rietx,
PyPI 1.11.0; depends on gemmi, numpy, scipy, spglib, pydantic and
`rietx-kernels`); `bridge.py` is the only module importing it. The comment in
Crystal-Builder's `pyproject.toml` describes rietx as MIT; PyPI metadata gives
no license field, so verify it in the rietx repository.

Crate: `xtalforge-powder`.

---

## Task 7.1 — Study rietx and fix the scope

- [ ] 7.1.1 Inventory the bridge
  Read `powder/bridge.py` and list every rietx entry point used, with inputs,
  outputs and meaning. Known from a first scan: `rietx.io.read_pattern`,
  `rietx.indexing.pick` / `pick_peaks_with_state` / `fit_peaks`,
  `rietx.schemas.indexing.{ObservedPeak, PeakList}`,
  `rietx.indexing.SearchSpec`, `rietx.index_pattern` (three indexing engines
  whose agreement is the confidence), `rietx.indexing.fom.predicted_lines`,
  `rietx.determine_extinction_symbol`, `rietx.schemas.structure.lebail_scaffold`,
  and the Pawley/Rietveld fits (`mode="pawley"`, `mode="rietveld"`, refinement
  plans).
- [ ] 7.1.2 Study the rietx sources
  Clone rietx and rietx-kernels into a scratch directory (never into this
  repository), read the mathematics of each entry point (profile functions,
  background, instrument presets, indexing engines, figures of merit, the
  least-squares engine, esd computation) and write
  `docs/algorithms/powder-scope.md`: per function, port, reimplement from the
  literature, or drop, with reasons.
- [ ] 7.1.3 License and attribution
  Record the rietx and rietx-kernels licenses in `THIRD_PARTY.md`; if porting
  is not allowed, reimplement from published methods only.
- [ ] 7.1.4 Test data and acceptance metrics
  Collect the measured patterns used by Crystal-Builder's powder tests and
  define acceptance: peak positions within a fraction of the step, indexing
  recovers the published cell, Pawley/Rietveld Rwp within a stated margin of
  the rietx results recorded as golden data, refined parameters within their
  esds.

## Task 7.2 — Pattern and radiation

- [x] 7.2.1 Measured pattern
  Port `powder/data.py`: three arrays (2θ, intensity, esd or √I), read from
  `.xy` (Task 3.5.3) and the other formats rietx reads that Crystal-Builder
  tests use; range selection, excluded regions.
- [ ] 7.2.2 Radiation and instrument
  `Radiation`: wavelength(s) with Kα2 ratio, polarization, geometry
  (Bragg-Brentano, Debye-Scherrer), zero shift, sample displacement;
  instrument presets equivalent to rietx's.

## Task 7.3 — Peaks

- [x] 7.3.1 Peak finding
  Equivalent of `rietx.indexing.pick_peaks_with_state`: a background envelope,
  candidate detection, each line fitted in its own window; fitting at named
  positions (`fit_peaks`).
- [ ] 7.3.2 Whole-range refinement
  Port `refine_peaks` from `powder/peaks.py` (604 lines; split): one
  least-squares fit of the whole range with a Chebyshev background of the
  requested order, every kept line's position, area, width and Lorentzian
  fraction free, its Kα2 placed by Bragg's law with the instrument ratio;
  esds for every parameter (TOPAS `xo_Is`).

## Task 7.4 — Least-squares engine

- [x] 7.4.1 Optimizer
  Levenberg-Marquardt with bounds, parameter constraints (equalities, linear
  relations), analytic derivatives where cheap and finite differences
  otherwise, convergence and divergence handling, covariance and esds.
  Decide between `argmin` and an in-tree implementation (pure Rust either
  way); record the decision.
- [x] 7.4.2 Profile functions
  Pseudo-Voigt, Thompson-Cox-Hastings, Caglioti widths, size and strain
  broadening, axial divergence asymmetry (Finger-Cox-Jephcoat), March-Dollase
  preferred orientation, as needed by the scope report.
- [x] 7.4.3 Background models
  Chebyshev polynomial of chosen order; fixed interpolated background.

## Task 7.5 — Indexing

- [x] 7.5.1 Indexing engines
  Port `powder/index.py` (398 lines) and reimplement the rietx search it
  calls: three engines that fail differently (to be identified in Task 7.1;
  dichotomy, TREOR-type and ITO-type methods are the usual set) and their
  agreement as confidence; search specification (crystal systems, volume
  range, impurity tolerance).
- [x] 7.5.2 Figures of merit and ranking
  Predicted lines per candidate cell, M20 and F_N figures of merit, ranking
  table (TOPAS `b_indy` and its `.ndx` table).
- [ ] 7.5.3 Free cell parameters per space group
  Port `powder/cell.py`: which of the six cell numbers a space group leaves
  free (cubic one length; monoclinic three lengths and the unique-axis angle),
  so a refinement never offers a = 4.59, b = 4.60 for a tetragonal cell.
- [ ] 7.5.4 Extinction symbol
  Equivalent of `rietx.determine_extinction_symbol`: candidate extinction
  symbols and space groups from the indexed reflections, using Phase 2 data.

## Task 7.6 — Pawley and Le Bail

- [x] 7.6.1 Le Bail scaffold
  Equivalent of `rietx.schemas.structure.lebail_scaffold`: a phase with a cell
  and space group and no atoms, its reflection list generated from Phase 2
  (systematic absences from the operations).
- [x] 7.6.2 Pawley fit
  Port `powder/pawley.py`: one cell and space group fitted to the whole
  pattern with every reflection intensity free (TOPAS `c_paw`), the test an
  indexed cell must pass. User settings: range, background order, which of
  zero shift, displacement, cell, size and strain broadening are free.
  Results: R-factors (Rp, Rwp, Rexp, χ²), refined cell with esds, reflection
  list. Document, as the source does, that a Bragg-Brentano cell also carries
  a systematic error of order 1e-4 no esd includes.
- [x] 7.6.3 Le Bail extraction
  Iterative intensity extraction as an alternative and a cross-check.

## Task 7.7 — Parameter sets

- [ ] 7.7.1 Persistent parameters
  Port `powder/parameters.py` (517 lines; split): a `ParameterSet` keeps
  every number (value, free or fixed, bounds, esd) between runs, so a Pawley
  fit's peak shape carries into the Rietveld step instead of restarting from
  a preset, and nothing is refitted where nobody can see it. Serialized in the
  run folder and the project.

## Task 7.8 — Rietveld

- [ ] 7.8.1 Structural model
  Port `powder/rietveld.py`: refine in McCusker's order (background, line
  positions, cell number by number, widths, then the structure: coordinates
  along each site's symmetry-allowed directions (Wyckoff free parameters,
  Task 2.4.3), displacements, occupancies, March-Dollase texture about an
  axis), or one of the predefined refinement plans. A refinement moves atoms;
  it never adds, removes or bonds them: the result is the input structure
  with new numbers on the same sites.
- [ ] 7.8.2 Structure to phase and back
  Port the phase-building logic of `bridge.py`: a phase is built from the
  sites (not by writing and reading a CIF), and the refined phase comes back
  onto the same sites with esds.
- [ ] 7.8.3 Validation
  Rietveld refinements of benchmark patterns match the rietx golden results
  within the acceptance metrics of Task 7.1.4.

## Task 7.9 — Rietveld with energies and the Pareto sweep

- [ ] 7.9.1 Combined objective
  Port `powder/energy.py` (572 lines; split): the objective
  `f = (1 - w) χ²/χ²₀ + w (E - E₀)/ΔE`, where χ²₀ and E₀ are the values at
  the start and ΔE is how far the energy falls when the atoms answer to the
  force field alone (one relaxation at the same cell, the w = 1 end), so both
  terms are of order one and a weight means the same thing on rutile and on a
  thousand-atom framework. Any Phase 4 calculator.
- [x] 7.9.2 Pareto sweep
  Port `powder/pareto.py`: refine at a series of weights from pattern-only to
  energy-only, keep the points no other point beats on both counts, and pick
  the knee where giving up fit starts to buy little energy.

## Task 7.10 — Automatic run and module steps

- [ ] 7.10.1 Automatic run
  Port `powder/auto.py`: peaks, indexing, then a Pawley fit of every leading
  cell in every leading space group, and a table ranked by fit quality (TOPAS
  `b_indy` into `c_paw` without a person copying the cell).
- [ ] 7.10.2 Module steps
  Port `modules/powder.py` (1893 lines; split by step): each step (peaks,
  index, pawley, rietveld, energy, pareto, auto) as a module action reading the
  previous step's output from the run folder, runnable headless
  (`xtalforge run pxrd.peaks --xy pattern.xy`), with logs and cancellation.
  Steps are not listed as generic menu modules, as in the source.

## Task 7.11 — Surfaces and docs

- [ ] 7.11.1 CLI
  `powder {peaks,index,pawley,rietveld,energy,pareto,auto}` with JSON results
  and run folders.
- [ ] 7.11.2 Python and WASM
  Cancellable refinements with progress; WASM supports small problems
  single-threaded and documents its limits.
- [ ] 7.11.3 Docs and demos
  Algorithm pages: profile functions, Levenberg-Marquardt with constraints,
  indexing engines and figures of merit, Pawley, Rietveld refinement order,
  energy-weighted objective, Pareto knee. Demo: load a pattern, fit peaks and
  index it in the browser.
- [ ] 7.11.4 Parity and release
  Parity matrix; changelog; MINOR release.
