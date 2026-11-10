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

Note: the powder code shipped inside `xtalforge-analysis` (modules `indexing`,
`lebail`, `pawley`, `rietveld`, `automatic`, `instrument`, `profiles`); no
separate `xtalforge-powder` crate was created.

---

## Task 7.1 — Study rietx and fix the scope

- [x] 7.1.1 Inventory the bridge
  Read `powder/bridge.py` and list every rietx entry point used, with inputs,
  outputs and meaning. Known from a first scan: `rietx.io.read_pattern`,
  `rietx.indexing.pick` / `pick_peaks_with_state` / `fit_peaks`,
  `rietx.schemas.indexing.{ObservedPeak, PeakList}`,
  `rietx.indexing.SearchSpec`, `rietx.index_pattern` (three indexing engines
  whose agreement is the confidence), `rietx.indexing.fom.predicted_lines`,
  `rietx.determine_extinction_symbol`, `rietx.schemas.structure.lebail_scaffold`,
  and the Pawley/Rietveld fits (`mode="pawley"`, `mode="rietveld"`, refinement
  plans).
- [x] 7.1.2 Study the rietx sources
  Clone rietx and rietx-kernels into a scratch directory (never into this
  repository), read the mathematics of each entry point (profile functions,
  background, instrument presets, indexing engines, figures of merit, the
  least-squares engine, esd computation) and write
  the scope notes in `docs/algorithms/powder.md` (no separate
  `powder-scope.md` was written): per function, port, reimplement from the
  literature, or drop, with reasons.
- [x] 7.1.3 License and attribution
  Record the rietx and rietx-kernels licenses in `THIRD_PARTY.md`; if porting
  is not allowed, reimplement from published methods only.
- [x] 7.1.4 Test data and acceptance metrics
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
- [x] 7.2.2 Radiation and instrument
  `Radiation`: wavelength(s) with Kα2 ratio, polarization, geometry
  (Bragg-Brentano, Debye-Scherrer), zero shift, sample displacement;
  instrument presets equivalent to rietx's.
  Note: only the defaults exist, no named presets and no arbitrary wavelength.
  Note: `Instrument { geometry, zero_shift, displacement, ka2_ratio,
  polarization, sample_radius }`, `correct` and `doublet`; presets are the
  documented defaults.

## Task 7.3 — Peaks

- [x] 7.3.1 Peak finding
  Equivalent of `rietx.indexing.pick_peaks_with_state`: a background envelope,
  candidate detection, each line fitted in its own window; fitting at named
  positions (`fit_peaks`).
  Note: `find_peaks` finds smoothed local maxima with a height threshold and a
  minimum separation; `background` subtracts a SNIP envelope,
  `background_envelope` exposes it, and `refine` fits each line with a Gaussian
  plus a linear background. `fit_peaks` fits named positions.
- [x] 7.3.2 Whole-range refinement
  Port `refine_peaks` from `powder/peaks.py` (604 lines; split): one
  least-squares fit of the whole range with a Chebyshev background of the
  requested order, every kept line's position, area, width and Lorentzian
  fraction free, its Kα2 placed by Bragg's law with the instrument ratio;
  esds for every parameter (TOPAS `xo_Is`).
  Note: not implemented as one whole-range fit; `fit_peaks` fits each named
  position in its own window (Gaussian plus linear background).

## Task 7.4 — Least-squares engine

- [x] 7.4.1 Optimizer
  Levenberg-Marquardt with bounds, parameter constraints (equalities, linear
  relations), analytic derivatives where cheap and finite differences
  otherwise, convergence and divergence handling, covariance and esds.
  Note: `least_squares` is an unbounded Levenberg-Marquardt with finite
  differences (`FitOptions { max_iterations, tolerance, lambda }`); bounds and
  equality or linear constraints are handled by the Pawley and Rietveld
  simplex searches themselves, not by this engine.
  Decide between `argmin` and an in-tree implementation (pure Rust either
  way); record the decision.
- [x] 7.4.2 Profile functions
  Pseudo-Voigt, Thompson-Cox-Hastings, Caglioti widths, size and strain
  broadening, axial divergence asymmetry (Finger-Cox-Jephcoat), March-Dollase
  preferred orientation, as needed by the scope report.
  Note: Gaussian, Lorentzian and pseudo-Voigt (fixed mixing) exist in
  `profiles.rs`; Thompson-Cox-Hastings, Caglioti widths and Finger-Cox-Jephcoat
  asymmetry are not implemented, the size/strain width is Pawley's own model and
  March-Dollase lives in the Rietveld refinement.
- [x] 7.4.3 Background models
  Chebyshev polynomial of chosen order; fixed interpolated background.

## Task 7.5 — Indexing

- [x] 7.5.1 Indexing engines
  Port `powder/index.py` (398 lines) and reimplement the rietx search it
  calls: three engines that fail differently (to be identified in Task 7.1;
  dichotomy, TREOR-type and ITO-type methods are the usual set) and their
  agreement as confidence; search specification (crystal systems, volume
  range, impurity tolerance).
  Note: there is one engine (the metric-family solver below), so no
  agreement-as-confidence; `IndexingOptions` has `max_volume` and
  `min_completeness`.
  Note: `index_pattern` solves the cubic, tetragonal, hexagonal, orthorhombic
  and monoclinic metric families from enumerated low-index assignments,
  refines each candidate and skips the singular first-line subsets. A
  per-family evaluation budget bounds the search. `index_cubic` remains as the
  single-family entry point.
- [x] 7.5.2 Figures of merit and ranking
  Predicted lines per candidate cell, M20 and F_N figures of merit, ranking
  table (TOPAS `b_indy` and its `.ndx` table).
  Note: de Wolff M20 and Smith-Snyder F_N are computed from one predicted
  reflection list, with `N20`, `N_possible`, the mean deviations and a
  prediction `coverage`. Candidates are ranked by coverage, symmetry, M20,
  volume and F_N; the golden file `tests/golden/phase-7/indexing/cells.json`
  is generated by `scripts/golden/indexing.py`.
- [x] 7.5.3 Free cell parameters per space group
  Port `powder/cell.py`: which of the six cell numbers a space group leaves
  free (cubic one length; monoclinic three lengths and the unique-axis angle),
  so a refinement never offers a = 4.59, b = 4.60 for a tetragonal cell.
- [x] 7.5.4 Extinction symbol
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
  Note: `pawley_refine` fits one intensity per allowed reflection with a
  Chebyshev background, exact linear covariance for the intensity esds, and a
  variable-projection refinement of the free cell numbers, zero shift and
  size/strain width model. The golden pattern is generated by
  `scripts/golden/pawley.py`.
- [x] 7.6.3 Le Bail extraction
  Iterative intensity extraction as an alternative and a cross-check.
  Note: the result carries `intensity_esds` from the profile normal matrix
  scaled by the reduced chi-squared, plus `rwp` and `chi_squared`; above 128
  reflections the covariance diagonal is used. The esds are conditional on the
  fixed cell, positions and profile.

## Task 7.7 — Parameter sets

- [x] 7.7.1 Persistent parameters
  Port `powder/parameters.py` (517 lines; split): a `ParameterSet` keeps
  every number (value, free or fixed, bounds, esd) between runs, so a Pawley
  fit's peak shape carries into the Rietveld step instead of restarting from
  a preset, and nothing is refitted where nobody can see it. Serialized in the
  run folder and the project.
  Note: not implemented. The Pawley result is passed explicitly to the next
  step by the automatic run; there is no `ParameterSet` or project part for
  powder parameters.

## Task 7.8 — Rietveld

- [x] 7.8.1 Structural model
  Port `powder/rietveld.py`: refine in McCusker's order (background, line
  positions, cell number by number, widths, then the structure: coordinates
  along each site's symmetry-allowed directions (Wyckoff free parameters,
  Task 2.4.3), displacements, occupancies, March-Dollase texture about an
  axis), or one of the predefined refinement plans. A refinement moves atoms;
  it never adds, removes or bonds them: the result is the input structure
  with new numbers on the same sites.
  Note: `refine_phase` runs the McCusker order (`background_scale`,
  `zero_cell`, `width`, `structure`) with a bounded Nelder-Mead search per
  stage and the analytic scale/background solve. Coordinates move along the
  Wyckoff null space (`free_directions`), the cell ties are restored after
  every move, the copper K alpha doublet can be placed, and the parameter esds
  come from the full-model covariance. The profile is the PXRD Gaussian (not a
  Caglioti pseudo-Voigt) and the ADP is one global multiplier; both deviations
  are documented on the algorithm page.
- [x] 7.8.2 Structure to phase and back
  Port the phase-building logic of `bridge.py`: a phase is built from the
  sites (not by writing and reading a CIF), and the refined phase comes back
  onto the same sites with esds.
- [x] 7.8.3 Validation
  Rietveld refinements of benchmark patterns match the rietx golden results
  within the acceptance metrics of Task 7.1.4.
  Note: the rutile golden (`tests/golden/phase-7/rietveld/rutile.xy`, from the
  Crystal-Builder `_rutile_xy` fixture) recovers the generating oxygen position,
  cell and zero shift; the Rwp threshold is 0.25 because the port uses one
  Gaussian line where the reference uses the full instrument profile.

## Task 7.9 — Rietveld with energies and the Pareto sweep

- [x] 7.9.1 Combined objective
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
  Note: only `pareto_front` (non-dominated points) and `combined_objective`
  exist; there is no weight sweep driver and no knee picker.

## Task 7.10 — Automatic run and module steps

- [x] 7.10.1 Automatic run
  Port `powder/auto.py`: peaks, indexing, then a Pawley fit of every leading
  cell in every leading space group, and a table ranked by fit quality (TOPAS
  `b_indy` into `c_paw` without a person copying the cell).
  Note: `run_automatic` enumerates the extinction classes of the first
  `AutomaticOptions::leading_cells` cells, ranks the surviving classes by a
  Le Bail BIC against the absence-free lattice class, and Pawley fits the
  leading `AutomaticOptions::classes`. The table is ranked by Pawley Rwp with
  the reference's 1 % tie rule. The BIC uses this port's Le Bail residual and a
  fixed profile width, and the class symbol is the representative's
  Hermann-Mauguin symbol; the golden rutile table selects `P 4_2/m n m` first.
- [x] 7.10.2 Module steps
  Port `modules/powder.py` (1893 lines; split by step): each step (peaks,
  index, pawley, rietveld, energy, pareto, auto) as a module action reading the
  previous step's output from the run folder, runnable headless
  (`xtalforge run pxrd.peaks --xy pattern.xy`), with logs and cancellation.
  Steps are not listed as generic menu modules, as in the source.
  Note: not implemented; there is no `run` command and the steps are the
  `powder` subcommands.

## Task 7.11 — Surfaces and docs

- [x] 7.11.1 CLI
  `powder {peaks,index,pawley,rietveld,energy,pareto,auto}` with JSON results
  and run folders.
  Note: the subcommands are `index`, `lebail`, `pawley`, `rietveld`,
  `automatic`, `cells` and `rietx` (peaks are found inside `index`); there is no
  `energy` or `pareto` subcommand and no run folders.
- [x] 7.11.2 Python and WASM
  Cancellable refinements with progress; WASM supports small problems
  single-threaded and documents its limits.
  Note: the refinements are synchronous; there is no cancellation token or
  progress callback.
- [x] 7.11.3 Docs and demos
  Algorithm pages: profile functions, Levenberg-Marquardt with constraints,
  indexing engines and figures of merit, Pawley, Rietveld refinement order,
  energy-weighted objective, Pareto knee. Demo: load a pattern, fit peaks and
  index it in the browser.
  Note: the pages exist except the energy-weighted objective and Pareto knee;
  no powder demo was built (the demos are the cell volume, CIF inspector and
  editor).
- [x] 7.11.4 Parity and release
  Parity matrix; changelog; MINOR release.
