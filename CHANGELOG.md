# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/) (pre-1.0: `0.0.x` until the first
essential release, then `0.x.y`).

## Unreleased

### Added

- Phase 8 Task 8.3: open extension registries and a Python plugin surface.
  `FormatRegistry`, `ScalarRegistry` and `StyleRegistry` join the existing
  calculator and module registries; `examples/plugin` is a buildable crate
  that registers a format, a calculator, a module, a scalar and a style.
  Python `Plugins` registers callables as modules, calculators and formats,
  with JSON in/out and typed errors. The draw-style enum stays closed; the
  open part is the settings preset, recorded in `docs/api/plugins.md`.
- Task 0.9: bundle an 18-complex subset of tmQM 2024 (MIT, one complex per
  metal from Sc to Ir) and validate the XYZ reader against the header
  stoichiometry, the perceived metal coordination and an XYZ round-trip. UFF
  builds, types and evaluates every metal with finite forces.
- Task 0.9: bundle an 18-structure subset of the CoRE MOF 2019 public archive
  (CC BY 4.0, Zenodo 14184621) and validate against it. Every file parses,
  expands to P1 and bonds; four frameworks are checked against their own
  `_chemical_formula_sum`, every declared `_cell_volume` is reproduced within
  0.05 A^3, CIF round-trips preserve sites and volume, and HKUST-1 and ZIF-8
  give bounded PXRD and porosity numbers.
- Phase 8 Task 8.2: `xtalforge mcp` serves the verb catalogue as MCP tools
  over stdio (JSON-RPC 2.0, protocol `2024-11-05`, no external MCP crate),
  with generated input schemas and an `instructions` payload for chat
  clients. Fourteen skill documents ship with the CLI (`xtalforge skill
  path|list|show|install`) and are published as the automation guide. Calls
  run under the standard size limits; cancellation is process-level and the
  window proxy is out of scope, both recorded in
  `docs/algorithms/agent.md` section 11.
- Phase 8 Task 8.1: the `xtalforge-agent` crate: one verb catalogue shared by
  the Rust, CLI, Python and WASM surfaces (22 verbs), a stateful session whose
  mutations use `Origin::Agent`, inspection and diagnostics (overlap,
  occupancy, charge, dangling and under-coordinated atoms, unusual
  coordination, symmetry), command-log change bursts, typed JSON answers and a
  pure-Rust orthographic renderer to SVG and PNG. CLI `capabilities`, `render`
  and `agent`; Python `AgentSession`, `agent_run`, `agent_capabilities`,
  `agent_verbs`; WASM `AgentSession`, `agent_run`, `agent_capabilities`,
  `agent_verbs`. The WASM size budget moves from 3648 to 3840 KiB to carry the
  new exported surface (3775 KiB).
- Phase 0 Task 0.9: external validation suite. EQeq is compared with the
  authors' program on the COD MOF-5 cell (424 charges within 1e-3 e), UFF with
  sci-form on methane and stretched water (within 1e-4 kcal/mol), and the
  DFTB+ parser with a recorded Si SCC output. `docs/algorithms/external-validation.md`
  records the bundled sources, their licences and the data sets still to add.
- Phase 5 Task 5.6: bundle the full RCSR 2019-06-01 snapshot in `data/rcsr/`
  and identify nets by coordination sequence and point symbol against all
  2,931 entries. `RcsrIndex::from_embedded()` returns `RcsrIdentification`
  with the unique name, an ambiguous set, or the nearest candidates; the
  2,909 keyed entries are validated against their canonical keys at index
  build time. The CGD reader now accepts both RCSR edge dialects.
- Phase 6 Task 6.8: the disordered-carbon recipe builder: convolution surface
  around a periodic net, split/collapse/flip remesh with tangential relaxation
  and projection, ribbon paths over the triangle graph, a three-coordinate
  dual and H/F/O edge termination. `CarbonBuild` reports density, ratios,
  ring census, mean C-C, pieces and periodicity; CLI `build carbon`, Python
  `build_carbon`, WASM `build_carbon_recipe_info`.
- Phase 7 Task 7.10.1: every leading-cell row now carries `space_groups`, the
  candidates whose general-position systematic absences every indexed
  reflection obeys. The screen enumerates the cell's crystal system, evaluates
  the structure factor of two generic points per reflection and ranks by
  group order; CLI, Python and WASM expose the candidates.
- Phase 7 Task 7.10.1: the automatic table now ranks extinction classes by a
  Le Bail BIC and Pawley-fits every leading (cell, class), as the reference
  does. A row carries `cell_rank`, `class_rank`, `class_symbol`, `delta_bic`
  and the groups that share the class; the 1 % Rwp tie rule orders the table.
- Phase 4 Task 4.9: validate the DFTB+ energy, force, stress and population
  parser against the recorded `si_scc_detailed.out` fixture. The xTB GFN
  methods remain available in-process through `sci-form` (Task 4.10).
- Phase 1 Task 1.12: `rings_explicit` and `ring_counts_explicit` census a
  builder's stored bonds without distance perception.

### Added

- Phase 4 Task 4.10: link the published MIT crate `sci-form` 2.0.0 behind the
  optional `sci-form` feature and implement `Calculator` for GFN0, GFN1 and
  GFN2-xTB. The methods are molecular only (a space-group structure is
  refused), GFN2 returns central-difference forces, and GFN0/GFN1 refuse
  `energy_forces`. The registry lists `gfn0`, `gfn1` and `gfn2`.

### Fixed

- Phase 5 Task 5.4: order the four crossings of a marching-tetrahedra 2-2 split around their quadrilateral. The fan taken in the generation order crossed itself, leaving non-manifold edges and a surface that was not watertight.
- Phase 7 Task 7.6.3: Le Bail extraction returns intensity esds from the
  profile normal matrix, plus Rwp and the reduced chi-squared; the CLI report
  includes them.
- Phase 7 Task 7.10.1: the automatic run fits a Pawley model of every leading
  indexing candidate and returns a table sorted by Pawley Rwp, with the cell
  esds; CLI, Python and WASM expose the table.
- Phase 7 Task 7.8: run the Rietveld refinement in McCusker order
  (analytic scale/background, zero shift and free cell numbers, width, then
  Wyckoff coordinates, isotropic ADPs, occupancies and March-Dollase texture)
  with a bounded direct search per stage, covariance-based parameter esds and
  the `stages` list. Coordinates move only along the Wyckoff null space
  (`free_directions`), the metric-family cell ties are restored after every
  move, and the copper K alpha doublet can be placed. A rutile golden pattern
  from Crystal-Builder validates the recovered position, cell and zero shift.
- Phase 7 Task 7.6.2: add the Pawley fit with one intensity per reflection,
  exact linear covariance and intensity esds, free cell numbers, zero shift and
  size/strain width model, R-factors and chi-squared; a Python-generated golden
  pattern validates the recovery.
- Phase 7 Task 7.10: the automatic run indexes the five metric families and
  adds a Pawley step of the indexed cell; Python `powder_pawley` and WASM
  `powder_pawley_info` expose the fit.
- Phase 7 Task 7.5.1/7.5.2: add general metric-family indexing (cubic,
  tetragonal, hexagonal, orthorhombic, monoclinic) with cell refinement, de
  Wolff M20, Smith-Snyder F_N and a coverage-ranked table; Python
  `powder_index`, WASM `powder_index_info` and `powder index --general`.
- Phase 7 Task 7.3.1: add the SNIP background envelope, per-line Gaussian plus
  linear background fitting and `fit_peaks`, with `--background` and
  `--refine-peaks` on the CLI.
- Phase 1 Task 1.15.4: map a non-P1 `_geom_bond` row from its site labels,
  operation code and image to P1 contacts (MIL53 matches the 36 golden
  contacts), collapsing symmetry-equivalent images by their minimum-image
  vector. Unmatched records annotate nothing instead of inventing an edge.
- Phase 1 Task 1.15.4: warn when a CIF operation loop's count disagrees with
  the named space group.
- Phase 3 Task 3.6.3: add LAMMPS text dump (element column or type map,
  orthogonal/triclinic, Cartesian/scaled coordinates) and DFTB+ `geo_end`
  readers to the trajectory dispatcher.
- Phase 1 Task 1.15.4: write non-P1 explicit bonds back as a validated
  P1-index CIF loop and compare the full symmetry-operation set with the named
  space group, not only its count.
- Phase 4 Task 4.7: choose overlapping disorder components from the declared
  CIF formula and cell charge, make the planner graph-aware, add hydrogen
  filling, the close-contact list, the prepared-structure undoable command
  path and their tests.
- Phase 4 Task 4.9.2/4.9.3: wire `DftbCalculator` and `XtbCalculator` to the
  runner, parse `detailed.out`/xTB stdout, validate force counts and support
  timeout, cancellation and scratch inspection.
- Phase 4 Task 4.9.1: add the scratch-folder process runner with wall-clock
  timeout, cancellation, capped stdout/stderr logs, relative input files and
  simulated-executable tests.
- Phase 4 Task 4.7: add duplicate merging, disorder-component choice, solvent
  removal, charge-balance reporting and the prepared-sample golden sanity
  pass; `PrepareReport` now returns the prepared structure.
- Phase 4 Task 4.8.3: add supercell force constants, shortest-image folding,
  the acoustic sum rule, a Monkhorst-Pack q-mesh, Gaussian phonon DOS and a
  phonopy spring-model golden comparison. An explicit `q_path` samples a band
  path.

## [0.0.3] - 2026-11-09

### Changed

- Phase 1 Task 1.11.1/1.11.5: perceive bonds with the reference defaults
  (`scale = 1.15`, metal-metal pairs forbidden) and validate the full golden
  corpus: every available bond set and coordination number matches
  Crystal-Builder after P1 atom mapping.
- Phase 1 Task 1.6/1.15: expand each site on its own duplicate grid, so
  disorder components at the same position both survive, and resolve a CIF
  Hall/H-M symbol into operations when the file omits the operation loop.
- Phase 1 Task 1.12: replace the chordless-cycle scan with the periodic
  Franzblau walk, expose `Ring::path`/`ring_paths` with lattice images on all
  four surfaces, and validate counts and paths against the golden corpus.
- Phase 1 Task 1.15.4: round-trip Uij (matched by label), disorder
  group/assembly, the oxidation-number charge fallback and the explicit-bond
  loop; split the CIF writer and structure extraction under `cif/`.
- Phase 3 Task 3.6.2/3.6.3: add lazy `XyzIndex` frame access and the VASP
  XDATCAR reader with XYZ/XDATCAR detection in the CLI.

## [0.0.2] - 2026-11-09

### Changed

- Phase 1 Task 1.10.1/1.10.2: replace the direct periodic pair scan with a
  linked-cell index in fractional space and exact image enumeration, add the
  10^3/10^4/10^5-atom benchmark target and compare it with Crystal-Builder.
- Phase 1 Task 1.10.2: P1 expansion now deduplicates through a fractional
  grid, removing the quadratic scan that dominated large structures.
- Phase 1 Task 1.7.1: supercell boundary duplicates use the same fractional
  grid, cutting the 2x2x2 MOF-5 build from 37.5 ms to 5.5 ms.

### Fixed

- Rename the Linux desktop entry to the reverse-DNS component ID so
  `appimagetool` finds the AppStream metadata; the AppImage build no longer
  reports the missing-metadata warning, and the `.deb`/`.rpm` payloads carry
  the matching desktop ID.

## [0.0.1] - 2026-11-08

### Changed

- Refresh the documentation against the implemented Phase 1-7 surface: the
  README status, guide, API pages, algorithm pages, validation table, benchmark
  notes, third-party attribution and roadmap current-status/gap list.
- Record the remaining documented subsets (neighbour acceleration, phonon
  q-mesh, preparation steps, general indexing and figures of merit, powder
  profile simplifications) in `docs/roadmap/index.md`.

### Fixed

- Upgrade `zip` to 9.4.2, clearing RUSTSEC-2026-0367, RUSTSEC-2026-0368,
  RUSTSEC-2026-0369 and RUSTSEC-2026-0429.
- Fix the rustdoc intra-doc links that failed `cargo doc -D warnings`.
- Expose C20, twists, ribbons, ZTC/Schwarz P surfaces and polymers in the CLI,
  Python and WASM, closing the Phase 6 surface-parity gaps; add regression
  tests for each surface.

### Added

- Phase 6 Task 6.1/6.5/6.7/6.8/6.9: SMILES 3D embedding, metal coordination,
  the sketch model and SMILES writer, commensurate twists, C20/C60 face
  spirals, gyroid/Schwarz P carbon surfaces and polymer chain packing.
- Phase 7 Task 7.1-7.11: the rietx inventory and acceptance metrics, the
  instrument model, whole-range and structural least-squares refinement,
  free cell parameters and extinction symbols, the combined objective, the
  automatic powder run, and the CLI/Python/WASM surfaces.
- Phase 6 Task 6.11-6.16: the `xtalforge-mof` crate with six topologies,
  eight blocks, framework assembly, rods, layers, merged nets, decoration,
  invented linkers and pore-space partition.
- Phase 6 Task 6.17: `mof` CLI commands, Python/WASM catalogues and builder,
  and the MOF API page.
- Phase 4 Task 4.3.3: native EEQ-style charges.
- Phase 4 Task 4.6: hydrogen placement.
- Phase 4 Task 4.7: the preparation pipeline and report.
- Phase 4 Task 4.9: engine probes plus DFTB+, xTB and ORCA input writers and
  output parsers.
- Phase 4 Task 4.10: the sci-form adapter contract and capability limits.
- Phase 4 Task 4.12: `ff` CLI commands, Python/WASM force-field surfaces and
  the algorithm pages.
- Phase 3 Task 3.7: the `.xtalproj` ZIP container with independent
  `format_version`, readable parts and a newer-version gate.
- Phase 3 Task 3.8.2: `symmetry`, `pxrd` and `uff` parameter presets shared by
  the CLI, Python and WASM schemas.
- Phase 3 Task 3.10.7: `SetCharges`, `Interpenetrate` and
  `MarkConnectionPoints` commands.
- Phase 3 Task 3.12.3: in-memory workspace archives for WASM and on-disk
  workspace binding for Python.
- Phase 3 Task 3.16: `formats`, `params`, `project`, `trajectory` and `script`
  CLI commands; Python `Document`/project/workspace bindings; WASM
  `WasmDocument` plus project and workspace exports; the editor demo.
- Phase 4 Task 4.8: Gamma-point Hessian, vibrational modes, IR spectrum,
  Gamma-point phonons, harmonic/ideal-gas thermochemistry and the
  quasi-harmonic free-energy scan.
- Phase 4 Task 4.3-4.5: QEq/EQeq charges, per-site charges, Ewald energy and
  forces, FIRE/L-BFGS optimization (P1/asymmetric, optional cell), equations of
  state and relaxed scans.
- Phase 4 Task 4.1/4.2: Calculator API, registry, run records and the native
  UFF/UFF4MOF engine with analytic gradients.
- Phase 5 Task 5.1-5.3: calculated PXRD with radiation, reflections,
  structure factors, LP and profiles; periodic distance grids and porosity.
- Phase 6 Task 6.5/6.6: graphene sheets/stacks and nanotube roll-up builders.
- Phase 3: format registry with capabilities/loss warnings, atomic writes,
  pymatgen JSON/CSSR/GEN/LAMMPS/SHELX/PDB, CGD nets, Grid3/cube/XY,
  trajectories, typed parameters, commands/document with undo, selection and
  tracking, workspace/modules, view data, colormaps and ASCII STL.
- Phase 3 Task 3.1-3.3: the format registry with capabilities and loss
  warnings, atomic writes, and pymatgen JSON, CSSR and DFTB+ GEN readers and
  writers.

- Phase 2 Task 2.8: layer-group detection for slabs, experimental collinear
  magnetic groups and molecular point groups.
- Phase 2 Task 2.10.2-2.10.4: Bravais classification, ASE-derived special
  points and band paths in the user's reciprocal basis.
- Phase 2 Task 2.3.5: change of setting by Hall number.

- Project bootstrap: `AGENTS.md` rules and versioning policy, roadmap (phases 0
  to 8, verified against the Crystal-Builder, sci-form and moyo sources),
  rule-check script, and VitePress documentation scaffold.
- Cargo workspace (Task 0.1.1) with `xtalforge-core`, `xtalforge-symmetry`,
  `xtalforge-io`, the facade `xtalforge`, `xtalforge-cli`, `xtalforge-python`
  and `xtalforge-wasm`; shared versions for `nalgebra`, `serde`, `thiserror`,
  `rayon`, `moyo` and `flate2`.
- Toolchain, formatting and lint configuration (Task 0.1.2):
  `rust-toolchain.toml`, `rustfmt.toml` and the workspace lints for unsafe
  code, missing docs, warnings and panic-prone patterns.
- Repository hygiene (Task 0.1.3): `.editorconfig` and `deny.toml` with the
  license allow list, advisory checks, duplicate-version ban and native-build
  bans.
- Facade crate (Task 0.1.4) with stable `core`, `symmetry` and `io` module
  names, a prelude, and the `parallel`, `sci-form`, `ml`, `external` and
  `mof-data` features documented in the crate and in the guide.
- Contributor guide (Task 0.1.5) and the Cargo features guide.
- Rule-check self-test (Task 0.2.1): positive and negative fixtures and a
  runner covering the 500-line boundary and the ban on warning-suppression
  attributes.
- Optional Git hook (Task 0.2.2): `scripts/install-hooks.sh` and a pre-commit
  hook running the rules, formatting and clippy on touched crates.
- CI workflow for the repository rules and their self-test (Task 0.2.4),
  plus blocking formatting and clippy jobs with warnings denied.
- Data bundle loader (Task 0.2.3): `xtalforge_core::data` with embedded gzip
  bundles, SHA-256 verification of the decompressed payload, once-only
  decompression and JSON decoding through `OnceLock`, and lazy external
  bundles for WASM and large databases.
- Data policy page and the `data/` directory convention, including the gzip
  fixture and checksum tests for the loader.
- moyo audit (Task 0.3): requirements matrix, gap decisions, eight-structure
  spglib comparison at four tolerances, wasm/pyo3 build check, and the report
  in `docs/algorithms/symmetry-backend.md`.
- sci-form reuse audit (Task 0.4): verified inventory, port quality gate,
  overlap decisions and the `sci-form` feature policy in
  `docs/algorithms/sci-form-reuse.md`.
- Hello-world slice (Task 0.5): `version()` and `cell_volume()` in the Rust
  API, CLI, Python wheel and WASM package; shared vectors, API pages, parity
  check, binding architecture and error mapping.
- Version and release tooling (Task 0.6): `scripts/sync-version.sh` with a
  self-test, agent-version check, pull request template and release workflow.
- Documentation site (Task 0.7): edit links, search, repository link, API and
  algorithm pages, four-surface code groups, KaTeX, and an offline cell-volume
  WASM demo.
- CI workflows (Task 0.8): Rust matrix, Python, WASM, docs, caching and path
  filters, plus the existing rules job.
- Golden-data oracle (Task 0.9): curated `tests/corpus/` with provenance,
  Phase 1 generator and data, Crystal-Builder test-data inventory,
  `xtalforge-testkit` and the regeneration policy.
- Dependency replacement map (Task 0.10): import survey, replacement
  decisions and license status for every Python dependency.
- Mandatory test layout (AGENTS.md section 2.2): Rust tests live in
  `crates/<crate>/tests/`, independent of the source files.
- Phase 1 Task 1.1: periodic table, deuterium and dummy species, Pyykko radii,
  symbol parsing, `element_info` in Rust, CLI, Python and WASM, the generated
  API table and the shared element vectors.
- Phase 1 Task 1.2: `Lattice` with matrix rows, standard orientation,
  parameters, volume, metric tensor, reciprocal cells, d-spacing and plane
  normals, coordinate conversions, wrapping, exact minimum-image search,
  transforms, strain, volume scaling and validation; `cell --params` in the
  CLI and `lattice_info`/`lattice_minimum_image` in Python and WASM.
- Phase 1 Task 1.3: the limits model with Standard/Generous/WarnOnly
  profiles, soft and hard verdicts, drawn/supercell/grid/phonon/Hessian
  estimates, `Error::LimitExceeded`, and `limits_bounds`/`limits_check` in
  Rust, Python and WASM.
- Phase 1 Task 1.4: `Site` (label, canonical element, occupancy, Uiso/Uij,
  charge, disorder group, props), `Structure`/`StructureBuilder`, molecule in
  a padded box, formula/mass/density/charge properties, the `info` CLI bridge
  and `structure_info` in Python and WASM.
- Phase 1 Task 1.5: `StructureDocument` with `SCHEMA_VERSION`, strict
  top-level fields, version checks and a versioned fixture.
- Phase 1 Task 1.6: `SymOp` parsing/printing/composition/inversion and
  `P1Cell`/`P1Atom` expansion with site/operation/lattice provenance and
  special-position deduplication.
- Phase 1 Task 1.7: diagonal and general supercells as P1, origin shifts,
  metadata/site-property preservation and the supercell limit check.
- Phase 1 Task 1.8: surface bases for arbitrary Miller planes, slab
  construction with layers, vacuum and P1 output, and plane angles.
- Phase 1 Task 1.9: finite-cell images, Cartesian point transforms, plane
  fitting/alignment and displacement ellipsoids.
- Phase 1 Task 1.10: periodic neighbor search with images, distances and
  deterministic ordering; `neighbors` and `neighbors_of`.
- Phase 1 Task 1.11.1/1.11.4: distance-based bond perception with covalent
  radii, dummy exclusion, bond images and the connectivity graph.

### Changed

- License is now `MIT OR Apache-2.0` (confirmed by the maintainer,
  Task 0.1.3), replacing the initial MIT-only `LICENSE` with `LICENSE-MIT` and
  `LICENSE-APACHE`.
- The English-only rule is retired (maintainer decision, 2026-11-08). The
  language checker, allow list, Spanish fixtures and commit-msg hook were
  removed; no language check runs in CI.
- Workspace lints deny warnings and missing docs. `scripts/check-rules.sh`
  rejects `allow`/`expect` attributes for dead code, unused lints and warnings.
- Commit discipline adopted (AGENTS.md section 2.1): one directory per commit,
  `feat(<scope>):` or `ref(<scope>):` subjects, and the completed tasks listed
  from hardest to easiest in the body. The optional pre-commit and commit-msg
  hooks enforce the directory scope and the message format.
- Workspace `unsafe_code` changed from `forbid` to `deny` so the pyo3 binding
  crate can relax it locally with a comment, as the roadmap allows.
- GitHub Actions workflows are committed but disabled at repository level by
  maintainer instruction; they have not run yet.
