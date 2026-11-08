# Phase 0 — Foundations

**Goal:** a working multi-crate workspace with the rules enforced by tooling,
CI on three operating systems, a documentation site, verified decisions on the
symmetry backend and on every Python dependency that must be replaced, a
concrete reuse plan for sci-form code, the binding architecture for the four
surfaces, and a golden-data oracle produced by Crystal-Builder.

**Exit criteria:**

- `cargo test --workspace` passes (the workspace may be almost empty).
- CI is green on Linux, macOS and Windows for Rust, Python, WASM, docs and rules.
- `npm run docs:build` works and the site is deployed.
- The hello-world slice (Task 0.5.7) works in all four surfaces.
- The reports of Tasks 0.3, 0.4 and 0.10 are written in `docs/algorithms/`.
- Golden data for the Phase 1 features exist in `tests/golden/`.

Nothing is published to crates.io, PyPI or npm in this phase except the name
reservations of Task 0.6.5.

---

## Task 0.1 — Repository and workspace skeleton

- [x] 0.1.1 Create the Cargo workspace
  Root `Cargo.toml` with `[workspace]`, `resolver = "3"`, and
  `[workspace.package]`: `version = "0.0.0"`, `edition = "2024"`,
  `rust-version = "1.89"`, `license`, `repository`, `homepage`,
  `documentation`, `authors`, `keywords`, `categories`. Shared versions in
  `[workspace.dependencies]`: `nalgebra = { version = "0.35", default-features = false, features = ["std"] }`,
  `serde` (derive), `serde_json`, `thiserror = "2"`, `rayon` (optional
  everywhere), `moyo = "0.22"`, `flate2`. Create only the crates needed by
  Phases 0 and 1: `xtalforge-core`, `xtalforge-symmetry`, `xtalforge-io`, the
  facade `xtalforge`, `xtalforge-cli`, `xtalforge-python`, `xtalforge-wasm`.
  Each uses `version.workspace = true` and the workspace lints.
  Acceptance: `cargo check --workspace` and `cargo test --workspace` pass on a
  clean checkout.
- [x] 0.1.2 Toolchain, formatting and lints
  `rust-toolchain.toml` pinned to stable with `rustfmt`, `clippy` and the
  `wasm32-unknown-unknown` target. `rustfmt.toml` (`max_width = 100`,
  `edition = "2024"`). `[workspace.lints]`: `unsafe_code = "forbid"` (binding
  crates may relax it locally with a comment explaining why),
  `missing_docs = "warn"` (raised to `deny` in Task 8.7.1), clippy
  `unwrap_used`, `expect_used`, `panic`, `indexing_slicing` as `warn` in
  library crates and `allow` in tests.
  Note: `unsafe_code` is `deny`, not `forbid`, so the pyo3 binding crate can
  relax it locally with a comment as the roadmap allows; it is the only
  exception.
- [x] 0.1.3 Repository hygiene and license
  `.gitignore` (done), `.editorconfig`, `LICENSE`, `THIRD_PARTY.md` (draft
  exists), `deny.toml` for `cargo-deny`: allowed licenses (MIT, Apache-2.0,
  BSD-2/3, ISC, Unicode-3.0, Zlib), advisories, bans on duplicate major
  versions of `nalgebra`, and a ban on crates that build C code in default
  features (`cc` reachable only through optional features). License:
  `MIT OR Apache-2.0`, confirmed by the maintainer on 2026-11-08 and recorded
  in `AGENTS.md` section 6; the texts are `LICENSE-MIT` and `LICENSE-APACHE`.
- [x] 0.1.4 Facade crate
  `crates/xtalforge` re-exports sub-crates under stable module names
  (`xtalforge::core`, `::symmetry`, `::io`, later `::ff`, `::analysis`, …) and
  a `prelude`. Feature flags: `parallel` (default, enables rayon in all
  crates), `sci-form` (quantum methods, Task 4.10), `ml` (Task 8.4),
  `external` (subprocess engines, Task 4.9), `mof-data` (embedded MOF bundle,
  Task 6.11.3). Document every feature in the crate docs and in
  `docs/guide/features.md`.
  Note: `parallel` is wired to optional Rayon in the domain crates; the four
  heavy features are declared and documented as reserved until their tasks
  land.
- [x] 0.1.5 Contributor guide
  `docs/guide/contributing.md`: the hard rules from `AGENTS.md`, the roadmap
  workflow, how to run every check locally, how to add a format, a command, a
  calculator, a binding, and how to write an algorithm page.

## Task 0.2 — Rule enforcement tooling

- [x] 0.2.1 Harden `scripts/check-rules.sh`
  The script fails on code files over 500 lines and on Rust warning
  suppression attributes. Add a self-test: `scripts/tests/rules/` with sample
  files that must pass (ordinary documentation, an exactly 500-line file) and
  must fail (a 501-line file, `#[allow(dead_code)]`,
  `#[expect(unused_variables)]`), and a runner
  `scripts/tests/run-rules-selftest.sh`.
  Note (2026-11-08): the language check was retired by the maintainer; its
  helper, allow list and Spanish fixtures were removed. Rule 2 is now the
  warning-suppression check. Negative fixtures use the `.fixture` suffix so the
  normal repository scan does not read them; the runner copies them to a
  temporary tree with their real extensions.
- [x] 0.2.2 Git hooks
  `scripts/install-hooks.sh` installs a pre-commit hook running
  `check-rules.sh`, `cargo fmt --check`, and `cargo clippy` on changed crates.
  Documented in the contributor guide; never mandatory for CI.
  Note: the language-specific commit-msg hook was removed with the language
  rule. A new commit-msg hook enforces the `feat`/`ref` subject and the task
  list from `AGENTS.md` section 2.1; the pre-commit hook rejects commits that
  span more than one directory.
- [x] 0.2.3 Large-data policy and loader
  Convention: any table with more than ~100 rows, any parameter set copied
  from a paper, and every database (RCSR, PORMAKE, fragments, fullerene
  isomers) lives in `data/` as JSON, compressed JSON (`.json.gz`) or a compact
  binary. Provide `xtalforge_core::data`: a loader for embedded compressed
  blobs (flate2, pure Rust backend `miniz_oxide`), parsed once with
  `std::sync::OnceLock`, with a SHA-256 checksum test per bundle so a corrupt
  or edited bundle is caught. Also support lazily loaded external bundles
  (needed for WASM and for large MOF data, Task 6.11.3).
  Note: implemented as `EmbeddedBundle`, `JsonBundle`, `load_external_gzip`
  and `load_external_json` in `xtalforge_core::data`, with a gzip fixture and
  checksum, digest, corruption and JSON tests. The policy is in
  `docs/guide/data.md`; embedded files live inside the owning crate for
  packaging, the workspace `data/` holds canonical and external bundles.
- [x] 0.2.4 CI job for rules
  A CI job that runs `check-rules.sh` and the self-test on every push and pull
  request; it is a required status check.
  Note: implemented in `.github/workflows/rules.yml` as `Rules / Repository
  rules`; enabling it as a required check is a repository setting that must be
  done once the remote exists.

## Task 0.3 — moyo capability audit (decision gate)

Decision already taken: use `moyo` **if** it lets us do everything required.
This task produces the evidence and the plan for every gap.

- [x] 0.3.1 Requirements matrix
  Read `Crystal-Builder/xtal/core/{spacegroup,symmetry,subgroups}.py`, the
  point-group part of `core/molecule.py`, and `analysis/kpath.py`, and list
  every symmetry function they use (they rely on gemmi tables and spglib).
  For each, record the moyo 0.22 API that covers it. Already verified on
  2026-11-08 by reading the moyo 0.22.0 sources:

  | Need | moyo 0.22 | Status |
  |---|---|---|
  | Space group detection with tolerance | `MoyoDataset::new(cell, symprec, angle_tolerance, setting, rotate_basis)` | public |
  | Number, Hall number, H-M symbol, Pearson symbol | `MoyoDataset` fields | public |
  | Operations, orbits, Wyckoff letters, site symmetry symbols | `MoyoDataset` fields | public |
  | Standardized and primitive standardized cells, transformations, origin shift, rotation | `std_cell`, `std_linear`, `std_origin_shift`, `std_rotation_matrix`, `prim_std_*`, `mapping_std_prim` | public |
  | Wyckoff tables per Hall number | `data::iter_wyckoff_positions_from_hall_number` | public |
  | Hall symbol database, settings, centering, arithmetic crystal class | `data::*` | public |
  | Group identification from operations | `identify` module | public |
  | Translationengleiche and klassengleiche subgroups | `subgroup::enumerate_*` | public |
  | Euclidean normalizer | `MoyoDataset::euclidean_normalizer` | public |
  | Layer groups, magnetic space groups | `MoyoLayerDataset`, magnetic datasets | public |
  | Niggli, Delaunay, Minkowski reduction | `math::{niggli,delaunay,minkowski}` | **crate-private** (`pub(super)`) |
  | Symmetrization / standardization internals | `symmetrize` module | **private**, only reachable through the dataset |
  | Brillouin zone, k-path | — | **absent** |
  | Point groups of molecules (non-periodic) | — | verify; likely absent |

- [x] 0.3.2 Gap decisions
  For every gap, write the decision: (a) open an upstream issue/PR to make the
  reductions public; (b) meanwhile port Niggli/Delaunay into
  `xtalforge-symmetry` with attribution (Task 2.6); (c) implement the
  Brillouin zone and k-path ourselves (Task 2.10), because Crystal-Builder gets
  them from ASE `Cell.bandpath`; (d) implement molecular point groups
  ourselves (Task 2.8.3). Verify also: whether `std_cell` positions are
  idealized (symmetrized) or only transformed, how partial occupancy must be
  encoded as atom types, and the behaviour for non-standard input settings.
  Note: `std_cell` refines both lattice and positions to the detected symmetry,
  not only transforms them; partial occupancy must be encoded as distinct
  species/atom types because moyo takes atomic numbers, not occupancies;
  `Setting::Spglib` reproduces spglib Hall numbers, while origin choice still
  changes some Wyckoff letters. Details in
  `docs/algorithms/symmetry-backend.md`.
- [x] 0.3.3 Behavioural comparison against spglib
  With the golden data of Task 0.9 (produced by spglib through
  Crystal-Builder), run moyo on the whole corpus and compare: space group
  number, Hall number, Wyckoff letters, number of operations, standardized
  lattice parameters. Repeat for `symprec` in {1e-5, 1e-3, 1e-2, 0.1}. Record
  every disagreement and its cause in the report.
  Note: `tests/golden/moyo-audit.json` holds eight structures x four symprec
  values; number, Hall number and operation count agree with
  `Setting::Spglib`. Si Wyckoff letters and four standardized-lattice cases
  are recorded as known origin/setting differences in the report.
- [x] 0.3.4 Build, target and size check
  Confirm that moyo compiles for `wasm32-unknown-unknown` (upstream ships
  `moyo-wasm`, so this is expected) and inside a pyo3 extension; confirm it
  pulls no native code; measure the compile time and the `.wasm` size it
  adds; check that its dependency versions (`nalgebra` 0.35, `serde`, `strum`
  0.28, `itertools`) do not duplicate ours.
  Note: wasm32 and pyo3 builds pass; nalgebra 0.35 is shared; release compile
  about 45 s; the hello-world WASM was 41 KiB because the unused moyo code is
  tree-shaken, and the Phase 1 element table brings it to about 190 KiB.
- [x] 0.3.5 Decision report
  `docs/algorithms/symmetry-backend.md`: verdict (adopt / adopt with gaps /
  reject), the matrix above with final statuses, the gap plan, links to
  upstream issues, and the isolation design: `xtalforge-symmetry` exposes our
  own types (`SpaceGroup`, `SymOp`, `SymmetryDataset`, `WyckoffSite`) so moyo
  never appears in the public API and could be replaced.

## Task 0.4 — sci-form reuse audit and port plan

Decision taken: reuse what we need from sci-form, because it is tested and in
production. Method: **port** (copy into our crates with attribution, adapted
to our types, the 500-line rule and our error policy) the crystal-related
code; depend on the `sci-form` crate itself, behind the optional `sci-form`
feature, only for heavy quantum chemistry (xTB, EHT, PM3), which is not worth
duplicating.

- [x] 0.4.1 Inventory with sizes and split plan
  Verify each row by reading the code, not just the file names; record the
  public functions, their tests, and how well they match the Crystal-Builder
  behaviour.

  | Source in `sci-form/src` | Lines | Plan |
  |---|---|---|
  | `materials/cif.rs` | 708 | Port into `xtalforge-io::cif`; split into tokenizer, loop/table parser, item resolver, writer |
  | `materials/cell.rs` | 1048 | Contains lattice, supercell, PXRD, symmetry expansion and porosity: split into five modules across `core`, `analysis`, `symmetry` |
  | `materials/space_group_table.rs` | 6750 | Do not port; replaced by moyo |
  | `materials/space_groups.rs` | 273 | Do not port; replaced by moyo; keep its tests as extra checks |
  | `materials/scattering_factors.rs` | 396 | Port into `xtalforge-analysis::pxrd`; coefficients to `data/` |
  | `materials/assembly.rs` | 829 | Study for the MOF builder (Task 6.12.5); the Crystal-Builder PORMAKE path is the reference |
  | `materials/sbu.rs`, `topology_data.rs` | 257 / 734 | Study; data to `data/` if used |
  | `materials/geometry_opt.rs`, `variable_cell.rs` | 661 / 501 | Port into `xtalforge-ff::optimize`; split |
  | `forcefield/uff.rs`, `params.rs`, `atom_typer.rs`, `minimizer.rs` | 513+ | Compare with Crystal-Builder `ff/uff/`; port the better-tested parts; parameters to `data/` |
  | `charges_eeq/*` | ~2800 | Port EEQ as an extra charge method; `eeqbc_data.rs` (1673) becomes a data file |
  | `smiles.rs`, `graph.rs`, `stereo.rs` | 875 / 627 / 1121 | Port into `xtalforge-build::smiles`; split each |
  | `rings/sssr.rs`, `rings/aromaticity.rs` | 568 / 219 | Port into `xtalforge-core::rings` (make periodic-aware); split `sssr.rs` |
  | `periodic.rs` | 399 | Port into `xtalforge-core::elements`; data to `data/` |
  | `surface/sasa.rs` | 335 | Port into `xtalforge-analysis` (accessible surface helpers) |
  | `eht/marching_cubes.rs` | — | Not used for periodic isosurfaces: Crystal-Builder deliberately uses marching tetrahedra (Task 5.4.1) |
  | `kmesh.rs` | 1094 | Study for k-point meshes (phonon DOS, Task 4.8.3) |
  | `fd_gradient.rs` | 86 | Port as the finite-difference gradient checker used in tests |
  | `xtb/`, `eht/`, `pm3/`, `conformer.rs`, `etkdg.rs` | large | Not ported; used through the optional `sci-form` feature |

- [x] 0.4.2 Quality gate for every ported file
  Checklist applied to each port and recorded in the pull request: read the
  file fully; remove `unwrap`/`expect`/panics on input; adapt to our types;
  translate any non-English text (none found in a first scan on 2026-11-08;
  accented author names such as Rappé are allowed); split under 500 lines;
  keep or recreate its tests; add the origin comment and the
  `THIRD_PARTY.md` entry.
- [x] 0.4.3 Overlap analysis with Crystal-Builder
  Where both projects implement the same thing (UFF, CIF, PXRD, porosity,
  rings, MOF assembly), compare behaviour on the corpus and decide which one
  is the reference for golden tests. Crystal-Builder is the default reference
  because xtalforge must reproduce it; sci-form wins only where it is proven
  more accurate, and the choice is documented.
- [x] 0.4.4 Policy for the `sci-form` feature
  Decide the version requirement (crates.io `sci-form = "1"`, not a path
  dependency, so published crates build), how its types are converted to
  ours, how the feature appears in Python wheels and in the CLI, and confirm it
  stays out of the default WASM build (size).
  Note: the audit, inventory, quality gate, overlap table and feature policy
  are recorded in `docs/algorithms/sci-form-reuse.md`.

## Task 0.5 — Binding architecture (Rust, CLI, Python, WASM at parity)

- [x] 0.5.1 Boundary design
  Every public operation is a function over `serde`-serializable data types
  defined once in Rust. Small values (cell parameters, options, results) cross
  the boundary as plain data; large objects (`Structure`, `Document`,
  trajectories, grids) are opaque handles with accessor methods, so data is
  not copied on every call. Arrays of coordinates cross as NumPy arrays
  (Python) and `Float64Array` (WASM). Write the design in
  `docs/api/bindings.md`.
- [x] 0.5.2 Error mapping
  Each crate's `Error` exposes `code()` (stable snake_case string),
  `message()`, and optional details (line number, atom indices). Mapped to a
  Python exception hierarchy rooted at `xtalforge.XtalforgeError` with one
  subclass per domain, to JS `Error` objects with a `code` property, and to CLI
  exit codes (2 usage, 3 input, 4 computation, 5 limits, 6 engine missing) plus
  a JSON error object under `--json`.
- [x] 0.5.3 CLI skeleton (`xtalforge-cli`)
  `clap` (derive) with subcommands grouped by domain. Global flags: `--json`,
  `--quiet`, `--threads N`, `--limits {standard,generous,warn-only}`. Reads a
  path or `-` for stdin, writes a path or stdout. Subcommands are added in each
  phase. Mirror the Crystal-Builder CLI verbs where they exist (`info`,
  `symmetry`, `convert`, `bonds`, `types`, `formats`, `modules`, `engines`,
  `run`, `capabilities`, `mcp`) so its users find the same names; record the
  mapping in the migration guide (Task 8.7.4).
- [x] 0.5.4 Python skeleton (`xtalforge-python`)
  pyo3 + maturin, abi3 wheels with minimum Python 3.9, no Python runtime
  dependencies except NumPy. Type stubs (`.pyi`) checked by `mypy --strict` in
  CI, `py.typed`. Module layout mirrors the Rust facade (`xtalforge.core`,
  `xtalforge.symmetry`, `xtalforge.io`, …). Optional pure-Python helpers
  converting to and from ASE `Atoms` and pymatgen `Structure` when those
  packages are installed (never required).
  Note: interop helpers are deferred to Phase 1, when the public `Structure`
  type they convert exists; the pyo3 module, stubs, `py.typed`, wheel and
  pytest suite are in place.
- [x] 0.5.5 WASM skeleton (`xtalforge-wasm`)
  wasm-bindgen with `wasm-pack` builds for `web`, `nodejs` and `bundler`;
  generated TypeScript definitions; no threads by default; an optional
  threaded build with `wasm-bindgen-rayon` documented; `console_error_panic_hook`
  only in debug builds. Size budget recorded and checked in CI.
  Note: web, nodejs and bundler builds work; the panic hook is debug-only; the
  threaded build is documented as a reserved opt-in; the 200 KiB budget and
  `scripts/check-wasm-size.sh` are in place.
- [x] 0.5.6 Parity matrix and shared test vectors
  `docs/api/parity.md` lists every public capability with four columns and a
  status. `scripts/check-parity.sh` compares the matrix with the actual
  exported symbols (rustdoc JSON, the Python module's `dir()`, the generated
  `.d.ts`, `xtalforge --help`) and fails on undocumented gaps. Shared JSON test
  vectors in `tests/vectors/<capability>/*.json` (input, options, expected
  output, tolerance) are executed by Rust tests, `pytest`, the WASM node tests
  and a CLI test runner; all four must agree.
- [x] 0.5.7 Hello-world slice
  Implement `version()` and `cell_volume(a, b, c, alpha, beta, gamma)` end to
  end in the four surfaces, with test vectors, API pages and parity entries, to
  prove the whole pipeline (build, test, docs, CI) before real features land.

## Task 0.6 — Versioning, release automation and agent memory

- [x] 0.6.1 `scripts/sync-version.sh`
  Reads `[workspace.package].version`; updates the `=X.Y.Z` pins between
  internal crates, `docs/package.json` and README install snippets; verifies
  that the maturin metadata and the generated `package.json` report the same
  version. `--check` mode only reports and exits non-zero on drift. Kept under
  500 lines; covered by a self-test.
- [x] 0.6.2 Changelog policy
  `CHANGELOG.md` follows Keep a Changelog with sections Added, Changed,
  Deprecated, Removed, Fixed, Security, and `Breaking` during `0.0.x`.
  `docs/changelog.md` includes the root file with VitePress
  `<!--@include: ../CHANGELOG.md-->` so there is a single source. Every pull
  request adds an `Unreleased` entry.
- [ ] 0.6.3 Release workflow
  Triggered by a tag `vX.Y.Z`: checks version sync and changelog; runs the full
  test matrix; publishes crates in dependency order (`cargo publish` with
  retries for index propagation); builds and uploads wheels (maturin-action,
  trusted publishing); publishes the npm package (provenance enabled); creates
  the GitHub release with notes from the changelog; deploys the docs. A
  `workflow_dispatch` dry-run mode must succeed before `v0.0.1`.
  Note: `.github/workflows/release.yml` is written; the dry run is blocked
  because GitHub Actions are disabled until the maintainer enables them.
- [x] 0.6.4 Agent memory discipline
  CI check: the version in `AGENTS.md` section 6 equals the workspace version
  on tagged commits. Pull request template with the checklist: roadmap ticked,
  parity matrix, docs, changelog, rules script.
- [ ] 0.6.5 Name reservation
  Reserve `xtalforge` on crates.io (publish `0.0.0` placeholder only if the
  maintainer agrees), PyPI and npm (`xtalforge-wasm`, optionally the
  `@xtalforge` scope). All were free on 2026-11-08.
  Note: `xtalforge`, `xtalforge-core` and `xtalforge-wasm` were checked free on
  2026-11-08 (crates.io search empty, PyPI 404, npm 404). No placeholder was
  published; reservation waits for explicit maintainer approval.

## Task 0.7 — Documentation site (VitePress)

- [x] 0.7.1 Site scaffold
  `docs/` already contains the VitePress configuration and the section stubs
  (Guide, API, Algorithms, Demos, Roadmap, Changelog). Run `npm install`,
  commit `package-lock.json`, verify `npm run docs:dev` and `docs:build`,
  enable local search, edit-on-GitHub links, last-updated dates, and the
  repository link once the remote exists.
  Note: lockfile committed, search/edit/last-updated/repository links enabled,
  `npm run docs:build` passes.
- [x] 0.7.2 API reference generation
  One page per domain, each showing the same operation in the four surfaces
  with VitePress code groups (tabs). Sources: rustdoc comments (summaries
  extracted from `cargo +nightly rustdoc -- --output-format json`, or
  hand-written pages checked against it), the `.pyi` stubs, the generated
  `.d.ts`, and `xtalforge <cmd> --help` captured by a script. A CI step fails
  when a public item has no page.
  Note: hand-written pages under `docs/api/`; `scripts/check-api-docs.sh` checks
  every public core item and `scripts/check-parity.sh` checks every matrix
  capability.
- [x] 0.7.3 Algorithms section
  Template for every algorithm page: problem, formulation (KaTeX math via
  `markdown-it-mathjax3`), references, pseudocode, complexity, numerical
  caveats, validation data and tolerances, differences from Crystal-Builder.
  Note: the template is documented in `docs/algorithms/index.md`; the Phase 0
  pages are the moyo, sci-form, dependency, golden-data and validation reports.
- [x] 0.7.4 Demos section
  Interactive demos importing the WASM package built from the workspace (Vite
  alias in dev, the npm package in production builds). A small shared demo
  component library (file input, structure summary table, a minimal 3D
  preview with plain WebGL or three.js, plot canvas). Demos must work offline
  after load.
  Note: `docs/demos/cell-volume.md` runs the WASM package in the browser; the
  VitePress build bundles the `.wasm` asset and the page works offline.
- [ ] 0.7.5 Deployment
  GitHub Pages workflow building on push to `main` and on tags; versioned
  documentation per MINOR release from `v0.1.0`; link checking.
  Note: the `docs` job of `.github/workflows/release.yml` builds and deploys to
  Pages on a version tag, after the verification and WASM jobs and reusing the
  built WASM artifact. The maintainer changed the policy to tag-only
  deployment (no `main` deployment). It has not run because GitHub Actions are
  disabled; versioned docs remain for the first MINOR release.

## Task 0.8 — Continuous integration

- [ ] 0.8.1 Rust matrix
  Linux, macOS (arm64), Windows; stable and MSRV 1.89; `fmt`, `clippy -D
  warnings`, tests with and without default features, doc tests, `cargo-deny`,
  `cargo hack --each-feature` on the facade crate.
- [ ] 0.8.2 Python CI
  maturin wheel build, `pytest` on Python 3.9 and the latest release, `mypy
  --strict` on the stubs, shared test vectors.
- [ ] 0.8.3 WASM CI
  `wasm-pack test --node` and a headless browser run; `.wasm` size budget
  check; TypeScript type check of the generated definitions.
- [ ] 0.8.4 Docs CI
  VitePress build, dead-link check, rules script, parity script.
  Note: `.github/workflows/ci.yml` builds the site with the WASM package and
  runs the parity script; the tag deployment lives in `release.yml` and reuses
  the `wasm` artifact.
- [ ] 0.8.5 Speed
  `Swatinem/rust-cache`, sharing build artifacts between jobs, path filters
  so documentation-only changes skip the heavy jobs; target under 15 minutes.
  Note: `.github/workflows/ci.yml` has the Rust/Python/WASM/docs matrix, caches
  and path filters. GitHub Actions are disabled by explicit maintainer
  instruction, so no job has run; the local equivalents pass. Enable Actions
  and run the workflows before ticking 0.8.1 to 0.8.5.

## Task 0.9 — Golden-data oracle from Crystal-Builder

Crystal-Builder is the behavioural reference; its sample structures and test
data become our oracle.

- [x] 0.9.1 Corpus selection and licensing
  Copy into `tests/corpus/` a curated set from
  `Crystal-Builder/resources/samples/`: `simple/` (small inorganic cells),
  `cod/` (21 MOF CIFs from the COD, CC0: Al-soc-MOF-1, Co-BDP, HKUST-1,
  MIL-100, MIL-101, MIL-53, MIL-88B, MOF-5, MOF-74, MOF-808, Mn-BTT, NU-1000,
  PCN-222, PCN-224, SIFSIX-1-Cu, SIFSIX-3-Ni, UiO-66, UiO-67, ZIF-8,
  cubic-EuHOTP, pbz-MOF-1), `prepared/`, `polymer/`, and the top-level CIFs.
  Read its `PROVENANCE.md` for every file and copy the provenance and license
  into our `tests/corpus/PROVENANCE.md`; skip any file whose license does not
  allow redistribution. Add hand-made edge cases: partial occupancy, disorder
  groups, special positions, a triclinic cell, a hexagonal setting, malformed
  CIFs.
- [x] 0.9.2 Golden-data generator
  `scripts/golden/` (Python, each file under 500 lines): runs in
  a virtual environment with Crystal-Builder installed (`pip install -e
  Crystal-Builder[dev]`, plus the extras each feature needs: `ase`, `build`
  for RDKit, phonopy where used). For each corpus entry and each feature it
  dumps JSON: cell parameters, P1 atoms, supercells, slabs, bonds, rings,
  functional groups, space group dataset, Wyckoff letters, k-path, PXRD
  pattern, UFF energy/forces/stress, charges, porosity numbers, topology.
  Each file records the Crystal-Builder commit hash, the versions of numpy,
  scipy, gemmi and spglib, and the options used. Output goes to
  `tests/golden/<feature>/<structure>.json`; features are generated phase by
  phase, Phase 1 first.
- [x] 0.9.3 Reuse of Crystal-Builder test data
  Inventory `Crystal-Builder/tests/data/` (`blocks`, `dftb` recorded outputs,
  `eqeq`, `rcsr`, `systre`, `zeopp_reference.json`, rod `.xyz` files,
  `pormake_upstream.json`) and map each item to the subtask that will use it.
  Copy only what is used, with provenance.
- [x] 0.9.4 Comparison library
  Rust test helpers in a `xtalforge-testkit` dev crate: float comparison with
  absolute and relative tolerances; atom-set equality up to lattice
  translation and permutation; bond-set equality; pattern comparison for
  PXRD (peak positions and normalized intensities); readable diff output.
  Tolerances per feature recorded in `docs/algorithms/validation.md`.
- [x] 0.9.5 Regeneration policy
  Golden data is regenerated only deliberately, with the reason in the commit
  message, and reviewed as a diff. CI never regenerates it.

## Task 0.10 — Dependency replacement map

Crystal-Builder leans on Python libraries. Each needs a pure Rust replacement
or an explicit decision before its phase starts. Write
`docs/algorithms/dependency-map.md` with one section per dependency.

- [x] 0.10.1 Survey the usages
  `grep` Crystal-Builder for imports of gemmi, spglib, scipy, ase, phonopy,
  rdkit, rietx, matplotlib, and list every function used, in which module,
  and for what.
- [x] 0.10.2 Decide each replacement
  Starting proposal to confirm or change:

  | Python dependency | Used for | Rust plan |
  |---|---|---|
  | gemmi | CIF parsing, space group tables | own CIF parser (sci-form port), moyo tables |
  | spglib | symmetry detection, standardization | moyo |
  | numpy | arrays | `nalgebra`, plain `Vec<[f64; 3]>`, `ndarray` only if needed |
  | scipy `cKDTree` | neighbor search | own cell list; `kiddo` k-d tree if needed |
  | scipy `optimize` | L-BFGS, minimizers, least squares | own L-BFGS/FIRE (sci-form port); Levenberg-Marquardt for Phase 7 |
  | scipy `spatial.Voronoi` | Brillouin zone | own Wigner-Seitz construction (half-space intersection) |
  | scipy `sparse`, `ndimage` | grids, labelling | own periodic flood fill / union-find |
  | ASE | band paths, reading foreign trajectories, ML calculators | own k-path (Task 2.10); own readers; ML via Task 8.4 |
  | phonopy | displacements, force constants, phonon DOS, thermal properties | own implementation (Task 4.8), validated against phonopy outputs |
  | RDKit | SMILES, embedding, sketcher chemistry | sci-form SMILES port; embedding via the `sci-form` feature |
  | rietx | peaks, indexing, Pawley, Rietveld | study then port natively (Phase 7) |
  | matplotlib | colormaps (sampled), plots | colormap tables as data; plots left to consumers |

- [x] 0.10.3 License check of every algorithm source
  For every replacement that ports code (moyo internals, rietx, phonopy
  algorithms, PORMAKE), record the license and the attribution requirements in
  `THIRD_PARTY.md` before porting. Reimplementing from papers is preferred
  where the license is unclear.
