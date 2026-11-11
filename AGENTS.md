# AGENTS.md — Project Memory and Rules

This file is the persistent memory of the repository. Every human and every AI
agent working here MUST read it first and keep it up to date. `CLAUDE.md` only
points here.

## 1. What this project is

`xtalforge` is a pure Rust library for building, manipulating, analysing and
exporting crystal structures. It ports everything the headless core of the
Python application **Crystal-Builder** does (`/home/lestad/github/Crystal-Builder`,
package `xtal/`, about 95 000 lines; the PySide6/VTK GUI in `xtalapp/` is out
of scope). It reuses proven code from **sci-form**
(`/home/lestad/github/sci-form`, by the same maintainer as this project) and
uses the **moyo** crate for symmetry.

It is delivered through four surfaces developed together, always at parity:

1. Rust API (the source of truth; facade crate `xtalforge`)
2. CLI (binary `xtalforge`)
3. Python bindings (pyo3 + maturin; PyPI package `xtalforge`)
4. WebAssembly bindings (wasm-bindgen; npm package `xtalforge-wasm`)

The names `xtalforge` (crates.io, PyPI, npm) and `xtalforge-wasm` (npm) were
free on 2026-11-08. Reserve them early (Task 0.6.5).

## 2. Hard rules (non-negotiable)

1. **No code file may exceed 500 lines.** This applies to every source file
   (`.rs`, `.py`, `.pyi`, `.ts`, `.js`, `.vue`, `.sh`, `.toml`, `.yml`),
   including tests, benches, examples, build scripts and generated code.
   Large tables are stored as data files (`.json`, `.json.gz`, `.bin`,
   `.cif`, `.cgd`, `.csv`) under `data/` and embedded with
   `include_bytes!`/`include_str!` or loaded lazily; they are never written as
   giant Rust literals. Data embedded into a published crate lives in that
   crate's `data/` directory (so `cargo package` includes it); the workspace
   `data/` holds canonical files, external bundles and repository-only assets.
   Markdown (`*.md`) has NO size limit. When a file
   approaches 400 lines, split it by responsibility before adding more.
   Enforced by `scripts/check-rules.sh` locally and in CI.
2. **Language is unrestricted.** The former English-only rule is retired
   (2026-11-08). No language check runs in CI; Spanish is allowed in code,
   docs, commit messages and data keys. Keep identifiers and stable error
   codes self-describing so the API remains usable across languages.
3. **Pure Rust by default.** Default features depend only on pure Rust crates.
   Anything needing a native library (ONNX Runtime, BLAS) or an external
   program (DFTB+, xtb, Zeo++, ORCA, Blender) is optional, behind a feature
   flag or a runtime engine probe, and documented.
4. **Surface parity.** A public capability is not done until it exists in the
   Rust API, the CLI, Python and WASM, is listed in `docs/api/parity.md`, and is
   documented in the VitePress site. A capability that cannot exist on a
   surface (filesystem, subprocesses or threads in WASM) is recorded as an
   explicit, justified exception in the parity matrix.
5. **No panics on user input.** Fallible functions return `Result<_, Error>`.
   `unwrap`, `expect`, `panic!` and unchecked indexing on external data are
   forbidden outside tests (clippy lints enforce `unwrap_used`, `expect_used`).
6. **Deterministic results.** Same input and options give the same output,
   independent of thread count and platform. Random algorithms take an
   explicit seed. Parallel reductions use a fixed order.
7. **Docs travel with code.** Every change updates, in the same commit, the
   relevant roadmap checkboxes, API pages, algorithm pages when behaviour
   changes, the parity matrix, and the `Unreleased` section of `CHANGELOG.md`.
8. **Attribution.** Code ported from sci-form or Crystal-Builder, and every
   third-party data set, is recorded in `THIRD_PARTY.md` with its license
   before it is merged. Each ported file starts with a comment naming its
   origin (`// Ported from Crystal-Builder xtal/core/slab.py (MIT).`).
9. **Golden-data validation.** Every ported algorithm has tests against golden
   data produced by the reference implementation (Task 0.9) with documented
   tolerances. A deliberate deviation is recorded in the algorithm page.
10. **Zero warnings.** Compiler, clippy and rustdoc warnings are errors.
    `[workspace.lints.rust]` denies warnings and missing docs. Never silence a
    warning with `#[allow(dead_code)]`, `#[allow(unused_*)]`,
    `#[allow(warnings)]`, `#[expect(...)]` or an equivalent tag; fix the cause.
    `scripts/check-rules.sh` rejects these attributes. The test-only clippy
    allowances, a binding-signature `too_many_arguments` allowance, and the
    documented `unsafe_code` relaxation in the pyo3 binding crate are the only
    exceptions.

### 2.1 Commit discipline (mandatory for every code agent)

- **Unit.** A commit contains only files that live in the same directory. A
  leaf directory (one with no child directories) is the minimum entity: all
  related files of that directory travel together. The repository root and
  intermediate directories are entities for the files they contain themselves;
  never mix those files with files from a child directory.
- **Subject.** `feat(<scope>): <summary>` for new behaviour or capabilities and
  `ref(<scope>): <summary>` for refactors, cleanups and maintenance. `<scope>`
  is the crate (`core`, `symmetry`, `io`, `cli`, `python`, `wasm`), `docs`, or
  the path being touched (`docs/guide`, `scripts`, `workspace`).
- **Body.** List the completed tasks from hardest to easiest, one `- ` bullet
  per task, with the roadmap number and a short result. A commit with no
  roadmap task lists its work items in the same order.
- **Splitting.** If a change spans directories, split it into one commit per
  directory instead of widening the commit.
- These rules apply to every agent and every change, including documentation,
  CI and small fixes. Amend or split the commit rather than bypassing them.
- `scripts/install-hooks.sh` installs optional hooks that enforce the
  one-directory rule (pre-commit), the subject format and the task list
  (commit-msg). Agents must not use `git commit --no-verify` unless the
  maintainer explicitly asks.

### 2.2 Test layout (mandatory for every code agent)

- Every Rust library crate keeps its tests in `crates/<crate>/tests/`, one
  file per functionality or module. Tests use the public API only.
- Do not put `#[cfg(test)]` unit-test modules in `src/`; the test list of a
  library is independent of its source files.
- Fixtures live in `crates/<crate>/tests/fixtures/` and are loaded with
  `include_bytes!` or `include_str!` from the test file.
- Every public operation has at least one test, including its error paths.
- Test-only clippy allowances may be declared at the top of each test file;
  they are the only `allow` attributes accepted in the repository.
- Integration tests for the CLI, Python and WASM live in the same
  `crates/<crate>/tests/` directory, next to the Rust tests of that surface.

## 3. Roadmap workflow

The roadmap lives in `docs/roadmap/`, one file per phase plus `index.md`:
`phase-00-foundations.md` … `phase-16-v2.md`.

Structure: **Phase → Task (`N.M`) → Subtask (`N.M.K`, with a checkbox)**.

```
## Task 1.2 — Lattice
- [ ] 1.2.1 Subtask title
  Description: what to build, the reference source, inputs and outputs,
  acceptance criteria, surfaces affected.
```

Rules for agents:

- Work on one subtask at a time; read its whole description and the reference
  source files it names before writing code.
- Mark a subtask `[x]` only when code, tests, all four surfaces (when the
  subtask is user-visible) and docs are done and merged.
- Never delete subtasks. If scope changes, strike the title with `~~text~~`
  and add a `Note:` line explaining why and where the work moved.
- New work discovered during a subtask becomes a new subtask with the next
  free number in the right task; numbers are never reused.
- Roadmap files may be as long as needed.

## 4. Repository layout (target)

```
xtalforge/
├─ AGENTS.md  CLAUDE.md  README.md  CHANGELOG.md  THIRD_PARTY.md
├─ LICENSE-MIT  LICENSE-APACHE
├─ Cargo.toml                 workspace; single source of truth for the version
├─ crates/
│  ├─ xtalforge-core/         elements, lattice, sites, structure, P1, supercell, slab,
│  │                          neighbors, bonding, rings, groups, measure, limits
│  ├─ xtalforge-symmetry/     moyo integration, space groups, Wyckoff, reduction,
│  │                          subgroups, point groups, Brillouin zone, k-path
│  ├─ xtalforge-io/           format registry and all readers/writers, project files
│  ├─ xtalforge-commands/     undoable commands, document, selection, tracking, workspace
│  ├─ xtalforge-view/         view settings, draw styles, colormaps, scene data, STL
│  ├─ xtalforge-ff/           Calculator trait, UFF/UFF4MOF, QEq/EQeq, Ewald,
│  │                          optimizers, scans, vibrations, phonons, thermo, engines
│  ├─ xtalforge-analysis/     PXRD, grids, porosity, voids, topology, RCSR, interpenetration
│  ├─ xtalforge-build/        molecules, substituents, fill, carbon, graphene, nanotubes,
│  │                          fullerenes, polymers
│  ├─ xtalforge-mof/          MOF builder (PORMAKE port) and its data
│  ├─ (xtalforge-powder)      planned; the powder code lives in xtalforge-analysis
│  ├─ xtalforge-modules/      module registry, jobs, run records, reports
│  ├─ xtalforge-agent/        session verbs, inspection, diagnostics, MCP server
│  ├─ xtalforge/              facade crate re-exporting everything (published)
│  ├─ xtalforge-cli/          binary `xtalforge`
│  ├─ xtalforge-python/       pyo3 bindings
│  └─ xtalforge-wasm/         wasm-bindgen bindings
├─ data/                      data files embedded or loaded by the crates
├─ tests/                     corpus, golden data, shared test vectors
├─ scripts/                   check-rules.sh, sync-version.sh, golden-data tools
└─ docs/                      VitePress site (guide, api, algorithms, demos, roadmap, changelog)
```

Crates are created when their phase starts; an empty crate is not published.

## 5. Versioning rules

All crates and all four surfaces share **one version** (lockstep). The single
source of truth is `[workspace.package] version` in the root `Cargo.toml`.
Tags are `vX.Y.Z`. Releases happen when a group of tasks is complete, never on
a timer.

### 5.1 Stages and when each number changes

| Stage | Versions | Enters when | Bump rules |
|---|---|---|---|
| Bootstrap | `0.0.0` (never published) | now | no releases |
| Early | `0.0.1`, `0.0.2`, … | Phase 1 complete **and** Task 2.1 merged → tag `v0.0.1` | Every release bumps PATCH, whatever it contains. Breaking changes are allowed but listed under `### Breaking` in the changelog. |
| Pre-stable | `0.1.0` … `0.y.z` | Phase 4 complete (force fields and optimizers usable) → tag `v0.1.0` | MINOR (`0.y`) for any breaking change on any surface, and when a phase completes. PATCH (`0.y.z`) for backward-compatible features, fixes, docs, data updates. |
| Stable | `1.0.0` and later | Phase 8 complete (API freeze, Task 8.8) | Strict SemVer: MAJOR for breaking changes on any surface, MINOR for new backward-compatible capabilities, PATCH for fixes only. |
| Stable 2 | `2.0.0` and later | Phase 16 complete (Task 16.4) | Same SemVer rules; 1.x items removed in 2.0 must have been deprecated in a 1.x release first. |

Release candidates (`X.Y.Z-rc.N`) are used only before `1.0.0` and before any
later MAJOR release.

### 5.2 What counts as public API (a change here can be "breaking")

- Public items of the facade crate `xtalforge` (names, signatures, types).
- CLI subcommand names, flags, exit codes and `--json` output shape.
- Python module, class and function names, keyword arguments, exception types.
- JS/TS exported names and TypeScript types of `xtalforge-wasm`.
- Error codes (`Error::code()` strings).
- The serialized structure schema (`schema_version`) and the project file
  format (`format_version`).

Independent version numbers that do NOT follow the crate version:
`schema_version` (integer, bumped on any change to serialized types),
`format_version` of project files (integer), and data bundles
(`data-YYYY.MM.DD`, for example the MOF and RCSR bundles).

### 5.3 Locations that carry the version

| Location | How it is kept in sync |
|---|---|
| `Cargo.toml` `[workspace.package].version` | edited by hand (source of truth) |
| Every `crates/*/Cargo.toml` | `version.workspace = true`; internal deps pinned `version = "=X.Y.Z"` by the script |
| `crates/xtalforge-python/pyproject.toml` | `dynamic = ["version"]`; maturin reads Cargo |
| `crates/xtalforge-wasm` generated `package.json` | generated by `wasm-pack`; the script verifies it |
| `docs/package.json` | rewritten by the script |
| `CHANGELOG.md` and `docs/changelog.md` | new section per release (`docs/changelog.md` includes the root file) |
| `README.md` install snippets | rewritten by the script |
| This file, section 6 | edited by hand at every release |

`scripts/sync-version.sh` (Task 0.6.1) writes the derived locations;
`scripts/sync-version.sh --check` fails CI on any drift.

### 5.4 Release checklist (copy into the release PR)

**There is no GitHub CI.** The maintainer runs the full local CI on the exact
commit that will be tagged; the tag-driven `release.yml` only builds and
publishes, and runs no tests. Never create a tag before this gate passes.

1. All roadmap checkboxes for the shipped scope are ticked.
2. Run the **local CI gate** from the repository root on the commit to tag:

   ```bash
   scripts/check-rules.sh
   scripts/sync-version.sh --check
   scripts/check-agent-version.sh
   scripts/check-parity.sh
   scripts/check-api-docs.sh
   cargo fmt --check
   cargo clippy --workspace --all-targets --all-features -- -D warnings
   RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
   cargo deny check
   cargo check -p xtalforge-wasm --target wasm32-unknown-unknown
   (cd crates/xtalforge-wasm && wasm-pack build --release --target web)
   scripts/check-wasm-size.sh
   (cd crates/xtalforge-wasm && wasm-pack test --node --release)
   scripts/check-surfaces.sh
   (cd crates/xtalforge-python && maturin build --release)
   python -m mypy --strict crates/xtalforge-python/python/xtalforge
   (cd docs && npm ci && npm run docs:build)
   ```

   The Python commands use the virtual environment that has the freshly built
   wheel installed. The docs link check is the same one the release job runs.
3. Decide the new version with section 5.1; bump the root `Cargo.toml`; run
   `scripts/sync-version.sh`.
4. Move `Unreleased` in `CHANGELOG.md` into the new version section with the
   date; update the parity matrix and section 6 below.
5. Commit and push everything; the working tree must be clean.
6. Tag the commit that passed the local gate: `git tag -a vX.Y.Z -m
   "xtalforge X.Y.Z"` and `git push origin vX.Y.Z`.
7. The `Release` workflow builds the three CLIs, creates the GitHub release,
   publishes PyPI/npm/crates.io from `ubuntu-24.04`, and deploys the docs in
   that order. It does not run tests.

## 6. Current state (update at every release)

- Version: `1.1.0` (published 2026-11-11 on crates.io, PyPI, npm, the GitHub
  release and GitHub Pages; `1.0.0` followed the `1.0.0-rc.1` stabilization
  period on 2026-11-10; `0.0.3` was the last pre-1.0 release, and `v0.0.1`
  was the first on 2026-11-09). `1.1.0` is backward compatible: the 3D
  playground, the audit fixes and additive WASM/JSON fields.
- Remote: https://github.com/jigonzalez930209/xtalforge (public). GitHub
  Actions are enabled (2026-11-08); only the tag-driven release workflow is
  kept. CI runs locally before each tag (section 5.4).
- Current phase: Phase 9 — Core features of a standard crystal
  library (`docs/roadmap/phase-09-core-features.md`), the first step of the
  road to `2.0.0` (Phases 9-16, planned 2026-11-11; see `docs/roadmap/index.md`).
  Phases 0-8 are complete. Phase 8 — Agent and
  distribution (`docs/roadmap/phase-08-distribution.md`) delivered tasks
  8.1-8.7:
  the `xtalforge-agent` crate (verbs, session, diagnostics, change bursts,
  SVG/PNG renderer, MCP over stdio, skill documents), the open extension
  registries and `examples/plugin`, the Python ML bridge with the model cache,
  distribution metadata and package checks, benchmarks, robustness/limits
  tests, tutorials, the migration guide and the algorithm coverage page.
  Tasks 8.8.1 and 8.8.2 are done (stability policy, freeze test and parity
  closure). Task 8.8.3 published `v1.0.0-rc.1` and then the stable `v1.0.0`
  on 2026-11-10 after the fixes-only stabilization period: 13 crates on
  crates.io, 7 platform wheels plus an sdist on PyPI, `xtalforge-wasm` on npm,
  the CLI/hosted assets on the GitHub release and the docs on Pages. Task 8.9
  added the 3D playground and the end-to-end audit battery (`tests/e2e/`,
  `scripts/e2e.sh`): 68 Playwright checks against gemmi, spglib, ASE and NumPy
  references plus deterministic screenshot goldens. The 2026-11-11 functional
  audit (Tasks 8.9.6-8.9.7) fixed PXRD form factors, twisted bilayers,
  polymers, the SMILES embedder and disordered carbon, added the playground
  tour, and also fixed the supercell symmetry expansion and the Hall-only
  `is_p1` check. The remaining coverage notes are recorded in
  `docs/guide/audit.md`; Task 8.9.5 closed with `scripts/check-surfaces.sh`
  and the refinement, porosity and MOF browser checks. `v1.1.0` published
  Task 8.9 on 2026-11-11. Phases 1-7 are
  complete and released at `0.0.3`; the remaining reference gaps (held cell
  quantities and pressure in the optimizer, the PXRD/Rietveld profile
  simplifications, bond carrying in slabs, the curated MOF subset and others)
  are listed in `docs/roadmap/index.md` under "Known subsets". `0.0.2` carried the linked-cell neighbor search, grid-based
  P1/supercell deduplication and the Linux desktop/AppStream fix; `0.0.3`
  adds the full-corpus bond validation, periodic Franzblau rings, the CIF
  round-trip fields and lazy XYZ/XDATCAR trajectories.
  Phase 2 covers moyo integration, the space-group catalogue, the
  detect/asymmetrize/P1/set workflow, Wyckoff assignments, standardized and
  primitive cells, lattice reduction, subgroup enumeration and lowering,
  tolerance diagnostics, the Brillouin zone, layer/magnetic/molecular point
  groups, band paths and change of setting. Phase 3 covers the format
  registry, atomic writes, the extra formats, grids, trajectories,
  commands/document, selection/tracking, workspace/modules, view data, STL,
  project files, parameter presets and the editor demo. Phase 4 covers the
  Calculator API, native UFF/UFF4MOF, charges/Ewald, FIRE/L-BFGS, scans,
  vibrations/thermochemistry, EEQ, the preparation pipeline, engine
  writers/probes and the quantum adapter contract. Phase 5 covers PXRD,
  porosity, isosurfaces, net quotient graphs, the RCSR subset, merge/overlay,
  peaks, indexing, Le Bail/Pawley and the Rietveld phase fit. Phase 6 covers
  graphene, nanotubes, C60/C20, SMILES parsing/embedding, sketches, fragments,
  the net/structure and MOF builders, twists, ribbons, carbon surfaces and
  polymers. Phase 7 covers measured patterns/instrument, whole-range and
  structural refinement, cells/extinction, the automatic run and the powder
  surfaces. Phase 2-7 releases follow the `0.0.x` line.
- Phase 0 blockers: Tasks 0.6.3 (release workflow), 0.6.5 (name reservation)
  and 0.7.5 (docs deployment) completed with the `v0.0.1` release on
  2026-11-09. Tasks 0.8.1 to 0.8.5 (CI) remain open by maintainer decision on
  2026-11-08 because the runner cost was not wanted; the checks are run
  locally. The names were free on 2026-11-08 and are now reserved by the
  publication.
- Testing: integration tests live in `crates/<crate>/tests/` (AGENTS.md
  section 2.2). Shared vectors are in `tests/vectors/`, the curated corpus in
  `tests/corpus/`, the Crystal-Builder golden data in `tests/golden/`, and the
  external data sets in `tests/external/` (CoRE MOF 2019 and tmQM 2024
  subsets, with licenses and citations in their READMEs). The agent tests add
  the verb catalogue, session, diagnostics, MCP, skill, renderer, MCP JSON-RPC
  and the ML bridge golden data.
- Data policy: embedded bundles live inside the owning crate
  (`crates/<crate>/data/`) so `cargo package` carries them; the workspace
  `data/` holds canonical files, external bundles and repository-only assets.
  Loader: `xtalforge_core::data` (Task 0.2.3).
- Workspace members (Task 0.1.1, extended by later phases): `xtalforge-core`,
  `xtalforge-symmetry`, `xtalforge-io`, `xtalforge-commands`, `xtalforge-view`,
  `xtalforge-ff`, `xtalforge-analysis`, `xtalforge-build`, `xtalforge-mof`,
  `xtalforge-modules`, `xtalforge-agent`, the facade `xtalforge`,
  `xtalforge-cli`, `xtalforge-python`, `xtalforge-wasm`, the test crate
  `xtalforge-testkit` and the unpublished example crate
  `xtalforge-plugin-example` under `examples/plugin`. Rust edition 2024, MSRV 1.89, license
  `MIT OR Apache-2.0` (`LICENSE-MIT`, `LICENSE-APACHE`).
- Decisions taken (2026-11-08): name `xtalforge`; symmetry through `moyo`
  (0.22, MIT OR Apache-2.0) behind our own types, subject to the audit in
  Task 0.3; reuse proven sci-form code by porting it (Task 0.4); the four
  surfaces are developed in parallel; VitePress docs from day one; license
  `MIT OR Apache-2.0` (confirmed by the maintainer on 2026-11-08, Task 0.1.3);
  the English-only rule is retired (2026-11-08), warnings are denied and
  warning-suppression attributes are forbidden (rule 10); commits are
  directory-scoped with `feat(scope):`/`ref(scope):` subjects and a
  hardest-to-easiest task list in the body (section 2.1).
- Facts verified on 2026-11-08: moyo's Niggli, Delaunay and Minkowski
  reductions exist but are crate-private (`pub(super)`); its standardization
  and symmetrization module is private but exposed through `MoyoDataset`;
  subgroup enumeration and Wyckoff tables are public; a WASM binding of moyo
  exists upstream (`moyo-wasm`). PORMAKE (vendored in Crystal-Builder) is MIT.
  COD sample CIFs are CC0. rietx (used by Crystal-Builder for powder
  refinement) is a separate Python project (https://github.com/yue-here/rietx).
- Known risks: see `docs/roadmap/index.md`.

## 7. Source map (where the reference behaviour lives)

| Topic | Reference implementation |
|---|---|
| Everything Crystal-Builder does | `/home/lestad/github/Crystal-Builder/xtal/` |
| Its design principles | `/home/lestad/github/Crystal-Builder/docs/PLAN.md` |
| CIF, cell, PXRD, porosity, MOF assembly, UFF, SMILES, rings, EEQ | `/home/lestad/github/sci-form/src/{materials,forcefield,rings,charges_eeq,smiles.rs,graph.rs,periodic.rs}` |
| Symmetry | `moyo` crate (https://github.com/spglib/moyo, docs https://docs.rs/moyo) |
| Powder refinement | rietx (https://github.com/yue-here/rietx) as used by `xtal/powder/bridge.py` |
| Test corpus | `Crystal-Builder/resources/samples/` (see its `PROVENANCE.md`) and `Crystal-Builder/tests/data/` |

## 8. Conventions

- Rust edition 2024. MSRV 1.89 (required by nalgebra 0.35; recorded as
  `rust-version` in `Cargo.toml` and tested in CI).
- Linear algebra: `nalgebra` 0.35 (the version used by moyo and sci-form).
- Units: Angstrom and degrees at every API boundary; radians internally; eV,
  eV/Angstrom, GPa at the boundary for energies, forces and stress.
- Fractional coordinates are `[f64; 3]`; wrapping into `[0, 1)` happens only
  where a function documents it.
- Errors: one `Error` enum per crate with `thiserror`, a stable `code()` string,
  converted at the bindings (Python exception classes, JS `Error.code`, CLI
  exit code and JSON error).
- Binding boundary: plain `serde` data structs, or opaque handles for large
  objects. JSON is the interchange format for the CLI and WASM.
- Long operations take a cancellation token and a progress callback.
- Commit messages: one directory per commit, `feat(<scope>): <summary>` or
  `ref(<scope>): <summary>`, with the completed tasks listed from hardest to
  easiest in the body (section 2.1).
