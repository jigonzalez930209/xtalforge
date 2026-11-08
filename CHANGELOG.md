# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/) (pre-1.0: `0.0.x` until the first
essential release, then `0.x.y`).

## Unreleased

### Added

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
