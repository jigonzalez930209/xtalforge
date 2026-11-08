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
