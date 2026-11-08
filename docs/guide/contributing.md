# Contributing

Read `AGENTS.md` at the repository root first: it holds the project rules, the
versioning policy and the current state. This page explains how to turn a
roadmap subtask into a reviewed change.

## Hard rules

1. No code file exceeds 500 lines. Large tables are data files under `data/`.
2. Language is unrestricted. The former English-only rule is retired
   (2026-11-08); no language check runs in CI.
3. The default build is pure Rust. Native libraries, subprocesses and ML
   runtimes are optional features, documented in the crate and in the guide.
4. Surface parity: a public capability ships in the Rust API, the CLI, Python
   and WASM, appears in `docs/api/parity.md` and has a VitePress page.
5. No panics on user input: return `Result` with a stable error code.
6. Deterministic results: same input and options, same output on every platform
   and thread count. Random algorithms take an explicit seed.
7. Docs travel with code: roadmap checkboxes, API pages, algorithm pages,
   parity matrix and the `Unreleased` changelog entry land in the same commit.
8. Ported code and third-party data are recorded in `THIRD_PARTY.md` before the
   port is merged, and every ported file names its origin in a header comment.
9. Every ported algorithm is validated against Crystal-Builder or sci-form
   golden data with documented tolerances.
10. Zero warnings. Compiler, clippy and rustdoc warnings are errors; do not
    silence them with `#[allow(dead_code)]`, `#[allow(unused_*)]`,
    `#[allow(warnings)]` or similar attributes. The rule checker rejects them.

## Roadmap workflow

- Work on one subtask at a time. Read its full description and the reference
  sources it names before writing code.
- Tick `[x]` only when code, tests, surfaces and docs are complete.
- New work discovered during a subtask becomes a new numbered subtask in the
  right task; numbers are never reused.
- If scope changes, strike the old title with `~~text~~` and add a `Note:` line
  that says why and where the work moved. Never delete a subtask.

## Run every check locally

```bash
scripts/check-rules.sh                 # 500-line limit and no warning suppressions
scripts/tests/run-rules-selftest.sh    # the rule checker itself
cargo fmt --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo check -p xtalforge-wasm --target wasm32-unknown-unknown
cargo deny check
cd docs && npm install && npm run docs:build
```

The `Rules` workflow runs the rule check and its self-test, `cargo fmt
--check` and clippy with warnings denied on every push and pull request; enable
both jobs as required status checks in the repository branch protection
settings. CI (Task 0.8) adds the full Rust, Python and WASM matrix.

## Local Git hooks

```bash
scripts/install-hooks.sh
```

The pre-commit hook runs the rule check, `cargo fmt --check`, rejects commits
whose staged files span more than one directory, and runs clippy on the crates
touched by the commit. The commit-msg hook enforces the `feat`/`ref` subject
and the task list in the body. Both are optional (`git commit --no-verify`
skips them) and never replace CI; agents must not bypass them unless the
maintainer explicitly asks. The installer never overwrites a hook it does not
manage.

## Adding a file format

1. Create the reader and writer under `crates/xtalforge-io/src/<format>/`, each
   file under 500 lines, with a module-level doc comment naming the format
   specification and the license of any sample data.
2. Register the format in the registry (Phase 3) and in the format table of the
   API page and the parity matrix.
3. Add corpus files under `tests/corpus/` with provenance in
   `tests/corpus/PROVENANCE.md`; add golden output from Crystal-Builder under
   `tests/golden/` when the format conversion is behaviourally significant.
4. Expose the operation in all four surfaces. Python stubs and the generated
   TypeScript definitions are part of the change.
5. Document the format in the guide and add a changelog entry.

## Adding a command

1. Every mutation is an undoable command in `xtalforge-commands` with a
   documented inverse; it never mutates the document directly.
2. Add the Rust test for apply/undo/redo and a shared JSON vector under
   `tests/vectors/` so the CLI, Python and WASM runners execute the same case.
3. Add the CLI verb (mirroring the Crystal-Builder name when it exists), the
   Python method and the WASM export; update the parity matrix.

## Adding a calculator

1. Implement the `Calculator` trait in `xtalforge-ff`: energy, forces and
   optional stress, with units Angstrom, eV and eV/Angstrom at the boundary.
2. Declare the elements and the system sizes it supports, and return a typed
   error when the structure is outside them.
3. Validate energy, forces and stress against the reference implementation and
   record the tolerances in `docs/algorithms/validation.md`.
4. Add the algorithm page: formulation, references, complexity, numerical
   caveats and differences from Crystal-Builder.

## Adding a binding

1. Define the data type once in Rust with `serde`. Small values cross the
   boundary as plain data; large objects cross as opaque handles.
2. Map the crate error to the surface: a Python exception subclass, a JS
   `Error` with `code`, or a CLI exit code plus the JSON error object.
3. Add the Python stub and the generated TypeScript declaration to the API page.
4. Run the shared vectors on all four surfaces before ticking the subtask.

## Writing an algorithm page

Every page in `docs/algorithms/` follows the same template: problem statement,
formulation with KaTeX math, references, pseudocode, complexity, numerical
caveats, validation data and tolerances, and differences from the reference
implementation. The page is written with the first implementation, never after
the release.

## Commits

A commit contains files from one directory only. A leaf directory (one with no
child directories) is the minimum unit: all related files of that directory
travel together. If a change spans directories, split it into one commit per
directory.

Subjects follow `feat(<scope>): <summary>` for new behaviour and
`ref(<scope>): <summary>` for refactors and maintenance. `<scope>` is the
crate, `docs`, or the path being touched. The body lists the completed tasks
from hardest to easiest, one bullet each:

```
feat(core): add embedded data bundles

- 0.2.3: add the gzip loader with SHA-256 verification and lazy external data
- 0.2.1: extend the rule self-test to warning suppressions
```

These rules are mandatory for every agent and are recorded in `AGENTS.md`
section 2.1.
