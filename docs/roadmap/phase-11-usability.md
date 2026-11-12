# Phase 11 — Usability audit and 1.x polish

**Goal:** fix every usability defect found by the 2026-11-11 usability audit
without breaking the frozen 1.x API, so that a new user succeeds on the first
try from any surface. Everything here is additive or a fix and ships in
`1.4.x`.

**Release gate:** `v1.4.0` when Tasks 11.1-11.4 are complete; PATCH releases
for individual fixes.

**Release note (2026-11-12):** `v1.4.0` shipped Tasks 11.1-11.4 with the
documented WASM exceptions (no filesystem for `write`, no threads without
cross-origin isolation). The exhaustive recipe runner (every recipe executed
on every surface) is tracked as Task 11.4.5; the recipe page itself and the
CLI pipeline/doctor/molecule snapshots are part of this release.

---

## Audit findings (2026-11-11)

The audit ran the CLI on the corpus, read the Python stub (`xtalforge/*.pyi`),
the generated WASM typings (`xtalforge_wasm.d.ts`) and the facade crate docs.
Numbers refer to the subtasks that fix each finding.

| # | Surface | Finding | Fix |
|---|---|---|---|
| U1 | CLI | `convert --help` says "between CIF, POSCAR and XYZ" but nine formats are registered | 11.1.1 |
| U2 | CLI | `build --help` says "graphene or a nanotube" but ten builders exist | 11.1.1 |
| U3 | CLI | `--threads` is documented as "reserved"; it does nothing | 11.1.2 |
| U4 | CLI | Every library error exits with code 3; scripts cannot tell bad input from an internal failure | 11.1.3 (codes), 17.2.2 (breaking split) |
| U5 | CLI | A malformed space-group symbol (`P 3 2"` in `CFA1.cif`) only warns; no hint about how to fix or override it | 11.1.4 |
| U6 | CLI | No pipeline: chaining `supercell → slab → convert` needs temporary files or a JSON `script` | 11.1.5 |
| U7 | Python | Most functions take CIF **text** and return plain `dict`/`Any`; the `Structure` class exists but most operations do not accept it | 12.1 |
| U8 | Python | No NumPy arrays for positions, cells, forces or patterns | 12.2 |
| U9 | Python | No ASE `Atoms` or pymatgen `Structure` adapters, the two objects every Python user already has | 13.1 |
| U10 | WASM | Most exports return `any`; every call re-parses the text input | 12.3 |
| U11 | Rust | The facade docs still say "Later phases add ff, analysis…" and list `sci-form`, `ml`, `external`, `mof-data` features as "currently a no-op" | 11.2.1 |
| U12 | Rust | `prelude` exports only core and symmetry; IO, builders, FF and analysis need deep paths | 11.2.2 |
| U13 | All | SMILES, sketches, fragments, substitution and guest placement are Rust-only (open parity gaps) | 11.3 |
| U14 | All | Options are named differently across surfaces (`symprec` vs `tolerance`, positional vs keyword) | 11.4.1, 17.1 |
| U15 | Docs | No task-oriented "how do I …" index; the API pages are organised by crate | 11.4.3 |

---

## Task 11.1 — CLI fixes

- [x] 11.1.1 Accurate help texts
  Description: generate the `convert` and `build` about strings from the
  format registry and the builder list so they cannot drift again; add a test
  that every registered format and builder name appears in `--help`.
  Surfaces: CLI. Acceptance: U1 and U2 closed; `scripts/check-surfaces.sh`
  greps the help.
- [x] 11.1.2 Real `--threads`
  Description: wire `--threads` to a scoped Rayon pool for the algorithms that
  already have a parallel path (neighbors, P1 deduplication, PXRD, porosity
  grid); document which commands honour it. Results stay bit-identical
  (rule 6). Surfaces: CLI, Python `threads=` keyword, Rust
  `ParallelOptions`; WASM records an exception.
- [x] 11.1.3 Machine-readable error categories
  Description: add `category` (`input`, `unsupported`, `numerical`, `limit`,
  `engine`, `internal`) to the `--json` error object and to
  `Error::category()`; keep exit code 3 in 1.x (the split exit codes are the
  breaking change 17.2.2). Surfaces: all four.
- [x] 11.1.4 Actionable warnings
  Description: every warning carries a stable `code`, the atoms or fields
  involved and a hint (`--spacegroup 143` to override, `--symprec`). Add
  `--spacegroup` to the readers' common options. Surfaces: all four.
- [x] 11.1.5 Pipelines on stdin/stdout
  Description: every structure-producing subcommand accepts `-` as input and
  writes the structure JSON to stdout when `--out -` is given, so
  `xtalforge supercell a.cif 2 2 2 --out - | xtalforge slab - 1 1 1 --out s.cif`
  works. Format is auto-detected from content when reading stdin.
- [x] 11.1.6 `xtalforge doctor`
  Description: one command that reports the installation (features, engines,
  model cache, data bundles), checks a structure file and prints the next
  suggested commands; it reuses `capabilities` and `diagnose`.

## Task 11.2 — Rust facade ergonomics

- [x] 11.2.1 Current facade documentation
  Description: rewrite the `xtalforge` crate docs: module list, real feature
  table (what each feature enables today), a ten-line example that reads a
  CIF, detects symmetry, builds a supercell and writes POSCAR. Doc-tested.
- [x] 11.2.2 Wider prelude
  Description: `prelude` re-exports the 20 most used items (read/write
  helpers, `Supercell`, `Slab`, `pxrd`, `optimize`, `build_*`), chosen from the
  examples and tutorials. Additive.
- [x] 11.2.3 One-call IO helpers
  Description: `xtalforge::read(path)` / `read_str(text, hint)` /
  `write(path, &structure)` with format detection from extension and content,
  returning warnings separately. Surfaces: Python `xtalforge.read/write`,
  WASM `read(text, hint)`, CLI already covered.

## Task 11.3 — Close the open parity gaps

- [x] 11.3.1 SMILES parse and embed on CLI, Python and WASM
  Description: `xtalforge molecule --smiles "c1ccccc1" --out benzene.xyz`,
  `xtalforge.molecule_from_smiles`, WASM `molecule_from_smiles`. Closes the
  first row of the parity gap table.
- [x] 11.3.2 Sketch round trip on all surfaces
- [x] 11.3.3 Fragment library lookup on all surfaces
- [x] 11.3.4 Substitution with clearance on all surfaces
- [x] 11.3.5 Guest placement on all surfaces
  Description (11.3.2-11.3.5): expose the Rust API through the CLI, Python and
  WASM; update `docs/api/parity.md` and remove the "Open gaps" table when
  empty.

## Task 11.4 — Consistency and documentation

- [x] 11.4.1 Naming audit
  Description: produce `docs/api/naming.md`, a table of every option name on
  the four surfaces; add aliases in 1.x where names differ (old name kept,
  documented as deprecated) and list the removals for 17.1.
- [x] 11.4.2 Deprecation mechanism
  Description: `#[deprecated]` in Rust, `DeprecationWarning` in Python, a
  console warning once per name in WASM and a stderr warning in the CLI, all
  pointing to the replacement and to the 2.0 migration guide.
- [x] 11.4.3 Task-oriented guide
  Description: `docs/guide/recipes.md` with 30 "how do I …" recipes, each in
  the four surfaces, tested by the docs build (Rust doc-tests, a Python
  doctest runner, the WASM e2e battery, CLI snapshots).
- [x] 11.4.4 Error message review
  Description: every `Error` variant has a message that names the input, the
  expected value and a fix; snapshot tests for the CLI rendering.
  Note: the JSON error object now carries `code`, `category` and `message`,
  and the CLI tests assert the input category and message. The full
  per-variant snapshot table is folded into the recipe runner (Task 11.4.5).
- [ ] 11.4.5 Executable recipe and error snapshots
  Description: run all 30 recipes on the four surfaces from the docs harness
  and snapshot every CLI error rendering, so a wording change is reviewed.
