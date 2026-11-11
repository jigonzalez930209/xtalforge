# Phase 17 — 2.0 API redesign and freeze

**Goal:** remove what Phases 11-16 deprecated, make the few breaking changes
that cannot be done additively, and freeze the 2.0 API. This is the only
phase allowed to break the 1.x API (AGENTS.md section 5.1, MAJOR bump).

**Release gate:** `v2.0.0-rc.N`, a fixes-only stabilization period, then
**v2.0.0**.

---

## Task 17.1 — Remove deprecated surface

- [ ] 17.1.1 Remove every item deprecated in 1.x (text-in/dict-out Python
  functions, `any`-returning WASM exports, option aliases from 11.4.1).
- [ ] 17.1.2 Unified option names and order across the four surfaces,
  enforced by a generated naming test from `docs/api/naming.md`.
- [ ] 17.1.3 One options struct per operation family with `#[non_exhaustive]`
  and builders, so future fields stop being breaking changes.

## Task 17.2 — Errors and exit codes

- [ ] 17.2.1 One public error type in the facade (`xtalforge::Error`) with
  `code()`, `category()` and source chaining; crate errors convert into it.
- [ ] 17.2.2 CLI exit codes by category (2 usage, 3 input, 4 unsupported,
  5 numerical, 6 limit, 7 engine, 70 internal); Python exception subclasses
  per category; JS `Error.category`.

## Task 17.3 — Data formats

- [ ] 17.3.1 `schema_version` 2 for the serialized structure (properties,
  provenance, disorder, magnetic moments as first-class fields) with a reader
  for version 1.
- [ ] 17.3.2 Project `format_version` bump with a 1.x importer.

## Task 17.4 — Freeze and release

- [ ] 17.4.1 `docs/guide/migration-2.md`: every removed or renamed item with
  its replacement, in the four surfaces, plus an automatic `xtalforge migrate`
  for workflow and project files.
- [ ] 17.4.2 Freeze test extended to the 2.0 surface (Rust public items, CLI
  `--help` and `--json` shapes, `.pyi`, `.d.ts`).
- [ ] 17.4.3 Release `v2.0.0-rc.1`, collect feedback, fixes only, then
  `v2.0.0` with the release checklist (AGENTS.md section 5.4).
