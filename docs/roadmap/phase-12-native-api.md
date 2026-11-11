# Phase 12 — Native, object-oriented bindings

**Goal:** make Python and JavaScript feel native: objects in, objects out,
arrays without copies where possible, full type information, no re-parsing.
Phase 12 is additive in `1.x`; the text/dict functions it supersedes are
deprecated (Task 11.4.2) and removed in 2.0 (Task 17.1).

**Release gate:** `v1.5.0`.

---

## Task 12.1 — Python object API

- [ ] 12.1.1 Every operation accepts `Structure`
  Description: every function that takes CIF text also accepts a
  `Structure`, and every structure-producing function returns a `Structure`
  (with `to_dict()` for the old shape). Methods mirror the free functions:
  `s.supercell(2, 2, 2).slab((1, 1, 1), layers=2, vacuum=10).write("s.cif")`.
- [ ] 12.1.2 Typed result classes
  Description: replace `dict`/`Any` results with frozen classes
  (`SymmetryDataset`, `PxrdPattern`, `OptimizeResult`, `PorosityReport`, …)
  generated from the serde types, with `to_dict()`, `__repr__` and rich
  `_repr_html_` for notebooks.
- [ ] 12.1.3 `mypy --strict` without `Any`
  Description: the `.pyi` stubs contain no `Any` in public signatures; a test
  greps for it.
- [ ] 12.1.4 Context managers and progress
  Description: long operations accept `progress=` callbacks and honour
  `KeyboardInterrupt` through the cancellation token; `tqdm` integration when
  installed.

## Task 12.2 — NumPy

- [ ] 12.2.1 Arrays out
  Description: `positions`, `frac_positions`, `cell`, `numbers`, `forces`,
  `stress`, PXRD `two_theta`/`intensity`, grids and meshes as `numpy.ndarray`
  (copy-free where layout allows). Feature `numpy` on by default for wheels.
- [ ] 12.2.2 Arrays in
  Description: `Structure.from_arrays(cell, numbers, frac=… | cart=…)`,
  `set_positions(array)`; shape and dtype validated with `InputError`.
- [ ] 12.2.3 Batch APIs
  Description: `xtalforge.batch.read(paths)`, `batch.detect_symmetry`,
  `batch.pxrd` releasing the GIL and running in parallel, returning lists in
  input order (rule 6).

## Task 12.3 — JavaScript/TypeScript object API

- [ ] 12.3.1 Structure handle
  Description: a `Structure` class owning the parsed Rust value
  (`Structure.fromText(text, hint)`), with methods for every operation, so a
  structure is parsed once; explicit `free()` plus `FinalizationRegistry`.
- [ ] 12.3.2 Typed results
  Description: generate TypeScript interfaces from the serde types
  (`tsify` or a schema-driven generator); no exported function returns `any`.
  A test fails on `any` in `xtalforge_wasm.d.ts`.
- [ ] 12.3.3 Typed arrays
  Description: positions, scene buffers, patterns and meshes as
  `Float64Array`/`Float32Array`/`Uint32Array` views for zero-copy rendering.
- [ ] 12.3.4 Web worker and Node entry points
  Description: an optional worker wrapper (`xtalforge-wasm/worker`) that runs
  heavy operations off the main thread with progress and cancellation; a Node
  entry with file helpers.
- [ ] 12.3.5 Optional threads
  Description: a second WASM build with `wasm-bindgen-rayon` for browsers with
  cross-origin isolation; same results as the single-threaded build.

## Task 12.4 — Surfaces and docs

- [ ] 12.4.1 Parity matrix, API pages, tutorials and the playground move to
  the object API; benchmarks show the parse-once gain.
