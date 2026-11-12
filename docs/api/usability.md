# Usability helpers (Phase 11)

These helpers close the audit findings without changing the 1.x API.

## One-call IO

- `read(path)` reads a structure from disk, detects the format from the
  extension or content and returns a `ReadResult { structure, warnings }`.
- `read_str(text, hint)` is the text variant; `hint` may be empty for
  content sniffing.
- `write(path, &structure)` chooses the writer from the path extension and
  uses the atomic writer.
- Python exposes `read(path)` and `write(path, text, format_hint="")`; WASM
  exposes `read(text, hint)` (no filesystem).

## Threads and errors

- `set_threads(threads)` installs the process-wide Rayon pool; the Rust
  entry point is `set_global_threads(threads)`. The CLI
  `--threads` uses a scoped `ParallelOptions` pool, the Rust
  `ParallelOptions::install` does the same for a single call. WASM records
  the browser-thread exception.
- `error_category(code)` maps a stable code to `input`, `unsupported`,
  `numerical`, `limit`, `engine` or `internal`. Every crate error also exposes
  `Error::category()`, and the CLI/WASM error objects carry `category`.
- `StructureWarning` is the structured warning type (`code`, `message`,
  `atoms`, `fields`, `hint`); `Structure::warnings()` collects it from the
  reader metadata and every `StructureSummary` includes it.
- `deprecation_warning(name, replacement)` emits a `DeprecationWarning` in
  Python and a `console.warn` in WASM; the Rust/CLI helper is
  `xtalforge_core::warn_deprecated`. See `docs/guide/deprecation.md`.

## Molecules, sketches, fragments, guests

- `molecule_from_smiles(smiles, seed=7, iterations=1500)` parses and embeds a
  molecule in 3D. CLI `xtalforge molecule --smiles …`.
- `sketch_roundtrip(smiles)` parses a sketch and writes it back as the
  canonical depth-first SMILES. CLI `xtalforge sketch`.
- `fragment_library()` lists the embedded fragments, `fragment_names()` is
  the name list, and `fragment_structure(name)` embeds one. CLI
  `xtalforge fragment [NAME]`.
- `substitute_structure(text, index, element)` replaces one asymmetric site
  after the clearance rules. CLI `xtalforge substitute`.
- `place_guest_structure(host, guest, offset, clearance)` places a guest
  molecule after a periodic clearance check. CLI `xtalforge guest`.
- WASM names are the same with `_info` for object-returning results
  (`fragment_info`, `substitute_info`, `place_guest_info`).

The option spellings on every surface are in [`docs/api/naming.md`](./naming.md).
