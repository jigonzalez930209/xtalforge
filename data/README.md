# Data directory

Canonical tables, parameter sets and databases for xtalforge. The policy is in
[`docs/guide/data.md`](../docs/guide/data.md).

- Use `.json` for tables that should stay readable in review, `.json.gz` for
  large tables, and `.bin` for compact numeric arrays with a documented layout.
- Compress with `gzip -9 -n` so the file is reproducible.
- Record the origin and license in `THIRD_PARTY.md` before adding a file.
- Embedded bundles are declared with `EmbeddedBundle` and a SHA-256 digest of
  the decompressed payload; external bundles are loaded with
  `load_external_gzip` or `load_external_json`.
- Data embedded into a published crate lives in that crate's `data/` directory
  so `cargo package` includes it. This directory holds canonical files,
  external bundles and repository-only assets.

No bundle is present yet; Phase 1 adds the element table (Task 1.1.1).
