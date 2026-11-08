# Projects and workspaces

## `.xtalproj`

A project is a ZIP of small readable JSON parts plus the structure as CIF. The
format version is independent of the crate version (`FORMAT_VERSION = 1`).

| Part | Contents |
|---|---|
| `structure.cif` | The structure, always present |
| `format.json` | `{"format_version": 1}` |
| `cell.json` | Bare cell when there are no atoms |
| `bonds.json` | Perceived graph, hand-drawn and suppressed bonds |
| `sites.json` | Per-site properties a CIF cannot hold |
| `view.json` | View settings (Task 3.14) |
| `session.json` | Selection, measurements, provenance |

```rust
use xtalforge::io::project::{ProjectParts, read_project, write_project};

let parts = ProjectParts::from_structure(&structure, None, None, None, None, None)?;
let bytes = write_project(&parts)?;
let back = read_project(&bytes)?;
```

The writer goes through `zip` with deflate and reads every older version; a
newer `format_version` is refused rather than misread. In WASM the project is
a `Uint8Array` (`project_pack` / `project_info`); the CLI has
`project {pack,info,extract}`.

## Workspace

`xtalforge-modules::Workspace` stores one directory per run:
`runs/<id>/run.json` with `RunRecord { id, module, input, status, outputs,
data, created }`. The browser has no disk, so `archive_records` serializes the
records into a ZIP and `read_archived_records` reads it back; Python and the
CLI use the on-disk workspace, WASM uses the in-memory archive
(`workspace_archive` / `workspace_read_archive`).

## Typed parameters

`xtalforge-core::preset("symmetry" | "pxrd" | "uff")` is the single definition
of a parameter set. The CLI (`params <name>`), Python (`param_schema`) and WASM
(`param_schema`) print the same schema, including defaults, bounds, units and
help text. `ParamSet::to_json`/`from_json` persist a set with revalidation.
