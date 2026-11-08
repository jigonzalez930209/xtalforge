# Structure schema

The serialized structure is wrapped in a versioned envelope so future changes
can migrate instead of silently misreading old files.

```json
{
  "schema_version": 1,
  "structure": {
    "name": "NaCl",
    "lattice": [[5.64, 0.0, 0.0], [0.0, 5.64, 0.0], [0.0, 0.0, 5.64]],
    "hall": "P 1",
    "sites": [],
    "metadata": {}
  }
}
```

## API

```rust
use xtalforge::core::{StructureDocument, SCHEMA_VERSION};

let document = StructureDocument::from_json(text)?;
assert_eq!(document.schema_version, SCHEMA_VERSION);
let text = document.to_json()?;
```

`StructureDocument::new(structure)` wraps a structure in the current version;
`Structure::to_document()` does the same. `from_json` rejects unknown top-level
fields and any version other than `SCHEMA_VERSION`. The site `props` map and
the structure `metadata` map stay open for extensions.

## Compatibility

A fixture for every schema version lives in `tests/schema/`; the current
fixture is `structure-v1.json`. Changing a serialized type bumps
`SCHEMA_VERSION` and adds migration code plus a new fixture. The crate version
is independent (AGENTS.md section 5.2).

The JSON Schema generation with `schemars` is part of the release tooling; the
Rust fixture tests are the compatibility gate until it lands.
