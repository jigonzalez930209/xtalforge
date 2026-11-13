# Provenance

`xtalforge_core::Provenance` records the xtalforge version, operation, options,
input SHA-256 digests, data bundle versions, seed and an optional timestamp.
`with_input_bytes` hashes the input; `attach` stores the record under the
`provenance` structure metadata key; `from_structure` reads it back.

`xtalforge_io::provenance(text)` extracts a record from CIF tags/comments,
extxyz `xtalforge_provenance=` comment pairs, project files or a structure
JSON document. The writers embed it automatically when the structure carries
the metadata.

```bash
xtalforge provenance result.cif
```

The command prints JSON with `--json` or a short version/digest summary.

Parity key: `provenance`.

The schema constants are `PROVENANCE_SCHEMA_VERSION` and
`PROVENANCE_METADATA_KEY`; `ProvenanceInput` carries an optional name and a
SHA-256 digest, and `sha256_hex` hashes bytes for callers that need the same
digest without building a record.
