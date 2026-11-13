# Embedded symmetry data

`kpaths.json.gz` holds the special points and band-path segments generated
from ASE 3.29 for the 14 Bravais types. The published crate decompresses it,
verifies its SHA-256 checksum and parses it on first use (Task 15.3.2):

```bash
gzip -dc kpaths.json.gz | sha256sum
```

The manifest constant lives in `../src/kpath.rs`. Regenerate the table from
ASE when the reference changes; the JSON schema is `{ "lattices": { "<label>":
{ "special_points": {...}, "path": [[label, label], ...] } } }`.
