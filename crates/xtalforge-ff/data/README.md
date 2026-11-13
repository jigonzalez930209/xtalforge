# Embedded force-field data

`uff.json.gz` holds the 218 UFF/UFF4MOF rows ported from Crystal-Builder
`xtal/ff/uff/params.py` (Rappé 1992, Addicoat 2014, Coupry 2016). The
published crate decompresses it on first use, after verifying its SHA-256
checksum (Task 15.3.2):

```bash
gzip -dc uff.json.gz | sha256sum
```

The manifest constant lives in `../src/uff/params.rs`; the provenance entry is
in `THIRD_PARTY.md` at the repository root.

`eeq2019.json` stays embedded uncompressed: it is 5 KiB, so gzip would save
about 3 KiB and the table is small enough to parse directly on first use.
