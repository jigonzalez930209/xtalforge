# CIF parsing

## Pipeline

1. `decode_text` strips a UTF-8 BOM, reads UTF-8 and falls back to Latin-1.
2. The tokenizer walks lines, removes comments outside quotes, keeps quoted
   values together and collects semicolon-delimited text blocks.
3. The parser builds data blocks, single items and loops. Loop values must be
   a multiple of the tag count.
4. `structure_from_block` reads the cell, Hall/H-M symbol, symmetry-operation
   loop, atom sites, occupancy, Uiso/Uij and charge, and builds a `Structure`.

## Symmetry priority

The operation loop is the setting when present. Without a loop, the reader
resolves the Hall symbol (or H-M symbol) through
`xtalforge_symmetry::SpaceGroup::parse` and installs its operations, recording
`symmetry_source = "Hall symbol"`. Only an unrecognised symbol falls back to
the identity and records a `symmetry_warning`. This closes the common CIF case
where the space group is named but the operation loop is omitted (graphene is
the corpus example). Inconsistency warnings beyond that are still limited.

## Numbers

`1.234(5)` is read as `1.234`; `?` and `.` become missing values. Malformed
numbers or truncated loops return `Error::InvalidInput`, never a panic. A
fuzz-style regression test feeds empty, truncated and mismatched inputs.

## Writer

`write_cif` writes the cell, Hall symbol, operation loop, sites, occupancy,
Uiso and charge. The round-trip test reads the written P1 structure back.
