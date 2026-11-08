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

The operation loop is the setting when present; the Hall symbol is stored for
the space-group identity; an H-M symbol alone cannot fix the operations, so the
reader records a `symmetry_warning` and the P1 expansion uses the identity.
This is the first pass of the Crystal-Builder priority rule; inconsistency
warnings and the full Hall lookup arrive with moyo in Task 2.1.

## Numbers

`1.234(5)` is read as `1.234`; `?` and `.` become missing values. Malformed
numbers or truncated loops return `Error::InvalidInput`, never a panic. A
fuzz-style regression test feeds empty, truncated and mismatched inputs.

## Writer

`write_cif` writes the cell, Hall symbol, operation loop, sites, occupancy,
Uiso and charge. The round-trip test reads the written P1 structure back.
