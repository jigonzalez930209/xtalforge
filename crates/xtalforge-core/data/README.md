# Element data

`elements.json.gz` holds the periodic table, Pyykko radii and isotope masses
used by `xtalforge-core::elements` (Task 1.1).

- Source: Crystal-Builder `core/elements.py`, `core/pyykko.py`,
  `ff/isotopes.py` and `ff/vibrations.py`.
- Generator: `scripts/generate-elements-data.py` in a virtual environment
  with Crystal-Builder installed.
- Checksum: `elements.rs` records the SHA-256 of the decompressed payload and
  `EmbeddedBundle` verifies it on first use.
- Compression: gzip with a fixed zero timestamp, so the file is reproducible.

The table contains H..Og (Z = 1..118), deuterium as its own species and the
dummy connection point `X`. Masses and radii are copied from gemmi through
Crystal-Builder; the Pyykko radii are in Angstrom after dividing the published
picometres by 100.
