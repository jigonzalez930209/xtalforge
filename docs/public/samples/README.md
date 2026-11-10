# Playground samples

These CIF files are byte-identical copies of the curated corpus:

| File | Corpus source | License |
|---|---|---|
| `Si.cif` | `tests/corpus/simple/Si.cif` | CC0 (COD) |
| `NaCl.cif` | `tests/corpus/simple/NaCl.cif` | CC0 (COD) |
| `diamond.cif` | `tests/corpus/simple/diamond.cif` | CC0 (COD) |
| `graphene.cif` | `tests/corpus/simple/graphene.cif` | CC0 (COD) |
| `TiO2.cif` | `tests/corpus/simple/TiO2.cif` | CC0 (COD) |
| `MFI.cif` | `tests/corpus/simple/MFI.cif` | CC0 (COD) |
| `MOF-5.cif` | `tests/corpus/prepared/MOF-5.cif` | project model derived from CC0 data |
| `Fe.cif` | `tests/corpus/simple/Fe.cif` | CC0 (COD) |

The full provenance and citations are in `tests/corpus/PROVENANCE.md`.
`tests/e2e/reference/generate.py --check-samples` fails when a copy drifts
from its corpus source.
