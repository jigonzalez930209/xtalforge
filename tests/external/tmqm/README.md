# tmQM subset

Eighteen transition-metal complexes copied unmodified from the tmQM 2024
release (`tmQM_X1.xyz.gz`, GitHub `uiocompcat/tmQM`). One complex per metal
across the 3d, 4d and 5d series, chosen with at most 70 atoms per record. Each
file keeps the original XYZ header: CSD code, total charge, multiplicity,
stoichiometry, metal coordination number (MND) and the CSD release.

## Licence and citation

The tmQM repository declares the **MIT** licence (`LICENSE`, Copyright (c)
2024 Bastian Bjerkem Skjelstad). Cite:

> B. B. Skjelstad, A. M. Ravn, U. Norby, A. M. A. A. Alagbe, A. N.
> Kristiansen, J. M. A. H. van der Ven, M. F. A. M. van der Ven, and
> D. Balcells. "tmQM 2024: 108k Transition Metal Complexes" (2024 release).
> The original dataset: B. Meyer, B. Sawatlon, S. Heinen, O. A. von Lilienfeld,
> C. Corminboeuf. "Machine learning meets quantum chemistry: the tmQM
> dataset." *J. Chem. Inf. Model.* **2021**, 61, 2585-2596.
> https://doi.org/10.1021/acs.jcim.0c01041

The MIT licence covers the repository; the structures are experimental
geometries from the CSD, relaxed with GFN2-xTB. They are redistributed here
with attribution under that licence.

## Files

| File | Metal | Atoms | Charge | MND |
|---|---|---|---|---|
| `ANAXOE_tmQM.xyz` | Nb | 31 | 0 | 6 |
| `CEYJEZ_tmQM.xyz` | W | 67 | 0 | 6 |
| `CILXOO_tmQM.xyz` | Mo | 39 | 0 | 6 |
| `CPMINC_tmQM.xyz` | Cr | 50 | 0 | 6 |
| `DAQHAI_tmQM.xyz` | Fe | 59 | 0 | 5 |
| `DASFIQ_tmQM.xyz` | Sc | 45 | -1 | 8 |
| `HEWTEM_tmQM.xyz` | Rh | 62 | 0 | 6 |
| `IPOSAK_tmQM.xyz` | Mn | 40 | 0 | 6 |
| `LAXLAB_tmQM.xyz` | Ir | 55 | 0 | 8 |
| `PITRAM_tmQM.xyz` | Zr | 63 | 0 | 12 |
| `QOSJIR_tmQM.xyz` | Ru | 58 | 0 | 9 |
| `REJFUL_tmQM.xyz` | V | 43 | -1 | 4 |
| `UMAHUO_tmQM.xyz` | Hf | 69 | 0 | 6 |
| `UTANIP_tmQM.xyz` | Re | 45 | +1 | 6 |
| `VUDZEC_tmQM.xyz` | Y | 56 | 0 | 7 |
| `WUGVEZ_tmQM.xyz` | Ta | 59 | 0 | 5 |
| `XARWUJ_tmQM.xyz` | Os | 31 | 0 | 6 |
| `ZOYPEL_tmQM.xyz` | Ti | 66 | 0 | 8 |

## Tests

- `crates/xtalforge-io/tests/external_tmqm.rs`: XYZ atom counts and element
  counts against the header stoichiometry, metal coordination from the
  perceived bonds, and an XYZ round-trip.
- `crates/xtalforge-ff/tests/external_tmqm.rs`: UFF/UFF4MOF builds for every
  metal and returns a finite energy and forces.
