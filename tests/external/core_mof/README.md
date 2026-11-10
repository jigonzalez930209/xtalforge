# CoRE MOF 2019 subset

Eighteen structures copied from the public CoRE MOF 2019 v2 archive
(`CoREMOF2019_public_v2_20241119.zip`, Zenodo record 14184621), directory
`CR/ASR/` (charged/relaxed, all solvent removed) and `CR/FSR/` where noted.
They exercise CIF parsing, P1 expansion, cell volume, density, formula,
bonding and porosity on real framework data with DDEC6/PACMAN charges.

## Licence and citation

The archive is **CC BY 4.0** per the Zenodo record. Attribution:

> Chung, Y. G.; Haldoupis, E.; Bucior, B. J.; Haranczyk, M.; Lee, S.;
> Zhang, H.; Vogiatzis, K. D.; Milisavljevic, M.; Ling, S.; Camp, J. S.;
> Slater, B.; Siepmann, J. I.; Sholl, D. S.; Snurr, R. Q. "Advances,
> Updates, and Analytics for the Computation-Ready, Experimental Metal-Organic
> Framework Database: CoRE MOF 2019." *J. Chem. Eng. Data* **2019**, 64,
> 5985-5998. https://doi.org/10.1021/acs.jced.9b00835
>
> CoRE MOF 2019 public v2 archive. Zenodo.
> https://doi.org/10.5281/zenodo.14184621

The files are unmodified copies; the original `.cif` names are kept.

## Files

| File | Metal(s) | Sites | Density (g/cm3) | Why it is here |
|---|---|---|---|---|
| `FIQCEN_clean_pacman.cif` | Cu | 156 | 0.879 | HKUST-1, Cu paddlewheel, charged |
| `OFERUN_clean_pacman.cif` | Zn | 138 | 0.925 | ZIF-8, tetrahedral Zn, imidazolate |
| `RUBTAK03_auto_pacman.cif` | Zr | 456 | 1.238 | UiO-66 family, Zr6 cluster |
| `BEKYOG_clean_pacman.cif` | Cr, Ag | 82 | 2.089 | mixed-metal, heavy scatterer |
| `ELOFIU_charged_pacman.cif` | Ca | 78 | 1.186 | charged framework, light metal |
| `IRISAD01_charged_pacman.cif` | Bi | 36 | 2.123 | heavy metal, small cell |
| `ACOGAB_clean_pacman.cif` | Cd | 174 | 1.155 | Cd, Br termination |
| `ABAVIJ_clean_pacman.cif` | Co | 108 | 1.525 | Co, N-donor linker |
| `AHOKAJ_clean_pacman.cif` | Fe | 36 | 2.089 | Fe phosphonate, no hydrogen |
| `AVIPAX_clean_pacman.cif` | Mg | 192 | 0.964 | alkaline earth |
| `AGUWUV_clean_pacman.cif` | Ni | 51 | 1.512 | Ni, small cell |
| `HAPBEI_clean_pacman.cif` | Ti, Cu | 444 | 1.183 | mixed-metal, iodine |
| `ALIPEQ_clean_pacman.cif` | Zr, Cd | 78 | 1.685 | mixed-metal cluster |
| `ABEXEM_clean_pacman.cif` | La | 144 | 1.761 | lanthanide |
| `CELZEA_clean_pacman.cif` | Mo | 52 | 1.630 | Mo, phosphine |
| `HUCBAJ_charged_pacman.cif` | Li, In | 212 | 1.200 | alkali + In, charged |
| `ADATOQ_clean_pacman.cif` | Ag | 692 | 0.632 | large cell, low density |
| `FALZOI_clean_pacman.cif` | Na, Al | 138 | 0.961 | aluminophosphate |

Densities are the values this library computes from the CIF cell and sites;
they are not copied from the archive.

## Tests

- `crates/xtalforge-io/tests/external_core_mof.rs`: every file parses, the
  summary is finite, the P1 expansion and bond graph succeed, and a CIF
  round-trip preserves the site count and volume. Four frameworks are checked
  against the file's own `_chemical_formula_sum` and `_cell_volume`.
- `crates/xtalforge-analysis/tests/external_core_mof.rs`: calculated PXRD and
  distance-grid porosity on HKUST-1 and ZIF-8 are finite and physically
  bounded.

Nothing else from the archive is redistributed. The full archive is 81 MB and
must be downloaded by the user.
