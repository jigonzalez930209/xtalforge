# QMOF subset

Fourteen structures from the QMOF Database (Figshare 13147324, v18), `relaxed_structures/` subset. They are unmodified DFT-relaxed CIFs (PBE-D3(BJ), VASP) with the original `qmof-*.cif` names.

## Licence and citation

The Figshare record declares **CC BY 4.0**. Cite:

> A. S. Rosen, S. M. Iyer, D. Ray, Z. Yao, A. Aspuru-Guzik, L. Gagliardi, J. M. Notestein, R. Q. Snurr. "Machine Learning the Quantum-Chemical Properties of Metal-Organic Frameworks for Accelerated Materials Discovery." *Matter* **2021**, 4, 1578-1596. https://doi.org/10.1016/j.matt.2021.02.015
>
> QMOF Database, Figshare. https://doi.org/10.6084/m9.figshare.13147324.v18

The full archive is 374 MB; only these 14 CIFs are redistributed.

## Files

| File | Metal(s) | P1 atoms | Volume (A^3) | Density (g/cm3) | Space group |
|---|---|---|---|---|---|
| `qmof-ea5d153.cif` | Sc | 26 | 294.295 | 2.031 | 167 |
| `qmof-9247cfe.cif` | Ti | 58 | 635.372 | 2.323 | 2 |
| `qmof-c3fa563.cif` | V | 24 | 236.591 | 1.854 | 23 |
| `qmof-7a9d846.cif` | Cr | 29 | 270.856 | 2.183 | 2 |
| `qmof-832a1ad.cif` | Mn | 17 | 220.982 | 1.888 | 12 |
| `qmof-b6b9d01.cif` | Fe/Pt | 20 | 390.466 | 2.003 | 10 |
| `qmof-89d2f92.cif` | Co | 23 | 326.182 | 2.214 | 12 |
| `qmof-421f50d.cif` | Ni | 26 | 283.396 | 2.142 | 15 |
| `qmof-40f119f.cif` | Cu | 17 | 210.892 | 2.046 | 12 |
| `qmof-36ade61.cif` | Zn | 26 | 330.782 | 3.065 | 51 |
| `qmof-7aa58f3.cif` | Zr/Ca/K | 52 | 573.296 | 2.044 | 82 |
| `qmof-e1456f8.cif` | Mo | 28 | 338.006 | 2.182 | 5 |
| `qmof-0331d6b.cif` | Ru | 52 | 634.438 | 2.184 | 2 |
| `qmof-5d3a22f.cif` | Ag | 20 | 193.351 | 3.091 | 2 |

Volumes, densities and atom counts are the QMOF tabulated values used as external checks; the formula column is the QMOF cell formula.

## Tests

`crates/xtalforge-io/tests/external_qmof.rs` checks the parse against these values and a CIF round-trip.
