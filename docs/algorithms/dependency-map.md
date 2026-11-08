# Dependency replacement map

Task 0.10: every Python dependency Crystal-Builder uses for the headless core
must have a pure Rust replacement or a documented decision before its phase
starts. The survey below was run on Crystal-Builder commit
`78d5476fd9e851fda00cde57b7ab8da3195dcf54`.

## Survey (Task 0.10.1)

| Dependency | Modules that import it | Used for |
|---|---|---|
| `gemmi` | `core/elements.py`, `core/spacegroup.py`, `io/cif_reader.py`, `io/shelx.py`, `analysis/pxrd.py`, `analysis/rcsr.py`, `analysis/netsearch.py`, `powder/bridge.py` | CIF and SHELX parsing, space-group tables, scattering factors, net search |
| `spglib` | `core/symmetry.py`, `core/supercell.py`, `core/subgroups.py` | symmetry detection, standard cells, subgroup reference |
| `numpy` | most of `xtal/` | arrays and vectorised geometry |
| `scipy` | `core/neighbors.py`, `analysis/{grid,kpath,voids}.py`, `build/{clearance,fill}.py`, `carbon/*`, `fullerene/geometry.py`, `ff/ewald.py`, `polymer/*`, `mof/pormake/*` | k-d trees, Voronoi, sparse grids, L-BFGS, special functions |
| `ase` | `mof/pormake/{building_block,framework,local_structure,neighbor_list,topology,utils}.py` | `Atoms`, neighbour lists and I/O inside the vendored PORMAKE |
| `phonopy` | `ff/phonons.py` (lazy import; `find_spec` probe) | symmetry-reduced displacements, phonon mesh and thermodynamics |
| `rdkit` | `build/chem.py` (lazy import in one module) | SMILES, embedding and sketcher chemistry |
| `rietx` | `powder/bridge.py` | peak fitting, indexing, Pawley and Rietveld |
| `matplotlib` | `view/colormaps.py` (sampled tables, no import at draw time) | colormap samples |

## Decisions (Task 0.10.2)

| Dependency | Rust plan | Phase |
|---|---|---|
| gemmi | Own CIF parser (sci-form port); moyo tables for space groups; own scattering-factor data | 1, 2, 5 |
| spglib | moyo for detection and standardization; own subgroup enumeration only if moyo's API is insufficient | 2 |
| numpy | `nalgebra` 0.35 and plain `Vec<[f64; 3]>`; `ndarray` only if a grid operation needs it | 1+ |
| scipy `cKDTree` | Own periodic cell list; `kiddo` k-d tree only if benchmarks require it | 1 |
| scipy `optimize` | Own L-BFGS/FIRE (sci-form port) and Levenberg-Marquardt for powder | 4, 7 |
| scipy `spatial.Voronoi` | Own Wigner-Seitz construction by half-space intersection | 2 |
| scipy `sparse`, `ndimage` | Own periodic flood fill and union-find | 5 |
| ASE | Own k-path (Task 2.10), own readers, ML calculators behind the `ml` feature | 2, 3, 8 |
| phonopy | Own displacements, force constants, phonon DOS and thermodynamics, validated against phonopy outputs | 4 |
| RDKit | sci-form SMILES port; embedding through the optional `sci-form` feature | 6 |
| rietx | Study, then port the numerical kernels natively after the license check | 7 |
| matplotlib | Colormap tables as data; plots are left to consumers | 3 |

No Python package becomes a runtime dependency of the Rust library. The Python
bindings expose the same Rust engine and need NumPy only for arrays.

## License check (Task 0.10.3)

| Package | License | Use in xtalforge | Status |
|---|---|---|---|
| gemmi | MPL-2.0 | Golden-data generator only; no code ported | verified, not redistributed |
| spglib | BSD-3-Clause | Golden-data generator and audit comparison only | verified, not redistributed |
| numpy | BSD-3-Clause | Golden-data generator only | verified, not redistributed |
| scipy | BSD-3-Clause | Golden-data generator only; algorithms reimplemented from papers | verified, not redistributed |
| ASE | LGPL-2.1-or-later | Reference for k-paths and C60 coordinates; the tables are redistributed as derived data, no code ported | verified, derived data |
| phonopy | BSD-3-Clause | Reference outputs only; the Gamma-point phonon path is native | verified, not redistributed |
| RDKit | BSD-3-Clause | Replaced by sci-form's MIT code | verified, not ported |
| rietx | MIT (stated in Crystal-Builder's `pyproject.toml`; no license field on PyPI) | Inventory only; no code or data copied | verified from the note, not redistributed |
| matplotlib | Matplotlib license (PSF-based and BSD-compatible) | Colormap anchor colours written by hand | verified, not redistributed |

Ported code and data are recorded in `THIRD_PARTY.md` before they are merged.
Where a license is unclear, the algorithm is reimplemented from the paper
instead of copied.
