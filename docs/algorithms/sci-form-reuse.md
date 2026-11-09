# sci-form reuse audit and port plan

Source: `sci-form` commit `6b78c7aa427519cb7f67d26c781a3f8e01117552`,
crate version `1.0.0`, MIT, MSRV 1.89. The line counts below were checked
against the source on 2026-11-08; the roadmap table is still accurate.

## Inventory

| Source | Lines | Public surface checked | Tests | Plan |
|---|---|---|---|---|
| `materials/cif.rs` | 708 | reader, loop/table parser, writer | 9 | Port into `xtalforge-io::cif`; split into tokenizer, loop parser, item resolver and writer. |
| `materials/cell.rs` | 1048 | lattice, supercell, PXRD, symmetry expansion, porosity helpers | 14 | Split across `core`, `analysis` and `symmetry`; do not port as one file. |
| `materials/space_group_table.rs` | 6750 | static tables and lookups | 2 | Do not port; replaced by moyo. Keep tests as extra checks. |
| `materials/space_groups.rs` | 273 | operation helpers | 6 | Do not port; replaced by moyo. |
| `materials/scattering_factors.rs` | 396 | Cromer-Mann coefficients | 3 | Port into `xtalforge-analysis::pxrd`; coefficients to `data/`. |
| `materials/assembly.rs` | 829 | SBU assembly and topology helpers | 19 | Study for the MOF builder; the PORMAKE path in Crystal-Builder is the reference. |
| `materials/sbu.rs` | 257 | SBU records | 2 | Study; data to `data/` if used. |
| `materials/topology_data.rs` | 734 | topology tables | 1 | Study; data to `data/` if used. |
| `materials/geometry_opt.rs` | 661 | geometry optimizer | 7 | Port into `xtalforge-ff::optimize`; split. |
| `materials/variable_cell.rs` | 501 | variable-cell relaxation | 5 | Port into `xtalforge-ff::optimize`; split. |
| `forcefield/uff.rs` | 513 | UFF terms and calculator | 8 | Compare with Crystal-Builder `ff/uff/`; port the better-tested parts. |
| `forcefield/params.rs` | 391 | UFF parameters | 3 | Parameters to `data/`; loader tested. |
| `forcefield/atom_typer.rs` | 247 | UFF atom typing | 4 | Compare with Crystal-Builder typer; port the reference behaviour. |
| `forcefield/minimizer.rs` | 222 | L-BFGS wrapper | 3 | Port into `xtalforge-ff::optimize` with cancellation and progress. |
| `charges_eeq/*` | 2928 | EEQ/EEQBC charges and data | 12 | Port EEQ as an extra charge method; `eeqbc_data.rs` becomes a data file. |
| `smiles.rs` | 875 | SMILES parser | 11 | Port into `xtalforge-build::smiles`; split. |
| `graph.rs` | 627 | molecular graph | 9 | Port into `xtalforge-build::smiles`; split. |
| `stereo.rs` | 1121 | stereochemistry | 10 | Port into `xtalforge-build::smiles`; split. |
| `rings/sssr.rs` | 568 | SSSR ring perception | 7 | Port into `xtalforge-core::rings`; make periodic-aware and split. |
| `rings/aromaticity.rs` | 219 | aromaticity | 4 | Port into `xtalforge-core::rings`. |
| `periodic.rs` | 399 | element data | 5 | Port into `xtalforge-core::elements`; data to `data/`. |
| `surface/sasa.rs` | 335 | accessible surface | 4 | Port into `xtalforge-analysis` as an optional helper. |
| `kmesh.rs` | 1094 | k-point meshes | 17 | Study for phonon DOS (Task 4.8.3). |
| `fd_gradient.rs` | 86 | finite-difference gradient check | 2 | Port as a test helper. |
| `xtb/`, `eht/`, `pm3/`, `conformer.rs`, `etkdg.rs` | large | quantum and embedding engines | many | Do not port; use through the optional `sci-form` feature. |

## Quality gate for every port (Task 0.4.2)

1. Read the complete source file and list its public functions.
2. Remove `unwrap`, `expect` and panics on input; return typed errors.
3. Adapt to our types, units and error codes.
4. Translate any non-English text (none found in the first scan on
   2026-11-08; accented author names such as Rappé are allowed).
5. Split under 500 lines with a written split plan.
6. Keep or recreate the tests; add the origin comment and the
   `THIRD_PARTY.md` entry.

## Overlap with Crystal-Builder (Task 0.4.3)

| Area | Crystal-Builder | sci-form | Reference for golden data |
|---|---|---|---|
| CIF | gemmi-based reader/writer | own parser | Crystal-Builder behaviour; sci-form code base |
| PXRD | `analysis/pxrd.py` and `materials/scattering_factors.rs` | `materials/cell.rs` helpers | Crystal-Builder output |
| UFF | `ff/uff/` | `forcefield/uff.rs` | Crystal-Builder unless sci-form is proven more accurate on the corpus |
| Porosity | grid and Zeo++ | cell helpers | Crystal-Builder |
| Rings | periodic `core/rings.py` | SSSR in `rings/` | Crystal-Builder periodic rings; sci-form algorithm adapted |
| MOF assembly | vendored PORMAKE | `materials/assembly.rs` | PORMAKE in Crystal-Builder |

Crystal-Builder is the default reference because xtalforge must reproduce it.
A sci-form implementation replaces it only where a corpus comparison proves it
more accurate, and the choice is recorded in the algorithm page.

## `sci-form` feature policy (Task 0.4.4)

- Dependency: `sci-form = "1"` from crates.io, optional, never a path
  dependency, so published crates build.
- Scope: xTB, EHT and PM3 quantum methods only. Crystal-related code is ported
  into our crates, not pulled through this dependency.
- Boundary: sci-form types are converted to ours in one private module; no
  sci-form type appears in the public API.
- Surfaces: the feature is available in the Rust API and CLI; Python wheels may
  include it as an opt-in build; it is disabled in the default WASM build for
  size and because the engines need threads and external data.
- Default: off everywhere; the default build stays pure Rust.
- Validation: the `sci-form` feature links the published crate and the quantum
  tests check finite water energies, conservative GFN2 forces and the refusal
  paths. A cross-check against a real `xtb` binary is not run here because the
  binary is not installed; the external-engine tests elsewhere use simulated
  executables.
