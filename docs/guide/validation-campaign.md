# Validation campaign

This page records the **2026-11-17/18 full validation campaign**: what was
executed, on which machine, against which independent implementations and data
sets, with the measured results and the gaps that remain. It complements the
[End-to-end audit](./audit) (browser/WASM battery) and the
[external validation](../algorithms/external-validation) inventory.

## Environment

| Item | Value |
|---|---|
| Host | Intel Xeon E5-2680 v4, 28 threads, 14.99 GiB RAM |
| OS | Linux 7.0.0-43-generic x86_64 |
| Rust | 1.97.1 (edition 2024, MSRV 1.89) |
| Python | 3.14.4 |
| xtalforge | 1.14.1 working tree (the fixes listed below are not released yet) |
| gemmi | 0.7.5 |
| spglib | 2.8.0 |
| ASE | 3.29.0 |
| pymatgen | 2026.10.24 |
| mendeleev | 1.3.0 |
| NumPy / SciPy | 2.5.5 / 1.18.1 |

Parallel work was run with `xargs -P 28` where the driver is per file. The
saturation test reached **2731 % CPU** (about 27.3 of 28 threads) with a
maximum resident set of 23 MB per process; the full QMOF run peaked at about
**8 GiB / 14 GiB (57 %)**. Memory stayed below the 80 % ceiling throughout.

## Surface suites

Run on the same working tree:

| Surface | Command | Result |
|---|---|---|
| Rust | `cargo test --workspace --all-features` | **1,456 passed, 0 failed** (303 test targets) |
| Browser/WASM | `scripts/e2e.sh --skip-build` | **68 passed, 0 failed** |
| Python | `pytest crates/xtalforge-python/tests` | **110 passed, 4 skipped** |
| WASM | `wasm-pack test --node --release` | **58 passed, 0 failed** (14 suites) |
| CLI/Python/WASM | `scripts/check-surfaces.sh` | **97 fields, 0 differences** |
| Split WASM | `scripts/build-wasm-packages.sh core symmetry analysis build full` | all five load and exercise their tier |
| WASM size | `scripts/check-wasm-size.sh` | full **24,104 KiB** (budget 25,600); core 2,051; symmetry 2,930; analysis 3,229; build 2,469 |
| Dependency audit | `cargo deny check` | advisories, bans, licences and sources OK |
| Benchmarks | `scripts/check-bench.sh --update` then `--reps 1` | baseline regenerated on this CPU; all seven rows within 5 % and the gate passes |
| Docs | `npm run docs:build` | clean |
| Static checks | clippy `-D warnings`, rustdoc `-D warnings`, fmt, rules, mypy strict, stubs | clean |

## Repository corpus

The 72-CIF corpus in `tests/corpus/` (simple, edge, COD, prepared and polymer
cases) was compared file by file with the CLI as the library entry point.

| Check | Reference | Scope | Result |
|---|---|---|---|
| Cell, volume, density, atom count | gemmi | 72 CIF | volume exact; 69 densities exact; MIL-100/101 differ by `1e-7` (atomic-weight rounding); UiO-66 clamps occupancies above one |
| Space group | spglib | 70 CIF | 70/70 identical number and Hall number |
| P1 positions | gemmi expansion | 71 CIF | maximum periodic distance `0.021 A` (MIL-100 special-position rounding); 70 within `2e-3 A` |
| Periodic neighbours | ASE | 71 CIF x 2 cutoffs | 142/142 identical pair counts; maximum distance difference `2.1e-5 A` |
| PXRD profiles | pymatgen `XRDCalculator` | 71 CIF | correlation `>= 0.9974`, median `0.99994`, maximum normalized difference `0.038` |
| Standardized and primitive cell | spglib | 70 CIF | lengths and angles exact; volume `<= 2.8e-7` relative; atom counts 70/70; pymatgen `StructureMatcher` confirms the same crystal for MFI, MOF-74, MIL-53, ZIF-8 and triclinic |
| Niggli reduction | spglib | 71 CIF | `<= 5e-7 A`, `<= 7e-7 deg`, `<= 7e-8` relative volume |
| CIF round trip | gemmi | 71 CIF | cell parameters and atom counts exact |
| Element table | mendeleev | 118 elements | Pyykkö radii 118/118; Bondi vdW 117/118 within `0.011 A` (Br 1.85 vs 1.83) |
| Magnetic space groups | spglib axial data | 3 collinear cells | UNI, operation and time-reversal counts 3/3 |

The experimental check compares the density from the measured COD cell with
the handbook value. Si, diamond, Fe, TiO2, quartz, graphite and CaF2 agree
within `0.12 %`; NaCl and ZnO differ by about 1–2 % because the archived COD
cell differs from the handbook cell, not because the expansion loses atoms.

## Full public archives

These archives are **not redistributed**; they were downloaded for the
campaign and processed in place.

| Archive | Version / licence | Scope | Result |
|---|---|---|---|
| CoRE MOF 2019 | public v2 2024-11-19, CC BY 4.0, Zenodo 14184621 | 19,081 CIF | 19,081 parsed; atom count, volume and density identical to gemmi; 18,848 space groups compared, 9 (`0.05 %`) differ at `symprec = 1e-3` and all 9 converge to the spglib answer at `1e-2` |
| QMOF | v18, CC BY 4.0, Figshare 13147324 | 20,372 relaxed CIF | 20,372 parsed; atom count, volume and density identical to gemmi; 19,921 (`97.8 %`) space groups agree with spglib; 15,401 match the tabulated QMOF group (spglib matches 15,158) |
| tmQM | 2024 release X1+X2+X3, MIT, GitHub `uiocompcat/tmQM` | 108,541 XYZ, 7,217,465 atoms | 0 parse errors and 0 stoichiometry mismatches |

The nine CoRE differences and the QMOF differences are tolerance-boundary
cases in flexible frameworks: at `1e-4 A` both searches return the lower
group and at `1e-2 A` both return the higher one. They are not missing atoms
or broken cells.

## External programs

| Program | Version | Case | Result |
|---|---|---|---|
| DFTB+ | 25.1 (conda-forge `nompi`, mio-1-1) | methane | generated input runs; energy `-87.774977538 eV`, max force `0.015703 eV/A`, matching `detailed.out` |
| xTB | 6.8.0 binary | methane and water | in-process GFN2 energy within `6.1e-6 eV`; forces within `1e-5 eV/A` |
| phonopy | 4.10.0 | diamond Si, 3x3x3 UFF supercell | Γ optical mode `749.212 cm^-1` (phonopy) vs `749.213 cm^-1` (`ff vibrations`); acoustic modes zero |
| MOFid | upstream tool (GPL, built locally with Open Babel + Java) | ABAVIJ | external tool returns `rtl` topology and the pyridine-carboxylate linker; the xtalforge subset returns `unknown` and a different decomposition |
| RASPA2 | 2.0.50 (built from source) | methane box; Cu-BTC/CO2 | RASPA gives `11.21 mol/kg` at 323 K / 10 bar; the native GCMC gives `0.15` molecules/cell (about `0.062 mmol/g`) and is nearly pressure-independent |
| Zeo++ | 0.3 | six prepared MOFs | `D_i` agrees within the grid resolution; `D_f` and POAV are not Zeo++ values (see the analysis page) |

### DFTB+ input details

The real run exposed three input defects that the recorded-output tests could
not see:

1. the inline `GenFormat` block was missing the species line, the origin and
   the three cell vectors;
2. H and He requested a `p` shell in the standard `mio` sets;
3. forces were never requested.

The writer now emits a standalone `geometry.gen` with the correct GEN header,
species, fractional coordinates, origin and cell; the HSD includes it with
`<<< "geometry.gen"`; H/He use `s`; and the `Analysis` block sets
`PrintForces = Yes`. The regression test
`dftb_writer_contains_the_required_blocks` checks the generated GEN file and
HSD.

### Phonopy atom-order caveat

The first phonopy comparison gave `697.5 cm^-1` because the displaced
supercell was written as a CIF and the CIF reader reordered the atoms before
the forces were returned to phonopy. Writing the displaced supercell as a
POSCAR preserves the order and the two implementations agree at
`749.212`/`749.213 cm^-1`. This is a useful warning for any external
force-constant workflow: **the force array must follow phonopy's supercell
order exactly**.

## Defects found and corrected

1. **Workspace manifest.** Cargo 1.97 rejected a member that disabled default
   features while the workspace entry left them implicit. The facade
   workspace dependency now sets `default-features = false`, and the CLI and
   Python crates opt back in.
2. **Standardized/primitive cell orientation.** moyo's column-major basis was
   read as row-major, shearing every non-orthogonal standardized cell. Quartz
   reported `a = 5.492 A, b = 4.254 A, gamma = 116.6` and `standardize`
   produced 9 sites that expanded to 45 atoms. All 70 parseable corpus cells
   now match spglib; a quartz regression test pins the hexagonal basis.
3. **Reduction orientation.** The same transpose affected Niggli, Delaunay and
   Minkowski reduction. Quartz and the full corpus now match spglib within
   `5e-7 A` / `7e-7 deg`.
4. **DFTB+ writer.** The three input defects listed above; fixed and covered by
   the external writer test and a real DFTB+ 25.1 run.
5. **PXRD audit reference.** The browser reference omitted Debye-Waller
   factors, so TiO2 (anisotropic `U_ij`) differed by 2.1 percentage points. The
   reference now applies the CIF `U_ij` in the crystal basis.
6. **Indexing test.** The monoclinic assertion rejected the equivalent
   `beta -> 180 - beta` setting; it now accepts both.
7. **Python stubs and clippy.** The engine override annotation uses
   `dict[str, object]` instead of `Any`; a `manual_range_contains` lint was
   corrected.
8. **cargo-deny.** The unmaintained `chrono` advisory is a transitive
   dependency of optional `arrow/parquet`; it is documented with
   `RUSTSEC-2026-0543` and the rest of the audit passes.

## Known gaps

- **MOFid parity**: the subset does not reproduce the external `rtl` /
  pyridine-carboxylate decomposition on ABAVIJ. It remains a documented
  subset, not a faithful port.
- **Native GCMC**: it is a Henry/Langmuir-level screening heuristic, not a
  detailed-balance GCMC. It is nearly pressure-independent and should not be
  used for production isotherms; RASPA remains the reference.
- **RASPA export**: `raspa_input` writes a skeleton (`simulation.input`,
  framework CIF, flat force-field list, mixing rules). It does not yet write
  `pseudo_atoms.def`, `force_field_mixing_rules.def` or molecule definitions,
  so it is not directly runnable.
- **QMOF/CoRE space groups**: differences are tolerance-boundary cases in
  flexible frameworks. The tabulated QMOF groups themselves match xtalforge
  for 75.6 % and spglib for 74.4 %, so neither tool is a definitive reference
  for those cells.
- **Benchmark baseline**: `benches/compare/baseline.json` now records this
  Xeon; timings on another CPU are advisory (the gate says so explicitly).
- **Licences**: CoRE MOF, QMOF, tmQM, MOFid, RASPA2 and the DFTB+ parameter
  sets are not redistributed with xtalforge; only the numeric/derived results
  and the small subsets already in `tests/external/` are committed.

## Where the evidence lives

- Campaign summary: this page.
- Browser/WASM battery: [`docs/guide/audit.md`](./audit).
- Dataset inventory and licences:
  [`docs/algorithms/external-validation.md`](../algorithms/external-validation).
- Algorithm-level tolerances: [`docs/algorithms/validation.md`](../algorithms/validation).
- DFTB+ regression: `crates/xtalforge-ff/tests/external.rs`.
- Quartz lattice regressions:
  `crates/xtalforge-symmetry/tests/standardize_lattice.rs`.
- Magnetic UNI reference:
  `crates/xtalforge-symmetry/tests/fixtures/magnetic-spglib.json` and
  `crates/xtalforge-symmetry/tests/magnetic_golden.rs`.
- PXRD reference: `tests/e2e/reference/generate.py` and
  `tests/e2e/reference/reference.json`.
- Benchmark baseline: `benches/compare/baseline.json`.
