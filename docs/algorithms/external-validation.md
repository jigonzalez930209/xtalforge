# External validation

The suite has two layers: golden data produced by the reference implementations
that travel with the repository, and external data sets that were produced by a
different program or measured independently. This page records what is used,
under which licence, and what it can and cannot prove.

## Bundled and exercised in CI

| Source | Licence | What it validates | Test |
|---|---|---|---|
| Crystal-Builder corpus (`tests/corpus/`, 72 CIF) | COD CC0 and project models | CIF parsing, P1 expansion, bonds, rings, cell contents | `crates/xtalforge-io/tests/golden_*.rs` |
| COD subset (`tests/corpus/cod/`, 21 CIF) | CC0 | symmetry settings, occupancies, disorder, framework cells | Phase 1 golden corpus |
| RCSR 2019-06-01 (`data/rcsr/`, 2,931 nets) | RCSR terms, bundled with authorisation | net fingerprints and canonical keys | `crates/xtalforge-analysis/tests/rcsr.rs`; index build validates 2,909 keyed entries |
| Crystal-Builder + rietx rutile pattern (`tests/golden/phase-7/rietveld`) | MIT / derived | PXRD, Pawley and Rietveld refinement | `rietveld_golden.rs` |
| Crystal-Builder extinction classes (`tests/golden/phase-7/extinction`) | MIT / derived | class partition and the leading rutile class | `extinction.rs`, `automatic.rs` |
| EQeq authors' program, MOF-5 (`tests/fixtures/eqeq`) | numbers only, GPL program not copied | 424 EQeq charges on the COD MOF-5 cell | `charges.rs`, maximum difference `1.6e-4 e` |
| sci-form UFF (`sci-form`, MIT) | MIT | UFF total energy on methane and stretched water | `external_validation.rs`, within `1e-4` kcal/mol |
| DFTB+ recorded Si SCC output (`tests/fixtures/dftb`) | MIT / recorded output | energy, forces, stress parser | `external.rs`, energy `-69.2989 eV` |
| Analytic spring phonon model (`tests/golden/phase-4/phonons`) | project | q-mesh folding and dispersion | `phonons.rs` |
| CoRE MOF 2019 subset (`tests/external/core_mof/`, 18 CIF) | CC BY 4.0 (Zenodo 14184621) | CIF parsing, P1 expansion, cell volume, formula, bonding, porosity and PXRD on charged and mixed-metal frameworks | `crates/xtalforge-io/tests/external_core_mof.rs`, `crates/xtalforge-analysis/tests/external_core_mof.rs` |
| tmQM 2024 subset (`tests/external/tmqm/`, 18 XYZ) | MIT repository licence; citation in that directory's README | XYZ parsing and stoichiometry, perceived metal coordination, UFF/UFF4MOF types and forces for Sc-Ir | `crates/xtalforge-io/tests/external_tmqm.rs`, `crates/xtalforge-ff/tests/external_tmqm.rs` |
| QMOF subset (`tests/external/qmof/`, 14 CIF) | CC BY 4.0 (Figshare 13147324 v18) | CIF parsing, cell volume, density, formula, round-trip and space-group recovery on DFT-relaxed MOFs | `crates/xtalforge-io/tests/external_qmof.rs` |

The CoRE MOF, tmQM and QMOF subsets are the **production data sets** in the
table: the files are unmodified structures from the public archives, with the
archives' own formula, volume, density and stoichiometry used as external
checks. They are 256 KB, 76 KB and 100 KB; the full archives are not
redistributed. QMOF adds a space-group check: four P1 files recover their
tabulated groups at `symprec = 0.01 A`.

The EQeq and UFF comparisons are the two that use a **different program** on
the same coordinates. The rutile tables use Crystal-Builder's powder bridge,
which is the reference behaviour the port follows; the RCSR index uses the
reference canonical keys.

## Considered but not bundled

These are the data sets a wider validation campaign should add. None is in the
repository yet; the licence and the download size must be checked before it is.

| Source | What it would validate | Status |
|---|---|---|
| CoRE MOF 2019/2024 (full archive) | framework porosity, MOF assembly, topology | CC BY 4.0; the 81 MB archive is downloadable and the 18-CIF subset is bundled; the rest is not redistributed |
| tmQM (full release) | transition-metal geometries | MIT repository licence; the 65 MB archive is downloadable and the 18-complex subset is bundled |
| QMOF (full archive) | quantum MOF geometries and properties | CC BY 4.0; the 374 MB archive is downloadable and the 14-CIF subset is bundled, including the space-group checks |
| Materials Project | structures, phonons, elastic constants | terms page returned 403 and an API key is required; licence not verified, not bundled |
| NIST SRM 640/660 | PXRD line positions and instrument profiles | certified reference material; redistribution terms not verified, not bundled |
| RRUFF | mineral XRD/Raman references | non-commercial terms; not bundled |
| AMCSD | mineral structures and powder patterns | public-domain status not verified; not bundled |
| Phonopy example cells | phonon dispersions and DOS | BSD-3-Clause; not bundled: no force-constant reference to compare against, the analytic spring golden covers the implementation |
| DFTB+ test suite | HSD/output coverage | LGPL-3.0 code and CC BY-SA 4.0 docs; not copied, one recorded SCC output covers the parser |
| xTB benchmark sets | GFN0/GFN1/GFN2 energies and geometries | licences vary by set; the in-process methods are compared with sci-form instead |

## What the current suite does not prove

- UFF is compared numerically with sci-form on two small molecules. On CoRE
  MOF and tmQM only the build, typing and finiteness are checked: no published
  UFF energy exists for those structures. The tmQM subset does cover one
  complex for each of 18 metals, Sc through Ir.
- EQeq is compared on one framework and uses no external reference on tmQM
  (the reference program has no ionisation data for several of those
  environments); the parameter table's coverage beyond Zn/O/C/H is checked
  only for neutrality and sign.
- PXRD uses the reference's rutile pattern; a measured instrument profile
  (NIST or a laboratory scan) is not included.
- Phonons use an analytic spring model and a phonopy golden file from Phase 4;
  a full periodic force-constant comparison is not in CI.
- The external engine **launchers** are tested with simulated executables; a
  real DFTB+/xTB installation is not exercised, although the xTB methods are
  available in-process through `sci-form`.
