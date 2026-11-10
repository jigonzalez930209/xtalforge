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

The EQeq and UFF comparisons are the two that use a **different program** on
the same coordinates. The rutile tables use Crystal-Builder's powder bridge,
which is the reference behaviour the port follows; the RCSR index uses the
reference canonical keys.

## Considered but not bundled

These are the data sets a wider validation campaign should add. None is in the
repository yet; the licence and the download size must be checked before it is.

| Source | What it would validate | Status |
|---|---|---|
| CoRE MOF 2019/2024 | framework porosity, MOF assembly, topology | licence and size to verify; no bundle |
| QMOF | quantum MOF geometries and properties | CC BY 4.0; not downloaded |
| tmQM | transition-metal complex geometries | CC BY 4.0; not downloaded |
| Materials Project | structures, phonons, elastic constants | CC BY 4.0; API key and download policy |
| NIST SRM 640/660 | PXRD line positions and instrument profiles | certified reference material, redistribution to verify |
| RRUFF | mineral XRD/Raman references | non-commercial terms to verify |
| AMCSD | mineral structures and powder patterns | public-domain status to verify |
| Phonopy example cells | phonon dispersions and DOS | BSD; small, a good next addition |
| DFTB+ test suite | HSD/output coverage | licence to verify; one recorded output already covers the parser |
| xTB benchmark sets | GFN0/GFN1/GFN2 energies and geometries | licences vary by set; the in-process methods are compared with sci-form |

## What the current suite does not prove

- UFF is compared on two small molecules. A metal-organic framework would
  exercise the UFF4MOF metal types and periodic van der Waals lists.
- EQeq is compared on one framework. The parameter table's elements beyond
  Zn/O/C/H are only checked for neutrality and sign.
- PXRD uses the reference's rutile pattern; a measured instrument profile
  (NIST or a laboratory scan) is not included.
- Phonons use an analytic spring model and a phonopy golden file from Phase 4;
  a full periodic force-constant comparison is not in CI.
- The external engine **launchers** are tested with simulated executables; a
  real DFTB+/xTB installation is not exercised, although the xTB methods are
  available in-process through `sci-form`.
