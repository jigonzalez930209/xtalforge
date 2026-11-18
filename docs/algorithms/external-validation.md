# External validation

The [2026-11-17 validation campaign](/guide/validation-campaign) records the
full archive runs, the real DFTB+/phonopy/MOFid/RASPA cases and the measured
gaps.

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

## Validation campaign (2026-11-17)

These references were executed during the full campaign and are recorded in the
[validation campaign guide](/guide/validation-campaign); they are not part of
the committed CI suite because the programs or archives are not redistributed.

| Source | Licence | What it validates | Result |
|---|---|---|---|
| Full CoRE MOF 2019 v2 | CC BY 4.0 | 19,081 CIF: parsing, formula, cell, density, space group | 0 parse/geometry/density differences; 9/18,848 space-group tolerance-boundary cases |
| Full QMOF v18 | CC BY 4.0 | 20,372 DFT-relaxed CIF: parsing, formula, cell, density, space group | 0 geometry/density differences; 97.8 % space-group agreement with spglib |
| Full tmQM X1-X3 | MIT | 108,541 XYZ, 7.2 M atoms: parsing and stoichiometry | 0 parse errors, 0 formula mismatches |
| Real DFTB+ 25.1 (`mio-1-1`) | DFTB+ open source | generated HSD/geometry, energy and forces | methane `-87.774977538 eV`, max force `0.015703 eV/A`, matching `detailed.out` |
| Real phonopy 4.10 + UFF | BSD-3-Clause | periodic force constants and Gamma frequencies | diamond Si optical mode `749.212` vs `749.213 cm^-1` |
| MOFid (GPL) | GPL | topology and MOFid/MOFkey on ABAVIJ | external tool `rtl`; xtalforge subset `unknown` (documented gap) |
| RASPA2 2.0.50 | RASPA licence | Cu-BTC/CO2 isotherm | `11.21 mol/kg` at 323 K/10 bar; native GCMC `0.15` molecules/cell (screening only) |

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
| CoRE MOF 2019/2024 (full archive) | framework porosity, MOF assembly, topology | CC BY 4.0; 19,081 CIFs validated on 2026-11-17 (0 parse/geometry/density differences) but not redistributed |
| tmQM (full release) | transition-metal geometries | MIT repository licence; 108,541 complexes / 7.2 M atoms validated on 2026-11-17, not redistributed |
| QMOF (full archive) | quantum MOF geometries and properties | CC BY 4.0; 20,372 CIFs validated on 2026-11-17, not redistributed |
| Materials Project | structures, phonons, elastic constants | terms page returned 403 and an API key is required; licence not verified, not bundled |
| NIST SRM 640/660 | PXRD line positions and instrument profiles | certified reference material; redistribution terms not verified, not bundled |
| RRUFF | mineral XRD/Raman references | non-commercial terms; not bundled |
| AMCSD | mineral structures and powder patterns | public-domain status not verified; not bundled |
| Phonopy example cells | phonon dispersions and DOS | BSD-3-Clause; a real phonopy 4.10 run with xtalforge UFF forces was executed on 2026-11-17; no external force-constant reference is bundled |
| DFTB+ test suite | HSD/output coverage | LGPL-3.0 code and CC BY-SA 4.0 docs; not copied; a real DFTB+ 25.1 methane run covers the generated input and parser |
| MOFid | MOFid/MOFkey identity | GPL tool built and run on 2026-11-17; xtalforge's documented subset still differs (see the campaign notes) |
| RASPA2 | GCMC isotherms and Widom insertion | RASPA 2.0.50 built and run on 2026-11-17; native GCMC is not quantitatively comparable yet |
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
  the 2026-11-17 run used xtalforge UFF forces with phonopy, not an external
  force-constant reference.
- The external engine **launchers** were exercised with real xTB 6.8 and
  DFTB+ 25.1 on this machine in the 2026-11-17 campaign; the committed CI
  suite still uses simulated executables for speed.
- The native GCMC is a screening heuristic, not a detailed-balance GCMC: on
  HKUST-1 at 323 K and 10 bar it gives `0.15` molecules per cell while RASPA
  gives `11.2 mol/kg`, and the loading barely changes between 1 and 100 bar.
  The RASPA export is a skeleton (force field, mixing rules and framework CIF)
  and is not a complete runnable input set.
- MOFid's external run on ABAVIJ gives `rtl` topology and a pyridine
  carboxylate linker; xtalforge's documented subset reports `unknown` and does
  not reproduce that decomposition.
