# Third-party attribution

Code or data may only be ported or bundled after its license has been verified
and recorded here. Each ported file also starts with a comment naming its
origin. Status: `verified` (license read in the source), `to verify`,
`blocked` (not usable until resolved), `planned` (license fine, not yet used).

| Source | License | What is used | Status |
|---|---|---|---|
| sci-form (https://github.com/jigonzalez930209/sci-form, same maintainer) | MIT | Crystal-related modules ported (roadmap Task 0.4); optional dependency for quantum methods | verified, planned |
| Crystal-Builder (https://github.com/JulesOpp/Crystal-Builder) | MIT, "Copyright (c) 2026 Jules Oppenheim" | Behavioural reference, algorithms ported from `xtal/`, test corpus, golden data | verified, planned |
| moyo (https://github.com/spglib/moyo) | MIT OR Apache-2.0 | Dependency for symmetry; Niggli, Delaunay and Minkowski reductions possibly ported (Task 2.6.2) | verified (crates.io metadata), planned |
| PORMAKE 0.2.3 (https://github.com/Sangwon91/PORMAKE), as vendored in Crystal-Builder | MIT, "Copyright (c) 2022 Sangwon" | MOF builder algorithms (Task 6.12) | verified, planned |
| PORMAKE database (topologies, building blocks) | to verify (may derive from RCSR and other databases) | MOF data bundle (Task 6.11) | to verify |
| RCSR net snapshots shipped by Crystal-Builder (`RCSRnets-2019-06-01.cgd.gz`, `rcsr-2019-06-01.json.gz`, `rcsr-transitivity-2026-09-24.json.gz`) | to verify (RCSR terms of use) | Net identification and search (Tasks 5.6, 5.7) | to verify |
| COD sample CIFs in `Crystal-Builder/resources/samples/cod/` | CC0 (per its `PROVENANCE.md`) | Test corpus (Task 0.9.1) | verified, planned |
| Other Crystal-Builder samples (`simple/`, `prepared/`, `polymer/`, top level) | per-file, see its `PROVENANCE.md` | Test corpus | to verify per file |
| Zeo++ radii (`Crystal-Builder/tests/data/zeo_radii.cc`) | to verify (Zeo++ license) | Radii set for porosity (Task 5.2.2) | to verify |
| rietx and rietx-kernels (https://github.com/yue-here/rietx) | stated as MIT in Crystal-Builder's `pyproject.toml`; no license field on PyPI | Reference for powder refinement (Phase 7) | to verify |
| UFF and UFF4MOF parameters (Rappé et al. 1992; Addicoat et al. 2014; Coupry et al. 2016) | published parameters, cited in the data file | Force field (Task 4.2.1) | planned |
| Pyykkö and Atsumi covalent radii (2009) | published data, cited in the data file | Element data (Task 1.1.3) | planned |
