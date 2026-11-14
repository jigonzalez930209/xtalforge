# Crystal-Builder test-data inventory

Task 0.9.3 maps the data shipped in Crystal-Builder's `tests/data/` to the
xtalforge subtask that uses it. Only the items that a phase actually needs are
copied, with provenance in `tests/corpus/PROVENANCE.md` or `THIRD_PARTY.md`.
Items marked "not copied" below are represented by hand-written fixtures or by
the documented subset instead.

| Item | Contents | Used by |
|---|---|---|
| `blocks/*.xyz` | MOF connection-point fragments (U4Ph, UOMe, UPh) | MOF builder, Task 6.12 |
| `dftb/*` | DFTB+ outputs: Hessian, modes, DOS, band lines, MD, relaxation | External engine parsers, Task 4.9; the Si SCC output is also in `crates/xtalforge-ff/tests/fixtures/dftb/` |
| `eqeq/mof5_cod_authors_program.json` | reference EQeq charges for MOF-5 | Charge methods, Task 4.4; copied to `crates/xtalforge-ff/tests/fixtures/eqeq/` |
| `eqeq/nist_asd_ionization.csv.gz`, `pubchem_periodic_table.csv` | ionization data and periodic-table parameters | QEq/EQeq parameters, Task 4.4 |
| `hkust1_wide_nodes.xyz` | wide-node MOF test fragment | MOF builder, Task 6.12 |
| `pormake_upstream.json` | recorded PORMAKE output | MOF builder golden data, Task 6.12 |
| `rcsr/2dall.txt.gz`, `rcsr/3dall.txt.gz` | RCSR net snapshots | Topology and RCSR search, Tasks 5.6 and 5.7 |
| `rod_*.xyz` | rod SBU fragments for polymer/MOF assembly | Builders, Tasks 6.11 and 6.12 |
| `systre/mof5_net.cgd`, `systre/mof5_net.systre.txt` | net and Systre output for MOF-5 | Topology, Task 5.6 |
| `systre/rutile_net.cgd`, `systre/rutile_net.systre.txt` | net and Systre output for rutile | Topology, Task 5.6 |
| `zeopp_reference.json` | Zeo++ porosity reference numbers | Porosity, Task 5.2 |
| `zeo_radii.cc` | Zeo++ radii table | Porosity radii, Task 5.2.2 |

## License status

The full RCSR snapshots are bundled in `data/rcsr/` with the maintainer's
authorisation (2026-11-09) and their checksums; the licensing note is in
`THIRD_PARTY.md`. The Zeo++ radii table and the PORMAKE recording are not
copied. `crates/xtalforge-mof/data/` keeps the hand-written illustrative
subset and, behind `mof-data`, the converted PORMAKE 0.2.3 net and block
bundle with its license in `THIRD_PARTY.md`. The DFTB+ and EQeq files are
reference outputs, not third-party data sets; their origin is recorded when
they are copied.

## Copy rule

A phase copies an item only when its tests use it. The copy lands with the
subtask that introduces the test and names its source in the provenance file.
The inventory is reviewed when a new phase starts.

## Refresh (Task 17.4.2)

`scripts/golden/indexing.py` and `scripts/golden/rcsr_index.sh` are
deterministic. The 1.10 refresh re-ran the indexing generator with no diff;
the RCSR index rebuild verifies every fingerprint against the bundled
snapshot. The Crystal-Builder generator (`scripts/golden/generate.py`) still
needs a Python environment with the reference package, as recorded in its
header, so the Phase 1 JSON files are unchanged.
