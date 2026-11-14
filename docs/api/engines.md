# External engines

## Availability, not existence

An external engine is optional and found at run time. The `external` Cargo
feature is reserved and a no-op (the probes and input writers are always
compiled, see [Cargo features](../guide/features)); a
probe checks that the binary actually runs (`--version` with a short read and
a timeout) instead of only checking that a file exists. `EngineStatus` is
`Available { path, version }` or `Missing { reason }`; `xtalforge engines`
prints the registry including the reason.

On WASM the probe always reports `Missing` with "not available in the browser";
this is the recorded parity exception for every external engine.

## Process runner

`run_in_scratch(&RunSpec, &AtomicBool)` writes the request's relative input
files into a private scratch folder, launches the program there with stdin
closed, and redirects stdout/stderr to log files so a chatty engine cannot
deadlock the poll loop. It waits with a wall-clock `timeout`, checks the
cancellation flag every 10 ms, kills the child on either condition, and returns
`RunOutcome { status: Success | Failed | TimedOut | Cancelled, exit_code,
stdout, stderr, stdout_truncated, stderr_truncated, scratch, duration }`.
Scratch is removed by default; `keep_scratch` leaves it for inspection. Input
names must be relative and may not escape the folder. The tests use generated
shell scripts for success, timeout, cancellation, log caps and cleanup; no
external engine binary is needed.

## DFTB+

`write_dftb_hsd(structure, &DftbOptions)` writes a `dftb_in.hsd` with
`Geometry = GenFormat`, the `Hamiltonian = DFTB` block, the Slater-Koster
folder naming and the k-point mesh. `parse_dftb_output` reads `detailed.out`
(or `results.tag`) for the total energy, forces and the stress tensor.
`DftbCalculator::new(options)` writes the input into a scratch folder, launches
`dftb+` through `run_in_scratch`, reads `detailed.out`, checks the force count
and exposes the result through `Calculator`. `with_program`, `with_timeout`,
`keep_scratch` and `cancel_handle` configure the run. The parser is validated
against the recorded `si_scc_detailed.out` fixture, not only against a
simulated executable.

sci-form carries the SCC-DFTB Fock/gamma building blocks
(`scf::fock_matrix`) but not a Slater-Koster parameter set or a complete DFTB
driver, so DFTB+ remains an external program; the xTB GFN methods are the ones
available in-process (Task 4.10).

## xTB

`write_xtb_input` writes the XYZ block plus `$chrg`/`$uhf`; `parse_xtb_output`
reads `TOTAL ENERGY`, the charge and the forces when present.
`XtbCalculator` maps GFN1-xTB, GFN2-xTB and GFN-FF to the corresponding `xtb`
arguments, runs it through the same scratch runner and parses stdout; an
unsupported method is refused before launch.

## ORCA

ORCA is molecular. `write_orca_input` makes every molecule whole across the
cell faces using the bond graph, then writes `! method basis` and a `* xyz`
block with the charge and multiplicity. `orca_catalogue` lists the supported
functionals and bases. Only inputs are written; ORCA is never launched.

## Quantum adapter

`quantum::adapter(method)` is the in-process path for sci-form GFN0/GFN1/GFN2,
behind the optional `sci-form` feature (the published MIT crate `sci-form`
2.0.0). Without the feature it returns a clear typed error telling the user to
rebuild with `--features sci-form`. The adapter implements `Calculator`, so
`energy` and `energy_forces` have the same units and error contract as UFF.

The methods are **molecular**. sci-form takes element numbers and Cartesian
positions with no lattice, so the adapter refuses a structure that carries
space-group symmetry and documents that the P1 cell of a
`Structure::molecule` is evaluated as an isolated cluster. GFN2 returns forces
(central-difference gradient); GFN0 and GFN1 return an energy only and refuse
`energy_forces` instead of reporting zero forces. None returns an analytic
stress. Occupancies other than one and dummy `X` atoms are refused, and the
methods solve the neutral molecule (GFN0 is passed charge `0` because the
GFN1/GFN2 entry points take no charge).

The registry lists `gfn0`, `gfn1` and `gfn2`; they are `Available` only when
the feature is enabled. The cross-check against a real `xtb` binary is not run
here because the binary is not installed; the tests check finite water
energies, conservative GFN2 forces and the refusal paths.

## Machine-learning potentials

The decision is recorded in [the ML backend page](../algorithms/ml-backend.md):
no native backend in Phases 1-7, no silent model, and the Python bridge is the
only supported path. The `ml` feature stays reserved and a no-op.

## Input sets (Phase 18)

`xtalforge_engines` turns a structure into a reproducible set of engine files.
`list_presets()` returns the 20 registered names, `preset(name)` returns a
typed writer and `write_input_set(structure, &InputSetRequest { preset,
overrides })` merges the nested JSON overrides over the versioned preset data
and returns `InputSetSummary { preset, files }`. The presets are stored under
`crates/xtalforge-engines/data/input-sets/` with a provenance block (template
version, source, reference and license). Unknown override keys are
`Error::UnknownOverride` unless the branch is declared open or the request
sets `allow_unknown`.

| Family | Presets | Files |
|---|---|---|
| VASP | `vasp/{relax,static,band,dos,mof}` | `INCAR`, `KPOINTS`, `POSCAR`, `POTCAR.spec` (names only, never the licensed POTCAR) |
| Quantum ESPRESSO | `qe/{relax,vc-relax,scf,nscf,bands}` | `pw.in` with SSSP pseudopotential names and the symmetry k-path for bands |
| CP2K | `cp2k/{geo-opt,mof-geo-opt,single-point}` | `cp2k.inp` with MOLOPT/GTH names, cell and coordinates, DFT-D3 |
| LAMMPS | `lammps/{uff,uff4mof}` | `data.lammps` plus the UFF/UFF4MOF input script |
| GULP | `gulp/{relax,lattice-relax}` | `gulp.gin` |
| RASPA | `raspa/widom` | the Phase 10 RASPA export, byte-identical |
| ORCA | `orca/{single-point,pbe0-def2-tzvp}` | the Phase 4 ORCA input, byte-identical |

The same structure, preset and overrides give byte-identical files on every
surface. VASP `MAGMOM` and `LDAU` read the `magmom`, `ldau_u`, `ldau_orbital`
and `ldau_j` site properties (or structure metadata), with the documented
settings fallback.

## Result model and parsers

`CalculationResult` carries the final structure, `trajectory`, energies,
forces, stress, magnetic moments, charges, convergence flags, timings,
warnings, scalars, Fermi level, bands, DOS, orbital projections and the
producing engine. Parsers map malformed or truncated files to
`Error::Parse` (`engine_parse`):

- `parse_vasp_vasprun`, `parse_vasp_outcar`, `parse_vasp_oszicar`,
  `parse_vasp_eigenval`, `parse_vasp_doscar`, `parse_vasp_procar`;
- `parse_qe`, `parse_qe_xml`, `parse_qe_bands`, `parse_qe_projwfc`;
- `parse_cp2k`, `parse_lammps`, `parse_lammps_dump`, `parse_gulp`;
- `bader_charges`/`bader_table` and `ddec_charges`/`ddec_table` map charges
  onto site order.

The fixtures under `crates/xtalforge-engines/tests/fixtures/engines/` are
hand-written representative outputs (repository license) with provenance in
their README; they are not copied third-party data.

## Electronic structure

`band_structure(&ElectronicStructure)` returns labeled bands, `band_gap`
returns the direct/indirect gap and edges, `effective_mass`/`band_edge_masses`
fit parabolic band edges, and `density_of_states`/`smear_dos` build total,
element and orbital DOS with Gaussian smearing and the d-band centre. The WASM
and Python surfaces expose JSON-shaped results; the CLI adds `engines bands`
and `engines dos`.

## Jobs and schedulers

`JobOptions`/`scheduler_script`/`write_scheduler_script` write deterministic
SLURM and PBS scripts; `submit`/`poll`/`collect` drive them behind a runtime
probe. `ComputeStep`/`run_job` write an input set, run the configured engine
and parse the result through `ParserRegistry`; process launch is behind the
`external` feature, and a dry run writes inputs only. `JobCheckpoint` records
the typed state for restart/resume. `detect_failure`/`failure_catalogue` map
SCF non-convergence, `ZBRENT`, walltime, out-of-memory and missing
pseudopotential output to a corrective action. The CLI is
`xtalforge engines script`; WASM/Python only write scripts because processes
are host-side (recorded exception).

## Surfaces

CLI: `xtalforge engines presets|inputs|parse|bands|dos|script|failures`.
Python: `engine_presets`, `engine_input_set`, `engine_parse`, `engine_bands`,
`engine_dos`, `engine_script`, `engine_failures`.
WASM: the same names with `_info`, except `engine_script_info` returns text.

Parity keys: `engine_inputs`, `engine_outputs`, `engine_electronic`, `engine_jobs`.
