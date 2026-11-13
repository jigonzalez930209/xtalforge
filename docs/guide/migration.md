# Migration guide from Crystal-Builder

The headless core of Crystal-Builder (`xtal/`, about 95,000 lines) is the
behavioural reference. xtalforge ports it to Rust and exposes the same work on
four surfaces. This page maps each module and command and records the
differences that will bite a script.

## Modules

| Crystal-Builder (`xtal/`) | xtalforge | Notes |
|---|---|---|
| `core/elements.py`, `core/pyykko.py` | `xtalforge-core` (`elements`, `data`) | same element table, one JSON bundle |
| `core/lattice.py`, `core/structure.py` | `xtalforge-core` (`lattice`, `structure`) | fractional `[f64; 3]`; Angstrom/degrees at the API |
| `core/p1.py` | `xtalforge-core::p1` | provenance kept on every expanded atom |
| `core/supercell.py`, `core/slab.py` | `xtalforge-core` and `xtalforge-commands` | commands are undoable documents |
| `core/neighbors.py`, `core/bonds.py`, `core/rings.py` | `xtalforge-core` | linked-cell neighbors; periodic rings keep images |
| `core/groups.py`, `core/measure.py` | `xtalforge-core` | group finder and measures |
| `io/` | `xtalforge-io` | format registry, project zip, trajectories |
| `symmetry/` | `xtalforge-symmetry` | moyo backend; our own dataset types |
| `ff/` | `xtalforge-ff` | UFF/UFF4MOF, QEq/EQeq/EEQ, Ewald, FIRE/L-BFGS, vibrations, engines |
| `ff/ase_engine.py`, `ff/mace`, `ff/orb`, `ff/mattersim` | Python `MlCalculator` | Python bridge, not native Rust |
| `builder/` (or `build/`) | `xtalforge-build` | molecules, carbon, polymers, hydrogens, fragments |
| `mof/` (PORMAKE port) | `xtalforge-mof` | curated topology/block subset, not the full database |
| `analysis/` (PXRD, grids, porosity, nets) | `xtalforge-analysis` | RCSR index embedded |
| `powder/` (indexing, Le Bail, Pawley, Rietveld) | `xtalforge-analysis` (modules `indexing`, `lebail`, `pawley`, `rietveld`) and `xtalforge powder ...` | there is no separate powder crate; refinement ranks are documented per page |
| `modules/`, `workspace/` | `xtalforge-modules` | on-disk runs and reports |
| `agent/` | `xtalforge-agent` | verbs, session, diagnostics, MCP, skill |
| `plugins.py` | `xtalforge` registries and Python `Plugins` | no global loader |
| `cli.py` | `xtalforge-cli` | same operations, `--json` on every command |
| `xtalapp/` (PySide6 + VTK) | **out of scope** | no GUI; `xtalforge-view` + the agent renderer |

## Ecosystem migration guides

- [Migrating from ASE](./from-ase)
- [Migrating from pymatgen](./from-pymatgen)
- [Migrating from spglib](./from-spglib)

## CLI commands

| Crystal-Builder | xtalforge | Notes |
|---|---|---|
| `info` | `xtalforge info` | summary; `--json` is the structured form |
| `render` | `xtalforge render` | SVG or PNG, orthographic; no VTK, no labels |
| `capabilities` | `xtalforge capabilities` | runtime features, engines, bundles, formats, verbs |
| `skill` | `xtalforge skill path|list|show|install` | embedded documents |
| `mcp` | `xtalforge mcp` | JSON-RPC over stdio; no window proxy |
| `symmetry` | `xtalforge symmetry ...` | detect, standardize, primitive, subgroups, k-path |
| `convert` | `xtalforge convert` | loss warnings become errors with `--strict` |
| `bonds` | `xtalforge bonds` | graph with periodic images |
| `types` | `xtalforge ff types` | UFF atom types and confidence |
| `formats` | `xtalforge formats` | registry and capabilities |
| `modules` | `xtalforge-modules` Rust API | no CLI command; the registry, jobs and run records are library-only |
| `engines` | `xtalforge ff engines` | probe result per engine |
| `run` | `xtalforge script` | JSON command list, undoable |
| builders | `xtalforge build ...`, `xtalforge mof ...` | names differ; `--help` lists the options |
| `agent` | `xtalforge agent script.json` | one process, one session |

## Agent verbs

Crystal-Builder's agent calls Python methods; xtalforge calls verbs by name.
The workflow maps one to one:

| Crystal-Builder agent | xtalforge verb |
|---|---|
| `session.open` | `load` |
| `session.info` | `info` |
| `session.warnings`/diagnostics | `diagnostics` |
| `session.supercell` | `supercell` |
| `session.wrap` | `wrap` |
| `session.add_hydrogens` | `hydrogens` |
| `session.optimize` | `optimize` |
| `session.symmetry` | `symmetry` |
| `session.pxrd` | `pxrd` |
| `session.porosity` | `porosity` |
| `session.uff` | `uff_energy` |
| `session.render` | `render` |
| `session.undo`/`redo` | `undo`/`redo` |
| change log | `history`, `changes` |
| `capabilities` | `capabilities` |

## Behavioural differences

- **Units.** Angstrom and degrees at every boundary; eV, eV/Angstrom and GPa
  for energies, forces and stress. Crystal-Builder's UFF internals are
  kcal/mol; the conversion happens inside the UFF model.
- **Errors.** Every error has a stable `code()`; the CLI prints it in JSON and
  uses exit code 3, Python raises `XtalforgeError` with the code as the first
  argument, WASM throws an object with `code`.
- **No panics on user input.** A malformed file or an impossible option is a
  typed error. A panic is a bug.
- **Determinism.** The same input and options give the same output independent
  of thread count. Randomised builders take an explicit seed.
- **Pictures.** The renderer is a small orthographic rasteriser, not VTK. It
  has no transparency, ambient occlusion or text labels; it is meant for an
  assistant to see the model.
- **No window.** The GUI, the window proxy and the `mcp.json` discovery are
  not ported. The MCP server is the transport.
- **Curated data.** The MOF builder ships a small illustrative topology/block
  subset. The full PORMAKE database is an external bundle, not part of the
  crate.
- **Powder ranking.** The automatic table ranks extinction classes by a Le
  Bail BIC and Pawley-fits every leading (cell, class); the 1 % Rwp tie rule
  orders it. Read `docs/algorithms/` for the exact procedure.
- **Python-only ML.** MACE, ORB-v3 and MatterSim run through ASE in Python;
  the Rust `ml` feature is a declared no-op.
- **MCP cancellation.** Verbs are synchronous; a client that must stop one
  terminates the server. The size limits bound the work.

## Moving a script

1. Replace `xtal` imports with the Rust crate, the Python package or the CLI.
2. Replace exceptions with the stable codes; do not parse the message.
3. Replace VTK screenshots with `render`.
4. Replace `session` method calls with verb names and JSON arguments.
5. Check every unit at the boundary once; the internal units are not public.
6. Run the same structure through `info` and `diagnostics` first; they catch
   most migration mistakes (wrong cell, lost symmetry, missing hydrogens).
