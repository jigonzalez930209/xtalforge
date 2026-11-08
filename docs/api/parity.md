# Surface parity matrix

Every public capability must exist in all four surfaces, or have a recorded and
justified exception. Status values: `todo`, `wip`, `done`, `n/a` (with the
reason in the exceptions table). `scripts/check-parity.sh` (Task 0.5.6)
compares this table with the symbols each surface actually exports.

| Capability | Phase | Rust | CLI | Python | WASM | Notes |
|---|---|---|---|---|---|---|
| `version` | 0 | done | done | done | done | Hello-world slice (Task 0.5.7) |
| `cell_volume` | 0 | done | done | done | done | Hello-world slice (Task 0.5.7) |
| `element_info` | 1 | done | done | done | done | Element data (Task 1.1.4) |
| `lattice_info` | 1 | done | done | done | done | Lattice (Task 1.2.5) |
| `lattice_minimum_image` | 1 | done | n/a | done | done | Rust, Python and WASM; the CLI reports geometry through `cell` |
| `limits_bounds` | 1 | done | n/a | done | done | Limits (Task 1.3.3); the CLI has the global `--limits` flag |
| `limits_check` | 1 | done | n/a | done | done | Limits (Task 1.3.3); the CLI has the global `--limits` flag |
| `structure_info` | 1 | done | done | done | done | Structure summary (Task 1.4.5); JSON bridge until CIF lands |
| `read_cif` | 1 | done | done | done | done | CIF reader (Task 1.15); CLI `info` and `convert` |
| `read_structure` | 1 | done | done | done | done | Format dispatch (Task 1.16.3); empty hint auto-detects |
| `supercell` | 1 | done | done | done | done | Supercells (Task 1.7); WASM `supercell_info` |
| `slab` | 1 | done | done | done | done | Slabs (Task 1.8); WASM `slab_info` |
| `bonds` | 1 | done | done | done | done | Bond graph (Task 1.11) |
| `rings` | 1 | done | done | done | done | Primitive rings (Task 1.12) |
| `groups` | 1 | done | done | done | done | Functional groups (Task 1.13) |
| `neighbors` | 1 | done | done | done | done | Periodic neighbors (Task 1.10) |
| `distance` | 1 | done | done | done | done | Measures (Task 1.14) |
| `angle` | 1 | done | done | done | done | Measures (Task 1.14) |
| `torsion` | 1 | done | done | done | done | Measures (Task 1.14) |

## Planned exceptions

| Capability | Surface | Reason | Alternative offered |
|---|---|---|---|
| Reading and writing paths, atomic writes | WASM | No filesystem in the browser | String and `Uint8Array` APIs |
| On-disk workspace of runs | WASM | No filesystem | In-memory workspace exported as a zip |
| External engines (DFTB+, xtb, Zeo++) and process probes | WASM | No subprocesses in the browser | Native engines (UFF, QEq, EQeq, grid porosity) |
| Multi-threading | WASM | Threads need cross-origin isolation | Optional threaded build |
| MCP server over stdio | Python, WASM | It is a process, served by the CLI | `xtalforge mcp` |
| Embedded MOF and RCSR bundles | WASM | Package size | Bundles fetched on demand |
