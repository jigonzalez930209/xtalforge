# Agent session and diagnostics

![Agent session](/diagrams/agent-session.svg)

Reference: Crystal-Builder `xtal/agent/` (MIT; `capabilities.py` 227 lines,
`session.py` 1834, `inspect.py` 352, `diagnostics.py` 253, `changes.py` 197,
`answers.py` 210, `render.py` 137). This page records how the port works and
where it deliberately differs.

## 1. Problem

An AI assistant needs a small, stable vocabulary of operations over a
document, answers it can print and pass on, and a machine-readable statement of
what the installation can do. Two failure modes are worse than a missing
feature: an assistant that quotes a capability the build does not have, and an
assistant whose edits a user cannot review or undo as a unit.

## 2. Catalogue as data

`verbs::VERBS` is one static table. Each entry declares:

- the name the surfaces call;
- a one-line summary;
- the JSON type, required flag and default of every parameter;
- the answer kind;
- a complete JSON example;
- whether it mutates the document.

The table is compiled into the crate. The CLI, Python and WebAssembly surfaces
read it; there is no second list. `capabilities` adds the facts that only a
running installation knows: compiled features, external binaries that answer
`--version` on `PATH`, embedded data bundles and the registered formats.

## 3. Session

A `Session` owns a map from name to `Document`. `run(verb, arguments)` does:

1. verify the name is in the catalogue (`unknown_verb`);
2. parse the arguments and report the first bad one (`invalid_argument`);
3. resolve the document, falling back to the active one;
4. call the operation;
5. wrap the result in an `Answer`.

Every mutation becomes a `Command` with `Origin::Agent`, applied through
`Document::apply`; the session never edits a structure in place. This is what
makes `undo`, `history` and `changes` see the assistant exactly like a person.

## 4. Diagnostics

`diagnose(structure)` expands the P1 cell once and runs six checks. Each
finding has a stable `code`, a `Severity` (`info`, `warning`, `error`), a
message, the atom indices and an optional verb that would fix it.

| Code | Severity | Rule |
|---|---|---|
| `empty_structure` | error | no atoms in the P1 cell |
| `overlapping_atoms` | error | periodic pair distance `< 0.4 A` |
| `partial_occupancy` | info | any site occupancy `< 0.999` |
| `charge_imbalance` | warning | `abs(total_charge) > 1e-6 e` |
| `dangling_atoms` | warning | non-metal atom with no perceived bond |
| `missing_hydrogens` | info | C/N/O atom with fewer bonds than its valence |
| `unusual_coordination` | warning | more than 12 bonds on one atom |
| `no_symmetry` | info | P1 with more than three atoms |
| `symmetry_failed` | info | the symmetry search did not finish |

The 0.4 A overlap threshold is the reference's; the valence table is H/D 1,
C 4, N 4, O 2, halogens 1, S 6, P 5, B 4. The valence check is deliberately
conservative: a carbon in a carbanion or a hypervalent sulfur produces an
`info`, not an error. Symmetry findings are `info` because a deliberately
asymmetric model is valid.

Findings are sorted by severity (error, warning, info) and then by code, so
two runs on the same structure produce byte-identical reports.

## 5. Change bursts

The command log records `(description, origin)` per edit. `changes_for` walks
the log and groups consecutive `Origin::Agent` records. A burst is
`(start, len, descriptions, undo_steps)`. `undo_steps` is `Some(len)` only when
the burst ends the log, because an interleaved user edit makes a single undo
insufficient. This is the reference behaviour, expressed over the typed log
instead of a Python object.

## 6. Pictures

The renderer takes the scene data of Task 3.14.4 (atoms, bonds, cell edges,
radii and colours) and projects it orthographically:

1. camera direction `f`, up hint; right `r = normalize(up x f)` and
   `u = normalize(f x r)`;
2. centre and half-extents of every visible sphere and cell edge in `(r, u)`;
3. `scale = min((W - 2m) / 2h_r, (H - 2m) / 2h_u)`;
4. screen `x = W/2 + (p - c).r * scale`, `y = H/2 - (p - c).u * scale`,
   depth `(p - c).f` (larger is closer).

SVG writes the primitives sorted by depth, cell edges first. PNG rasterises
with a depth buffer: sphere pixels take the near surface
`z + sqrt(R^2 - d^2)`; bond and cell segments take the interpolated depth
along the segment. One fixed light shades spheres and capsules. The PNG
encoder writes the signature, `IHDR`, one zlib-compressed `IDAT` and `IEND`;
`flate2` supplies deflate and the CRC-32 is computed locally. No C library,
image crate or GPU is involved.

The reference uses VTK and off-screen rendering. This port is intentionally
smaller: no transparency, no ambient occlusion, no text labels. It is good
enough to see a folded linker, a missing cluster or a wrong cell.

## 7. Complexity

| Operation | Cost |
|---|---|
| Load and summary | one P1 expansion, `O(N)` |
| Diagnostics | `O(N^2)` periodic pairs plus one symmetry search |
| Supercell, wrap, hydrogens | one command, `O(N)` to `O(N log N)` |
| PXRD | `O(reflections x profile)` |
| Porosity | one grid of `prod(ceil(L_i / s))` points |
| UFF energy | one model build plus the neighbour list |
| Render PNG | `O(W H)` plus the primitives' bounding boxes |

## 8. Validation

- `crates/xtalforge-agent/tests/verbs.rs`: catalogue uniqueness, every example
  only names declared parameters, required parameters have no default.
- `crates/xtalforge-agent/tests/session.rs`: load/info/diagnostics on NaCl,
  unknown verb/document errors, supercell plus undo/redo plus change bursts,
  agent origin in the history, save round trip, Fm-3m detection, diamond
  six-rings, PXRD and porosity finiteness, UFF energy and a water optimization
  that does not raise the energy, close.
- `crates/xtalforge-agent/tests/render.rs`: SVG size and atom circles, PNG
  signature and `IHDR` dimensions, the no-visible-atom error, dimension checks
  and the CLI path.
- `crates/xtalforge-wasm/tests/agent.rs`: catalogue and capabilities JSON, a
  session that loads CIF text and answers `info`, and the coded unknown-verb
  error.

## 9. MCP transport

The server is a minimal JSON-RPC 2.0 implementation: newline-delimited
messages on stdio, methods `initialize`, `notifications/initialized`, `ping`,
`tools/list`, `tools/call` and `notifications/cancelled`, protocol version
`2024-11-05`. A full MCP crate was considered and rejected: the surface is one
request/response pair per tool, the schema is generated from `VERBS`, and the
message set fits in one file (Task 8.2.1 records the choice). Malformed JSON is
`-32700`, an unknown method is `-32601`, a bad `tools/call` parameter is
`-32602`, and a tool failure is a successful response with `isError: true` and
the stable code in the text payload.

`initialize` carries `instructions` (Task 8.2.2), so a client that cannot read
the skill files still learns the workflow, the answer shape and the limit
behaviour.

## 10. Skill documents

Fourteen markdown files are embedded with `include_str!`: the main skill, the
MCP setup, three worked prompts (MOF build, porosity report, deposited CIF)
and nine references (API, calculations, carbon, diagnostics, MOF, polymer,
porosity, preparation, workspace). `xtalforge skill show` prints one,
`install` writes the tree, and the published form is the automation guide. The
files are not generated at runtime: a release ships the text that was
reviewed.

## 11. Safety

Every call runs in the session's standard `Limits` profile. The supercell
check runs before the command allocates, the porosity grid check before the
distance grid, and the render check before the raster buffer. A refusal is
`Error::LimitExceeded` with the verdict's sentence.

Verbs are synchronous. `notifications/cancelled` is accepted but can only take
effect between calls; the server cannot preempt a running verb without
threads and a way to stop the native computation. A client that must interrupt
a call terminates the `xtalforge mcp` process. This is a recorded deviation
from the reference, which runs the same synchronous verbs in a Python process.
Timeouts are therefore enforced by the size limits (bounded work) and by
`optimize`'s step limit, not by a wall-clock timer that would report a result
the computation already produced.

The window proxy and discovery of Task 8.2.4 are out of scope: there is no
window. The reference's `mcp.json` shape (`url`, `pid`) is documented in the
automation guide for a future front end, but no code reads it.

## 12. Differences from Crystal-Builder

- The reference's window proxy and desktop discovery are out of scope; there
  is no window to proxy. The MCP transport is Task 8.2.
- Pictures use a small orthographic rasteriser instead of VTK; see section 6.
- The verb catalogue is shared with the bindings rather than generated from
  the CLI at runtime; both express the same operations.
- Diagnostics return atom indices in the P1 cell, not object identities.
