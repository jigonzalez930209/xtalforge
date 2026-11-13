# How to drive xtalforge

You are working with a crystal-structure library. The same operations exist in
Rust, the CLI, Python and WebAssembly, and an MCP server exposes them as tools.
This skill tells you the order to call them and what the answers mean.

## Before the first call

Call `capabilities`. It reports:

- compiled features (`sci-form`, `ml`, `external`);
- external engines and whether they answer on `PATH` (DFTB+, xTB, ORCA);
- embedded data bundles (the RCSR net index);
- the registered formats and the complete verb list.

Never promise an engine or a format that `capabilities` does not list.

## Documents

`load` reads a structure from `path` or from `text` and gives it a name.
Every later verb defaults to the active document; pass `document` to address
another one. `list_documents` and `close` manage the set. `save` writes a
registered format; the file extension chooses it when `format` is omitted.

All edits are undoable. `history` shows the command log with origins,
`changes` groups the assistant's edits into bursts, and `undo`/`redo` move
through the history. Prefer one macro edit over many small edits so the burst
is a useful review unit.

## Order of work

1. **Inspect.** `info` returns the summary, cell, composition, site table and
   diagnostics in one answer. `diagnostics` alone when you only need problems.
2. **Fix the model.** `hydrogens` for a deposited X-ray model without H;
   `wrap` to put sites in the cell; `supercell` for a larger cell.
3. **Understand it.** `symmetry` (space group, Wyckoff letters, standardized
   cell), `bonds` (perceived graph with images), `rings` (primitive rings).
4. **Measure.** `pxrd` (calculated pattern), `porosity` (distance grid),
   `uff_energy` (energy, forces, stress), `optimize` (relax with UFF).
5. **Show it.** `render` writes SVG or PNG. SVG is smaller and searchable;
   PNG is what most chat clients display.

## Answers

Every answer is one JSON object with a `kind` field: `capabilities`,
`documents`, `loaded`, `saved`, `info`, `diagnostics`, `structure`, `symmetry`,
`bonds`, `rings`, `pxrd`, `porosity`, `energy`, `optimization`, `picture`,
`history`, `changes`, `workflow_schema`, `workflow_validation` or
`workflow_run`.

Numbers are Angstrom, degrees, eV, eV/Angstrom and GPa at the boundary. PXRD
intensities are normalized to 100; the untouched profile values are not part of
the answer.

## Errors

Errors carry a stable `code`:

| Code | What to do |
|---|---|
| `unknown_verb` | Re-read `capabilities`; the name is wrong |
| `invalid_argument` | Fix the named argument; the message says why |
| `unknown_document` | Use `list_documents`; the name is not loaded |
| `no_document` | Call `load` first |
| `io` | The path is not readable or writable |
| domain codes | The operation rejected the geometry or the options |

An MCP error result has `isError: true`; the text payload is the JSON object
with `code` and `message`.

## Diagnostics

Read `diagnostics` before editing. `overlapping_atoms` and `empty_structure`
are errors and usually mean a bad file or a bad conversion. `charge_imbalance`,
`dangling_atoms` and `unusual_coordination` are warnings: check whether the
model is supposed to be ionic, whether the bond rules fit the chemistry, and
whether a metal really has that many neighbours. `partial_occupancy`,
`missing_hydrogens` and `no_symmetry` are information.

## Limits and cost

The session runs under the standard size limits. `supercell` and `porosity`
check the estimate before allocating and refuse with `LimitExceeded`.
Porosity cost grows with `1/spacing^3`; 0.4 A is a good default, 1.0 A is a
quick look. `optimize` defaults to 200 steps; raise it only when the answer
says the run did not converge.

## Files in this skill

- `mcp.md`: client configuration.
- `prompts/`: worked sequences for a MOF build, a porosity report and a
  deposited CIF.
- `references/`: one page per area, with the verbs, options and caveats.
- `references/workflows.md`: the declarative workflow document and its runner.
