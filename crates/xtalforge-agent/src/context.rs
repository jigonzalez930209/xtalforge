//! Context for chat clients (Task 8.2.2).
//!
//! The `initialize` result carries an `instructions` string: what a client
//! that cannot read the skill files from disk is told before it calls a tool.
//! It names the workflow, the answer shape and the limits, and points at
//! `capabilities` instead of quoting an engine from memory.

use crate::verbs;

/// Instructions sent to an MCP client.
#[must_use]
pub fn instructions() -> String {
    format!(
        r#"xtalforge {version}: build, manipulate, analyse and export crystal structures.

Start with the `capabilities` tool. It reports what this installation can do
(features, external engines on PATH, embedded data bundles) and the verbs that
are actually available. Do not assume an engine is installed.

Typical workflow:
1. `load` a structure from a path or from CIF/POSCAR/XYZ text and give it a
   document name; later verbs default to the active document.
2. `info` returns the summary, cell, composition, site table and diagnostics
   in one answer. Read `diagnostics` before editing: `overlapping_atoms` and
   `empty_structure` are errors, the others are warnings or information.
3. Edit with `supercell`, `wrap`, `hydrogens` or `optimize`; every edit goes
   through the undo history with origin `agent`.
4. Analyse with `symmetry`, `bonds`, `rings`, `pxrd`, `porosity` and
   `uff_energy`.
5. `render` draws SVG or PNG; `save` writes a file in a registered format.
6. `changes` lists the assistant's edits in bursts; `undo` and `redo` move
   through the history.

Every tool takes one JSON object and returns one JSON object with a `kind`
field. Errors carry a stable `code` (`unknown_verb`, `invalid_argument`,
`unknown_document`, `no_document`, `io`, or a domain code). An MCP error
result sets `isError: true` and puts the code in the text payload.

Safety: the session runs under the standard size limits and checks them before
a supercell or a distance grid allocates. Verbs are synchronous; the server
accepts `notifications/cancelled` but can only apply it between calls. To stop
a long-running verb, terminate the `xtalforge mcp` process.

There are {verbs} tools. Engine binaries (DFTB+, xTB, ORCA) are optional; the
native UFF, QEq/EQeq, PXRD, porosity and net analysis work without them.

Full guide: docs/guide/automation.md. Skill files: `xtalforge skill path`."#,
        version = xtalforge_core::version(),
        verbs = verbs::VERBS.len(),
    )
}
