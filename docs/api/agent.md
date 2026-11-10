# Agent API

The agent layer is what an AI assistant drives: a stateful session, one verb
catalogue and typed JSON answers. The same catalogue is exported by the Rust
crate, the CLI (`xtalforge capabilities` and `xtalforge agent`), the Python
package (`AgentSession`) and the WebAssembly package (`AgentSession`).

- [Algorithm page](../algorithms/agent)
- Verb catalogue: [`xtalforge::agent::verbs`](https://docs.rs/xtalforge-agent)
- Parity and known exceptions: [surface parity](./parity)

## Rust

```rust
use serde_json::json;
use xtalforge::agent::{Answer, Session};

let mut session = Session::new();
session.run(
    "load",
    &json!({"path": "MOF-5.cif", "name": "mof5"}),
)?;
let Answer::Info { inspection, diagnostics, .. } =
    session.run("info", &json!({}))?
else {
    unreachable!()
};
println!("{}: {} problems", inspection.summary.formula, diagnostics.diagnostics.len());
# Ok::<(), xtalforge::agent::Error>(())
```

The batch form is `agent_run` / `run_verbs`: it executes a JSON array of
`{verb, arguments}` objects against a fresh session and returns one answer
each. The catalogue itself is `agent_verbs` (`xtalforge::agent::catalogue`)
and the runtime report is `agent_capabilities`
(`xtalforge::agent::capability_report`).

`Session::run` takes the verb name and a JSON object. Every missing or
mistyped argument is an [`xtalforge::agent::Error::Argument`], never a panic.
The answer is one of the `Answer` variants; `Answer::to_json` serializes it.

## CLI

```bash
xtalforge capabilities
xtalforge render MOF-5.cif --out mof5.png --format png --width 1200 --height 900
xtalforge --json agent session.json
```

An agent script is one JSON array (or an object with a `verbs` array):

```json
{
  "verbs": [
    {"verb": "load", "arguments": {"path": "MOF-5.cif", "name": "mof5"}},
    {"verb": "info", "arguments": {}},
    {"verb": "render", "arguments": {"format": "png", "path": "mof5.png"}}
  ]
}
```

`--json` prints one array with one answer per call. Without it, each answer is
printed as `verb: {json}`.

## Python

```python
import json
import xtalforge

print(xtalforge.agent_capabilities())
session = xtalforge.AgentSession()
session.run("load", json.dumps({"path": "MOF-5.cif", "name": "mof5"}))
info = json.loads(session.run("info", "{}"))
print(info["inspection"]["summary"]["formula"])
structure_json = session.structure_json("mof5")
```

`run` takes and returns JSON strings so the Python surface has the same
arguments and answers as the Rust, CLI and WebAssembly ones. Errors are
`InputError` (bad arguments) or `XtalforgeError`; the first exception argument
is the stable code.

## WebAssembly

```js
import init, { AgentSession, agent_capabilities } from "xtalforge-wasm";
await init();
console.log(JSON.parse(agent_capabilities()));
const session = new AgentSession();
session.run("load", JSON.stringify({ text: cif, format: "cif", name: "mof5" }));
const info = JSON.parse(session.run("info", "{}"));
```

The browser has no filesystem: `load` works from `text`, `save` and the `path`
form of `render` are not available, and every in-memory verb works. The
exception is recorded in the parity matrix.

## MCP server

```bash
xtalforge mcp
```

The server speaks JSON-RPC 2.0 over stdio (protocol `2024-11-05`), one tool
per verb. `initialize` returns the protocol, the server name, the tool
capability and an `instructions` string for clients that cannot read files
from disk. `tools/list` builds each `inputSchema` from the same catalogue;
`tools/call` returns the answer as text with `isError` set on failure. The full
message set, error codes and lifecycle are in the
[automation guide](../guide/automation#start-the-mcp-server).

## Skill documents

```bash
xtalforge skill path       # default install directory
xtalforge skill list       # every embedded document
xtalforge skill show SKILL # one document
xtalforge skill install    # write all of them
```

`xtalforge skill install <dir>` writes the main skill, the MCP setup, three
worked prompts and nine reference pages. They are embedded in the binary; no
download is needed.

## Verbs

`xtalforge capabilities` lists the verbs. They are grouped as:

| Group | Verbs |
|---|---|
| Session | `capabilities`, `list_documents`, `load`, `save`, `close` |
| Inspection | `info`, `diagnostics` |
| Structure edits | `supercell`, `wrap`, `hydrogens`, `optimize` |
| Analysis | `symmetry`, `bonds`, `rings`, `pxrd`, `porosity` |
| Force field | `uff_energy` |
| Pictures | `render` |
| History | `undo`, `redo`, `history`, `changes` |

Every mutating verb goes through `xtalforge_commands::Document` with
`Origin::Agent`. `changes` groups consecutive agent edits into bursts and
reports how many `undo` calls remove the last one.

## Errors

| Code | Meaning |
|---|---|
| `unknown_verb` | The name is not in the catalogue |
| `invalid_argument` | A required argument is missing or has the wrong JSON type |
| `unknown_document` | The session has no document with that name |
| `no_document` | The verb needs a document and none is active |
| `io` | File input or output failed |
| core, symmetry, analysis and force-field codes | The underlying operation failed |

## Limits

The session runs under the standard `Limits` profile. PXRD, porosity, UFF and
optimization are the expensive verbs; a client should call `capabilities`
first and avoid `optimize` on a large framework without a step limit.
