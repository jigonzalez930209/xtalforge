# Automation guide

This guide is the published form of the skill an AI assistant reads: the MCP
server, the verb catalogue, the workflow and the safety rules. The same text
is shipped with the CLI (`xtalforge skill show SKILL`); `xtalforge skill
install` writes every file into `~/.local/share/xtalforge/skill`.

- [Agent API](../api/agent)
- [Agent algorithms](../algorithms/agent)
- [Surface parity](../api/parity)

## What an assistant can do

The agent layer exposes the library as **verbs**. A verb takes one JSON object
and returns one JSON object with a `kind` field. The verbs cover documents,
inspection, diagnostics, structural edits, symmetry, bonding, rings, PXRD,
porosity, UFF energy and optimization, pictures and history.

The catalogue is data, not prose: every surface reads
`xtalforge::agent::verbs::VERBS`, so the CLI, Python, WebAssembly and MCP
cannot drift from the Rust API.

## Start the MCP server

```bash
xtalforge mcp
```

It is JSON-RPC 2.0 over stdio, newline-delimited, protocol `2024-11-05`. A
client configuration is:

```json
{
  "mcpServers": {
    "xtalforge": {
      "command": "xtalforge",
      "args": ["mcp"]
    }
  }
}
```

`initialize` returns `instructions` for clients that cannot read files from
disk. The tools are named after the verbs and their `inputSchema` is generated
from the catalogue.

## First call

Call `capabilities` before anything else:

```json
{"verb": "capabilities", "arguments": {}}
```

It reports the version, compiled features, external engines that answered on
`PATH`, embedded data bundles, registered formats and the verb list. An
assistant must not claim an engine or a format that this report does not list.

## A session

```json
{"verbs": [
  {"verb": "load", "arguments": {"path": "MOF-5.cif", "name": "mof5"}},
  {"verb": "info", "arguments": {}},
  {"verb": "symmetry", "arguments": {"symprec": 0.001}},
  {"verb": "porosity", "arguments": {"spacing": 0.4}},
  {"verb": "render", "arguments": {"format": "png", "path": "mof5.png"}}
]}
```

Run it with `xtalforge --json agent session.json`. In Python:

```python
import json
import xtalforge

session = xtalforge.AgentSession()
session.run("load", json.dumps({"path": "MOF-5.cif", "name": "mof5"}))
info = json.loads(session.run("info", "{}"))
print(info["inspection"]["summary"]["formula"])
```

In WebAssembly the same class exists; `load` takes `text` because the browser
has no filesystem.

## Inspection and diagnostics

`info` is the one call that says what the structure is: summary, cell,
composition, site table and diagnostics. `diagnostics` alone returns the
findings with severity, P1 atom indices and a suggested verb. The codes and
thresholds are in [the algorithm page](../algorithms/agent#4-diagnostics).

Read the errors first. `overlapping_atoms` and `empty_structure` mean the
model or the conversion is wrong; measuring it is wasted work.

## Editing

Every edit is a command with `Origin::Agent`. `changes` groups consecutive
agent edits into bursts; when the burst ends the log, `undo_steps` says how
many `undo` calls remove it. Prefer one meaningful edit over many small ones
so the review unit is useful.

## Pictures

`render` draws the scene to SVG or PNG with a small orthographic renderer. It
is enough to see a folded linker, a misplaced cluster or a wrong cell. It is
not a publication renderer: no transparency, no ambient occlusion, no labels.
The SVG is searchable and small; the PNG is what most chat clients display.

## Safety

- The session runs under the standard `Limits` profile. `supercell`,
  `porosity` and `render` check their estimates before allocating and refuse
  with `limit_exceeded`.
- `optimize` is bounded by `max_steps` (default 200).
- Verbs are synchronous. The MCP server accepts `notifications/cancelled`, but
  a running verb cannot be preempted; a client that must stop one terminates
  the process. This is recorded as a deliberate deviation from the reference.
- The server has the file access of the user that starts it. Run it as the
  account that should own the files.

## Window proxy is out of scope

Crystal-Builder can forward MCP calls to a running GUI window through an
`mcp.json` discovery file. xtalforge has no window, so there is nothing to
proxy; there is no `proxy`/`discovery` surface. If a front end is built on
xtalforge and wants the same behaviour, the reference format is a JSON object
with the server URL and the process id:

```json
{"url": "http://127.0.0.1:8765/mcp", "pid": 12345}
```

The port is chosen by the front end; the file lives next to the front end's
configuration. xtalforge does not read it.
