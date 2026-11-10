# Reference: API and surfaces

| Surface | Entry point | Session |
|---|---|---|
| Rust | `xtalforge::agent` | `Session::run(verb, &json!({...}))` |
| CLI | `xtalforge capabilities`, `xtalforge agent script.json` | one process |
| Python | `xtalforge.AgentSession` | `run(verb, json_string)` |
| WASM | `AgentSession` | `run(verb, json_string)` |
| MCP | `xtalforge mcp` | one stdio connection |

The verb catalogue is `xtalforge::agent::catalogue()` in Rust, `agent_verbs()`
in Python and WASM, and `xtalforge capabilities` on the CLI. The runtime report
is `capability_report()`, `agent_capabilities()` or `xtalforge capabilities
--json`.

Arguments are one JSON object. Answers are one JSON object with `kind`. Both
are stable public API (see `AGENTS.md` section 5.2): renaming a verb, an
argument or a `kind` is a breaking change.

## Rust example

```rust
use serde_json::json;
use xtalforge::agent::Session;

let mut session = Session::new();
session.run("load", &json!({"path": "MOF-5.cif", "name": "mof5"}))?;
let answer = session.run("symmetry", &json!({}))?;
println!("{}", answer.to_json()?);
# Ok::<(), xtalforge::agent::Error>(())
```

## Error codes

`unknown_verb`, `invalid_argument`, `unknown_document`, `no_document`, `io`,
`json`, plus the codes of the core, force-field and analysis crates. The first
argument of a Python exception is the code; a WASM error is an object with
`code` and `message`; an MCP error is a text payload with both.

## Filesystem

The browser has no filesystem. `load` from `text` works, `save` and the `path`
form of `render` do not. Python and Rust read paths. The MCP server runs as a
process and has the same file access as the user.
