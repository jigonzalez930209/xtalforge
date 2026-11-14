# Workflows and the `run` command

`xtalforge_agent::workflow` provides `Workflow`/`Step` (model version 1),
`Workflow::from_json`/`from_yaml`/`from_text`/`validate`,
`schema_json`/`schema_text`,
`expand_inputs`, `run_workflow` with `RunOptions`/`WorkflowReport`, and the
`WorkflowCheckpoint` resume format. The same model backs the CLI
`xtalforge run <workflow> --inputs ... --jobs ...`, the agent verbs
`workflow_schema`, `workflow_validate` and `workflow_run`, and the MCP tools
with the same names.

The workflow steps are the agent verbs, so a workflow and an agent session
share one vocabulary. The JSON Schema is returned by `workflow schema` and is
also available at `docs/public/workflow.schema.json`.

## YAML

YAML is a second syntax for the same model; JSON stays canonical. The parser
is a small pure-Rust subset: block mappings and sequences, compact
`- name: value` entries, flow lists/maps of simple values, quoted and plain
scalars, comments, `null`/booleans/numbers and an optional leading `---`.
Anchors, aliases, tags, merge keys, block scalars, multiline plain scalars and
multi-document files are rejected with a line and column. `.yaml`/`.yml`
paths select YAML; other paths try JSON first and fall back to YAML. The JSON
Schema marks YAML as accepted.

## In-memory execution

`run_workflow_memory(&workflow, &files)` takes a `BTreeMap<String, String>` of
input paths to text and returns `MemoryWorkflowReport { report, outputs }`
without touching a filesystem. The WASM export `run_workflow_memory` wraps it,
so the browser runs the same steps and returns written documents as a map.
Checkpoints are not produced because there is no disk to resume from.

Parity key: `workflow`.
