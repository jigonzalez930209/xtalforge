# Reference: workflows

A workflow is a versioned JSON pipeline applied to every input structure. It is
the shortest way to run the same operations over thousands of files, and the
CLI (`xtalforge run`), the Rust API (`xtalforge_agent::workflow`) and the MCP
tools use the same runner.

## Shape

```json
{
  "workflow_version": 1,
  "name": "pxrd-batch",
  "description": "Supercell, pattern and energy per input",
  "steps": [
    {"kind": "read", "path": "{{input}}"},
    {"kind": "transform", "name": "cell", "verb": "supercell",
     "arguments": {"repeats": [2, 2, 2]}},
    {"kind": "analyse", "name": "pattern", "verb": "pxrd",
     "arguments": {"two_theta_max": 60.0}},
    {"kind": "compute", "name": "energy", "verb": "uff_energy"},
    {"kind": "write", "path": "out/{{stem}}.cif"}
  ]
}
```

Rules:

- `workflow_version` is `1`; other versions are a typed `workflow_version`
  error.
- The first step is always `read`; its path defaults to `{{input}}`.
- `transform` runs a mutating verb, `analyse` and `compute` run a read-only
  verb. The lists come from `capabilities`; session plumbing (`load`, `save`,
  `render`, `undo`, `history`, ...) is not a step.
- `write` needs a path; `read` and `write` accept a `format`.
- Step labels default to `step-1`, `step-2`, ... and must be unique.
- Placeholders: `{{input}}`, `{{stem}}`, `{{ext}}`, `{{dir}}` and `{{index}}`
  (1-based position in the sorted input list).

The published schema is `xtalforge run --schema` and the `workflow_schema`
tool. Validation returns the first problem with its JSON path, for example
`steps[1].verb`.

## Running

```bash
xtalforge run wf.json --inputs "structures/*.cif" --jobs 8 --json
xtalforge run wf.json --inputs structures/ --checkpoint run.checkpoint.json
xtalforge run wf.json --validate
```

Inputs may be files, directories (their direct files) or patterns with `*`,
`?` and `**`. The runner sorts and deduplicates them, gives every item a fresh
session and keeps the report in that order. A failed item stops only itself.

The checkpoint records the workflow fingerprint and every completed item. Run
the same command again to skip the `ok` items and retry the failed ones; use
`--fresh` to start over. A checkpoint from another workflow is refused.

Over MCP, call `workflow_validate`, then `workflow_run` with `inputs` (an
array), `jobs`, an optional `checkpoint` and `resume`. The answer kind is
`workflow_run` and contains the same report as the CLI.

## Limitations

The document is JSON. A `.yaml` file loads when it contains JSON (JSON is a
subset of YAML); other YAML syntax is not supported by this build. Python and
WebAssembly do not expose the runner yet, so a batch with filesystem access
runs through the Rust API, the CLI or the MCP server.
