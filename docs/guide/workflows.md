# Declarative workflows

A workflow is a versioned JSON (or YAML) document that runs the same steps
over every input, with deterministic output order and a checkpoint after each
item. JSON is canonical; the YAML subset is documented in
[the API page](../api/workflows.md).

```json
{
  "workflow_version": 1,
  "name": "inspect-core-mofs",
  "steps": [
    {"kind": "read", "path": "{{input}}"},
    {"kind": "analyse", "name": "info", "verb": "info"},
    {"kind": "analyse", "name": "symmetry", "verb": "symmetry",
     "arguments": {"symprec": 0.001}},
    {"kind": "write", "name": "copy", "path": "{{dir}}/out/{{stem}}.cif"}
  ]
}
```

```bash
xtalforge run inspect.json \
  --inputs 'tests/external/core_mof/*.cif' \
  --jobs 4 --checkpoint work/checkpoint.json --out work/report.json
```

`{{input}}` is the current path, `{{stem}}` the file name without extension,
`{{ext}}` the extension, `{{dir}}` its directory and `{{index}}` the sorted
input index. A step is `read`, `transform`, `analyse`, `compute` or `write`;
the verb names and arguments are the agent verbs documented in
[Automation](../guide/automation.md).

## Reproducibility

Inputs are expanded, sorted and deduplicated before execution. Each item gets a
fresh session, so one failure cannot leak into the next. The report lists every
item in input order and the process exits with code 3 if any item failed. A
checkpoint stores the completed item fingerprints; rerunning the same command
skips them and reports the resumed count.

## From Python

```python
import xtalforge

workflow = xtalforge.Workflow(open("inspect.json").read())
print(workflow.validate())
report = workflow.run(["a.cif", "b.cif"], jobs=2)
```

`Workflow.schema()` returns the published JSON Schema. The WASM surface runs
the same document with `run_workflow_memory`, taking input documents as a
path-to-text map and returning the written documents as a map; no filesystem
is needed.
