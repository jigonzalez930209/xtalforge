# Workflows and the `run` command

`xtalforge_agent::workflow` provides `Workflow`/`Step` (model version 1),
`Workflow::from_json`/`validate`, `schema_json`/`schema_text`,
`expand_inputs`, `run_workflow` with `RunOptions`/`WorkflowReport`, and the
`WorkflowCheckpoint` resume format. The same model backs the CLI
`xtalforge run <workflow> --inputs ... --jobs ...`, the agent verbs
`workflow_schema`, `workflow_validate` and `workflow_run`, and the MCP tools
with the same names.

The workflow steps are the agent verbs, so a workflow and an agent session
share one vocabulary. The JSON Schema is returned by `workflow schema` and is
also available at `docs/public/workflow.schema.json`.

Parity key: `workflow`.
