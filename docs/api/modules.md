# Modules and workspace

`Module` is the registry contract: `info()` returns `ModuleInfo { name,
description, version, parameters }`. `ModuleRegistry` registers, unregisters,
gets and lists modules; duplicate or empty names are typed errors.

`Workspace::open(root)` creates the directory. Each run lives in
`runs/<id>/run.json` as a `RunRecord { id, module, input, status, outputs, data,
created }`; `write_run`, `read_run` and `runs` are the on-disk API. The browser
has no filesystem and will export the records as JSON/zip instead (Task 3.12.3).

## Jobs

`JobQueue` is a FIFO of `Job { id, module, parameters, status }`.
`submit` returns the id; `start_next`, `finish`, `fail` and `cancel` move the
state through `Queued`, `Running`, `Done`, `Failed` and `Cancelled`.
