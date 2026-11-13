# Modules and workspace

`Module` is the registry contract: `info()` returns `ModuleInfo { name,
description, version, parameters }`. `ModuleRegistry` registers, unregisters,
gets and lists modules; duplicate or empty names are typed errors.

`Workspace::open(root)` creates the directory. Each run lives in
`runs/<id>/run.json` as a `RunRecord { id, module, input, status, outputs, data,
created }`; `write_run`, `read_run` and `runs` are the on-disk API. The browser
has no filesystem, so `archive_records` serializes the records into a ZIP and
`read_archived_records` reads it back (WASM `workspace_archive` and
`workspace_read_archive`).

`RunStatus` is `Queued`, `Running`, `Done`, `Failed` or `Cancelled`;
`Report { run_id, title, summary, sections, references }` with
`ReportSection { heading, body }` turns a run into Markdown with `to_markdown`;
`TestModule`/`test_module()` is the trivial module used by the registry tests;
`JobStatus` is the job-queue state.

## Jobs

`JobQueue` is a FIFO of `Job { id, module, parameters, status }`.
`submit` returns the id; `start_next`, `finish`, `fail` and `cancel` move the
state through `Queued`, `Running`, `Done`, `Failed` and `Cancelled`.

## Reports and listings

`Report::to_markdown`, `to_html` and `to_json` render the same record; HTML
escapes every dynamic string and shows section bodies inside `<pre>` so the
data model never becomes executable markup. `Report::from_json` round trips the
serde form.

`ModuleRegistry::listing()` returns sorted `ModuleListing { name, description,
version, parameters, available }`; `listing_json()` serializes it. The
`xtalforge modules` command prints the built-in registry (currently `version`)
as text or JSON.

Parity key: `report_renderers`.
