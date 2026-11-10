# Modules, jobs, runs and reports

Reference: Crystal-Builder `xtal/modules/` and its workspace (MIT). The modules
crate is the bookkeeping layer around a calculation, not a calculator.

## Module registry

A `Module` declares a `ModuleInfo` (name, description, version, parameter
names) and a `run(parameters_json) -> result_json` method. The registry is a
`BTreeMap`, so listing is deterministic and a duplicate name is a typed error.
The default `run` refuses, which makes an unimplemented module explicit.

Modules are stateless by contract. A module that needs a structure receives it
through its parameters or through the job that supplies the input; the registry
does not hold documents.

## Jobs

`JobQueue` is a deterministic FIFO queue:

- `Job` has an id, module name, parameters and status;
- `JobStatus` is `Queued`, `Running`, `Done` or `Failed`;
- the queue keeps the order of submission and reports one job at a time;
- a failure carries a message and does not stop the queue.

The queue has no threads: a caller decides when to run the next job. This keeps
a run reproducible and keeps the WASM surface possible (a browser has no
threads).

## Workspace

`Workspace` is the on-disk layout:

```
runs/<id>/run.json      the RunRecord
runs/<id>/report.md     the rendered report, when present
runs/<id>/...           module-specific files
```

`RunRecord` stores the id, module, input summary, status, timestamps and the
files it wrote. `archive_records` packs the records into a zip;
`read_archived_records` reads them back. The CLI reads and writes `run.json`
directly; the archive is the portable form.

A record is written before the module runs and updated after it, so an
interrupted run leaves evidence instead of nothing. The status is the source of
truth; a report without a matching record is not trusted.

## Reports

A `Report` is a title and a list of `ReportSection` (heading + Markdown body).
`to_markdown` is the only renderer; the structured form serializes to JSON for
the bindings. A report contains no number that was not first written to the run
folder, so deleting the report never loses a result.

## Validation

- `crates/xtalforge-modules/tests/` covers registration duplicates, the FIFO
  order, failure propagation, the workspace layout, archive round-trips and
  Markdown/JSON report round-trips.
- The agent and the CLI use the same records; the parity matrix marks the
  in-memory workspace as `n/a` for the CLI because it reads `run.json` directly.

## Differences from Crystal-Builder

- No background thread pool; the caller drives the queue.
- A run record is written before execution, not only after it.
- The report is presentation only; the reference may compute while rendering.
