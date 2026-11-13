# Reports and module registry

A run report is presentation over results that were already computed. The
Markdown renderer is the original surface; Phase 14 adds HTML and JSON.

## Data model

`Report` keeps the run id, title, summary, ordered sections and references.
`to_markdown`, `to_html` and `to_json` are pure renderers: they never change a
number and they escape user text.

- Markdown: headings and plain paragraphs, one trailing newline.
- HTML: a standalone fragment with `<h1>`, `<h2>`, `<ul>` and HTML-escaped
  text. No script, style or external asset is emitted.
- JSON: the serde representation of the same struct.

Escaping rules are the standard five (`&`, `<`, `>`, `"`, `'`); a section body
is treated as pre-escaped Markdown for HTML only after escaping, never as raw
HTML.

## Module registry listing

`ModuleRegistry::listing()` returns a stable, sorted list of
`ModuleInfo { name, description, parameters, result }`. The CLI `modules`
command prints that list as text or JSON; the agent capability report uses the
same data. Registration order never affects output.

## Validation

- The three renderers agree on the title, summary and section headings.
- HTML escapes `<script>` and attribute-breaking quotes.
- JSON round trips through `serde_json`.
- The registry listing is sorted and duplicate names are rejected.
