# Reference: documents and history

A session owns documents by name. `load` creates one, `close` removes it, and
`list_documents` shows the set with the active name.

## Commands and origin

Every edit is a `xtalforge_commands::Command` recorded with an
`Origin`:

- `user`: a person in a front end;
- `script`: the CLI or a script;
- `agent`: this layer.

The session always uses `agent`. The history is the document's undo stack; an
edit that is undone leaves the log, and a new edit after an undo discards the
redo branch.

## Change bursts

`changes` groups consecutive `agent` records. A burst has:

- `start`: position in the history;
- `len`: number of edits;
- `descriptions`: one line per edit;
- `undo_steps`: how many `undo` calls remove the burst, but only when it ends
  the log (a later user edit makes one undo insufficient).

Use one session and one document per task so the burst is a meaningful review
unit. At the end of a task, report the burst or undo it.

## Saving

`save` writes the active document in a registered format and marks the
document clean. It does not close it. A save does not create a history entry;
the structure did not change.
