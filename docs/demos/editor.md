# Structure editor demo

Edit a structure in the browser through the WebAssembly bindings. The
`WasmDocument` class keeps an undo/redo timeline, the command history and the
dirty flag, and writes the edited structure back as CIF. Nothing leaves the
page.

<script setup>
import DocumentEditor from '../.vitepress/components/DocumentEditor.vue'
</script>

<DocumentEditor />

The same `Document` exists in Rust (`xtalforge-commands`), Python
(`xtalforge.Document`) and the CLI (`xtalforge script`), so an edit recorded in
one place has identical semantics everywhere.
