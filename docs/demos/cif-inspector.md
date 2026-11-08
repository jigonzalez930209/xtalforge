# CIF inspector demo

Paste or edit a CIF and inspect it entirely in the browser through the
WebAssembly package: cell parameters, formula, density, the perceived bond
list, primitive-ring statistics, and downloadable supercells and slabs. The
demo never uploads the file.

<script setup>
import CifInspector from '../.vitepress/components/CifInspector.vue'
</script>

<CifInspector />

The default structure is silicon. `read_cif`, `bonds` and `rings` are the same
functions the Rust API, the CLI and Python expose; `supercell_text` and
`slab_text` build the result and write it as CIF for the download button.
