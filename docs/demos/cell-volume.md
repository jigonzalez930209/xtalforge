# Cell volume demo

The same `cell_volume` function as the Rust API, the CLI and Python, running in
the browser through the WebAssembly package. The demo works offline after the
page loads; it never calls a server.

<script setup>
import CellVolume from '../.vitepress/components/CellVolume.vue'
</script>

<CellVolume />

The volume follows

$$
V = abc\sqrt{1 - \cos^2\alpha - \cos^2\beta - \cos^2\gamma
  + 2\cos\alpha\cos\beta\cos\gamma}
$$

with lengths in Angstrom and angles in degrees. An invalid cell returns the
stable `invalid_input` code from the WASM module.
