<script setup lang="ts">
import type { SceneStats } from '../playground/wasm'

const props = defineProps<{
  summary: any
  scene: SceneStats | null
  label: string
  text: string
  busy: boolean
}>()

const emit = defineEmits<{
  (event: 'download', kind: string): void
}>()

async function copyCif() {
  try {
    await navigator.clipboard.writeText(props.text)
  } catch {
    // Clippy-free fallback is intentionally silent: the textarea still has
    // the CIF when the browser denies clipboard access.
  }
}
</script>

<template>
  <section class="tab-section">
    <h3>{{ label }}</h3>
    <dl class="summary-grid">
      <div><dt>Source length</dt><dd>{{ text.length }} chars</dd></div>
      <div><dt>Scene atoms</dt><dd>{{ scene?.atoms ?? 0 }}</dd></div>
      <div><dt>Scene bonds</dt><dd>{{ scene?.bonds ?? 0 }}</dd></div>
      <div><dt>Cell edges</dt><dd>{{ scene?.cellEdges ?? 0 }}</dd></div>
      <div v-if="summary"><dt>Formula</dt><dd>{{ summary.formula }}</dd></div>
      <div v-if="summary"><dt>Volume</dt><dd>{{ Number(summary.volume).toFixed(2) }} Å³</dd></div>
    </dl>
    <div class="row">
      <button :disabled="busy" @click="emit('download', 'cif')">CIF</button>
      <button :disabled="busy" @click="emit('download', 'poscar')">POSCAR</button>
      <button :disabled="busy" @click="emit('download', 'xyz')">XYZ</button>
      <button :disabled="busy" @click="emit('download', 'stl')">STL</button>
      <button :disabled="busy" @click="emit('download', 'png')">PNG</button>
      <button :disabled="busy" @click="emit('download', 'scene')">Scene JSON</button>
      <button :disabled="busy" @click="copyCif">Copy CIF</button>
    </div>
    <p class="hint">
      CIF, POSCAR and XYZ use the same writers as the CLI; STL is the printable
      mesh from the view crate; Scene JSON is the geometry currently drawn.
    </p>
  </section>
</template>
