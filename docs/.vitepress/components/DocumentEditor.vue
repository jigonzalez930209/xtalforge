<script setup>
import { computed, onMounted, ref } from 'vue'

const SILICON = `data_Si
_cell_length_a 5.430941
_cell_length_b 5.430941
_cell_length_c 5.430941
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_symmetry_equiv_pos_as_xyz
'x,y,z'
'-x,-y,z'
'-x,y,-z'
'x,-y,-z'
'y,x,-z'
'-y,-x,-z'
'y,-x,z'
'-y,x,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 Si 0 0 0
`

const wasm = ref(null)
const version = ref('')
const text = ref(SILICON)
const document = ref(null)
const summary = ref(null)
const history = ref([])
const canUndo = ref(false)
const canRedo = ref(false)
const dirty = ref(false)
const delta = ref(0.05)
const element = ref('Ge')
const error = ref('')
const ready = computed(() => wasm.value !== null)

onMounted(async () => {
  const module = await import('xtalforge-wasm')
  await module.default()
  wasm.value = module
  version.value = module.version()
  reset()
})

function sync() {
  if (!document.value) return
  summary.value = document.value.summary()
  history.value = document.value.history
  canUndo.value = document.value.can_undo
  canRedo.value = document.value.can_redo
  dirty.value = document.value.dirty
}

function reset() {
  if (!wasm.value) return
  try {
    document.value = new wasm.value.WasmDocument(text.value)
    sync()
    error.value = ''
  } catch (caught) {
    document.value = null
    error.value = `${caught.code ?? 'error'}: ${caught.message}`
  }
}

function move() {
  try {
    document.value.move_sites([0], Number(delta.value), 0, 0)
    sync()
    error.value = ''
  } catch (caught) {
    error.value = `${caught.code ?? 'error'}: ${caught.message}`
  }
}

function setElement() {
  try {
    document.value.set_element(0, element.value)
    sync()
    error.value = ''
  } catch (caught) {
    error.value = `${caught.code ?? 'error'}: ${caught.message}`
  }
}

function undo() {
  document.value.undo()
  sync()
}

function redo() {
  document.value.redo()
  sync()
}

function download() {
  const contents = document.value.to_cif()
  const blob = new Blob([contents], { type: 'chemical/x-cif' })
  const link = document.createElement('a')
  link.href = URL.createObjectURL(blob)
  link.download = 'edited.cif'
  link.click()
  URL.revokeObjectURL(link.href)
}
</script>

<template>
  <div class="document-editor">
    <p v-if="ready" class="status">WASM {{ version }}</p>
    <p v-else class="status">Loading the WebAssembly module…</p>
    <textarea v-model="text" spellcheck="false" rows="12" aria-label="CIF text" />
    <button :disabled="!ready" @click="reset">Load structure</button>

    <p v-if="error" class="error">{{ error }}</p>

    <section v-if="summary" class="summary">
      <h3>{{ summary.name }} — {{ summary.formula }}</h3>
      <p>
        {{ summary.n_sites }} sites, volume {{ summary.volume.toFixed(3) }} Å³,
        dirty: {{ dirty }}
      </p>
      <div class="row">
        <button :disabled="!canUndo" @click="undo">Undo</button>
        <button :disabled="!canRedo" @click="redo">Redo</button>
        <button @click="download">Download CIF</button>
      </div>
      <div class="row">
        <label>Δx (fractional)<input v-model.number="delta" type="number" step="0.01" /></label>
        <button @click="move">Move atom 0</button>
        <label>Element<input v-model="element" type="text" /></label>
        <button @click="setElement">Set atom 0</button>
      </div>
      <details>
        <summary>Command history ({{ history.length }})</summary>
        <ol>
          <li v-for="(entry, index) in history" :key="index">{{ entry }}</li>
        </ol>
      </details>
    </section>
  </div>
</template>

<style scoped>
.document-editor {
  border: 1px solid var(--vp-c-divider);
  border-radius: 8px;
  padding: 1rem;
}
.status {
  color: var(--vp-c-text-2);
  font-size: 0.9rem;
}
textarea {
  width: 100%;
  font-family: var(--vp-font-family-mono);
  font-size: 0.8rem;
  border: 1px solid var(--vp-c-divider);
  border-radius: 4px;
  padding: 0.5rem;
  background: var(--vp-c-bg);
  color: var(--vp-c-text-1);
}
button {
  border: 1px solid var(--vp-c-brand-1);
  border-radius: 4px;
  padding: 0.3rem 0.6rem;
  color: var(--vp-c-brand-1);
}
button:disabled {
  opacity: 0.4;
}
.row {
  display: flex;
  flex-wrap: wrap;
  gap: 0.75rem;
  align-items: end;
  margin-top: 0.6rem;
}
.row label {
  display: flex;
  flex-direction: column;
  gap: 0.2rem;
  font-size: 0.85rem;
}
.row input {
  width: 6rem;
  border: 1px solid var(--vp-c-divider);
  border-radius: 4px;
  padding: 0.3rem;
  background: var(--vp-c-bg);
  color: var(--vp-c-text-1);
}
.error {
  color: var(--vp-c-danger-1);
}
</style>
