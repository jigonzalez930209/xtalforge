<script setup>
import { computed, onMounted, ref } from 'vue'

const SILICON = `data_Si
_cell_length_a 5.430941
_cell_length_b 5.430941
_cell_length_c 5.430941
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_Hall 'F 4d 2 3 -1d'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 Si 0 0 0
`

const wasm = ref(null)
const wasmVersion = ref('')
const text = ref(SILICON)
const summary = ref(null)
const bonds = ref([])
const ringCounts = ref([])
const symmetry = ref(null)
const wyckoff = ref([])
const zone = ref(null)
const error = ref('')
const na = ref(2)
const nb = ref(2)
const nc = ref(2)
const h = ref(1)
const k = ref(1)
const l = ref(1)
const layers = ref(2)
const vacuum = ref(10)

const ready = computed(() => wasm.value !== null)
const ringSummary = computed(() =>
  ringCounts.value.map((entry) => `${entry[1]} x ${entry[0]}-ring`).join(', ') || 'none',
)

onMounted(async () => {
  const module = await import('xtalforge-wasm')
  await module.default()
  wasm.value = module
  wasmVersion.value = module.version()
  inspect()
})

function inspect() {
  if (!wasm.value) return
  try {
    summary.value = wasm.value.read_cif(text.value)
    bonds.value = wasm.value.bonds(text.value)
    const rings = wasm.value.rings(text.value, 8)
    const counts = new Map()
    for (const ring of rings) {
      counts.set(ring.size, (counts.get(ring.size) ?? 0) + 1)
    }
    ringCounts.value = [...counts.entries()].sort((a, b) => a[0] - b[0])
    symmetry.value = wasm.value.detect_symmetry(text.value, 0.001)
    wyckoff.value = wasm.value.wyckoff_info(text.value, 0.001)
    zone.value = wasm.value.brillouin_zone(text.value)
    error.value = ''
  } catch (caught) {
    summary.value = null
    error.value = `${caught.code ?? 'error'}: ${caught.message}`
  }
}

function download(name, contents) {
  const blob = new Blob([contents], { type: 'chemical/x-cif' })
  const link = document.createElement('a')
  link.href = URL.createObjectURL(blob)
  link.download = name
  link.click()
  URL.revokeObjectURL(link.href)
}

function downloadSupercell() {
  try {
    const contents = wasm.value.supercell_text(
      text.value,
      Number(na.value),
      Number(nb.value),
      Number(nc.value),
    )
    download(`${summary.value?.name ?? 'structure'}-${na.value}x${nb.value}x${nc.value}.cif`, contents)
    error.value = ''
  } catch (caught) {
    error.value = `${caught.code ?? 'error'}: ${caught.message}`
  }
}

function downloadSlab() {
  try {
    const contents = wasm.value.slab_text(
      text.value,
      Number(h.value),
      Number(k.value),
      Number(l.value),
      Number(layers.value),
      Number(vacuum.value),
    )
    download(`${summary.value?.name ?? 'structure'}-${h.value}${k.value}${l.value}.cif`, contents)
    error.value = ''
  } catch (caught) {
    error.value = `${caught.code ?? 'error'}: ${caught.message}`
  }
}

</script>

<template>
  <div class="cif-inspector">
    <p v-if="ready" class="status">WASM {{ wasmVersion }}</p>
    <p v-else class="status">Loading the WebAssembly module…</p>
    <textarea v-model="text" spellcheck="false" rows="14" aria-label="CIF text" />
    <button :disabled="!ready" @click="inspect">Inspect</button>

    <p v-if="error" class="error">{{ error }}</p>

    <section v-if="summary" class="summary">
      <h3>{{ summary.name }} — {{ summary.formula }}</h3>
      <dl>
        <div><dt>Cell volume</dt><dd>{{ summary.volume.toFixed(3) }} Å³</dd></div>
        <div><dt>Density</dt><dd>{{ summary.density.toFixed(4) }} g/cm³</dd></div>
        <div><dt>Asymmetric sites</dt><dd>{{ summary.n_sites }}</dd></div>
        <div><dt>P1 atoms</dt><dd>{{ summary.n_atoms }}</dd></div>
        <div><dt>Bonds</dt><dd>{{ bonds.length }}</dd></div>
        <div><dt>Rings (≤ 8)</dt><dd>{{ ringSummary }}</dd></div>
        <div><dt>Space group</dt><dd>{{ symmetry.number }} {{ symmetry.hm_symbol }}</dd></div>
        <div><dt>Hall number</dt><dd>{{ symmetry.hall_number }}</dd></div>
        <div><dt>Wyckoff</dt><dd>{{ wyckoff.map((item) => item.letter).join(' ') }}</dd></div>
        <div><dt>Brillouin zone</dt><dd>{{ zone.faces.length }} faces / {{ zone.vertices.length }} vertices</dd></div>
      </dl>
    </section>

    <section v-if="summary" class="builder">
      <h4>Supercell</h4>
      <div class="row">
        <label>a<input v-model.number="na" type="number" min="1" /></label>
        <label>b<input v-model.number="nb" type="number" min="1" /></label>
        <label>c<input v-model.number="nc" type="number" min="1" /></label>
        <button @click="downloadSupercell">Download CIF</button>
      </div>
      <h4>Slab</h4>
      <div class="row">
        <label>h<input v-model.number="h" type="number" /></label>
        <label>k<input v-model.number="k" type="number" /></label>
        <label>l<input v-model.number="l" type="number" /></label>
        <label>layers<input v-model.number="layers" type="number" min="1" /></label>
        <label>vacuum (Å)<input v-model.number="vacuum" type="number" min="0" /></label>
        <button @click="downloadSlab">Download CIF</button>
      </div>
    </section>

    <details>
      <summary>Bond list</summary>
      <table>
        <thead><tr><th>i</th><th>j</th><th>image</th><th>distance (Å)</th></tr></thead>
        <tbody>
          <tr v-for="(bond, index) in bonds" :key="index">
            <td>{{ bond.i }}</td>
            <td>{{ bond.j }}</td>
            <td>{{ bond.image.join(',') }}</td>
            <td>{{ bond.distance.toFixed(4) }}</td>
          </tr>
        </tbody>
      </table>
    </details>
  </div>
</template>

<style scoped>
.cif-inspector {
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
  margin-top: 0.75rem;
  border: 1px solid var(--vp-c-brand-1);
  border-radius: 4px;
  padding: 0.35rem 0.7rem;
  color: var(--vp-c-brand-1);
}
.summary dl {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
  gap: 0.5rem;
}
.summary dt {
  color: var(--vp-c-text-2);
  font-size: 0.8rem;
}
.summary dd {
  font-family: var(--vp-font-family-mono);
}
.builder {
  margin-top: 1rem;
}
.row {
  display: flex;
  flex-wrap: wrap;
  gap: 0.75rem;
  align-items: end;
}
.row label {
  display: flex;
  flex-direction: column;
  gap: 0.2rem;
  font-size: 0.85rem;
}
.row input {
  width: 5rem;
  border: 1px solid var(--vp-c-divider);
  border-radius: 4px;
  padding: 0.3rem;
  background: var(--vp-c-bg);
  color: var(--vp-c-text-1);
}
table {
  width: 100%;
  font-family: var(--vp-font-family-mono);
  font-size: 0.8rem;
  border-collapse: collapse;
}
th,
td {
  border-bottom: 1px solid var(--vp-c-divider);
  padding: 0.25rem;
  text-align: left;
}
.error {
  color: var(--vp-c-danger-1);
}
</style>
