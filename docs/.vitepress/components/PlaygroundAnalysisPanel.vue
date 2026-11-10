<script setup lang="ts">
import { ref } from 'vue'
import type { Inspection } from '../playground/wasm'

interface Measure {
  kind: string
  value: number
  indices: number[]
}

const props = defineProps<{
  wasm: any
  text: string
  inspection: Inspection | null
  selection: number[]
  measure: Measure | null
  busy: boolean
}>()

const emit = defineEmits<{
  (event: 'apply', text: string, name: string): void
  (event: 'clear-selection'): void
}>()

const cutoff = ref(3.0)
const neighbors = ref<any[] | null>(null)
const spacing = ref(0.6)
const probe = ref(1.2)
const porosity = ref<any | null>(null)
const groups = ref<any[] | null>(null)
const hydrogenDistance = ref(1.0)
const hydrogenSummary = ref<any | null>(null)
const error = ref('')

function run(action: () => void) {
  try {
    action()
    error.value = ''
  } catch (caught) {
    const value = caught as { code?: string; message?: string }
    error.value = value.code && value.message ? `${value.code}: ${value.message}` : String(caught)
  }
}

function findNeighbors() {
  run(() => {
    neighbors.value = props.wasm.neighbors(props.text, Number(cutoff.value))
  })
}

function findGroups() {
  run(() => {
    groups.value = props.wasm.groups(props.text)
  })
}

function runPorosity() {
  run(() => {
    porosity.value = props.wasm.porosity_info(
      props.text,
      Number(spacing.value),
      Number(probe.value),
    )
  })
}

function addHydrogens() {
  run(() => {
    const filled = props.wasm.hydrogens_text(props.text, Number(hydrogenDistance.value))
    hydrogenSummary.value = props.wasm.read_cif(filled)
    emit('apply', filled, 'with hydrogens')
  })
}
</script>

<template>
  <section class="tab-section">
    <details open>
      <summary>Bonds and rings</summary>
      <dl class="summary-grid">
        <div><dt>Perceived bonds</dt><dd>{{ inspection?.bondsCount ?? 0 }}</dd></div>
        <div><dt>Ring summary</dt><dd>{{ inspection?.ringStats.map((entry: any) => `${entry.count} × ${entry.size}`).join(', ') || 'none' }}</dd></div>
        <div><dt>Total atoms</dt><dd>{{ inspection?.summary?.n_atoms ?? 0 }}</dd></div>
        <div><dt>Asymmetric sites</dt><dd>{{ inspection?.summary?.n_sites ?? 0 }}</dd></div>
      </dl>
      <button :disabled="busy" @click="findGroups">Find functional groups</button>
      <p v-if="groups" class="mono">
        {{ groups.length }} groups:
        {{ groups.slice(0, 12).map((group: any) => `${group.kind}(${group.atoms.join(',')})`).join(' ') }}
      </p>
    </details>

    <details>
      <summary>Periodic neighbours</summary>
      <div class="row">
        <label>cutoff (Å)<input v-model.number="cutoff" type="number" min="0.5" step="0.1" /></label>
        <button :disabled="busy" @click="findNeighbors">List neighbours</button>
      </div>
      <template v-if="neighbors">
        <p class="hint">{{ neighbors.length }} contacts (pairs with images).</p>
        <table>
          <thead><tr><th>i</th><th>j</th><th>image</th><th>distance</th></tr></thead>
          <tbody>
            <tr v-for="(item, index) in neighbors.slice(0, 10)" :key="index">
              <td>{{ item.i }}</td><td>{{ item.j }}</td><td>{{ item.image.join(',') }}</td>
              <td>{{ Number(item.distance).toFixed(4) }}</td>
            </tr>
          </tbody>
        </table>
      </template>
    </details>

    <details>
      <summary>Porosity</summary>
      <div class="row">
        <label>spacing (Å)<input v-model.number="spacing" type="number" min="0.2" step="0.1" /></label>
        <label>probe (Å)<input v-model.number="probe" type="number" min="0" step="0.1" /></label>
        <button :disabled="busy" @click="runPorosity">Run distance grid</button>
      </div>
      <dl v-if="porosity" class="summary-grid">
        <div><dt>Cell volume</dt><dd>{{ Number(porosity.cell_volume).toFixed(2) }} Å³</dd></div>
        <div><dt>Accessible</dt><dd>{{ Number(porosity.accessible_volume).toFixed(2) }} Å³</dd></div>
        <div><dt>Pore volume</dt><dd>{{ Number(porosity.pore_volume).toFixed(2) }} Å³</dd></div>
        <div><dt>Surface area</dt><dd>{{ Number(porosity.surface_area).toFixed(2) }} Å²</dd></div>
        <div><dt>Largest included</dt><dd>{{ Number(porosity.largest_included_diameter).toFixed(2) }} Å</dd></div>
      </dl>
      <p class="hint">A coarse grid is fast; refine the spacing for production numbers.</p>
    </details>

    <details>
      <summary>Hydrogens</summary>
      <div class="row">
        <label>X-H distance<input v-model.number="hydrogenDistance" type="number" min="0.5" step="0.05" /></label>
        <button :disabled="busy" @click="addHydrogens">Add hydrogens</button>
      </div>
      <p v-if="hydrogenSummary" class="mono">New formula: {{ hydrogenSummary.formula }} · {{ hydrogenSummary.n_atoms }} atoms</p>
    </details>

    <details open>
      <summary>Selection and measures</summary>
      <p class="hint">Click one, two, three or four atoms in the canvas.</p>
      <p class="mono">Selected: {{ selection.length ? selection.join(', ') : 'none' }}</p>
      <p v-if="measure" class="mono">{{ measure.kind }} = {{ measure.value.toFixed(4) }}</p>
      <button :disabled="busy || selection.length === 0" @click="emit('clear-selection')">Clear selection</button>
    </details>
    <p v-if="error" class="error">{{ error }}</p>
  </section>
</template>
