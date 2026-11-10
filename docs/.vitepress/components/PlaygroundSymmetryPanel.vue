<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import type { Inspection } from '../playground/wasm'

const props = defineProps<{
  wasm: any
  text: string
  inspection: Inspection | null
  zone: any | null
  brillouin: boolean
  busy: boolean
}>()

const emit = defineEmits<{
  (event: 'toggle-brillouin'): void
}>()

const kpath = ref<any | null>(null)
const ladder = ref<any[]>([])
const diagnostics = ref<any | null>(null)
const symbol = ref('')
const error = ref('')

watch(
  () => props.inspection?.symmetry?.hm_symbol,
  (value) => {
    if (value && !symbol.value) symbol.value = value
  },
  { immediate: true },
)

const wyckoffRows = computed(() => {
  const seen = new Set<string>()
  const rows: any[] = []
  for (const entry of props.inspection?.wyckoff ?? []) {
    const key = `${entry.label}|${entry.letter}|${entry.multiplicity}`
    if (seen.has(key)) continue
    seen.add(key)
    rows.push(entry)
  }
  return rows
})

function run(action: () => unknown) {
  try {
    action()
    error.value = ''
  } catch (caught) {
    const value = caught as { code?: string; message?: string }
    error.value = value.code && value.message ? `${value.code}: ${value.message}` : String(caught)
  }
}

function showKpath() {
  run(() => {
    kpath.value = props.wasm.kpath_info(props.text)
  })
}

function showLadder() {
  run(() => {
    ladder.value = props.wasm.tolerance_ladder_info(props.text)
  })
}

function diagnose() {
  run(() => {
    diagnostics.value = props.wasm.diagnose_info(props.text, symbol.value)
  })
}
</script>

<template>
  <section class="tab-section">
    <template v-if="inspection?.symmetry">
      <h3>Space group {{ inspection.symmetry.number }} — {{ inspection.symmetry.hm_symbol }}</h3>
      <dl class="summary-grid">
        <div><dt>Hall number</dt><dd>{{ inspection.symmetry.hall_number }}</dd></div>
        <div><dt>Pearson</dt><dd>{{ inspection.symmetry.pearson_symbol }}</dd></div>
        <div><dt>Operations</dt><dd>{{ inspection.symmetry.operations.length }}</dd></div>
        <div><dt>Orbits</dt><dd>{{ inspection.symmetry.orbits.length }}</dd></div>
        <div><dt>Wyckoff sites</dt><dd>{{ wyckoffRows.length }}</dd></div>
      </dl>
    </template>

    <details open>
      <summary>Wyckoff positions</summary>
      <table v-if="wyckoffRows.length">
        <thead><tr><th>site</th><th>letter</th><th>mult.</th><th>symmetry</th><th>coordinates</th></tr></thead>
        <tbody>
          <tr v-for="row in wyckoffRows" :key="`${row.label}-${row.letter}-${row.multiplicity}`">
            <td>{{ row.label }}</td>
            <td>{{ row.letter }}</td>
            <td>{{ row.multiplicity }}</td>
            <td>{{ row.site_symmetry }}</td>
            <td>{{ row.coordinates }}</td>
          </tr>
        </tbody>
      </table>
      <p v-else class="hint">No Wyckoff assignment for this structure.</p>
    </details>

    <details>
      <summary>Brillouin zone</summary>
      <p v-if="zone" class="hint">{{ zone.vertices.length }} vertices · {{ zone.faces.length }} faces · {{ zone.edges.length }} edges</p>
      <button :disabled="busy" @click="emit('toggle-brillouin')">
        {{ brillouin ? 'Back to the structure' : 'Show Brillouin zone' }}
      </button>
      <p class="hint">The zone is drawn in the same canvas; orbit and zoom still work.</p>
    </details>

    <details>
      <summary>Band path</summary>
      <button :disabled="busy" @click="showKpath">Compute k-path</button>
      <template v-if="kpath">
        <p class="hint">{{ kpath.bravais }} · {{ kpath.path.length }} segments</p>
        <p class="mono">{{ kpath.path.map((pair: string[]) => pair.join('→')).join('  ') }}</p>
        <table>
          <thead><tr><th>label</th><th>x</th><th>y</th><th>z</th><th>distance</th></tr></thead>
          <tbody>
            <tr v-for="vertex in kpath.vertices" :key="vertex.label">
              <td>{{ vertex.label }}</td>
              <td v-for="(value, index) in vertex.cartesian" :key="index">{{ Number(value).toFixed(3) }}</td>
              <td>{{ Number(vertex.distance).toFixed(3) }}</td>
            </tr>
          </tbody>
        </table>
      </template>
    </details>

    <details>
      <summary>Tolerance ladder</summary>
      <button :disabled="busy" @click="showLadder">Run ladder</button>
      <table v-if="ladder.length">
        <thead><tr><th>min</th><th>max</th><th>number</th><th>Hall</th><th>symbol</th></tr></thead>
        <tbody>
          <tr v-for="row in ladder" :key="`${row.min}-${row.max}`">
            <td>{{ Number(row.min).toFixed(4) }}</td>
            <td>{{ Number(row.max).toFixed(4) }}</td>
            <td>{{ row.number }}</td>
            <td>{{ row.hall_number }}</td>
            <td>{{ row.hm_symbol }}</td>
          </tr>
        </tbody>
      </table>
    </details>

    <details>
      <summary>Diagnose against a group</summary>
      <div class="row">
        <label>symbol<input v-model="symbol" class="wide" placeholder="F d -3 m" /></label>
        <button :disabled="busy" @click="diagnose">Diagnose</button>
      </div>
      <template v-if="diagnostics">
        <p class="mono">max {{ Number(diagnostics.max_displacement).toFixed(4) }} Å · mean {{ Number(diagnostics.mean_displacement).toFixed(4) }} Å</p>
        <table>
          <thead><tr><th>atom</th><th>element</th><th>displacement</th></tr></thead>
          <tbody>
            <tr v-for="atom in diagnostics.atoms.slice(0, 10)" :key="atom.atom">
              <td>{{ atom.atom }}</td><td>{{ atom.element }}</td><td>{{ Number(atom.displacement).toFixed(4) }}</td>
            </tr>
          </tbody>
        </table>
      </template>
    </details>
    <p v-if="error" class="error">{{ error }}</p>
  </section>
</template>
