<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import type { Inspection, SceneStats } from '../playground/wasm'

const props = defineProps<{
  wasm: any
  text: string
  inspection: Inspection | null
  scene: SceneStats | null
  busy: boolean
  label: string
}>()

const emit = defineEmits<{
  (event: 'operation', name: string, params?: Record<string, any>): void
  (event: 'apply', text: string, name: string): void
  (event: 'restore'): void
}>()

const na = ref(2)
const nb = ref(2)
const nc = ref(2)
const h = ref(1)
const k = ref(0)
const l = ref(0)
const layers = ref(2)
const vacuum = ref(10)
const nx = ref(4)
const ny = ref(4)
const tubeN = ref(5)
const tubeM = ref(5)
const tubeLength = ref(2)
const ribbonWidth = ref(6)
const ribbonLength = ref(6)
const twistN = ref(1)
const twistM = ref(1)
const twistAngle = ref(21.79)
const twistLayers = ref(2)
const ztcSamples = ref(16)
const ztcCell = ref(12)
const ztcIso = ref(0)
const carbonNet = ref('dia')
const carbonDensity = ref(0.42)
const carbonSeed = ref(0)
const polymerSmiles = ref('C(c1ccccc1)C')
const polymerMonomers = ref(4)
const polymerDensity = ref(1.0)
const topology = ref('pcu')
const node = ref('Zn4O')
const linker = ref('benzene')
const topologies = ref<string[]>([])
const nodes = ref<string[]>([])
const linkers = ref<string[]>([])
const draft = ref(props.text)
const showText = ref(false)

watch(
  () => props.text,
  (value) => {
    draft.value = value
  },
)
watch(
  () => props.wasm,
  (wasm) => {
    if (!wasm) return
    const blocks = wasm.mof_blocks_info()
    topologies.value = wasm.mof_topologies_info().map((item: any) => item.name)
    nodes.value = blocks
      .filter((item: any) => item.kind === 'Node')
      .map((item: any) => item.name)
    linkers.value = blocks
      .filter((item: any) => item.kind === 'Linker')
      .map((item: any) => item.name)
    topology.value = topologies.value[0] ?? topology.value
    node.value = nodes.value[0] ?? node.value
    linker.value = linkers.value[0] ?? linker.value
  },
  { immediate: true },
)

const summary = computed(() => props.inspection?.summary)
const ringText = computed(
  () =>
    props.inspection?.ringStats
      .map((entry) => `${entry.count} × ${entry.size}-ring`)
      .join(', ') || 'none',
)
</script>

<template>
  <section class="tab-section">
    <template v-if="summary">
      <h3>{{ label }} — {{ summary.formula }}</h3>
      <dl class="summary-grid">
        <div><dt>Cell volume</dt><dd>{{ summary.volume.toFixed(3) }} Å³</dd></div>
        <div><dt>Density</dt><dd>{{ summary.density.toFixed(4) }} g/cm³</dd></div>
        <div><dt>Sites</dt><dd>{{ summary.n_sites }} → {{ summary.n_atoms }} P1 atoms</dd></div>
        <div><dt>Scene</dt><dd>{{ scene?.atoms ?? 0 }} atoms · {{ scene?.bonds ?? 0 }} bonds</dd></div>
        <div><dt>Bonds</dt><dd>{{ inspection?.bondsCount ?? 0 }}</dd></div>
        <div><dt>Rings ≤ 8</dt><dd>{{ ringText }}</dd></div>
      </dl>
      <div class="row">
        <button :disabled="busy" @click="emit('restore')">Restore input</button>
        <button :disabled="busy" @click="emit('operation', 'p1')">P1</button>
        <button :disabled="busy" @click="emit('operation', 'asymmetrize')">Asymmetric unit</button>
        <button :disabled="busy" @click="emit('operation', 'primitive')">Primitive</button>
        <button :disabled="busy" @click="emit('operation', 'standardize')">Standardize</button>
        <button :disabled="busy" @click="emit('operation', 'reduce')">Niggli</button>
      </div>
    </template>

    <details open>
      <summary>Supercell and slab</summary>
      <div class="row">
        <label>a<input v-model.number="na" type="number" min="1" max="8" /></label>
        <label>b<input v-model.number="nb" type="number" min="1" max="8" /></label>
        <label>c<input v-model.number="nc" type="number" min="1" max="8" /></label>
        <button :disabled="busy" @click="emit('operation', 'supercell', { na, nb, nc })">Supercell</button>
      </div>
      <div class="row">
        <label>h<input v-model.number="h" type="number" /></label>
        <label>k<input v-model.number="k" type="number" /></label>
        <label>l<input v-model.number="l" type="number" /></label>
        <label>layers<input v-model.number="layers" type="number" min="1" /></label>
        <label>vacuum<input v-model.number="vacuum" type="number" min="0" /></label>
        <button :disabled="busy" @click="emit('operation', 'slab', { h, k, l, layers, vacuum })">Slab</button>
      </div>
    </details>

    <details>
      <summary>Carbon builders</summary>
      <div class="row">
        <label>nx<input v-model.number="nx" type="number" min="1" /></label>
        <label>ny<input v-model.number="ny" type="number" min="1" /></label>
        <button :disabled="busy" @click="emit('operation', 'graphene', { nx, ny })">Graphene</button>
        <button :disabled="busy" @click="emit('operation', 'c60')">C60</button>
        <button :disabled="busy" @click="emit('operation', 'c20')">C20</button>
      </div>
      <div class="row">
        <label>n<input v-model.number="tubeN" type="number" /></label>
        <label>m<input v-model.number="tubeM" type="number" /></label>
        <label>length<input v-model.number="tubeLength" type="number" min="1" /></label>
        <button :disabled="busy" @click="emit('operation', 'nanotube', { n: tubeN, m: tubeM, length: tubeLength })">Nanotube</button>
      </div>
      <div class="row">
        <label>width<input v-model.number="ribbonWidth" type="number" min="1" /></label>
        <label>length<input v-model.number="ribbonLength" type="number" min="1" /></label>
        <button :disabled="busy" @click="emit('operation', 'ribbon', { width: ribbonWidth, length: ribbonLength })">Nanoribbon</button>
      </div>
      <div class="row">
        <label>angle°<input v-model.number="twistAngle" type="number" step="0.01" /></label>
        <label>repeat a<input v-model.number="twistN" type="number" min="1" /></label>
        <label>repeat b<input v-model.number="twistM" type="number" min="1" /></label>
        <label>layers<input v-model.number="twistLayers" type="number" min="1" /></label>
        <button :disabled="busy" @click="emit('operation', 'twist', { n: twistN, m: twistM, theta: twistAngle, layers: twistLayers })">Twist</button>
      </div>
      <div class="row">
        <label>samples<input v-model.number="ztcSamples" type="number" min="4" /></label>
        <label>cell<input v-model.number="ztcCell" type="number" step="0.5" /></label>
        <label>iso<input v-model.number="ztcIso" type="number" step="0.1" /></label>
        <button :disabled="busy" @click="emit('operation', 'ztc', { samples: ztcSamples, cell: ztcCell, isovalue: ztcIso })">ZTC</button>
        <button :disabled="busy" @click="emit('operation', 'schwarzite', { samples: ztcSamples, cell: ztcCell, isovalue: ztcIso })">Schwarz P</button>
      </div>
      <p class="hint">Commensurate moiré cells: 21.79°, 13.17°, 9.43°, 7.34°, 6.01°…; other angles use the closest cell.</p>
      <div class="row">
        <label>net<select v-model="carbonNet"><option>dia</option><option>pcu</option></select></label>
        <label>density<input v-model.number="carbonDensity" type="number" step="0.02" /></label>
        <label>seed<input v-model.number="carbonSeed" type="number" min="0" /></label>
        <button :disabled="busy" @click="emit('operation', 'carbon', { net: carbonNet, repeat: 1, density: carbonDensity, seed: carbonSeed })">Disordered carbon</button>
      </div>
      <div class="row">
        <label>SMILES<input v-model="polymerSmiles" class="wide" /></label>
        <label>monomers<input v-model.number="polymerMonomers" type="number" min="1" /></label>
        <label>density<input v-model.number="polymerDensity" type="number" step="0.1" /></label>
        <button :disabled="busy" @click="emit('operation', 'polymer', { smiles: polymerSmiles, monomers: polymerMonomers, density: polymerDensity })">Polymer</button>
      </div>
    </details>

    <details>
      <summary>MOF builder</summary>
      <div class="row">
        <label>topology<select v-model="topology"><option v-for="item in topologies" :key="item">{{ item }}</option></select></label>
        <label>node<select v-model="node"><option v-for="item in nodes" :key="item">{{ item }}</option></select></label>
        <label>linker<select v-model="linker"><option v-for="item in linkers" :key="item">{{ item }}</option></select></label>
        <button :disabled="busy" @click="emit('operation', 'mof', { topology: topology, node: node, linker: linker })">Build MOF</button>
      </div>
      <p class="hint">The catalogue is read from the embedded MOF data; names come from the same list the CLI exposes.</p>
    </details>

    <details>
      <summary>Input text</summary>
      <button @click="showText = !showText">{{ showText ? 'Hide' : 'Show' }} CIF / POSCAR / XYZ</button>
      <textarea v-if="showText" v-model="draft" rows="12" spellcheck="false" aria-label="Structure text" />
      <button :disabled="busy || !showText" @click="emit('apply', draft, 'Pasted structure')">Apply text</button>
    </details>
  </section>
</template>
