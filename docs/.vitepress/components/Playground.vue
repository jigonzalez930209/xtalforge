<script setup lang="ts">
import { computed, onBeforeUnmount, onMounted, ref } from 'vue'
import { SceneRenderer } from '../playground/scene'
import { SAMPLES, loadSampleText } from '../playground/samples'
import {
  errorText,
  inspect,
  loadWasm,
  sceneStats,
  type Inspection,
  type SceneStats,
} from '../playground/wasm'

const host = ref<HTMLElement | null>(null)
const ready = ref(false)
const busy = ref(false)
const error = ref('')
const text = ref('')
const baseText = ref('')
const label = ref('Silicon')
const baseLabel = ref('Silicon')
const selected = ref('Si')
const inspection = ref<Inspection | null>(null)
const scene = ref<SceneStats | null>(null)
const pxrd = ref<any | null>(null)
const style = ref('ball-stick')
const scale = ref(1)
const showCell = ref(true)
const showAxes = ref(true)
const autoRotate = ref(false)
const showText = ref(false)
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
const topology = ref('pcu')
const node = ref('Zn4O')
const linker = ref('benzene')
const topologies = ref<string[]>([])
const nodes = ref<string[]>([])
const linkers = ref<string[]>([])
let renderer: SceneRenderer | null = null
let wasm: any = null

const summary = computed(() => inspection.value?.summary)
const symmetry = computed(() => inspection.value?.symmetry)
const ringText = computed(
  () =>
    inspection.value?.ringStats
      .map((entry) => `${entry.count} × ${entry.size}-ring`)
      .join(', ') || 'none',
)
const wyckoffText = computed(
  () => inspection.value?.wyckoff.map((entry: any) => entry.letter).join(' ') || '—',
)
const pxrdPeaks = computed(() => (pxrd.value?.reflections ?? []).slice(0, 12))

function applyVisibility() {
  renderer?.setCellVisible(showCell.value)
  renderer?.setAxesVisible(showAxes.value)
  renderer?.setAutoRotate(autoRotate.value)
}

function renderScene() {
  if (!wasm || !renderer || !text.value) return
  try {
    const stats = sceneStats(wasm, text.value, style.value, Number(scale.value))
    scene.value = stats
    renderer.setScene(stats.data)
    applyVisibility()
    error.value = ''
  } catch (caught) {
    error.value = errorText(caught)
  }
}

async function loadText(next: string, name: string, keepBase = false) {
  busy.value = true
  error.value = ''
  try {
    text.value = next
    if (!keepBase) {
      baseText.value = next
      baseLabel.value = name
    }
    label.value = name
    inspection.value = inspect(wasm, next)
    if (inspection.value.error) throw new Error(inspection.value.error)
    pxrd.value = null
    renderScene()
  } catch (caught) {
    error.value = errorText(caught)
    inspection.value = null
    scene.value = null
  } finally {
    busy.value = false
  }
}

async function loadSample(id: string) {
  const sample = SAMPLES.find((item) => item.id === id)
  if (!sample) return
  selected.value = id
  busy.value = true
  try {
    const contents = await loadSampleText(sample)
    await loadText(contents, sample.label)
  } catch (caught) {
    error.value = errorText(caught)
  } finally {
    busy.value = false
  }
}

async function restore() {
  if (baseText.value) await loadText(baseText.value, baseLabel.value, true)
}

async function run(name: string, params: Record<string, number> = {}) {
  if (!wasm || !text.value) return
  busy.value = true
  error.value = ''
  try {
    let result = ''
    if (name === 'supercell') {
      result = wasm.supercell_text(text.value, params.na ?? na.value, params.nb ?? nb.value, params.nc ?? nc.value)
    } else if (name === 'slab') {
      result = wasm.slab_text(text.value, params.h ?? h.value, params.k ?? k.value, params.l ?? l.value, params.layers ?? layers.value, params.vacuum ?? vacuum.value)
    } else if (name === 'p1') {
      result = wasm.to_p1_text(text.value)
    } else if (name === 'asymmetrize') {
      result = wasm.asymmetrize_text(text.value, params.tolerance ?? 0.001)
    } else if (name === 'primitive') {
      result = wasm.primitive_text(text.value, params.tolerance ?? 0.001)
    } else if (name === 'standardize') {
      result = wasm.standardize_text(text.value, params.tolerance ?? 0.001)
    } else if (name === 'reduce') {
      result = wasm.reduce_text(text.value, 'niggli')
    } else if (name === 'graphene') {
      result = wasm.build_graphene_text(params.nx ?? nx.value, params.ny ?? ny.value, 1, 3.35, 15)
    } else if (name === 'nanotube') {
      result = wasm.build_nanotube_text(params.n ?? tubeN.value, params.m ?? tubeM.value, params.length ?? tubeLength.value, 10)
    } else if (name === 'c60') {
      result = wasm.build_c60_text()
    } else if (name === 'c20') {
      result = wasm.build_c20_text()
    } else if (name === 'mof') {
      result = wasm.build_mof_text(params.topology ?? topology.value, params.node ?? node.value, params.linker ?? linker.value)
    } else {
      throw new Error(`unknown operation ${name}`)
    }
    await loadText(result, name, true)
  } catch (caught) {
    error.value = errorText(caught)
  } finally {
    busy.value = false
  }
}

function calculatePxrd() {
  if (!wasm || !text.value) return
  try {
    pxrd.value = wasm.pxrd_info(text.value, 'cu-ka', 60, 0.02)
    error.value = ''
  } catch (caught) {
    error.value = errorText(caught)
  }
}

function download(name: string, contents: string, type: string) {
  const blob = new Blob([contents], { type })
  const link = document.createElement('a')
  link.href = URL.createObjectURL(blob)
  link.download = name
  link.click()
  URL.revokeObjectURL(link.href)
}

function downloadCif() {
  download(`${summary.value?.name ?? 'structure'}.cif`, text.value, 'chemical/x-cif')
}

function downloadPng() {
  if (!renderer) return
  const link = document.createElement('a')
  link.href = renderer.snapshot()
  link.download = `${summary.value?.name ?? 'structure'}.png`
  link.click()
}

function wyckoffOrbitMultiplicities(): number[] {
  const seen = new Set<string>()
  const values: number[] = []
  for (const entry of inspection.value?.wyckoff ?? []) {
    const key = `${entry.label}|${entry.letter}`
    if (seen.has(key)) continue
    seen.add(key)
    values.push(Number(entry.multiplicity))
  }
  return values
}

function state() {
  return {
    label: label.value,
    text: text.value,
    summary: summary.value,
    bonds: inspection.value?.bondsCount ?? 0,
    ringStats: inspection.value?.ringStats ?? [],
    symmetry: symmetry.value,
    wyckoff: inspection.value?.wyckoff.map((entry: any) => entry.letter) ?? [],
    wyckoffOrbits: wyckoffOrbitMultiplicities(),
    scene: scene.value
      ? {
          atoms: scene.value.atoms,
          bonds: scene.value.bonds,
          cell_edges: scene.value.cellEdges,
        }
      : null,
    atomPositions:
      scene.value?.data.atoms.map((atom) => ({
        element: atom.element,
        position: atom.position,
      })) ?? [],
    bondDistances:
      inspection.value?.bonds.map((bond: any) => Number(bond.distance)) ?? [],
    cellEdges: scene.value?.data.cell_edges ?? [],
    projected: renderer?.projectAtoms() ?? [],
    pxrd: pxrd.value,
    base: baseText.value,
  }
}

onMounted(async () => {
  if (!host.value) return
  renderer = new SceneRenderer(host.value)
  try {
    wasm = await loadWasm()
    const blocks = wasm.mof_blocks_info()
    topologies.value = wasm.mof_topologies_info().map((item: any) => item.name)
    nodes.value = blocks.filter((item: any) => item.kind === 'Node').map((item: any) => item.name)
    linkers.value = blocks.filter((item: any) => item.kind === 'Linker').map((item: any) => item.name)
    topology.value = topologies.value[0] ?? 'pcu'
    node.value = nodes.value[0] ?? 'Zn4O'
    linker.value = linkers.value[0] ?? 'benzene'
    ready.value = true
    const params = new URLSearchParams(window.location.search)
    await loadSample(params.get('sample') ?? 'Si')
  } catch (caught) {
    error.value = errorText(caught)
  }
  ;(window as any).__xf = {
    ready: true,
    loadText: async (contents: string, name = 'custom') => {
      await loadText(contents, name)
      return state()
    },
    loadSample: async (id: string) => {
      await loadSample(id)
      return state()
    },
    operation: async (name: string, params: Record<string, number> = {}) => {
      await run(name, params)
      return state()
    },
    setView: (options: { azimuth?: number; elevation?: number; zoom?: number }) => renderer?.setView(options),
    setStyle: (next: string, nextScale = 1) => {
      style.value = next
      scale.value = nextScale
      renderScene()
      return state()
    },
    setVisibility: (options: { cell?: boolean; axes?: boolean; rotate?: boolean }) => {
      if (options.cell !== undefined) showCell.value = options.cell
      if (options.axes !== undefined) showAxes.value = options.axes
      if (options.rotate !== undefined) autoRotate.value = options.rotate
      applyVisibility()
    },
    pxrd: (radiation = 'cu-ka') => {
      pxrd.value = wasm.pxrd_info(text.value, radiation, 60, 0.02)
      return pxrd.value
    },
    state,
    screenshot: () => renderer?.snapshot() ?? '',
  }
})

onBeforeUnmount(() => {
  renderer?.dispose()
  renderer = null
  delete (window as any).__xf
})
</script>

<template>
  <div class="playground">
    <div class="toolbar">
      <label>
        Sample
        <select v-model="selected" :disabled="!ready || busy" @change="loadSample(selected)">
          <option v-for="sample in SAMPLES" :key="sample.id" :value="sample.id">{{ sample.label }}</option>
        </select>
      </label>
      <button :disabled="!ready || busy" @click="restore">Restore input</button>
      <button :disabled="!ready || !text" @click="downloadCif">Download CIF</button>
      <span class="status">{{ ready ? (busy ? 'Working…' : label) : 'Loading the WebAssembly module…' }}</span>
    </div>

    <p v-if="error" class="error">{{ error }}</p>

    <div class="layout">
      <div class="viewer">
        <div ref="host" class="canvas-host" />
        <div class="view-controls">
          <label>
            Style
            <select v-model="style" :disabled="!ready" @change="renderScene">
              <option value="ball-stick">Ball and stick</option>
              <option value="space-fill">Space fill</option>
              <option value="wireframe">Wireframe</option>
              <option value="sticks">Sticks</option>
            </select>
          </label>
          <label>
            Scale
            <input v-model.number="scale" type="range" min="0.2" max="2.5" step="0.1" @change="renderScene" />
          </label>
          <label class="check"><input v-model="showCell" type="checkbox" @change="applyVisibility" /> Cell</label>
          <label class="check"><input v-model="showAxes" type="checkbox" @change="applyVisibility" /> Axes</label>
          <label class="check"><input v-model="autoRotate" type="checkbox" @change="applyVisibility" /> Rotate</label>
          <button :disabled="!ready" @click="renderer?.resetView()">Reset view</button>
          <button :disabled="!ready" @click="downloadPng">PNG</button>
        </div>
      </div>

      <aside class="panel">
        <section v-if="summary" class="summary">
          <h3>{{ summary.formula }}</h3>
          <dl>
            <div><dt>Cell volume</dt><dd>{{ summary.volume.toFixed(3) }} Å³</dd></div>
            <div><dt>Density</dt><dd>{{ summary.density.toFixed(4) }} g/cm³</dd></div>
            <div><dt>Sites</dt><dd>{{ summary.n_sites }} → {{ summary.n_atoms }} P1 atoms</dd></div>
            <div><dt>Scene</dt><dd>{{ scene?.atoms ?? 0 }} atoms · {{ scene?.bonds ?? 0 }} bonds</dd></div>
            <div><dt>Bonds</dt><dd>{{ inspection?.bondsCount ?? 0 }}</dd></div>
            <div><dt>Rings ≤ 8</dt><dd>{{ ringText }}</dd></div>
            <div><dt>Space group</dt><dd>{{ symmetry ? `${symmetry.number} ${symmetry.hm_symbol}` : '—' }}</dd></div>
            <div><dt>Hall</dt><dd>{{ symmetry?.hall_number ?? '—' }}</dd></div>
            <div><dt>Wyckoff</dt><dd>{{ wyckoffText }}</dd></div>
          </dl>
        </section>

        <details open>
          <summary>Transforms</summary>
          <div class="row">
            <label>a<input v-model.number="na" type="number" min="1" max="8" /></label>
            <label>b<input v-model.number="nb" type="number" min="1" max="8" /></label>
            <label>c<input v-model.number="nc" type="number" min="1" max="8" /></label>
            <button :disabled="busy" @click="run('supercell')">Supercell</button>
          </div>
          <div class="row">
            <label>h<input v-model.number="h" type="number" /></label>
            <label>k<input v-model.number="k" type="number" /></label>
            <label>l<input v-model.number="l" type="number" /></label>
            <label>layers<input v-model.number="layers" type="number" min="1" /></label>
            <label>vacuum<input v-model.number="vacuum" type="number" min="0" /></label>
            <button :disabled="busy" @click="run('slab')">Slab</button>
          </div>
          <div class="row buttons">
            <button :disabled="busy" @click="run('p1')">P1</button>
            <button :disabled="busy" @click="run('asymmetrize')">Asymmetric unit</button>
            <button :disabled="busy" @click="run('primitive')">Primitive</button>
            <button :disabled="busy" @click="run('standardize')">Standardize</button>
            <button :disabled="busy" @click="run('reduce')">Niggli</button>
          </div>
        </details>

        <details>
          <summary>Builders</summary>
          <div class="row">
            <label>nx<input v-model.number="nx" type="number" min="1" /></label>
            <label>ny<input v-model.number="ny" type="number" min="1" /></label>
            <button :disabled="busy" @click="run('graphene')">Graphene</button>
            <button :disabled="busy" @click="run('c60')">C60</button>
            <button :disabled="busy" @click="run('c20')">C20</button>
          </div>
          <div class="row">
            <label>n<input v-model.number="tubeN" type="number" /></label>
            <label>m<input v-model.number="tubeM" type="number" /></label>
            <label>length<input v-model.number="tubeLength" type="number" min="1" /></label>
            <button :disabled="busy" @click="run('nanotube')">Nanotube</button>
          </div>
          <div class="row">
            <label>topology<select v-model="topology"><option v-for="item in topologies" :key="item">{{ item }}</option></select></label>
            <label>node<select v-model="node"><option v-for="item in nodes" :key="item">{{ item }}</option></select></label>
            <label>linker<select v-model="linker"><option v-for="item in linkers" :key="item">{{ item }}</option></select></label>
            <button :disabled="busy" @click="run('mof')">MOF</button>
          </div>
        </details>

        <details>
          <summary>Powder pattern</summary>
          <button :disabled="busy" @click="calculatePxrd">Calculate Cu Kα PXRD</button>
          <table v-if="pxrdPeaks.length">
            <thead><tr><th>hkl</th><th>2θ</th><th>d</th><th>I</th></tr></thead>
            <tbody>
              <tr v-for="peak in pxrdPeaks" :key="peak.hkl.join(',')">
                <td>({{ peak.hkl.join(' ') }})</td>
                <td>{{ peak.two_theta.toFixed(3) }}</td>
                <td>{{ peak.d.toFixed(3) }}</td>
                <td>{{ peak.intensity.toFixed(1) }}</td>
              </tr>
            </tbody>
          </table>
        </details>

        <details>
          <summary>Edit the input</summary>
          <button @click="showText = !showText">{{ showText ? 'Hide' : 'Show' }} CIF / POSCAR / XYZ</button>
          <textarea v-if="showText" v-model="text" rows="12" spellcheck="false" aria-label="Structure text" />
          <button :disabled="busy || !showText" @click="loadText(text, 'Pasted structure')">Apply text</button>
        </details>
      </aside>
    </div>
  </div>
</template>

<style scoped src="../playground/playground.css"></style>
