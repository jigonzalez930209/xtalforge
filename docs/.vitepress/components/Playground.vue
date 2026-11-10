<script setup lang="ts">
import { computed, onBeforeUnmount, onMounted, ref, shallowRef } from 'vue'
import { SceneRenderer } from '../playground/scene'
import { SAMPLES, TABS, loadSampleText } from '../playground/samples'
import {
  errorText,
  inspect,
  loadWasm,
  sceneStats,
  type Inspection,
  type SceneStats,
} from '../playground/wasm'
import PlaygroundStructurePanel from './PlaygroundStructurePanel.vue'
import PlaygroundSymmetryPanel from './PlaygroundSymmetryPanel.vue'
import PlaygroundAnalysisPanel from './PlaygroundAnalysisPanel.vue'
import PlaygroundPowderPanel from './PlaygroundPowderPanel.vue'
import PlaygroundForcefieldPanel from './PlaygroundForcefieldPanel.vue'
import PlaygroundExportPanel from './PlaygroundExportPanel.vue'

type Measure = { kind: string; value: number; indices: number[] }
const host = ref<HTMLElement | null>(null)
const viewer = ref<HTMLElement | null>(null)
const wasmModule = shallowRef<any>(null)
const ready = ref(false)
const busy = ref(false)
const error = ref('')
const text = ref('')
const baseText = ref('')
const label = ref('Silicon')
const baseLabel = ref('Silicon')
const selected = ref('Si')
const tab = ref('structure')
const inspection = ref<Inspection | null>(null)
const scene = ref<SceneStats | null>(null)
const selection = ref<number[]>([])
const measure = ref<Measure | null>(null)
const zone = ref<any | null>(null)
const brillouin = ref(false)
const fullscreen = ref(false)
const style = ref('ball-stick')
const scale = ref(1)
const showCell = ref(true)
const showAxes = ref(true)
const autoRotate = ref(false)
let renderer: SceneRenderer | null = null

const summary = computed(() => inspection.value?.summary)

function applyVisibility() {
  renderer?.setCellVisible(showCell.value)
  renderer?.setAxesVisible(showAxes.value)
  renderer?.setAutoRotate(autoRotate.value)
}

function renderScene() {
  if (!wasmModule.value || !renderer || !text.value) return
  try {
    const stats = sceneStats(wasmModule.value, text.value, style.value, Number(scale.value))
    scene.value = stats
    renderer.setScene(stats.data)
    applyVisibility()
    error.value = ''
  } catch (caught) {
    error.value = errorText(caught)
  }
}

function resetInspection() {
  selection.value = []
  measure.value = null
  zone.value = null
  brillouin.value = false
  renderer?.setSelection([])
  renderer?.hideBrillouin()
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
    inspection.value = inspect(wasmModule.value, next)
    if (inspection.value.error) throw new Error(inspection.value.error)
    resetInspection()
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

async function run(name: string, params: Record<string, any> = {}) {
  if (!wasmModule.value || !text.value) return
  busy.value = true
  error.value = ''
  try {
    const wasm = wasmModule.value
    let result = ''
    if (name === 'supercell') {
      result = wasm.supercell_text(text.value, params.na ?? 2, params.nb ?? 2, params.nc ?? 2)
    } else if (name === 'slab') {
      result = wasm.slab_text(text.value, params.h ?? 1, params.k ?? 0, params.l ?? 0, params.layers ?? 2, params.vacuum ?? 10)
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
      result = wasm.build_graphene_text(params.nx ?? 4, params.ny ?? 4, 1, 3.35, 15)
    } else if (name === 'nanotube') {
      result = wasm.build_nanotube_text(params.n ?? 5, params.m ?? 5, params.length ?? 2, 10)
    } else if (name === 'c60') {
      result = wasm.build_c60_text()
    } else if (name === 'c20') {
      result = wasm.build_c20_text()
    } else if (name === 'ribbon') {
      result = wasm.build_ribbon_text(params.width ?? 6, params.length ?? 6, params.vacuum ?? 12)
    } else if (name === 'twist') {
      result = wasm.build_twist_text(params.n ?? 6, params.m ?? 0, params.theta ?? 21.8, params.layers ?? 2)
    } else if (name === 'ztc') {
      result = wasm.build_ztc_text(params.samples ?? 16, params.cell ?? 12, params.isovalue ?? 0)
    } else if (name === 'schwarzite') {
      result = wasm.build_schwarzite_text(params.samples ?? 16, params.cell ?? 12, params.isovalue ?? 0)
    } else if (name === 'polymer') {
      result = wasm.build_polymer_text(params.smiles ?? 'C=C', params.monomers ?? 8, params.density ?? 1.0)
    } else if (name === 'mof') {
      result = wasm.build_mof_text(params.topology ?? 'pcu', params.node ?? 'Zn4O', params.linker ?? 'benzene')
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

function toggleBrillouin() {
  if (!wasmModule.value || !renderer) return
  if (brillouin.value) {
    renderer.hideBrillouin()
    brillouin.value = false
    return
  }
  try {
    zone.value = wasmModule.value.brillouin_zone(text.value)
    renderer.showBrillouin(zone.value.vertices, zone.value.faces)
    brillouin.value = true
    error.value = ''
  } catch (caught) {
    error.value = errorText(caught)
  }
}

function updateMeasure() {
  renderer?.setSelection(selection.value)
  measure.value = null
  if (!wasmModule.value) return
  const wasm = wasmModule.value
  const [i, j, k, l] = selection.value
  try {
    if (selection.value.length === 2) {
      measure.value = {
        kind: 'distance (Å)',
        value: wasm.distance(text.value, i, j),
        indices: [...selection.value],
      }
    } else if (selection.value.length === 3) {
      measure.value = {
        kind: 'angle (deg)',
        value: wasm.angle(text.value, i, j, k),
        indices: [...selection.value],
      }
    } else if (selection.value.length === 4) {
      measure.value = {
        kind: 'torsion (deg)',
        value: wasm.torsion(text.value, i, j, k, l),
        indices: [...selection.value],
      }
    }
  } catch (caught) {
    error.value = errorText(caught)
  }
}

function onCanvasClick(event: MouseEvent) {
  if (!renderer || brillouin.value) return
  const index = renderer.pick(event.clientX, event.clientY)
  if (index === null) return
  const position = selection.value.indexOf(index)
  if (position >= 0) {
    selection.value.splice(position, 1)
  } else {
    if (selection.value.length >= 4) selection.value.shift()
    selection.value.push(index)
  }
  updateMeasure()
}

function selectAtoms(indices: number[]) {
  selection.value = indices.slice(0, 4)
  updateMeasure()
}

function download(name: string, contents: string, type: string) {
  const blob = new Blob([contents], { type })
  const link = document.createElement('a')
  link.href = URL.createObjectURL(blob)
  link.download = name
  link.click()
  URL.revokeObjectURL(link.href)
}

function exportAs(kind: string) {
  if (!wasmModule.value) return
  const wasm = wasmModule.value
  const base = summary.value?.name ?? 'structure'
  try {
    if (kind === 'cif') download(`${base}.cif`, text.value, 'chemical/x-cif')
    if (kind === 'poscar') download(`${base}.poscar`, wasm.convert_structure(text.value, '', 'poscar'), 'text/plain')
    if (kind === 'xyz') download(`${base}.xyz`, wasm.convert_structure(text.value, '', 'xyz'), 'text/plain')
    if (kind === 'stl') download(`${base}.stl`, wasm.mesh_stl(text.value, style.value, Number(scale.value)), 'model/stl')
    if (kind === 'scene') download(`${base}-scene.json`, JSON.stringify(scene.value?.data ?? {}, null, 2), 'application/json')
    if (kind === 'png' && renderer) {
      const link = document.createElement('a')
      link.href = renderer.snapshot()
      link.download = `${base}.png`
      link.click()
    }
    error.value = ''
  } catch (caught) {
    error.value = errorText(caught)
  }
}

async function toggleFullscreen() {
  try {
    if (document.fullscreenElement) await document.exitFullscreen()
    else await viewer.value?.requestFullscreen()
  } catch (caught) {
    error.value = errorText(caught)
  }
}

function state() {
  return {
    label: label.value,
    text: text.value,
    summary: summary.value,
    bonds: inspection.value?.bondsCount ?? 0,
    ringStats: inspection.value?.ringStats ?? [],
    symmetry: inspection.value?.symmetry,
    wyckoff: inspection.value?.wyckoff.map((entry: any) => entry.letter) ?? [],
    wyckoffOrbits: orbitMultiplicities(),
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
    selection: selection.value,
    measure: measure.value,
    brillouin: brillouin.value,
    zone: zone.value
      ? {
          vertices: zone.value.vertices.length,
          faces: zone.value.faces.length,
          edges: zone.value.edges.length,
        }
      : null,
    base: baseText.value,
  }
}

function orbitMultiplicities(): number[] {
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

onMounted(async () => {
  if (!host.value) return
  renderer = new SceneRenderer(host.value)
  document.addEventListener('fullscreenchange', () => {
    fullscreen.value = document.fullscreenElement === viewer.value
  })
  try {
    const module = await loadWasm()
    wasmModule.value = module
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
    operation: async (name: string, params: Record<string, any> = {}) => {
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
    selectAtoms,
    toggleBrillouin,
    brillouinZone: () => wasmModule.value.brillouin_zone(text.value),
    pxrd: (radiation = 'cu-ka') => {
      return wasmModule.value.pxrd_info(text.value, radiation, 60, 0.02)
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
      <span class="status">{{ ready ? (busy ? 'Working…' : label) : 'Loading the WebAssembly module…' }}</span>
      <button v-if="!fullscreen" :disabled="!ready" @click="toggleFullscreen">Fullscreen</button>
    </div>
    <p v-if="error" class="error">{{ error }}</p>

    <div class="layout">
      <div ref="viewer" class="viewer">
        <div ref="host" class="canvas-host" @click="onCanvasClick" />
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
          <button v-if="fullscreen" :disabled="!ready" @click="toggleFullscreen">Exit fullscreen</button>
        </div>
        <p class="hint canvas-hint">
          Click atoms to measure distance, angle or torsion. {{ brillouin ? 'Brillouin zone shown.' : '' }}
        </p>
      </div>
      <aside class="panel">
        <nav class="tabs">
          <button
            v-for="[id, title] in TABS"
            :key="id"
            :class="{ active: tab === id }"
            @click="tab = id"
          >{{ title }}</button>
        </nav>
        <PlaygroundStructurePanel
          v-if="tab === 'structure'"
          :wasm="wasmModule"
          :text="text"
          :inspection="inspection"
          :scene="scene"
          :busy="busy"
          :label="label"
          @operation="run"
          @apply="(contents, name) => loadText(contents, name)"
          @restore="restore"
        />
        <PlaygroundSymmetryPanel
          v-else-if="tab === 'symmetry'"
          :wasm="wasmModule"
          :text="text"
          :inspection="inspection"
          :zone="zone"
          :brillouin="brillouin"
          :busy="busy"
          @toggle-brillouin="toggleBrillouin"
        />
        <PlaygroundAnalysisPanel
          v-else-if="tab === 'analysis'"
          :wasm="wasmModule"
          :text="text"
          :inspection="inspection"
          :selection="selection"
          :measure="measure"
          :busy="busy"
          @apply="(contents, name) => loadText(contents, name)"
          @clear-selection="selectAtoms([])"
        />
        <PlaygroundPowderPanel
          v-else-if="tab === 'powder'"
          :wasm="wasmModule"
          :text="text"
          :busy="busy"
        />
        <PlaygroundForcefieldPanel
          v-else-if="tab === 'forcefield'"
          :wasm="wasmModule"
          :text="text"
          :busy="busy"
          @apply="(contents, name) => loadText(contents, name)"
        />
        <PlaygroundExportPanel
          v-else
          :summary="summary"
          :scene="scene"
          :label="label"
          :text="text"
          :busy="busy"
          @download="exportAs"
        />
      </aside>
    </div>
  </div>
</template>
<style src="../playground/playground.css"></style>
