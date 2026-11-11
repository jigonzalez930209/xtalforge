<script setup lang="ts">
import { computed, ref } from 'vue'

const props = defineProps<{ wasm: any; text: string; busy: boolean }>()

type Point = [number, number]
const patternText = ref('')
const patternName = ref('')
const method = ref('pawley')
const result = ref<any | null>(null)
const calculated = ref<number[]>([])
const error = ref('')
const running = ref(false)

const observed = computed<Point[]>(() =>
  patternText.value
    .split('\n')
    .map((line) => line.trim())
    .filter((line) => line && !line.startsWith('#'))
    .map((line) => line.split(/\s+/).map(Number))
    .filter((values) => values.length >= 2 && values.every(Number.isFinite))
    .map((values) => [values[0], values[1]] as Point),
)

function plain(value: any): any {
  return JSON.parse(JSON.stringify(value, (_, item) => (item instanceof Map ? Object.fromEntries(item) : item)))
}

async function useReference() {
  const base = import.meta.env.BASE_URL
  const response = await fetch(`${base}samples/rutile-cuka.xy`)
  patternText.value = await response.text()
  patternName.value = 'rutile-cuka.xy (synthetic rutile, Cu Kα1+Kα2)'
  result.value = null
}

async function onFile(event: Event) {
  const file = (event.target as HTMLInputElement).files?.[0]
  if (!file) return
  patternText.value = await file.text()
  patternName.value = file.name
  result.value = null
}

function run() {
  if (!props.wasm || !patternText.value) return
  running.value = true
  error.value = ''
  calculated.value = []
  try {
    const wasm = props.wasm
    if (method.value === 'lebail') result.value = plain(wasm.powder_lebail_info(props.text, patternText.value))
    if (method.value === 'pawley') {
      result.value = plain(wasm.powder_pawley_width_info(props.text, patternText.value, 'pseudo_voigt', true, true))
      calculated.value = result.value.calculated ?? []
    }
    if (method.value === 'rietveld') {
      result.value = plain(wasm.powder_refine_info(props.text, patternText.value, true, true, false, true, true, true, false, true))
    }
    if (method.value === 'automatic') result.value = plain(wasm.powder_automatic_info(props.text, patternText.value))
  } catch (caught: any) {
    error.value = caught?.code ? `${caught.code}: ${caught.message}` : String(caught?.message ?? caught)
    result.value = null
  } finally {
    running.value = false
  }
}

const rows = computed<[string, string][]>(() => {
  const value = result.value
  if (!value) return []
  const out: [string, string][] = []
  const number = (x: any, digits = 4) => (typeof x === 'number' ? x.toFixed(digits) : String(x))
  if (value.r_factor !== undefined) out.push(['R (Le Bail)', number(value.r_factor)])
  if (value.rp !== undefined) out.push(['Rp', number(value.rp)])
  if (value.rwp !== undefined) out.push(['Rwp', number(value.rwp)])
  if (value.chi_squared !== undefined) out.push(['χ²', number(value.chi_squared, 3)])
  if (Array.isArray(value.cell)) out.push(['cell a, b, c (Å)', value.cell.slice(0, 3).map((x: number) => x.toFixed(4)).join(', ')])
  const leading = value.leading_cells?.[0]
  if (leading) {
    out.push(['best cell', `${leading.system} ${leading.class_symbol}`])
    out.push(['best a, b, c (Å)', [leading.a, leading.b, leading.c].map((x: number) => Number(x).toFixed(4)).join(', ')])
  }
  if (value.zero_shift !== undefined) out.push(['zero shift (°)', number(value.zero_shift)])
  if (value.size !== undefined) out.push(['FWHM size, strain (°)', `${number(value.size)}, ${number(value.strain)}`])
  if (value.intensities) out.push(['reflections', String(value.intensities.length)])
  if (value.reflections) out.push(['reflections', String(value.reflections.length)])
  if (value.stages) out.push(['stages', value.stages.join(' → ')])
  const parameters = Array.isArray(value.parameters) ? value.parameters : Object.entries(value.parameters ?? {})
  for (const [name, x] of parameters.slice(0, 10)) out.push([name, number(x, 5)])
  if (value.pawley_rwp !== undefined) out.push(['best Pawley Rwp', number(value.pawley_rwp)])
  if (value.steps) out.push(['steps', value.steps.join(' → ')])
  return out
})

const plot = computed(() => {
  const points = observed.value
  if (!points.length) return null
  const [width, height] = [420, 150]
  const first = points[0][0]
  const last = points[points.length - 1][0]
  const top = Math.max(...points.map((p) => p[1]), ...calculated.value, 1)
  const x = (angle: number) => ((angle - first) / Math.max(last - first, 1e-9)) * width
  const y = (value: number) => height - (value / top) * (height - 10)
  const line = (values: number[]) => values.map((v, i) => `${x(points[i][0]).toFixed(1)},${y(v).toFixed(1)}`).join(' ')
  const obs = line(points.map((p) => p[1]))
  const calc = calculated.value.length === points.length ? line(calculated.value) : ''
  const diff = calc
    ? points.map((p, i) => `${x(p[0]).toFixed(1)},${(height + 18 - ((p[1] - calculated.value[i]) / top) * 40).toFixed(1)}`).join(' ')
    : ''
  return { obs, calc, diff, first, last }
})
</script>

<template>
  <details class="refine" open>
    <summary>Refine against a measured pattern</summary>
    <p class="hint">
      Le Bail, Pawley, Rietveld (Pawley cell → background/scale → zero/cell → width → structure)
      and the automatic run, on the current structure and a two-column 2θ/intensity pattern.
    </p>
    <div class="row">
      <button :disabled="busy" @click="useReference">Reference rutile pattern</button>
      <label>file<input type="file" accept=".xy,.xye,.dat,.txt,.csv" @change="onFile" /></label>
    </div>
    <p v-if="patternName" class="hint">{{ patternName }} · {{ observed.length }} points</p>
    <div class="row">
      <label>
        method
        <select v-model="method">
          <option value="lebail">Le Bail</option>
          <option value="pawley">Pawley</option>
          <option value="rietveld">Rietveld</option>
          <option value="automatic">Automatic</option>
        </select>
      </label>
      <button class="primary" :disabled="busy || running || !patternText" @click="run">Refine</button>
    </div>
    <p v-if="error" class="error">{{ error }}</p>
    <svg v-if="plot" class="pxrd-plot" viewBox="-6 -10 432 200" role="img" aria-label="Observed and calculated pattern">
      <line x1="0" y1="150" x2="420" y2="150" class="axis" />
      <polyline :points="plot.obs" class="observed" />
      <polyline v-if="plot.calc" :points="plot.calc" class="pattern" />
      <polyline v-if="plot.diff" :points="plot.diff" class="difference" />
      <text x="0" y="186" class="tick-label">{{ plot.first.toFixed(0) }}°</text>
      <text x="420" y="186" class="tick-label" text-anchor="end">{{ plot.last.toFixed(0) }}° 2θ</text>
    </svg>
    <table v-if="rows.length" class="refine-results">
      <tbody>
        <tr v-for="[name, value] in rows" :key="name"><th>{{ name }}</th><td>{{ value }}</td></tr>
      </tbody>
    </table>
  </details>
</template>
