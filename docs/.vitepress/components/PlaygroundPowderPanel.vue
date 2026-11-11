<script setup lang="ts">
import { computed, onMounted, ref, watch } from 'vue'

const props = defineProps<{
  wasm: any
  text: string
  busy: boolean
}>()

// Rendered with v-for: five or more static <option> nodes are hoisted into
// one static HTML block, which Vue cannot insert inside a <select>.
const RADIATIONS = [
  ['cu-ka', 'Cu Kα'],
  ['cu-ka1', 'Cu Kα1'],
  ['mo-ka', 'Mo Kα'],
  ['co-ka', 'Co Kα'],
  ['fe-ka', 'Fe Kα'],
  ['cr-ka', 'Cr Kα'],
  ['ag-ka', 'Ag Kα'],
] as const
const radiation = ref('cu-ka')
const maxAngle = ref(60)
const pattern = ref<any | null>(null)
const error = ref('')

const points = computed(() => (pattern.value?.points ?? []) as number[][])
const peaks = computed(() => (pattern.value?.reflections ?? []) as any[])
const maximum = computed(() =>
  points.value.reduce((value, point) => Math.max(value, Number(point[1])), 0) || 1,
)
const polyline = computed(() => {
  const width = 420
  const height = 150
  return points.value
    .map((point) => {
      const [x, y] = [Number(point[0]), Number(point[1])]
      const sx = ((x - 5) / Math.max(Number(maxAngle.value) - 5, 1)) * width
      const sy = height - (y / maximum.value) * (height - 10)
      return `${sx.toFixed(1)},${sy.toFixed(1)}`
    })
    .join(' ')
})
const strongest = computed(() =>
  peaks.value.reduce((value, peak) => Math.max(value, Number(peak.intensity)), 0) || 1,
)
const toX = (angle: number) => ((angle - 5) / Math.max(Number(maxAngle.value) - 5, 1)) * 420
// Sticks are integrated intensities, scaled to their own maximum so they
// read as I/I0 next to the profile (whose height also depends on the width).
const sticks = computed(() =>
  peaks.value.map((peak) => {
    const relative = Number(peak.intensity) / strongest.value
    return {
      x: toX(Number(peak.two_theta)),
      y: 150 - relative * 140,
      relative,
      label: `${peak.hkl.join('')}`,
    }
  }),
)
const labels = computed(() =>
  [...sticks.value]
    .sort((a, b) => b.relative - a.relative)
    .slice(0, 6)
    .filter((stick) => stick.relative > 0.05),
)
const ticks = computed(() => {
  const out: number[] = []
  for (let angle = 10; angle <= Number(maxAngle.value); angle += 10) out.push(angle)
  return out
})

function calculate() {
  try {
    pattern.value = props.wasm.pxrd_info(
      props.text,
      radiation.value,
      Number(maxAngle.value),
      0.02,
    )
    error.value = ''
  } catch (caught) {
    const value = caught as { code?: string; message?: string }
    error.value = value.code && value.message ? `${value.code}: ${value.message}` : String(caught)
  }
}

watch(
  () => props.text,
  () => {
    if (props.wasm && props.text) calculate()
  },
)
onMounted(() => {
  if (props.wasm && props.text) calculate()
})

function downloadCsv() {
  const lines = ['two_theta,intensity']
  for (const point of points.value) lines.push(`${point[0]},${point[1]}`)
  const blob = new Blob([`${lines.join('\n')}\n`], { type: 'text/csv' })
  const link = document.createElement('a')
  link.href = URL.createObjectURL(blob)
  link.download = 'pxrd.csv'
  link.click()
  URL.revokeObjectURL(link.href)
}
</script>

<template>
  <section class="tab-section">
    <div class="row">
      <label>
        radiation
        <select v-model="radiation">
          <option v-for="[value, name] in RADIATIONS" :key="value" :value="value">{{ name }}</option>
        </select>
      </label>
      <label>2θ max<input v-model.number="maxAngle" type="number" min="10" max="150" /></label>
      <button :disabled="busy" @click="calculate">Calculate</button>
      <button :disabled="!pattern" @click="downloadCsv">CSV</button>
    </div>
    <p v-if="error" class="error">{{ error }}</p>
    <template v-if="pattern">
      <p class="hint">{{ peaks.length }} reflections · {{ points.length }} pattern points</p>
      <svg class="pxrd-plot" viewBox="-6 -14 432 190" role="img" aria-label="Calculated powder pattern">
        <line x1="0" y1="150" x2="420" y2="150" class="axis" />
        <line x1="0" y1="0" x2="0" y2="150" class="axis" />
        <polyline :points="polyline" class="pattern" />
        <line
          v-for="stick in sticks"
          :key="`${stick.x}-${stick.label}`"
          :x1="stick.x"
          :x2="stick.x"
          :y1="150"
          :y2="stick.y"
          class="stick"
        />
        <text
          v-for="stick in labels"
          :key="`label-${stick.label}-${stick.x}`"
          :x="stick.x"
          :y="Math.max(stick.y - 3, -4)"
          class="peak-label"
          text-anchor="middle"
        >{{ stick.label }}</text>
        <g v-for="angle in ticks" :key="`tick-${angle}`">
          <line :x1="toX(angle)" :x2="toX(angle)" y1="150" y2="154" class="axis" />
          <text :x="toX(angle)" y="166" class="tick-label" text-anchor="middle">{{ angle }}</text>
        </g>
        <text x="420" y="176" class="tick-label" text-anchor="end">2θ (°)</text>
      </svg>
      <table>
        <thead><tr><th>hkl</th><th>2θ</th><th>d (Å)</th><th>mult.</th><th>I/I₀ (%)</th></tr></thead>
        <tbody>
          <tr v-for="peak in peaks.slice(0, 15)" :key="peak.hkl.join(',')">
            <td>({{ peak.hkl.join(' ') }})</td>
            <td>{{ Number(peak.two_theta).toFixed(4) }}</td>
            <td>{{ Number(peak.d).toFixed(4) }}</td>
            <td>{{ peak.multiplicity }}</td>
            <td>{{ (100 * Number(peak.intensity) / strongest).toFixed(1) }}</td>
          </tr>
        </tbody>
      </table>
    </template>
  </section>
</template>
