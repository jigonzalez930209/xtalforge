<script setup lang="ts">
import { computed, ref } from 'vue'

const props = defineProps<{
  wasm: any
  text: string
  busy: boolean
}>()

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
const sticks = computed(() =>
  peaks.value.map((peak) => {
    const x = ((Number(peak.two_theta) - 5) / Math.max(Number(maxAngle.value) - 5, 1)) * 420
    const y = 150 - (Number(peak.intensity) / maximum.value) * 140
    return { x, y, label: `(${peak.hkl.join(' ')})` }
  }),
)

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
          <option value="cu-ka">Cu Kα</option>
          <option value="cu-ka1">Cu Kα1</option>
          <option value="mo-ka">Mo Kα</option>
          <option value="co-ka">Co Kα</option>
          <option value="fe-ka">Fe Kα</option>
          <option value="cr-ka">Cr Kα</option>
          <option value="ag-ka">Ag Kα</option>
        </select>
      </label>
      <label>2θ max<input v-model.number="maxAngle" type="number" min="10" max="150" /></label>
      <button :disabled="busy" @click="calculate">Calculate</button>
      <button :disabled="!pattern" @click="downloadCsv">CSV</button>
    </div>
    <p v-if="error" class="error">{{ error }}</p>
    <template v-if="pattern">
      <p class="hint">{{ peaks.length }} reflections · {{ points.length }} pattern points</p>
      <svg class="pxrd-plot" viewBox="0 0 420 160" role="img" aria-label="Calculated powder pattern">
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
      </svg>
      <table>
        <thead><tr><th>hkl</th><th>2θ</th><th>d (Å)</th><th>multiplicity</th><th>intensity</th></tr></thead>
        <tbody>
          <tr v-for="peak in peaks.slice(0, 15)" :key="peak.hkl.join(',')">
            <td>({{ peak.hkl.join(' ') }})</td>
            <td>{{ Number(peak.two_theta).toFixed(4) }}</td>
            <td>{{ Number(peak.d).toFixed(4) }}</td>
            <td>{{ peak.multiplicity }}</td>
            <td>{{ Number(peak.intensity).toFixed(1) }}</td>
          </tr>
        </tbody>
      </table>
    </template>
  </section>
</template>
