<script setup>
import { onMounted, ref } from 'vue'

const ready = ref(false)
const wasmVersion = ref('')
const a = ref(3)
const b = ref(4)
const c = ref(5)
const alpha = ref(60)
const beta = ref(70)
const gamma = ref(80)
const result = ref('')
const error = ref('')
let calculate = null

onMounted(async () => {
  const module = await import('xtalforge-wasm')
  await module.default()
  calculate = module.cell_volume
  wasmVersion.value = module.version()
  ready.value = true
  run()
})

function run() {
  if (!calculate) return
  try {
    result.value = calculate(
      Number(a.value),
      Number(b.value),
      Number(c.value),
      Number(alpha.value),
      Number(beta.value),
      Number(gamma.value),
    ).toFixed(6)
    error.value = ''
  } catch (caught) {
    result.value = ''
    error.value = `${caught.code ?? 'error'}: ${caught.message}`
  }
}
</script>

<template>
  <div class="cell-volume-demo">
    <p v-if="ready" class="status">WASM {{ wasmVersion }}</p>
    <p v-else class="status">Loading the WebAssembly module…</p>
    <div class="grid">
      <label>a (Å)<input v-model.number="a" type="number" step="0.1" /></label>
      <label>b (Å)<input v-model.number="b" type="number" step="0.1" /></label>
      <label>c (Å)<input v-model.number="c" type="number" step="0.1" /></label>
      <label>α (°)<input v-model.number="alpha" type="number" step="1" /></label>
      <label>β (°)<input v-model.number="beta" type="number" step="1" /></label>
      <label>γ (°)<input v-model.number="gamma" type="number" step="1" /></label>
    </div>
    <button :disabled="!ready" @click="run">Calculate</button>
    <p v-if="result" class="result">V = {{ result }} Å³</p>
    <p v-if="error" class="error">{{ error }}</p>
  </div>
</template>

<style scoped>
.cell-volume-demo {
  border: 1px solid var(--vp-c-divider);
  border-radius: 8px;
  padding: 1rem;
}
.status {
  color: var(--vp-c-text-2);
  font-size: 0.9rem;
}
.grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
  gap: 0.75rem;
}
label {
  display: flex;
  flex-direction: column;
  gap: 0.25rem;
  font-size: 0.9rem;
}
input {
  border: 1px solid var(--vp-c-divider);
  border-radius: 4px;
  padding: 0.35rem;
  background: var(--vp-c-bg);
  color: var(--vp-c-text-1);
}
button {
  margin-top: 1rem;
  border: 1px solid var(--vp-c-brand-1);
  border-radius: 4px;
  padding: 0.4rem 0.8rem;
  color: var(--vp-c-brand-1);
}
.result {
  font-family: var(--vp-font-family-mono);
}
.error {
  color: var(--vp-c-danger-1);
}
</style>
