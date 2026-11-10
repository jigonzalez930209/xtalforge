<script setup lang="ts">
import { ref } from 'vue'

const props = defineProps<{
  wasm: any
  text: string
  busy: boolean
}>()

const emit = defineEmits<{
  (event: 'apply', text: string, name: string): void
}>()

const energy = ref<number | null>(null)
const method = ref('qeq')
const totalCharge = ref(0)
const charges = ref<any | null>(null)
const types = ref<any[] | null>(null)
const maxSteps = ref(200)
const forceTolerance = ref(0.05)
const optimizeCell = ref(false)
const optimizeSummary = ref<any | null>(null)
const vibrations = ref<any | null>(null)
const running = ref(false)
const error = ref('')

function run(action: () => void) {
  running.value = true
  try {
    action()
    error.value = ''
  } catch (caught) {
    const value = caught as { code?: string; message?: string }
    error.value = value.code && value.message ? `${value.code}: ${value.message}` : String(caught)
  } finally {
    running.value = false
  }
}

function calculateEnergy() {
  run(() => {
    energy.value = Number(props.wasm.ff_energy(props.text))
    optimizeSummary.value = null
  })
}

function calculateCharges() {
  run(() => {
    charges.value = props.wasm.ff_charges_info(
      props.text,
      method.value,
      Number(totalCharge.value),
    )
  })
}

function calculateTypes() {
  run(() => {
    types.value = props.wasm.ff_types_info(props.text)
  })
}

function optimize() {
  run(() => {
    const cif = props.wasm.ff_optimize_text(
      props.text,
      'lbfgs',
      Number(maxSteps.value),
      Number(forceTolerance.value),
      optimizeCell.value,
    )
    optimizeSummary.value = props.wasm.read_cif(cif)
    emit('apply', cif, 'UFF minimized')
  })
}

function calculateVibrations() {
  run(() => {
    vibrations.value = props.wasm.ff_vibrations_info(props.text, 0.01)
  })
}
</script>

<template>
  <section class="tab-section">
    <p class="hint">
      UFF and QEq run in the browser. Large frameworks can take a few seconds;
      the energy call is the cheapest first check.
    </p>
    <div class="row">
      <button :disabled="busy || running" @click="calculateEnergy">UFF energy</button>
      <button :disabled="busy || running" @click="calculateTypes">Atom types</button>
      <button :disabled="busy || running" @click="calculateVibrations">Vibrations</button>
    </div>
    <p v-if="energy !== null" class="mono">Energy = {{ energy.toFixed(6) }} eV</p>
    <p v-if="types" class="mono">{{ types.length }} typed atoms; first: {{ types.slice(0, 6).map((item: any) => `#${item.atom}:${item.type}`).join(' ') }}</p>
    <template v-if="vibrations">
      <p class="mono">
        {{ vibrations.frequencies.length }} modes · imaginary {{ vibrations.imaginary.length }} ·
        zero {{ vibrations.zero_modes.length }}
      </p>
      <p class="mono">
        range {{ Number(vibrations.frequencies[0]).toFixed(2) }} to
        {{ Number(vibrations.frequencies[vibrations.frequencies.length - 1]).toFixed(2) }} cm⁻¹
      </p>
    </template>

    <details>
      <summary>QEq / EQeq charges</summary>
      <div class="row">
        <label>
          method
          <select v-model="method">
            <option value="qeq">QEq</option>
            <option value="eqeq">EQeq</option>
          </select>
        </label>
        <label>total charge<input v-model.number="totalCharge" type="number" step="1" /></label>
        <button :disabled="busy || running" @click="calculateCharges">Calculate charges</button>
      </div>
      <template v-if="charges">
        <p class="mono">
          total {{ Number(charges.total).toFixed(4) }} e · {{ charges.iterations }} iterations ·
          range {{ Math.min(...charges.per_atom).toFixed(3) }} to {{ Math.max(...charges.per_atom).toFixed(3) }}
        </p>
        <p class="mono">{{ charges.per_atom.slice(0, 12).map((value: number) => Number(value).toFixed(3)).join(' ') }}</p>
      </template>
    </details>

    <details>
      <summary>UFF minimization</summary>
      <div class="row">
        <label>max steps<input v-model.number="maxSteps" type="number" min="1" max="2000" /></label>
        <label>force tol.<input v-model.number="forceTolerance" type="number" min="0.001" step="0.01" /></label>
        <label class="check"><input v-model="optimizeCell" type="checkbox" /> cell</label>
        <button :disabled="busy || running" @click="optimize">Minimize</button>
      </div>
      <p v-if="optimizeSummary" class="mono">
        New cell: {{ Number(optimizeSummary.volume).toFixed(3) }} Å³ · {{ optimizeSummary.n_atoms }} atoms
      </p>
      <p class="hint">The result replaces the current structure and can be undone with Restore input.</p>
    </details>
    <p v-if="error" class="error">{{ error }}</p>
  </section>
</template>
