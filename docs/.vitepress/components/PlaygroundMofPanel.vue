<script setup lang="ts">
import { computed, ref, watch } from 'vue'

const props = defineProps<{
  wasm: any
  text: string
  busy: boolean
}>()

const error = ref('')
const check = ref<any | null>(null)
const pores = ref<any | null>(null)
const identity = ref<any | null>(null)
const deconstruction = ref<any | null>(null)
const spacing = ref(0.8)
const probe = ref(1.2)

const summary = computed(() => {
  if (!check.value) return null
  return {
    clean: check.value.clean,
    errors: check.value.counts?.error ?? 0,
    warnings: check.value.counts?.warning ?? 0,
    findings: (check.value.findings ?? []).slice(0, 8),
  }
})

function run(action: () => void) {
  try {
    action()
    error.value = ''
  } catch (caught) {
    const value = caught as { code?: string; message?: string }
    error.value = value.code && value.message ? `${value.code}: ${value.message}` : String(caught)
  }
}

function parse(text: string) {
  return JSON.parse(text)
}

function runCheck() {
  run(() => {
    check.value = parse(props.wasm.mof_check(props.text, ''))
  })
}

function runPores() {
  run(() => {
    pores.value = parse(
      props.wasm.mof_pores(
        props.text,
        '',
        JSON.stringify({
          spacing: Number(spacing.value),
          probe_radius: Number(probe.value),
        }),
      ),
    )
  })
}

function runIdentity() {
  run(() => {
    identity.value = parse(props.wasm.mof_mofid(props.text, ''))
  })
}

function runDeconstruction() {
  run(() => {
    deconstruction.value = parse(props.wasm.mof_deconstruct(props.text, ''))
  })
}

watch(
  () => props.text,
  () => {
    check.value = null
    pores.value = null
    identity.value = null
    deconstruction.value = null
  },
)
</script>

<template>
  <section class="tab-section">
    <details open>
      <summary>MOF health</summary>
      <button :disabled="busy" @click="runCheck">Check framework</button>
      <template v-if="summary">
        <dl class="summary-grid">
          <div><dt>Clean</dt><dd>{{ summary.clean }}</dd></div>
          <div><dt>Errors</dt><dd>{{ summary.errors }}</dd></div>
          <div><dt>Warnings</dt><dd>{{ summary.warnings }}</dd></div>
        </dl>
        <ul class="mono">
          <li v-for="finding in summary.findings" :key="`${finding.code}-${finding.message}`">
            [{{ finding.severity }}] {{ finding.code }}: {{ finding.message }}
          </li>
        </ul>
      </template>
    </details>

    <details>
      <summary>Pores</summary>
      <div class="row">
        <label>spacing (Å)<input v-model.number="spacing" type="number" min="0.2" step="0.1" /></label>
        <label>probe (Å)<input v-model.number="probe" type="number" min="0" step="0.1" /></label>
        <button :disabled="busy" @click="runPores">Compute</button>
      </div>
      <dl v-if="pores" class="summary-grid">
        <div><dt>PLD</dt><dd>{{ Number(pores.descriptors?.pld ?? 0).toFixed(3) }} Å</dd></div>
        <div><dt>LCD</dt><dd>{{ Number(pores.descriptors?.lcd ?? 0).toFixed(3) }} Å</dd></div>
        <div><dt>ASA</dt><dd>{{ Number(pores.descriptors?.asa_m2_g ?? 0).toFixed(2) }} m²/g</dd></div>
        <div><dt>Void</dt><dd>{{ (Number(pores.descriptors?.void_fraction ?? 0) * 100).toFixed(1) }}%</dd></div>
        <div><dt>Channels</dt><dd>{{ pores.descriptors?.channel_dimension ?? 0 }}</dd></div>
      </dl>
    </details>

    <details>
      <summary>Identity and building units</summary>
      <div class="row">
        <button :disabled="busy" @click="runIdentity">MOFid</button>
        <button :disabled="busy" @click="runDeconstruction">Deconstruct</button>
      </div>
      <p v-if="identity" class="mono">MOFid: {{ identity.mofid }}</p>
      <p v-if="identity" class="mono">MOFkey: {{ identity.mofkey }}</p>
      <template v-if="deconstruction">
        <p class="hint">
          {{ deconstruction.units?.length ?? 0 }} units; topology
          {{ deconstruction.topology ?? 'unknown' }}; interpenetration
          {{ deconstruction.interpenetration ?? 0 }}
        </p>
        <ul class="mono">
          <li v-for="unit in (deconstruction.units ?? []).slice(0, 10)" :key="`${unit.kind}-${unit.atoms?.join(',')}`">
            {{ unit.kind }}: {{ unit.formula }} ({{ unit.atoms?.length ?? 0 }} atoms)
          </li>
        </ul>
      </template>
    </details>

    <p v-if="error" class="error">{{ error }}</p>
  </section>
</template>
