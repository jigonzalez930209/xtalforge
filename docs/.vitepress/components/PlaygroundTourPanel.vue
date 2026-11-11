<script setup lang="ts">
import { computed } from 'vue'
import { TOUR, type TourStep } from '../playground/tour'

const props = defineProps<{
  busy: boolean
  current: number
}>()

const emit = defineEmits<{
  (event: 'step', step: TourStep, index: number): void
}>()

const areas = computed(() => [...new Set(TOUR.map((step) => step.area))])
const next = computed(() => Math.min(props.current + 1, TOUR.length - 1))

function go(index: number) {
  const step = TOUR[index]
  if (step) emit('step', step, index)
}
</script>

<template>
  <section class="tab-section tour">
    <p class="hint">
      A guided walk through the library. Each step runs the real WebAssembly
      function, updates the 3D view and opens the tab where the result is read.
    </p>
    <div class="row">
      <button class="primary" :disabled="busy" @click="go(current < 0 ? 0 : next)">
        {{ current < 0 ? 'Start the tour' : current >= TOUR.length - 1 ? 'Last step' : 'Next step →' }}
      </button>
      <span class="hint">{{ current < 0 ? `${TOUR.length} steps` : `Step ${current + 1} of ${TOUR.length}` }}</span>
    </div>
    <div class="tour-areas">
      <span v-for="area in areas" :key="area" class="chip">{{ area }}</span>
    </div>
    <ol class="tour-steps">
      <li
        v-for="(step, index) in TOUR"
        :key="step.id"
        :class="{ active: index === current, done: index < current }"
      >
        <button :disabled="busy" :data-step="step.id" @click="go(index)">
          <strong>{{ step.title }}</strong>
          <span class="chip">{{ step.area }}</span>
        </button>
        <p v-if="index === current" class="tour-note">{{ step.note }}</p>
      </li>
    </ol>
  </section>
</template>
