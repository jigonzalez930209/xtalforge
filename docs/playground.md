# 3D playground

Build, transform and inspect a structure in the browser, then look at the
result in three dimensions. Every value and every atom position comes from the
WebAssembly build of the same Rust functions the API, the CLI and Python use;
nothing is recalculated in JavaScript and no data leaves the page.

<script setup>
import Playground from './.vitepress/components/Playground.vue'
</script>

<Playground />

## What can be tested

- **Load** the bundled samples (silicon, NaCl, diamond, graphene, rutile,
  MFI, MOF-5, iron) or paste any CIF, POSCAR or XYZ file.
- **Transform**: supercell, slab (h k l, layers, vacuum), P1 expansion,
  asymmetric unit, primitive cell, standardized cell and Niggli reduction.
- **Build**: graphene, nanotube `(n, m)`, C60, C20 and a MOF from the
  topology/node/linker catalogue.
- **Analyse**: formula, cell volume, density, P1 atom count, bonds, primitive
  rings up to eight atoms, detected space group, Hall number, Wyckoff letters
  and a Cu Kα powder pattern.
- **View**: ball-and-stick, space-fill, wireframe and sticks, with the cell
  edges, axes, orbit, zoom, rotation and PNG export.

The 3D frame is the [`scene_data`](/api/view) geometry: atom positions,
display radii, element colours, bond endpoints and the 12 cell-edge segments.

## Reproducible checks

The playground is also the end-to-end test surface. `tests/e2e/` drives it
with Playwright: each scenario compares the numbers against the reference
JSON generated with gemmi, ASE and spglib, then takes a deterministic
screenshot and compares it with the committed golden image. Run the full
battery from the repository root:

```bash
scripts/e2e.sh
```

The script builds the WASM package and the docs site, regenerates the
reference data, runs the visual and numeric tests, and leaves the captured
screenshots in `tests/e2e/test-results/` when a check fails.

## Direct links

`/playground?sample=Si` opens a sample. The `sample` query accepts any id from
the sample list. Automation can also call `window.__xf`:

```js
await window.__xf.ready
const state = await window.__xf.operation('supercell', { na: 2, nb: 2, nc: 2 })
window.__xf.setView({ azimuth: 0.9, elevation: 0.4, zoom: 1 })
```
