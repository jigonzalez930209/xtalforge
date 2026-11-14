---
aside: false
---

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

- **Tour**: the first tab walks through seventeen steps, from reading a CIF and
  detecting its space group to the Brillouin zone, measures, supercells,
  slabs, PXRD, nanotubes, twisted bilayers, C60, MOFs, disordered carbon,
  polymers, porosity, the force field and the exports. Each step runs the
  real function, sets the camera and opens the tab where the result is read.
- **Load** the bundled samples (silicon, NaCl, diamond, graphene, rutile,
  MFI, MOF-5, iron) or paste any CIF, POSCAR or XYZ file.
- **Transform**: supercell, slab (h k l, layers, vacuum), P1 expansion,
  asymmetric unit, primitive cell, standardized cell and Niggli reduction.
- **Build**: graphene, nanotube `(n, m)`, C60, C20, nanoribbon, commensurate
  twisted bilayer (21.79°, 13.17°, 9.43°…), ZTC, Schwarz P, disordered carbon
  on the dia/pcu nets, polymer from a SMILES repeat unit and a MOF from the
  topology/node/linker catalogue.
- **Analyse**: formula, cell volume, density, P1 atom count, bonds, primitive
  rings up to eight atoms, detected space group, Hall number, Wyckoff letters,
  tolerance ladder, diagnostics, Brillouin zone, band path, periodic
  neighbours, functional groups, porosity, hydrogens and atom-picked
  distance/angle/torsion.
- **Force field**: UFF energy, atom types, QEq/EQeq charges, vibrations and a
  UFF minimization that replaces the current structure.
- **MOF**: health check, MOFid/MOFkey, deconstruction into building units and
  the native PLD/LCD/ASA/void/channel descriptors.
- **Powder**: Cu Kα and other radiations, International Tables form factors,
  a plotted pattern with hkl labels, an I/I₀ peak table and CSV export, plus a
  refinement panel: Le Bail, Pawley, Rietveld and the automatic run against
  the bundled synthetic rutile pattern or your own 2θ/intensity file, with an
  observed/calculated/difference plot.
- **View**: ball-and-stick, space-fill, wireframe and sticks, with the cell
  edges, axes, orbit, zoom, rotation, atom selection, fullscreen and PNG
  export.
- **Export**: CIF, POSCAR, XYZ, STL, PNG, scene JSON and CSV.

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

<style>
body:has(.playground) .VPDoc .container,
body:has(.playground) .VPDoc .content,
body:has(.playground) .VPDoc .content-container {
  max-width: 1440px !important;
  width: 100%;
}
</style>
