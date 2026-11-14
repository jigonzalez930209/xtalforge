// Loads one split WASM package and exercises its tier (Task 15.3.1).
//
// Usage: node crates/xtalforge-wasm/tests/packages.mjs PKG_DIR TIER
//   PKG_DIR  directory produced by wasm-pack (pkg or pkg/<tier>)
//   TIER     core | symmetry | analysis | build | full
//
// The test fails when a package does not load, when an export of a lower tier
// is missing, or when an export of a higher tier leaks into a smaller package.
// The wasm-pack build script runs it after every build; it is a plain Node
// script because the four packages are four separate wasm-bindgen artifacts.
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import { pathToFileURL } from 'node:url';

const [directory, tier] = process.argv.slice(2);
if (!directory || !tier) {
  console.error('usage: packages.mjs PKG_DIR TIER');
  process.exit(2);
}

const dir = path.resolve(directory);
const base = path.join(dir, tier === 'full' ? 'xtalforge_wasm' : `xtalforge_wasm_${tier}`);
const module = await import(pathToFileURL(`${base}.js`).href);
module.initSync({ module: fs.readFileSync(`${base}_bg.wasm`) });

assert.equal(typeof module.version(), 'string', 'version() export');

const expected = {
  core: [
    'read_structure', 'convert_structure', 'supercell_info', 'slab_info', 'scene_buffers',
    'mesh_buffers', 'project_pack', 'deform', 'radial_distribution', 'Structure',
  ],
  symmetry: ['detect_symmetry', 'wyckoff_info', 'standardize_text', 'match_structures'],
  analysis: [
    'pxrd_info', 'porosity_info', 'powder_refine_info', 'bond_valence_sums', 'pxrd_buffers',
  ],
  build: [
    'build_graphene_info', 'build_carbon_text', 'molecule_from_smiles', 'hydrogens_info',
  ],
  full: [
    'ff_energy', 'mof_check', 'agent_capabilities', 'elastic_tensor', 'workspace_archive',
    // Every WASM symbol named in docs/guide/recipes.md (Task 17.1.5).
    'to_p1_info', 'bonds', 'distance', 'ff_optimize_text', 'ff_charges_info',
    'group_structures', 'rattle', 'interpolate', 'optimal_supercell', 'voronoi',
    'coordination', 'oxidation_states', 'unique_miller_indices', 'wulff',
    'adsorption_sites', 'mof_mofid', 'mof_pores', 'mof_widom', 'build_polymer_text',
    // Phase 18 engine bridge.
    'engine_presets_info', 'engine_input_set_info', 'engine_parse_info', 'engine_bands_info',
    'engine_dos_info', 'engine_script_info', 'engine_failures_info',
  ],
};
const required = {
  core: expected.core,
  symmetry: [...expected.core, ...expected.symmetry],
  analysis: [...expected.core, ...expected.analysis],
  build: [...expected.core, ...expected.build],
  full: [
    ...expected.core, ...expected.symmetry, ...expected.analysis, ...expected.build, ...expected.full,
  ],
};
const forbidden = {
  core: ['detect_symmetry', 'pxrd_info', 'build_graphene_info', 'ff_energy', 'mof_check'],
  symmetry: ['pxrd_info', 'build_graphene_info', 'ff_energy', 'mof_check'],
  analysis: ['detect_symmetry', 'build_graphene_info', 'ff_energy', 'mof_check'],
  build: ['detect_symmetry', 'pxrd_info', 'ff_energy', 'mof_check'],
  full: [],
};

for (const exportName of new Set(required[tier])) {
  assert.ok(exportName in module, `${tier} is missing ${exportName}`);
}
for (const exportName of forbidden[tier]) {
  assert.ok(!(exportName in module), `${tier} must not export ${exportName}`);
}

// The object handle carries the methods of its tier and no others.
const proto = module.Structure.prototype;
assert.equal(typeof proto.to_cif, 'function', 'Structure.to_cif');
assert.equal(typeof proto.supercell, 'function', 'Structure.supercell');
const has_symmetry = tier === 'symmetry' || tier === 'full';
const has_analysis = tier === 'analysis' || tier === 'full';
const has_build = tier === 'build' || tier === 'full';
assert.equal(typeof proto.detect_symmetry === 'function', has_symmetry, 'Structure.detect_symmetry');
assert.equal(typeof proto.pxrd === 'function', has_analysis, 'Structure.pxrd');
assert.equal(typeof proto.add_hydrogens === 'function', has_build, 'Structure.add_hydrogens');
assert.equal(typeof proto.ff_energy === 'function', tier === 'full', 'Structure.ff_energy');

const cif = fs.readFileSync(
  new URL('../../../docs/public/samples/Si.cif', import.meta.url),
  'utf8',
);
const summary = module.read_structure(cif, '');
assert.ok(summary.n_atoms > 0, 'read_structure returns atoms');
const structure = module.Structure.fromCif(cif);
assert.ok(structure.nAtoms > 0, 'Structure.fromCif');
structure.free();

if (tier === 'symmetry' || tier === 'full') {
  assert.equal(module.detect_symmetry(cif, 0.001).number, 227, 'silicon space group');
  assert.ok(module.kpath_info(cif).path.length > 0, 'band path (gzip k-path data)');
}
if (tier === 'analysis' || tier === 'full') {
  assert.ok(module.pxrd_info(cif, 'cu-ka', 60, 0.02).reflections.length > 0, 'PXRD');
  assert.ok(module.bond_valence_sums(cif, '', 3.2).length > 0, 'bond valence (gzip ICSD data)');
}
if (tier === 'build' || tier === 'full') {
  assert.ok(module.build_graphene_info(2, 2, 1, 3.35, 10).n_atoms > 0, 'graphene');
}
if (tier === 'full') {
  assert.ok(Number.isFinite(module.ff_energy(cif)), 'UFF energy');
  assert.ok(module.agent_capabilities().length > 0, 'agent capabilities');
}

console.log(`OK: ${tier} package loads and exercises its tier`);
