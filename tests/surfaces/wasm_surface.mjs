// Runs the shared surface checks through the WASM package under Node.
// Usage: node tests/surfaces/wasm_surface.mjs OUT.json (from the repository root).
import fs from 'fs';
const root = new URL('../../', import.meta.url).pathname;
const pkg = root + 'crates/xtalforge-wasm/pkg/';
const m = await import(pkg + 'xtalforge_wasm.js');
m.initSync({ module: fs.readFileSync(pkg + 'xtalforge_wasm_bg.wasm') });
const out = {};
for (const s of ['Si','NaCl','diamond','graphene','TiO2','MFI','MOF-5','Fe']) {
  const t = fs.readFileSync(`${root}docs/public/samples/${s}.cif`, 'utf8');
  const r = {};
  const sum = m.read_structure(t, ''); r.n_atoms = sum.n_atoms; r.volume = sum.volume; r.formula = sum.formula;
  r.sg = m.detect_symmetry(t, 0.001).number;
  r.wyckoff = m.wyckoff_info(t, 0.001).map(w => w.letter).sort().join('');
  r.bonds = m.bonds(t).length;
  if (sum.n_atoms <= 200) r.rings = m.rings(t, 8).length;
  const p = m.pxrd_info(t, 'cu-ka', 60, 0.02); r.pxrd = p.reflections.map(x => [x.hkl.join(''), +x.two_theta.toFixed(4), x.intensity]);
  if (sum.n_atoms <= 200) { r.energy = m.ff_energy(t); for (const meth of ['qeq','eqeq','eeq']) { try { r[meth] = m.ff_charges_info(t, meth, 0).per_atom; } catch (e) { r[meth] = 'ERR ' + (e.code || e.message); } } }
  out[s] = r;
}
fs.writeFileSync(process.argv[2], JSON.stringify(out));
