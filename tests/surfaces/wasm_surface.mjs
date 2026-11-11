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
// Porosity, MOF assembly and powder refinement on fixed inputs.
const tio2 = fs.readFileSync(root + 'tests/corpus/simple/TiO2.cif', 'utf8');
const rutile = fs.readFileSync(root + 'tests/golden/phase-7/rietveld/rutile.xy', 'utf8');
const plain = (value) => JSON.parse(JSON.stringify(value, (k, v) => (v instanceof Map ? Object.fromEntries(v) : v)));
const extra = {};
const porosity = plain(m.porosity_info(fs.readFileSync(root + 'docs/public/samples/MFI.cif', 'utf8'), 0.5, 1.2));
extra.porosity = ['accessible_volume', 'pore_volume', 'surface_area', 'largest_included_diameter'].map((k) => porosity[k]);
const mof = plain(m.build_mof_info('pcu', 'Zn4O', 'benzene'));
extra.mof = [mof.n_sites, Math.round(mof.volume * 1000) / 1000, (mof.warnings || []).length];
const lebail = plain(m.powder_lebail_info(tio2, rutile));
extra.lebail = [lebail.r_factor, lebail.rwp, lebail.intensities.length];
const pawley = plain(m.powder_pawley_info(tio2, rutile, 'pseudo_voigt'));
extra.pawley = [pawley.cell[0], pawley.cell[2], pawley.rwp];
const refined = plain(m.powder_refine_info(tio2, rutile, true, true, false, true, true, true, false, true));
const parameters = Array.isArray(refined.parameters) ? Object.fromEntries(refined.parameters) : refined.parameters;
extra.rietveld = [refined.rwp, refined.stages, parameters['cell.0']];
const automatic = plain(m.powder_automatic_info(tio2, rutile));
extra.automatic = [automatic.rwp, automatic.pawley_rwp];
out.extra = extra;
fs.writeFileSync(process.argv[2], JSON.stringify(out));
