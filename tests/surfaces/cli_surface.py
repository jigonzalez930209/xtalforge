"""Runs the shared surface checks through the release CLI.

Usage: python tests/surfaces/cli_surface.py OUT.json (from the repository root).
"""
import json,subprocess,sys
X='./target/release/xtalforge'
def run(*a):
    r=subprocess.run([X,'--json',*a],capture_output=True,text=True)
    if r.returncode: return 'ERR '+r.stderr.strip()[:80]
    return json.loads(r.stdout)
out={}
for s in ['Si','NaCl','diamond','graphene','TiO2','MFI','MOF-5','Fe']:
    f=f'docs/public/samples/{s}.cif'; r={}
    i=run('info',f); r['n_atoms']=i.get('n_atoms'); r['volume']=i.get('volume'); r['formula']=i.get('formula')
    r['sg']=run('symmetry','detect',f)['number']
    r['wyckoff']=''.join(sorted(a['letter'] for a in run('symmetry','wyckoff',f)))
    b=run('bonds',f); r['bonds']=len(b['bonds'] if isinstance(b,dict) else b)
    if r['n_atoms']<=200:
        g=run('rings','--max-size','8',f); r['rings']=len(g['rings'] if isinstance(g,dict) else g) if not isinstance(g,str) else g
    p=run('pxrd','--max-two-theta','60',f); r['pxrd']=[[''.join(map(str,q['hkl'])),round(q['two_theta'],4),q['intensity']] for q in p['reflections']]
    if r['n_atoms']<=200:
        r['energy']=run('ff','energy',f)['energy']
        for m in ['qeq','eqeq','eeq']:
            c=run('ff','charges','--method',m,f); r[m]=c['per_atom'] if isinstance(c,dict) else c
    out[s]=r
json.dump(out,open(sys.argv[1],'w'))

# Porosity, MOF assembly and powder refinement on fixed inputs.
TIO2 = "tests/corpus/simple/TiO2.cif"
RUTILE = "tests/golden/phase-7/rietveld/rutile.xy"
extra = {}
porosity = run("analysis", "porosity", "--spacing", "0.5", "--probe", "1.2", "docs/public/samples/MFI.cif")
extra["porosity"] = [porosity[k] for k in ("accessible_volume", "pore_volume", "surface_area", "largest_included_diameter")]
mof = run("mof", "build", "pcu", "Zn4O", "benzene")
lattice = mof["structure"]["lattice"]
volume = abs(
    lattice[0][0] * (lattice[1][1] * lattice[2][2] - lattice[1][2] * lattice[2][1])
    - lattice[0][1] * (lattice[1][0] * lattice[2][2] - lattice[1][2] * lattice[2][0])
    + lattice[0][2] * (lattice[1][0] * lattice[2][1] - lattice[1][1] * lattice[2][0])
)
extra["mof"] = [len(mof["structure"]["sites"]), round(volume, 3), len(mof["warnings"])]
lebail = run("powder", "lebail", TIO2, RUTILE)
extra["lebail"] = [lebail["r_factor"], lebail["rwp"], len(lebail["intensities"])]
pawley = run("powder", "pawley", TIO2, RUTILE)
extra["pawley"] = [pawley["cell"][0], pawley["cell"][2], pawley["rwp"]]
refined = run("powder", "rietveld", "--refine", TIO2, RUTILE)
extra["rietveld"] = [refined["rwp"], refined["stages"], dict(refined["parameters"])["cell.0"]]
automatic = run("powder", "automatic", TIO2, RUTILE)
extra["automatic"] = [automatic["rwp"], automatic["pawley_rwp"], automatic["best_cell"]]
out["extra"] = extra
json.dump(out, open(sys.argv[1], "w"))
