"""Runs the shared surface checks through the Python bindings.

Usage: python tests/surfaces/python_surface.py OUT.json (from the repository root).
"""
import json,sys,xtalforge as x
out={}
for s in ['Si','NaCl','diamond','graphene','TiO2','MFI','MOF-5','Fe']:
    t=open(f'docs/public/samples/{s}.cif').read(); r={}
    sm=x.read_structure(t,'cif'); sm=json.loads(sm) if isinstance(sm,str) else sm
    info=sm
    r['n_atoms']=info.get('n_atoms'); r['volume']=info.get('volume'); r['formula']=info.get('formula')
    d=x.detect_symmetry(t); d=json.loads(d) if isinstance(d,str) else d; r['sg']=d['number']
    w=x.wyckoff_assignments(t); w=json.loads(w) if isinstance(w,str) else w; r['wyckoff']=''.join(sorted(a['letter'] for a in w))
    b=x.bonds(t); b=json.loads(b) if isinstance(b,str) else b; r['bonds']=len(b)
    if (r['n_atoms'] or 0)<=200:
        rg=x.rings(t,8); rg=json.loads(rg) if isinstance(rg,str) else rg; r['rings']=len(rg)
    p=x.pxrd(t,'cu-ka',60.0,0.02); p=json.loads(p) if isinstance(p,str) else p
    r['pxrd']=[[''.join(map(str,q['hkl'])),round(q['two_theta'],4),q['intensity']] for q in p['reflections']]
    if (r['n_atoms'] or 0)<=200:
        r['energy']=x.ff_energy(t)
        for m in ['qeq','eqeq','eeq']:
            try:
                c=x.ff_charges(t,m); c=json.loads(c) if isinstance(c,str) else c; r[m]=c['per_atom'] if isinstance(c,dict) else c
            except Exception as e: r[m]='ERR '+str(e)
    out[s]=r
json.dump(out,open(sys.argv[1],'w'))

# Porosity, MOF assembly and powder refinement on fixed inputs.
tio2 = open("tests/corpus/simple/TiO2.cif").read()
rutile = open("tests/golden/phase-7/rietveld/rutile.xy").read()
extra = {}
porosity = x.porosity(open("docs/public/samples/MFI.cif").read(), 0.5, 1.2)
extra["porosity"] = [porosity[k] for k in ("accessible_volume", "pore_volume", "surface_area", "largest_included_diameter")]
mof = x.build_mof("pcu", "Zn4O", "benzene")
extra["mof"] = [mof["n_sites"], round(mof["volume"], 3), len(mof["warnings"])]
lebail = x.powder_lebail(tio2, rutile)
extra["lebail"] = [lebail["r_factor"], lebail["rwp"], len(lebail["intensities"])]
pawley = x.powder_pawley(tio2, rutile)
extra["pawley"] = [pawley["cell"][0], pawley["cell"][2], pawley["rwp"]]
refined = x.powder_refine(tio2, rutile, refine_positions=True, refine_adps=True, refine_zero=True, refine_cell=True, refine_width=True, doublet=True)
extra["rietveld"] = [refined["rwp"], refined["stages"], dict(refined["parameters"])["cell.0"]]
automatic = x.powder_automatic(tio2, rutile)
extra["automatic"] = [automatic["rwp"], automatic["pawley_rwp"]]
out["extra"] = extra
json.dump(out, open(sys.argv[1], "w"))
