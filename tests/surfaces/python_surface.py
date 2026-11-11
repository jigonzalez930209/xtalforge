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
