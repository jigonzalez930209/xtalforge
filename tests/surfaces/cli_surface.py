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
