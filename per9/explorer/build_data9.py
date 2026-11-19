#!/usr/bin/env python3
"""E95 viewer: explorer/data.json extended for per9_core.js / per9_live.html (E94 explorer data format, n = 9).
Keeps 'centers' (x, y, ks, real; same order: the u5 scripts read them) and adds label (even numerator of ks), angles; adds
n, coef (trackmp Curve().coef_int, 43 rows), y0, x0, bps, kinds, nodes (branch.json nodes_detour), punctures = [] (A3 not met).
No leadQ / leadPow (the exact leading coefficient is not used at n = 9). Backs the old file up once (data_pre_viewer.json)."""
import json, os, shutil, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(here, '..'))
import trackmp as tm
fn = os.path.join(here, 'data.json'); bk = os.path.join(here, 'data_pre_viewer.json')
if not os.path.exists(bk): shutil.copy(fn, bk)
D = json.load(open(fn)); B = json.load(open(os.path.join(here, '..', 'branch.json')))
cen = []
for c in D['centers']:
    ev = [k for k in c['ks'] if k % 2 == 0]; assert len(ev) == 1, c['ks']
    cen.append(dict(x=c['x'], y=c['y'], ks=c['ks'], real=c['real'], label=ev[0], angles=c['ks']))
C = tm.Curve()
coef = [[int(a) for a in row] for row in C.coef_int]
out = dict(n=9, coef=coef, y0=B['y0'], x0=[[complex(*v).real, complex(*v).imag] if isinstance(v, (list, tuple)) else [complex(v).real, complex(v).imag] for v in B['x0']],
           bps=B['bps'], kinds=B['kinds'], nodes=[[a, b, R] for a, b, R in B['nodes_detour']], centers=cen, punctures=[])
json.dump(out, open(fn, 'w'))
print('data.json: n', out['n'], '| coef rows', len(coef), 'max row length', max(len(r) for r in coef), '| x0', len(out['x0']), '| bps', len(out['bps']),
      '| nodes', out['nodes'], '| centres', len(cen), '(labels distinct:', len({c['label'] for c in cen}) == len(cen), ') | punctures', len(out['punctures']),
      '| bytes', os.path.getsize(fn))
