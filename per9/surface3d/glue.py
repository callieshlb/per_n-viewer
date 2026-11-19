#!/usr/bin/env python3
"""E95 full run F2: cut-system gluing (../cutglue.py, canaried on E94) -> DIR/edge_perms.json in E94's format, with no roots:
'R' holds None for every vertex and 'd' gives the number of sheets (vertex labels = slit labels of the cut system).
Asserts no mesh vertex within 1e-12 of a cut. usage: glue.py DIR. Raw values only."""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
from cutglue import CutSystem, edge_path
DIR = sys.argv[1]
t0 = time.time()
B = json.load(open(os.path.join(HERE, '..', 'branch.json')))
M = json.load(open(os.path.join(DIR, 'sphere_mesh.json')))
Y = np.array([complex(*v) for v in M['y']]); bp = M['bp']; nV = len(Y)
CS = CutSystem(B); d = CS.d
# degeneracy: distance of every vertex to every cut segment (excluding a branch vertex's own cut endpoint)
A, Bc = CS.A, CS.B; dmin = np.inf; worst = None
for s in range(0, nV, 4000):
    yy = Y[s:s+4000, None]; v = Bc[None, :] - A[None, :]
    t = np.clip(((yy - A[None, :]).conjugate()*v).real/np.abs(v)**2, 0, 1)
    dd = np.abs(A[None, :] + t*v - yy)
    for j in range(dd.shape[0]):
        if bp[s + j] >= 0: dd[j, bp[s + j]] = np.inf
    m = dd.min()
    if m < dmin: dmin = m; worst = s + int(np.unravel_index(dd.argmin(), dd.shape)[0])
print(f'vertices {nV} | cuts {len(CS.tails)} | min distance vertex-cut {dmin:.3e} (vertex {worst}) | tau product round y0 = id (asserted)', flush=True)
assert dmin > 1e-12, 'mesh vertex on a cut: STOP (rotate the mesh)'
edges = []; ncross = 0
for p, q in M['E']:
    if bp[p] >= 0 and bp[q] >= 0: raise SystemExit('branch-branch edge')
    if bp[q] >= 0 or bp[p] >= 0:
        a, c = (p, q) if bp[q] >= 0 else (q, p)
        P = CS.edge_perm([Y[a], Y[c]], end_b=bp[c]); g, ng = CS.groups(bp[c])
        edges.append(dict(a=a, c=c, m=[g[P[l]] for l in range(d)], toGroup=True))
    else:
        P = CS.edge_perm(edge_path(Y[p], Y[q], []))
        ncross += P != list(range(d))
        edges.append(dict(a=p, c=q, m=P))
group = {}
for v in range(nV):
    if bp[v] >= 0:
        g, ng = CS.groups(bp[v]); group[str(v)] = dict(g=g, ng=ng, perm=CS.perms[bp[v]])
json.dump(dict(R=[None]*nV, d=d, group=group, edges=edges), open(os.path.join(DIR, 'edge_perms.json'), 'w'))
print(f'edges {len(edges)} | normal edges with a non-identity map {ncross} | branch vertices {len(group)}, points over them '
      f'{sorted(set(v["ng"] for v in group.values()))} | {time.time()-t0:.1f} s | wrote {DIR}/edge_perms.json', flush=True)
