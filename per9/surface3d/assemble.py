#!/usr/bin/env python3
"""E93 mesh step 3 (E92 assemble.py, 5 -> d sheets): lift the sphere triangulation to the d = 10 sheets; check closed,
orientable, chi = 2 - 2g (A1 second route). usage: assemble.py DIR. Surface face d*t + i lies over sphere face t."""
import json
import os
import sys
import numpy as np
DIR = sys.argv[1] if len(sys.argv) > 1 else "."
os.chdir(DIR)
M = json.load(open('sphere_mesh.json'))
P = json.load(open('edge_perms.json'))
d = P.get('d') or len(next(r for r in P['R'] if r is not None))   # E95: glue.py writes 'd'
nV = len(M['y'])
isB = np.array(M['bp']) >= 0
emap = {}
for (p, q), e in zip(M['E'], P['edges']):
    emap[(p, q)] = e
ng = {int(k): v['ng'] for k, v in P['group'].items()}
def m_normal(n, m):          # root i at normal n -> root at normal m
    e = emap[(min(n, m), max(n, m))]
    if e['a'] == n:
        return e['m']
    inv = [0]*d
    for i, j in enumerate(e['m']):
        inv[j] = i
    return inv
def m_group(n, b):           # root i at normal n -> group at branch vertex b
    e = emap[(min(n, b), max(n, b))]
    assert e.get('toGroup') and e['a'] == n
    return e['m']
vid = {}
def X(v, i):
    k = (v, i)
    if k not in vid:
        vid[k] = len(vid)
    return vid[k]
tris, closure_fail, group_fail = [], 0, 0
for t in M['T']:
    t = list(t)
    bpos = [i for i in range(3) if isB[t[i]]]
    assert len(bpos) <= 1
    if bpos:
        r = (bpos[0]+1) % 3
        t = t[r:]+t[:r]          # rotate so the branch vertex is last (orientation kept)
    n0, n1, n2 = t
    for i in range(d):
        j = m_normal(n0, n1)[i]
        if not isB[n2]:
            k = m_normal(n1, n2)[j]
            if m_normal(n2, n0)[k] != i:
                closure_fail += 1
            tris.append((X(n0, i), X(n1, j), X(n2, k)))
        else:
            g = m_group(n0, n2)[i]
            if m_group(n1, n2)[j] != g:
                group_fail += 1
            tris.append((X(n0, i), X(n1, j), X(n2, ('g', g))))
tris = np.array(tris)
de = {}
for a, b, c in tris:
    for p, q in ((a, b), (b, c), (c, a)):
        de[(p, q)] = de.get((p, q), 0) + 1
und = {}
for (p, q), cnt in de.items():
    und[(min(p, q), max(p, q))] = und.get((min(p, q), max(p, q)), 0) + cnt
bad_count = sum(1 for v in und.values() if v != 2)
bad_orient = sum(1 for (p, q), cnt in de.items() if cnt != 1 or de.get((q, p), 0) != 1)
V, E, F = len(vid), len(und), len(tris)
exp_V = d*int((~isB).sum()) + sum(ng.values())
print(f'surface mesh: V={V} (expected {exp_V}) E={E} F={F} chi={V-E+F} genus={(2-(V-E+F))/2}')
print(f'triangle closure failures {closure_fail}, branch-group mismatches {group_fail}, '
      f'edges not in exactly 2 triangles {bad_count}, directed edges breaking orientation {bad_orient}')
inv = [None]*V
for k, i in vid.items():
    inv[i] = k
json.dump(dict(V=[[int(v), (i if isinstance(i, int) else -1-int(i[1]))] for v, i in inv], T=tris.tolist()),
          open('surface_mesh.json', 'w'))
