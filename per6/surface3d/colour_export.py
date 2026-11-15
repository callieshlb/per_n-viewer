#!/usr/bin/env python3
"""E92 3-D step 4: per surface vertex: sphere point, root x, sheet, dynamics colour; markers. Writes
surface_data.js (base64 typed arrays) for the 3-D page."""
import base64
import json
import numpy as np
M = json.load(open('sphere_mesh.json'))
P = json.load(open('edge_perms.json'))
S = json.load(open('surface_mesh.json'))
D = json.load(open('../explorer/data.json'))
Y = np.array([complex(*y) for y in M['y']])
r2 = np.abs(Y)**2
U = np.stack([2*Y.real, 2*Y.imag, r2-1], -1)/(r2+1)[:, None]
XV = np.array(S['V'])                     # (sphere vertex, root index | -1-group)
T = np.array(S['T'])
nX = len(XV)
x = np.full(nX, np.nan+0j)
for k, (v, i) in enumerate(XV):
    if i >= 0:
        x[k] = complex(*P['R'][v][i])
br = np.where(XV[:, 1] < 0)[0]
# branch-vertex points: x from the adjacent vertices (the colliding roots meet there)
nbr = {int(k): [] for k in br}
for a, b, c in T:
    for p, q in ((a, b), (b, c), (c, a), (b, a), (c, b), (a, c)):
        if p in nbr and XV[q, 1] >= 0:
            nbr[p].append(q)
for k in br:
    x[k] = np.mean(x[nbr[int(k)]])
y = Y[XV[:, 0]]
sheet = np.where(XV[:, 1] >= 0, XV[:, 1], -1)
# dynamics (same rules as the explorer: relative capture radius 1e-3, 600 iterates, period tol 1e-7)
def chord(a, b):
    return np.abs(a-b)/np.sqrt((1+np.abs(a)**2)*(1+np.abs(b)**2))
def chord_inf(a):
    return 1/np.sqrt(1+np.abs(a)**2)
with np.errstate(all='ignore'):
    w = (x+y-1)/(x-1)
    s, p = x+w, x*w
    def f(z):
        big = ~np.isfinite(z) | (np.abs(z) > 1e150)
        zz = np.where(big, 1.0, z)
        zz = np.where(np.abs(zz) < 1e-150, 1e-150, zz)
        out = 1 - s/zz + p/(zz*zz)
        return np.where(big, 1.0+0j, out)
    om = 2*p/s
    a4 = f(y)
    marked = [np.zeros_like(x), None, np.ones_like(x), y, a4, x]
    def dist(z, a):
        return chord_inf(z) if a is None else chord(z, a)
    scale = np.ones(nX)
    th = []
    for j, a in enumerate(marked):
        m = np.ones(nX)
        for i2, b in enumerate(marked):
            if i2 == j:
                continue
            dd = chord_inf(a) if b is None and a is not None else (chord_inf(b) if a is None and b is not None else (np.zeros(nX) if a is None else chord(a, b)))
            m = np.minimum(m, dd)
        th.append(1e-3*m)
    z = om.copy()
    done = np.zeros(nX, bool)
    phase = np.full(nX, -1)
    ctime = np.zeros(nX, int)
    for it in range(600):
        z = f(z)
        for j, a in enumerate(marked):
            hit = (~done) & (dist(z, a) < th[j])
            phase[hit] = (j-(it+1)) % 6
            ctime[hit] = it+1
            done |= hit
    tail = [z]
    for it in range(24):
        z = f(z); tail.append(z)
    per = np.zeros(nX, int)
    for pp in range(1, 13):
        hit = (~done) & (chord(tail[-1], tail[-1-pp]) < 1e-7)
        per[hit] = pp
        done |= hit
cap = np.array([[42, 120, 214], [235, 104, 52], [27, 175, 122], [237, 161, 0], [232, 123, 164], [127, 63, 191]], float)
col = np.zeros((nX, 3))
col[:] = [18, 18, 18]
m = phase >= 0
sh = np.clip(1 - 0.6*np.log1p(ctime[m])/np.log1p(80), 0.3, 1)[:, None]
col[m] = cap[phase[m]]*sh + 20*(1-sh)
col[per == 1] = [85, 85, 85]
col[per >= 2] = [163, 163, 163]
for k in br:     # branch points: average of neighbours
    col[k] = col[nbr[int(k)]].mean(0)
print('surface vertices', nX, '| captured', int(m.sum()), '| fixed-point', int((per == 1).sum()), '| other cycle', int((per >= 2).sum()),
      '| undecided', int((~done).sum()), '| non-finite x', int((~np.isfinite(x)).sum()))
# markers: nearest surface vertex with matching root
def marker(yy, xx):
    yy, xx = complex(*yy), complex(*xx)
    v = int(np.argmin(np.abs(U - (np.array([2*yy.real, 2*yy.imag, abs(yy)**2-1])/(abs(yy)**2+1))).sum(1)))
    cand = np.where(XV[:, 0] == v)[0]
    k = cand[np.argmin(chord(x[cand], xx))]
    return int(k), float(chord(x[k], xx))
cm = [dict(k=marker(c['y'], c['x'])[0], label=str(c['label']), sheet=c['sheet']) for c in D['centers']]
pm = [dict(k=marker(q['y'], q['x'])[0], label=q['name'], sheet=q['sheet']) for q in D['punctures']]
print('marker match (chordal x) max: centres', max(marker(c['y'], c['x'])[1] for c in D['centers']),
      '| punctures', max(marker(q['y'], q['x'])[1] for q in D['punctures']))
b64 = lambda a: base64.b64encode(a.tobytes()).decode()
out = dict(nS=len(Y), nX=nX, U=b64(U.astype(np.float32)), sv=b64(XV[:, 0].astype(np.uint32)),
           x=b64(np.stack([x.real, x.imag], -1).astype(np.float64)), y=b64(np.stack([y.real, y.imag], -1).astype(np.float64)),
           sheet=b64(sheet.astype(np.int8)), col=b64(col.astype(np.uint8)), T=b64(T.astype(np.uint32)), nT=len(T),
           centers=cm, punctures=pm)
open('surface_data.js', 'w').write('const SURF = ' + json.dumps(out) + ';\n')
import os
print('surface_data.js', round(os.path.getsize('surface_data.js')/1e6, 2), 'MB')
