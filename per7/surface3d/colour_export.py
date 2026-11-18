#!/usr/bin/env python3
"""E93 (E92 colour_export.py, n = 7): per surface vertex of the fine glued mesh: y = a_3, x = a_6 (root from edge_perms.json;
branch-vertex points: mean of the adjacent lifts), and the dynamics colour by the 2-D explorer's rule (fate of the free
critical point omega under f(z) = 1 - s/z + p/z^2; capture by the marked cycle {0, inf, 1, a_3, a_4, a_5, a_6} with relative
radius 1e-3, 600 iterates; phase (j - t) mod 7; else attracting cycle by period (tol 1e-7); else undecided).
Writes surface_data.json (base64 float64 y, x; uint8 rgb col; int8 root index per vertex, -1 at branch points)."""
import base64
import json
import numpy as np
M = json.load(open('sphere_mesh.json'))
P = json.load(open('edge_perms.json'))
S = json.load(open('surface_mesh.json'))
NM = 7
CAP = np.array([[42, 120, 214], [235, 104, 52], [27, 175, 122], [237, 161, 0], [232, 123, 164], [127, 63, 191], [160, 82, 45]], float)
Y = np.array([complex(*y) for y in M['y']])
XV = np.array(S['V']); T = np.array(S['T'])
nX = len(XV)
x = np.full(nX, np.nan + 0j)
for k, (v, i) in enumerate(XV):
    if i >= 0:
        x[k] = complex(*P['R'][v][i])
br = np.where(XV[:, 1] < 0)[0]
nbr = {int(k): [] for k in br}
for a, b, c in T:
    for p, q in ((a, b), (b, c), (c, a), (b, a), (c, b), (a, c)):
        if p in nbr and XV[q, 1] >= 0:
            nbr[p].append(q)
for k in br:
    x[k] = np.mean(x[nbr[int(k)]])
y = Y[XV[:, 0]]
def chord(a, b): return np.abs(a - b)/np.sqrt((1 + np.abs(a)**2)*(1 + np.abs(b)**2))
def chord_inf(a): return 1/np.sqrt(1 + np.abs(a)**2)
with np.errstate(all='ignore'):
    w = (x + y - 1)/(x - 1)
    s, p = x + w, x*w
    def f(z):
        big = ~np.isfinite(z) | (np.abs(z) > 1e150)
        zz = np.where(big, 1.0, z)
        zz = np.where(np.abs(zz) < 1e-150, 1e-150, zz)
        out = 1 - s/zz + p/(zz*zz)
        return np.where(big, 1.0 + 0j, out)
    om = 2*p/s
    marked = [np.zeros_like(x), None, np.ones_like(x), y]
    for k in range(4, NM - 1):
        marked.append(f(marked[-1]))
    marked.append(x)
    def dist(z, a):
        return chord_inf(z) if a is None else chord(z, a)
    th = []
    for j, a in enumerate(marked):
        m = np.ones(nX)
        for i2, b in enumerate(marked):
            if i2 == j: continue
            dd = chord_inf(a) if b is None else (chord_inf(b) if a is None else chord(a, b))
            m = np.minimum(m, dd)
        th.append(1e-3*m)
    z = om.copy()
    done = np.zeros(nX, bool); phase = np.full(nX, -1); ctime = np.zeros(nX, int)
    for it in range(600):
        z = f(z)
        for j, a in enumerate(marked):
            hit = (~done) & (dist(z, a) < th[j])
            phase[hit] = (j - (it + 1)) % NM; ctime[hit] = it + 1; done |= hit
    tail = [z]
    for it in range(24):
        z = f(z); tail.append(z)
    per = np.zeros(nX, int)
    for pp in range(1, 13):
        hit = (~done) & (chord(tail[-1], tail[-1 - pp]) < 1e-7)
        per[hit] = pp; done |= hit
col = np.zeros((nX, 3))
cap = phase >= 0
sh = np.clip(1 - 0.6*np.log1p(ctime)/np.log1p(80), 0.3, 1)
col[cap] = CAP[phase[cap]]*sh[cap, None] + 20*(1 - sh[cap, None])
col[per == 1] = 85; col[per > 1] = 163
col[~done] = 18
col = np.clip(col, 0, 255).astype(np.uint8)
print('vertices', nX, '| captured', int(cap.sum()), '| period 1', int((per == 1).sum()), '| other cycles', int((per > 1).sum()),
      '| undecided', int((~done).sum()), '| phase histogram', np.bincount(phase[cap], minlength=NM).tolist())
b64 = lambda a: base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()
json.dump(dict(y=b64(np.stack([y.real, y.imag], 1).astype(np.float64)), x=b64(np.stack([x.real, x.imag], 1).astype(np.float64)),
               col=b64(col), root=b64(XV[:, 1].clip(-1, 127).astype(np.int8))), open('surface_data.json', 'w'))
print('wrote surface_data.json')
