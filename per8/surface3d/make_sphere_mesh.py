#!/usr/bin/env python3
"""E94 mesh step 1 (E93 make_sphere_mesh.py with the refinement read from branch.json, as registered for n = 8).
usage: make_sphere_mesh.py fine|coarse [count]   (writes sphere_mesh.json in . (fine) or coarse/; 'count' only prints sizes)
- icosphere level 5 (fine) / 4 (coarse), randomly rotated (seed 7) as E93;
- hub grid (hexagonal, 0.012 fine / 0.024 coarse): grid points within 0.1 of a candidate whose nearest neighbour is < 0.05;
- fine grids (spacing nn/4 fine, nn/2 coarse) within 3 nn of every candidate with nn < 0.012;
- rings of 12 round every non-node branch point: r0 = min(0.003, 0.25 nn), radii r0 1.6^m (m <= 6 fine, 4 coarse), radius <= 0.4 nn,
  ring m rotated by 0.37 m; points within 0.5 r0 of a branch point removed; branch-vertex offset (edge_perms) 0.3 r0;
- nodes (../nodes.json): no vertex within 'excl', one ring of 12 at 'ring'; paths detour at 'detour' (edge_perms, meshloc)."""
import json
import os
import sys
import numpy as np
from scipy.spatial import ConvexHull
HERE = os.path.dirname(os.path.abspath(__file__))
mode = sys.argv[1] if len(sys.argv) > 1 else 'fine'
FINE = mode == 'fine'
out_dir = HERE if FINE else os.path.join(HERE, 'coarse')
B = json.load(open(os.path.join(HERE, '..', 'branch.json')))
NOD = json.load(open(os.path.join(HERE, '..', 'nodes.json')))
bps = np.array([complex(*b) for b in B['bps']]); kinds = B['kinds']
isbranch = np.array([k == 'branch' for k in kinds])
nn = np.array([np.abs(np.delete(bps, i) - b).min() for i, b in enumerate(bps)])
def ico(level):
    t = (1+5**0.5)/2
    V = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t), (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    V = [np.array(v)/np.linalg.norm(v) for v in V]
    F = [(0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),(1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),(3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),(4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1)]
    for _ in range(level):
        cache, nF = {}, []
        def mid(a, b):
            k = (min(a, b), max(a, b))
            if k not in cache:
                m = V[a]+V[b]; V.append(m/np.linalg.norm(m)); cache[k] = len(V)-1
            return cache[k]
        for a, b, c in F:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            nF += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        F = nF
    return np.array(V)
rng = np.random.default_rng(7)
Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
U = ico(5 if FINE else 4) @ Q.T
def to_y(u): return (u[:, 0] + 1j*u[:, 1])/(1 - u[:, 2])
def to_u(y):
    y = np.asarray(y); r2 = np.abs(y)**2
    return np.stack([2*y.real, 2*y.imag, r2-1], -1)/(r2+1)[:, None]
def hexgrid(h, x0, x1, y0_, y1_):
    xs = np.arange(x0, x1, h); ys = np.arange(y0_, y1_, h*np.sqrt(3)/2)
    X, Yg = np.meshgrid(xs, ys); X = X + (np.arange(len(ys)) % 2)[:, None]*h/2
    return (X + 1j*Yg).ravel()
counts = {'ico': len(U)}
pts = [to_y(U)]
s = 1 if FINE else 2
cl = bps[nn < 0.05]
g = hexgrid(0.012*s, cl.real.min() - 0.1, cl.real.max() + 0.1, cl.imag.min() - 0.1, cl.imag.max() + 0.1)
g = g[np.abs(g[:, None] - cl[None, :]).min(1) < 0.1]
pts.append(g); counts['hub'] = len(g)
nf = 0
for b, n_ in zip(bps, nn):
    if n_ < 0.012:
        h = n_/4*s; g = hexgrid(h, b.real - 3*n_, b.real + 3*n_, b.imag - 3*n_, b.imag + 3*n_)
        g = g[np.abs(g - b) < 3*n_]; pts.append(g); nf += len(g)
counts['fine grids'] = nf
r0 = np.minimum(0.003, 0.25*nn)
M = 7 if FINE else 5
ring = []
for i in np.where(isbranch)[0]:
    for m in range(M):
        r = r0[i]*1.6**m
        if r > 0.4*nn[i]: break
        ring += list(bps[i] + r*np.exp(1j*(2*np.pi*np.arange(12)/12 + 0.37*m)))
for q in NOD:
    c = complex(*q['c']); ring += list(c + q['ring']*np.exp(1j*(2*np.pi*np.arange(12)/12 + 0.11)))
pts.append(np.array(ring)); counts['rings'] = len(ring)
pts = np.concatenate(pts)
keep = np.ones(len(pts), bool)
for i in np.where(isbranch)[0]:
    keep &= np.abs(pts - bps[i]) > 0.5*r0[i]
for q in NOD:
    keep &= np.abs(pts - complex(*q['c'])) > q['excl'] - 1e-12
pts = pts[keep]
# E94: overlapping grids produce duplicates / near-duplicates (qhull leaves them unused: chi 2820 in the first build); keep one
# point per cluster of points closer than 1e-7, and none within 1e-7 of a branch point.
from scipy.spatial import cKDTree
kd = cKDTree(np.stack([pts.real, pts.imag], 1)); pr = kd.query_pairs(1e-7, output_type='ndarray')
drop = np.zeros(len(pts), bool); drop[pr[:, 1]] = True
pts = pts[~drop]
pts = pts[np.abs(pts[:, None] - bps[isbranch][None, :]).min(1) > 1e-7]
counts['removed near-duplicates'] = int(drop.sum())
allp = np.concatenate([pts, bps[isbranch]])
is_bp = np.concatenate([np.zeros(len(pts), int)-1, np.where(isbranch)[0]])
counts['total vertices'] = len(allp)
if len(sys.argv) > 2 and sys.argv[2] == 'count':
    print(mode, counts); sys.exit()
os.makedirs(out_dir, exist_ok=True)
Uall = to_u(allp)
hull = ConvexHull(Uall)
T = hull.simplices.copy()
c = Uall[T].mean(1); nrm = np.cross(Uall[T[:, 1]]-Uall[T[:, 0]], Uall[T[:, 2]]-Uall[T[:, 0]])
flip = (nrm*c).sum(1) < 0
T[flip] = T[flip][:, [0, 2, 1]]
unused = len(allp) - len(np.unique(T))
assert unused == 0, f'{unused} vertices unused by the hull'
Es = np.sort(np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]]), axis=1)
Eu = np.unique(Es, axis=0)
bb = int(((is_bp[Eu[:, 0]] >= 0) & (is_bp[Eu[:, 1]] >= 0)).sum())
el = np.linalg.norm(Uall[Eu[:, 0]] - Uall[Eu[:, 1]], axis=1)
dn = min(float(np.abs(allp - complex(*q['c'])).min()) for q in NOD)
print(f'{mode} sphere mesh: V {len(allp)} E {len(Eu)} F {len(T)} chi {len(allp) - len(Eu) + len(T)} | {counts} | branch-point vertices '
      f'{int((is_bp >= 0).sum())} | branch-branch edges {bb} | chordal edge length min {el.min():.2e} max {el.max():.2e} | min distance to a node {dn:.4f}')
eps = {int(v): float(0.3*r0[is_bp[v]]) for v in np.where(is_bp >= 0)[0]}
json.dump(dict(y=[[z.real, z.imag] for z in allp], bp=is_bp.tolist(), T=T.tolist(), E=Eu.tolist(), eps=eps),
          open(os.path.join(out_dir, 'sphere_mesh.json'), 'w'))
