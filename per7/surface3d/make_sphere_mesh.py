#!/usr/bin/env python3
"""E93 mesh step 1 (generalizes E92 make_sphere_mesh.py): triangulate the y-sphere (stereographic, north pole = y = inf).
usage: make_sphere_mesh.py fine|coarse  (writes sphere_mesh.json in . (fine) or coarse/ (coarse))
- icosphere (level 5 fine, 4 coarse), randomly rotated (seed 7): no vertex at a pole, no conjugation symmetry;
- hub grid (hexagonal, spacing 0.012 fine / 0.024 coarse) on the cluster of branch points Re y in [0.2, 1.1], |Im y| <= 0.25;
- finer grids (0.002 / 0.004) round the three tightest groups (nearest-neighbour distances 0.0047 .. 0.0117);
- rings of 12 points round every non-node branch point, radii r0 * 1.6^m (m = 0..6 fine, 0..4 coarse), r0 = 0.0015 when the
  nearest other candidate is closer than 0.01, else 0.003; rings with radius > 0.4 x (nearest-candidate distance) dropped;
  ring m rotated by 0.37 m;
- the node y = 0 (trivial monodromy, two 5-clusters tangent to order 4) is NOT a vertex: no vertex within 0.03 of it; rings
  of radius 0.03 * 1.6^m round it (radius <= 0.4 x 0.246);
- every other point within 0.5 r0 of a branch point removed; every non-node branch point is a vertex. Pole points (one root
  -> x = inf, trivial lasso) get no special treatment."""
import json
import os
import sys
import numpy as np
from scipy.spatial import ConvexHull
HERE = os.path.dirname(os.path.abspath(__file__))
mode = sys.argv[1] if len(sys.argv) > 1 else 'fine'
FINE = mode == 'fine'
out_dir = HERE if FINE else os.path.join(HERE, 'coarse')
os.makedirs(out_dir, exist_ok=True)
B = json.load(open(os.path.join(HERE, '..', 'branch.json')))
bps = np.array([complex(*b) for b in B['bps']]); kinds = B['kinds']
isbranch = np.array([k == 'branch' for k in kinds]); isnode = np.array([k == 'node' for k in kinds])
nodes = bps[isnode]
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
def to_y(u):
    return (u[:, 0] + 1j*u[:, 1])/(1 - u[:, 2])
def to_u(y):
    y = np.asarray(y); r2 = np.abs(y)**2
    return np.stack([2*y.real, 2*y.imag, r2-1], -1)/(r2+1)[:, None]
pts = list(to_y(U))
def hexgrid(h, x0, x1, y0_, y1_):
    out = []
    for i, yy in enumerate(np.arange(y0_, y1_, h*np.sqrt(3)/2)):
        for xx in np.arange(x0, x1, h):
            out.append(xx + (h/2 if i % 2 else 0) + 1j*yy)
    return out
s = 1 if FINE else 2
pts += hexgrid(0.012*s, 0.2, 1.1, -0.25, 0.25)
for box in ((0.48, 0.53, -0.025, 0.025), (0.58, 0.61, -0.02, 0.02), (0.28, 0.32, -0.02, 0.02)):
    pts += hexgrid(0.002*s, *box)
nn = np.array([np.abs(np.delete(bps, i) - b).min() for i, b in enumerate(bps)])
r0 = np.where(nn < 0.01, 0.0015, 0.003)
M = 7 if FINE else 5
nring = 0
for i in np.where(isbranch)[0]:
    for m in range(M):
        r = r0[i]*1.6**m
        if r > 0.4*nn[i]: break
        for a in range(12):
            pts.append(bps[i] + r*np.exp(1j*(2*np.pi*a/12 + 0.37*m))); nring += 1
for c in nodes:
    nnc = np.abs(bps[~isnode] - c).min()
    for m in range(M):
        r = 0.03*1.6**m
        if r > 0.4*nnc: break
        for a in range(12):
            pts.append(c + r*np.exp(1j*(2*np.pi*a/12 + 0.37*m))); nring += 1
pts = np.array(pts)
keep = np.ones(len(pts), bool)
for i in np.where(isbranch)[0]:
    keep &= np.abs(pts - bps[i]) > 0.5*r0[i]
for c in nodes:
    keep &= np.abs(pts - c) > 0.03 - 1e-12
pts = pts[keep]
allp = np.concatenate([pts, bps[isbranch]])
is_bp = np.concatenate([np.zeros(len(pts), int)-1, np.where(isbranch)[0]])
Uall = to_u(allp)
hull = ConvexHull(Uall)
T = hull.simplices.copy()
c = Uall[T].mean(1); nrm = np.cross(Uall[T[:, 1]]-Uall[T[:, 0]], Uall[T[:, 2]]-Uall[T[:, 0]])
flip = (nrm*c).sum(1) < 0
T[flip] = T[flip][:, [0, 2, 1]]
E = set()
for a, b_, c_ in T:
    for p, q in ((a, b_), (b_, c_), (c_, a)):
        E.add((min(p, q), max(p, q)))
V, Fn, En = len(allp), len(T), len(E)
# branch-branch edges are not allowed (the gluing needs a normal vertex on every edge at a branch point)
bb = sum(1 for p, q in E if is_bp[p] >= 0 and is_bp[q] >= 0)
# min edge length (chordal) and min distance of any vertex to the node
Ea = np.array(sorted(E)); el = np.linalg.norm(Uall[Ea[:, 0]] - Uall[Ea[:, 1]], axis=1)
print(f'{mode} sphere mesh: V {V} E {En} F {Fn} chi {V-En+Fn} | ring points {nring} | branch-point vertices {int((is_bp >= 0).sum())} '
      f'| branch-branch edges {bb} | chordal edge length min {el.min():.2e} max {el.max():.2e} '
      f'| min |y - node| {np.abs(allp[:, None] - nodes[None, :]).min():.4f}')
json.dump(dict(y=[[z.real, z.imag] for z in allp], bp=is_bp.tolist(), T=T.tolist(), E=[[int(p), int(q)] for p, q in sorted(E)]),
          open(os.path.join(out_dir, 'sphere_mesh.json'), 'w'))
