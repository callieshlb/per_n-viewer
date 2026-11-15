#!/usr/bin/env python3
"""E92 3-D step 1: triangulate the y-sphere (stereographic, north pole = y = inf), refined near the
branch points and the hub, with every branch point (except the node y = 0) as a vertex."""
import json
import numpy as np
from scipy.spatial import ConvexHull
B = json.load(open('../branch.json'))
bps = np.array([complex(*b) for b in B['bps']])
node = np.abs(bps) < 1e-9
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
U = ico(5) @ Q.T                                     # randomly rotated: no vertex at a pole
def to_y(u):
    return (u[:, 0] + 1j*u[:, 1])/(1 - u[:, 2])
def to_u(y):
    y = np.asarray(y); r2 = np.abs(y)**2
    return np.stack([2*y.real, 2*y.imag, r2-1], -1)/(r2+1)[:, None]
pts = list(to_y(U))
# hub refinement: |y - 0.3| < 0.75, spacing 0.012 (hexagonal)
h = 0.012
for i, yy in enumerate(np.arange(-0.75, 0.75, h*np.sqrt(3)/2)):
    for xx in np.arange(-0.75, 0.75, h):
        z = 0.3 + xx + (h/2 if i % 2 else 0) + 1j*yy
        if abs(z-0.3) < 0.75: pts.append(z)
# fine grids where the 0.361 slit, the 0.508 slit and the hub tree pass ~0.004 apart
for hh, (x0_, x1_, y0_, y1_) in ((0.0015, (0.34, 0.62, -0.13, 0.13)), (0.0007, (0.47, 0.53, -0.04, 0.04))):
    for i, yy in enumerate(np.arange(y0_, y1_, hh*np.sqrt(3)/2)):
        for xx in np.arange(x0_, x1_, hh):
            pts.append(xx + (hh/2 if i % 2 else 0) + 1j*yy)
# rings around every branch point, and around the node (which is not a vertex)
for b in bps[~node]:
    for m in range(7):
        r = 0.003*1.6**m
        for a in range(12):
            pts.append(b + r*np.exp(1j*(2*np.pi*a/12 + 0.37*m)))
pts = np.array(pts)
keep = np.ones(len(pts), bool)
for b in bps:
    keep &= np.abs(pts-b) > 0.0006
keep &= np.abs(pts) > 0.012          # the node y = 0 needs no vertices (trivial monodromy)
pts = pts[keep]
allp = np.concatenate([pts, bps[~node]])
is_bp = np.concatenate([np.zeros(len(pts), int)-1, np.where(~node)[0]])
Uall = to_u(allp)
hull = ConvexHull(Uall)
T = hull.simplices
# orient triangles outward
c = Uall[T].mean(1); n = np.cross(Uall[T[:, 1]]-Uall[T[:, 0]], Uall[T[:, 2]]-Uall[T[:, 0]])
flip = (n*c).sum(1) < 0
T[flip] = T[flip][:, [0, 2, 1]]
E = set()
for a, b_, c_ in T:
    for p, q in ((a, b_), (b_, c_), (c_, a)):
        E.add((min(p, q), max(p, q)))
V, Fn, En = len(allp), len(T), len(E)
print('sphere mesh: V', V, 'E', En, 'F', Fn, 'chi', V-En+Fn, '| branch-point vertices', int((is_bp >= 0).sum()),
      '| min |y| (node at 0 is not a vertex)', float(np.abs(allp).min()))
json.dump(dict(y=[[z.real, z.imag] for z in allp], bp=is_bp.tolist(), T=T.tolist(), E=[[int(p), int(q)] for p, q in sorted(E)]),
          open('sphere_mesh.json', 'w'))
