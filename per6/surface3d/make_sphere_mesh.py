#!/usr/bin/env python3
"""E99 (copied from E92 3-D step 1): triangulate the y-sphere (stereographic, north pole = y = inf), refined near the
branch points and the hub, with every branch point (except the node y = 0) as a vertex.
E99 addition: geometric rings about the 7 puncture centres (REGISTRATION.md table), E92 points inside each ring disk
removed except the rings of other branch points; centre vertices added at the unramified finite centres."""
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
pts = list(to_y(U)); src_lab = [-1]*len(pts)          # -1 ico/grid, b >= 0: ring of branch point b
# hub refinement: |y - 0.3| < 0.75, spacing 0.012 (hexagonal)
h = 0.012
for i, yy in enumerate(np.arange(-0.75, 0.75, h*np.sqrt(3)/2)):
    for xx in np.arange(-0.75, 0.75, h):
        z = 0.3 + xx + (h/2 if i % 2 else 0) + 1j*yy
        if abs(z-0.3) < 0.75: pts.append(z); src_lab.append(-1)
# fine grids where the 0.361 slit, the 0.508 slit and the hub tree pass ~0.004 apart
for hh, (x0_, x1_, y0_, y1_) in ((0.0015, (0.34, 0.62, -0.13, 0.13)), (0.0007, (0.47, 0.53, -0.04, 0.04))):
    for i, yy in enumerate(np.arange(y0_, y1_, hh*np.sqrt(3)/2)):
        for xx in np.arange(x0_, x1_, hh):
            pts.append(xx + (hh/2 if i % 2 else 0) + 1j*yy); src_lab.append(-1)
# rings around every branch point, and around the node (which is not a vertex)
for bi, b in zip(np.where(~node)[0], bps[~node]):
    for m in range(7):
        r = 0.003*1.6**m
        for a in range(12):
            pts.append(b + r*np.exp(1j*(2*np.pi*a/12 + 0.37*m))); src_lab.append(int(bi))
pts = np.array(pts); src_lab = np.array(src_lab)
keep = np.ones(len(pts), bool)
for b in bps:
    keep &= np.abs(pts-b) > 0.0006
keep &= np.abs(pts) > 0.012          # the node y = 0 needs no vertices (trivial monodromy)
pts = pts[keep]; src_lab = src_lab[keep]
n92 = len(pts)
# ---- E99: rings about the puncture centres (REGISTRATION.md). coordinate w = y - c, or w = 1/y for c = inf.
CEN = [  # name, centre (None = inf), e, r_out, N, r_in, extra inner ring (radius, n) or None, add centre vertex
    ('1/2', 0.5, 3, 3e-3, 28, 1e-9, None, False),
    ('1', 1.0, 3, 5e-2, 40, 1e-9, None, False),
    ('0', 0.0, 1, 1e-1, 72, 5e-4, (1.5e-4, 6), False),
    ('2/3', 2/3, 1, 4e-2, 56, 1e-4, None, True),
    ('(sqrt5-1)/2', (5**0.5 - 1)/2, 1, 1.2e-2, 56, 1e-4, None, True),
    ('-(sqrt5+1)/2', -(5**0.5 + 1)/2, 1, 2.5e-1, 48, 1e-4, None, True),
    ('inf', None, 1, 1e-1, 48, 1e-4, (3e-5, 6), False)]
def coord(y, c): return 1/y if c is None else y - c
def uncoord(w, c): return 1/w if c is None else c + w
newp, newlab = [], []
for name, c, e, r_out, N, r_in, extra, cv in CEN:
    w = coord(pts, c)
    inside = np.abs(w) < r_out
    own = np.zeros(len(pts), bool)
    if c is not None:                          # the centre's own E92 rings (centre is a branch point) are replaced too
        for bi in np.where(~node)[0]:
            if abs(bps[bi] - c) < 1e-12: own |= src_lab == bi
    drop = inside & ((src_lab < 0) | own)
    kept_in = inside & ~drop
    rho = (1 - 2*np.pi/(e*N))**e
    rings = []; k = 0
    while r_out*rho**k >= r_in*(1 - 1e-12):
        rings.append((r_out*rho**k, N, 0.123 + 0.37*k + np.pi*k/N)); k += 1
    if extra is not None: rings.append((extra[0], extra[1], 0.123 + 0.37*k))
    kw = w[kept_in]; ring_pts = []; ndrop = 0
    for r, n, ph in rings:
        q = r*np.exp(1j*(ph + 2*np.pi*np.arange(n)/n)); sp = 2*np.pi*r/n
        if len(kw):
            dmin = np.abs(q[:, None] - kw[None, :]).min(1); ndrop += int((dmin < 0.3*sp).sum()); q = q[dmin >= 0.3*sp]
        ring_pts.append(q)
    ring_pts = np.concatenate(ring_pts)
    if cv: ring_pts = np.concatenate([ring_pts, [0j]])
    print(f'centre {name}: e {e} r_out {r_out:g} N {N} rho {rho:.4f} rings {len(rings)} (radius {rings[0][0]:.3g} .. {rings[-1][0]:.3g}); '
          f'E92 points removed {int(drop.sum())}, kept inside (other rings) {int(kept_in.sum())}; new points {len(ring_pts)} '
          f'(dropped near kept {ndrop}){"; centre vertex added" if cv else ""}')
    pts = pts[~drop]; src_lab = src_lab[~drop]
    newp.append(uncoord(ring_pts, c))
pts = np.concatenate([pts] + newp)
assert np.all(np.isfinite(pts))
dmin_bp = min(float(np.abs(pts - b).min()) for b in bps[~node])
print(f'points: E92 {n92}, after removal {len(src_lab)}, with rings {len(pts)}; min distance of a non-branch point to a branch point {dmin_bp:.3e}')
allp = np.concatenate([pts, bps[~node]])
is_bp = np.concatenate([np.zeros(len(pts), int)-1, np.where(~node)[0]])
Uall = to_u(allp)
# E99: qhull cannot resolve the sphere's curvature at ring spacings << 1e-6 (sagitta ~ L^2/2 below round-off), so inside
# |y - c| < R_CUT of the two cone centres the points are triangulated locally (2-D Delaunay in w = (y - c)/R_CUT) and
# stitched to the hull along the innermost ring with radius >= R_CUT (a convex polygon: the hull's cap triangles on it,
# those with all three corners on that ring, are replaced).
from scipy.spatial import Delaunay
R_CUT = 1e-4
loc_sets = []
for c in (0.5, 1.0):
    w = allp - c; a = np.abs(w)
    inner = np.where(a < R_CUT)[0]
    rr = a[(a >= R_CUT) & (a < 1.0)]; r_ring = rr.min()
    ring = np.where(np.abs(a - r_ring) < 1e-9*r_ring)[0]
    loc_sets.append((c, inner, ring))
    print(f'local patch at {c}: inner points {len(inner)} (incl. branch vertex {bool((is_bp[inner] >= 0).any())}), stitching ring radius {r_ring:.4e} with {len(ring)} points')
excl = np.zeros(len(allp), bool)
for c, inner, ring in loc_sets: excl[inner] = True
gidx = np.where(~excl)[0]
hull = ConvexHull(Uall[gidx])
T = gidx[hull.simplices]
assert len(np.unique(T)) == len(gidx), ('hull dropped points', len(gidx) - len(np.unique(T)))
# orient triangles outward
c_ = Uall[T].mean(1); n_ = np.cross(Uall[T[:, 1]]-Uall[T[:, 0]], Uall[T[:, 2]]-Uall[T[:, 0]])
flip = (n_*c_).sum(1) < 0
T[flip] = T[flip][:, [0, 2, 1]]
def sarea(tri):                     # signed area in the y-plane (finite y only)
    p = allp[tri]; return ((p[:, 1] - p[:, 0]).conj()*(p[:, 2] - p[:, 0])).imag
for c, inner, ring in loc_sets:
    rs = set(ring.tolist())
    cap = np.array([all(int(v) in rs for v in t) for t in T])
    sgn = np.sign(np.median(sarea(T[np.isin(T, ring).any(1) & ~cap])))
    T = T[~cap]
    li = np.concatenate([ring, inner]); wl = (allp[li] - c)/R_CUT
    dl = Delaunay(np.stack([wl.real, wl.imag], 1))
    TL = li[dl.simplices]
    s = np.sign(sarea(TL)); TL[s != sgn] = TL[s != sgn][:, [0, 2, 1]]
    print(f'  patch {c}: removed {int(cap.sum())} cap triangles, added {len(TL)} local triangles (outward orientation sign in y {int(sgn)}); '
          f'min |signed area| / R_CUT^2 {np.abs(sarea(TL)).min()/R_CUT**2:.3e}')
    T = np.concatenate([T, TL])
dirc = {}
for a, b_, c_ in T:
    for p, q in ((a, b_), (b_, c_), (c_, a)): dirc[(p, q)] = dirc.get((p, q), 0) + 1
bad = sum(1 for (p, q), k in dirc.items() if k != 1 or dirc.get((q, p), 0) != 1)
print('sphere triangulation: directed edges not matched exactly once by the reverse edge:', bad)
E = set()
for a, b_, c_ in T:
    for p, q in ((a, b_), (b_, c_), (c_, a)):
        E.add((min(p, q), max(p, q)))
V, Fn, En = len(allp), len(T), len(E)
print('sphere mesh: V', V, 'E', En, 'F', Fn, 'chi', V-En+Fn, '| branch-point vertices', int((is_bp >= 0).sum()),
      '| min |y| (node at 0 is not a vertex)', float(np.abs(allp).min()))
json.dump(dict(y=[[z.real, z.imag] for z in allp], bp=is_bp.tolist(), T=T.tolist(), E=[[int(p), int(q)] for p, q in sorted(E)]),
          open('sphere_mesh.json', 'w'))
print('min |y| over vertices', float(np.abs(allp).min()), '| max |y|', float(np.abs(allp).max()))
