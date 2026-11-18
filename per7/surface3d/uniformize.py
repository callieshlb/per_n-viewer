#!/usr/bin/env python3
"""E93-U step 1 (E92 uniformize.py; genus and sheet count read from the data): discrete uniformization of compact Per_7(0)
(Bobenko-Pinkall-Springborn, hyperbolic target). Input: a closed triangle mesh of the genus-24 surface (surface_mesh.json + sphere_mesh.json in the working directory),
Euclidean edge lengths = chordal lengths of the y-sphere (pullback conformal structure).
Unknowns: u_i per vertex. Hyperbolic lengths: 2 sinh(l_ij/2) = exp((u_i+u_j)/2) * l0_ij.
Solve: angle sum = 2 pi at every vertex (damped Newton; Jacobian assembled from the hyperbolic-triangle angle derivatives
dA/da = sinh a/(sinh b sinh c sin A), dA/db = -dA/da cos C, dA/dc = -dA/da cos B, and dl/du_i = tanh(l/2)).
Usage: python3 uniformize.py [selftest] (run in a mesh directory). Writes uniform_u.npy, uniform_log.json."""
import os
import json
import sys
import time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
TAU = 2*np.pi
def angles(a, b, c):
    """angles opposite sides a, b, c of hyperbolic triangles (half-angle formula, stable for small triangles)."""
    s = (a + b + c)/2
    sa, sb, sc = s - a, s - b, s - c
    ok = (sa > 0) & (sb > 0) & (sc > 0)
    ss = np.sinh(s); ssa, ssb, ssc = np.sinh(np.maximum(sa, 0)), np.sinh(np.maximum(sb, 0)), np.sinh(np.maximum(sc, 0))
    with np.errstate(divide='ignore', invalid='ignore'):
        A = 2*np.arctan(np.sqrt(ssb*ssc/(ss*ssa)))
        B = 2*np.arctan(np.sqrt(ssc*ssa/(ss*ssb)))
        C = 2*np.arctan(np.sqrt(ssa*ssb/(ss*ssc)))
    return A, B, C, ok
def lengths(l0, u, E):
    return 2*np.arcsinh(0.5*l0*np.exp(0.5*(u[E[:, 0]] + u[E[:, 1]])))
def selftest():
    rng = np.random.default_rng(1)
    worst = 0
    for _ in range(200):
        a, b, c = rng.uniform(0.01, 2, 3)
        if not (a < b + c and b < a + c and c < a + b): continue
        A, B, C, _ = angles(np.array([a]), np.array([b]), np.array([c]))
        h = 1e-6
        for (da, db, dc), name in (((h, 0, 0), 'a'), ((0, h, 0), 'b'), ((0, 0, h), 'c')):
            Ap, _, _, _ = angles(np.array([a+da]), np.array([b+db]), np.array([c+dc]))
            Am, _, _, _ = angles(np.array([a-da]), np.array([b-db]), np.array([c-dc]))
            fd = (Ap - Am)/(2*h)
            dAda = np.sinh(a)/(np.sinh(b)*np.sinh(c)*np.sin(A))
            an = {'a': dAda, 'b': -dAda*np.cos(C), 'c': -dAda*np.cos(B)}[name]
            worst = max(worst, float(np.max(np.abs(fd - an)/np.maximum(1, np.abs(an)))))  # numpy 2: no float() of 1-element arrays
    print('selftest: max relative error of analytic angle derivatives vs finite differences', f'{worst:.2e}')
    # Jacobian of angle sums vs finite differences on a small random closed-ish patch is checked in solve() (symmetry)
    return worst
def setup(T, l0E=None, E=None):
    """edge ids of the side opposite each corner, consistent with edges(T)."""
    nF = len(T)
    Es = np.sort(np.concatenate([T[:, [1, 2]], T[:, [2, 0]], T[:, [0, 1]]]), axis=1)
    _, inv = np.unique(Es, axis=0, return_inverse=True)
    return inv.reshape(-1).reshape(3, nF).T
def edges(T):
    Es = np.sort(np.concatenate([T[:, [1, 2]], T[:, [2, 0]], T[:, [0, 1]]]), axis=1)
    return np.unique(Es, axis=0)
def residual(u, T, fe, l0, E, nV):
    l = lengths(l0, u, E); a, b, c = l[fe[:, 0]], l[fe[:, 1]], l[fe[:, 2]]
    A, B, C, ok = angles(a, b, c)
    S = np.bincount(T[:, 0], A, nV) + np.bincount(T[:, 1], B, nV) + np.bincount(T[:, 2], C, nV)
    return S - TAU, l, (a, b, c, A, B, C, ok)
def jacobian(T, fe, l, tri, E, nV):
    a, b, c, A, B, C, ok = tri
    th = np.tanh(l/2); ta, tb, tc = th[fe[:, 0]], th[fe[:, 1]], th[fe[:, 2]]
    sa_, sb_, sc_ = np.sinh(a), np.sinh(b), np.sinh(c)
    dA_a = sa_/(sb_*sc_*np.sin(A)); dB_b = sb_/(sc_*sa_*np.sin(B)); dC_c = sc_/(sa_*sb_*np.sin(C))
    dA = (dA_a, -dA_a*np.cos(C), -dA_a*np.cos(B))      # wrt a, b, c
    dB = (-dB_b*np.cos(C), dB_b, -dB_b*np.cos(A))
    dC = (-dC_c*np.cos(B), -dC_c*np.cos(A), dC_c)
    # side a = (j,k), b = (k,i), c = (i,j) for corners (i,j,k) = T[:,0..2]
    i, j, k = T[:, 0], T[:, 1], T[:, 2]
    rows, cols, vals = [], [], []
    for corner, d in ((i, dA), (j, dB), (k, dC)):
        da_, db_, dc_ = d
        # d(angle)/du_m = sum over sides containing m of d(angle)/d(side) * tanh(side/2)
        for m, v in ((i, db_*tb + dc_*tc), (j, da_*ta + dc_*tc), (k, da_*ta + db_*tb)):
            rows.append(corner); cols.append(m); vals.append(np.where(ok, v, 0.0))
    return sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nV, nV)).tocsc()
def solve(T, l0, E, nV, u0, tol=1e-11, maxit=60, log=print):
    fe = setup(T, l0, E)
    u = u0.copy(); hist = []
    F, l, tri = residual(u, T, fe, l0, E, nV)
    for it in range(maxit):
        nrm = float(np.abs(F).max()); bad = int((~tri[6]).sum())
        hist.append(dict(it=it, max_defect=nrm, l2=float(np.linalg.norm(F)), degenerate=bad))
        log(f'newton {it}: max |angle sum - 2pi| {nrm:.3e}  l2 {np.linalg.norm(F):.3e}  degenerate triangles {bad}')
        if nrm < tol: break
        t0 = time.time(); J = jacobian(T, fe, l, tri, E, nV)
        if it == 0:
            asym = abs(J - J.T).max(); log(f'jacobian symmetry |J - J^T| max {asym:.2e} (|J| max {abs(J).max():.2e})')
        du = splu((-J).tocsc(), permc_spec='COLAMD').solve(F)
        step = 1.0
        while True:
            un = u + step*du
            Fn, ln, trin = residual(un, T, fe, l0, E, nV)
            if trin[6].all() and np.linalg.norm(Fn) < (1 - 1e-4*step)*np.linalg.norm(F): break
            step /= 2
            if step < 1e-6: log('line search failed'); return u, hist, False
        u, F, l, tri = un, Fn, ln, trin
        log(f'   step {step:g}, max |du| {np.abs(step*du).max():.3e}, solve {time.time()-t0:.1f}s')
    return u, hist, bool(np.abs(F).max() < tol)
def load_mesh():
    S = json.load(open('surface_mesh.json')); M = json.load(open('sphere_mesh.json'))
    T = np.array(S['T']); SV = np.array(S['V'])
    Y = np.array([complex(*y) for y in M['y']]); r2 = np.abs(Y)**2
    Us = np.stack([2*Y.real, 2*Y.imag, r2 - 1], -1)/(r2 + 1)[:, None]
    P = Us[SV[:, 0]]
    E = edges(T)
    l0 = np.linalg.norm(P[E[:, 0]] - P[E[:, 1]], axis=1)
    return T, SV, Y, P, E, l0
if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'selftest':
        selftest(); sys.exit()
    T, SV, Y, P, E, l0 = load_mesh(); nV = len(SV)
    GEN = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'branch.json')))['genus']
    NSH = len(json.load(open('edge_perms.json'))['R'][next(i for i, r in enumerate(json.load(open('edge_perms.json'))['R']) if r)])
    AREA = 4*np.pi*(GEN - 1)
    print('genus', GEN, '| sheets', NSH, '| hyperbolic target area 4 pi (g - 1) =', f'{AREA:.12f}')
    print('mesh: V', nV, 'E', len(E), 'F', len(T), 'chi', nV - len(E) + len(T), '| edge length min %.3e max %.3e' % (l0.min(), l0.max()))
    # Euclidean angle sums of the input (cone angles: 2 pi * ramification at branch vertices, ~2 pi elsewhere)
    a, b, c = (np.linalg.norm(P[T[:, (j+1) % 3]] - P[T[:, (j+2) % 3]], axis=1) for j in range(3))
    def eang(a, b, c): return np.arccos(np.clip((b*b + c*c - a*a)/(2*b*c), -1, 1))
    S0 = np.bincount(T[:, 0], eang(a, b, c), nV) + np.bincount(T[:, 1], eang(b, c, a), nV) + np.bincount(T[:, 2], eang(c, a, b), nV)
    ram = np.round(S0/TAU).astype(int)
    print('input cone angles / 2pi: counts', dict(zip(*[x.tolist() for x in np.unique(ram, return_counts=True)])),
          '| vertices with ramification >= 2:', int((ram >= 2).sum()), '| total ramification', int((ram - 1).clip(0).sum()))
    area0 = float(np.sum(0.5*np.linalg.norm(np.cross(P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]]), axis=1)))
    print(f'input Euclidean area {area0:.6f} ({NSH} x 4 pi = {NSH*4*np.pi:.6f}; hyperbolic target {AREA:.6f})')
    u0 = np.full(nV, -np.log(2.0))
    # start with the known cone behaviour: on the ramified sheet near a vertex of cone angle 2 pi e, the hyperbolic length
    # scale is ~ |y - y_b|^(-(e-1)/e); cut off at chordal radius R0 (continuous there)
    import scipy.sparse.csgraph as csg
    # E93 change (run 1 failed: 52 degenerate triangles at the start, u0 up to 13.3, because the corrections of the ~25 branch
    # points within 0.3 of each other in the Per_7 cluster were summed): radius per branch point R_b = min(0.3, 0.45 x chordal
    # distance to the nearest other candidate point), so the disks (in y, hence on the surface) do not overlap.
    Bj = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'branch.json')))
    cand = np.array([complex(*b) for b in Bj['bps']])
    def chordal(a, b): return np.abs(a - b)*2/np.sqrt((1 + np.abs(a)**2)*(1 + np.abs(b)**2))
    Ad = sp.coo_matrix((l0, (E[:, 0], E[:, 1])), shape=(nV, nV)); Ad = (Ad + Ad.T).tocsr()
    rv = np.where(ram >= 2)[0]
    yb = Y[SV[rv, 0]]
    Rb = np.array([min(0.3, 0.45*np.sort(chordal(cand, yy))[1]) for yy in yb])
    dist = csg.dijkstra(Ad, indices=rv, limit=float(Rb.max()))
    for n_, b in enumerate(rv):
        e = ram[b]; d = dist[n_]; m = d < Rb[n_]
        rmin = 0.5*l0[(E[:, 0] == b) | (E[:, 1] == b)].min()
        u0[m] += -((e - 1)/e)*np.log(np.maximum(d[m], rmin)/Rb[n_])
    print('cone-correction radii (chordal): min', f'{Rb.min():.4f}', 'max', f'{Rb.max():.4f}')
    print('initial guess: cone correction at', len(rv), 'ramified vertices, u0 range', f'{u0.min():.3f} .. {u0.max():.3f}')
    u, hist, conv = solve(T, l0, E, nV, u0)
    F, l, tri = residual(u, T, setup(T, l0, E), l0, E, nV)
    a, b, c, A, B, C, ok = tri
    area = float(np.sum(np.pi - A - B - C))
    print(f'U1: converged {conv}, max |angle sum - 2pi| {np.abs(F).max():.3e}')
    print(f'U2: hyperbolic area {area:.12f}  4 pi (g - 1) = {AREA:.12f}  difference {area - AREA:.3e}')
    print(f'u range {u.min():.4f} .. {u.max():.4f}; at ramified vertices: {np.round(u[ram >= 2], 3).tolist()}')
    np.save('uniform_u.npy', u)
    json.dump(dict(hist=hist, converged=conv, max_defect=float(np.abs(F).max()), area=area, ramified=np.where(ram >= 2)[0].tolist()),
              open('uniform_log.json', 'w'), indent=1)
