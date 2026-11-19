#!/usr/bin/env python3
"""E95 V3 prerequisite: x (= a_8) at every surface vertex of one mesh (the cut-system gluing never computed vertex roots; the
viewer's charts store (y, x)). Harness, not mathematics; canaried below.
1. roots: the 42 roots at every normal sphere vertex (flint certified, escalating precision; trackmp.Curve.roots fallback, counted),
   3 workers, checkpointed per block (vx_roots/).
2. labels: x[v, l] = root carrying slit label l at v. Gauge at v0 (normal vertex nearest y0, straight segment y0 -> v0 crossing no
   cut): X0 continued (512 bits), X0 index j -> label h_at(v0)[j] (cut-system sector map). BFS over normal-normal sphere edges with
   the glued label maps m (edge_perms.json, label at a -> label at c): x[w, m[l]] = the root at w nearest x[u, l], accepted when
   every one of the 42 nearest distances is < 1/10 of its runner-up; otherwise the 42 roots are tracked along the edge polyline
   (cutglue.edge_path, as glue.py) at 512 bits (counted).
3. checks (raw): (a) every NON-tree normal edge re-predicts x[w] from x[u] (same rule) -> disagreements; (b) 200 random vertices:
   continuation y0 -> v (meshloc.Locator.roots_X0) with h_at(v) -> labels compared with the BFS labels.
4. branch-vertex surface vertices (V entry [b, -1 - g]): x = the root at b nearest the mean of x over the labels mapped to group g
   from one normal neighbour.
Writes vertex_x.npz (x per surface vertex, y per surface vertex). usage: vertex_x.py mesh_dir"""
import json, os, sys, time
import numpy as np
from collections import deque
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here); sys.path.insert(0, os.path.join(here, '..'))
d = os.path.abspath(sys.argv[1]); os.chdir(d)
import options_probe_lib as L
from meshloc import CURVE as C, CS, X0, Y0, NODES, Locator
from cutglue import edge_path
from gmpy2 import mpc
M = json.load(open('sphere_mesh.json')); Y = np.array([complex(*v) for v in M['y']]); bp = np.array(M['bp']); nS = len(Y)
EP = json.load(open('edge_perms.json')); n = EP.get('d') or 42
normal = np.where(bp < 0)[0]
RDISC = 0.1                                  # E95: inside |y| < 0.1 (the y = 0 disc of the gluing) roots collapse; x left NaN there
outer = set(int(v) for v in normal if abs(Y[v]) >= RDISC)
os.makedirs('vx_roots', exist_ok=True)
BLK = 2000
def rootblock(b):
    fn = f'vx_roots/b{b:04d}.npy'
    if os.path.exists(fn): return 0
    vs = normal[b*BLK:(b + 1)*BLK]; out = np.zeros((len(vs), n), complex); fb = 0
    for k, v in enumerate(vs):
        try:
            out[k] = [complex(r.mid()) for r in L.froots(C, complex(Y[v]))]
        except ValueError:
            out[k] = [complex(r) for r in C.roots(mpc(complex(Y[v])))]; fb += 1
    np.save(fn + '.tmp.npy', out); os.replace(fn + '.tmp.npy', fn); return fb

def transport_task(args):
    """x at w by label from x at u (E95 run 6): nearest root after a first-order prediction; roots whose nearest / runner-up
    ratio is >= 0.1 (and their neighbours within 5x their runner-up distance) are tracked along the edge polyline at 512 bits;
    every root is then assigned to its nearest root at w and the assignment must be a bijection (else all 42 are tracked)."""
    ya, yc, xu, Rw, m = args
    det = 0
    dd = np.abs(Rw[None, :] - xu[:, None]); o = np.sort(dd, 1)
    bad = np.where(~(o[:, 0] < 0.1*o[:, 1]))[0]
    out = np.empty(n, complex)
    nn = np.argmin(dd, 1)
    if len(bad) == 0 and len(set(nn)) == n:
        for l in range(n): out[m[l]] = Rw[nn[l]]
        return out, False, 0
    def track(idx):
        pts = edge_path(ya, yc, NODES)
        xt = C.track_path([mpc(complex(xu[i])) for i in idx], pts)
        d_ = 0
        if xt is None:
            mref = CS.edge_perm(pts)
            for sgn in (1, -1, 2, -2):
                pth = [ya, (ya + yc)/2 + 0.3*sgn*1j*(yc - ya), yc]
                if CS.edge_perm(pth) != mref: continue
                xt = C.track_path([mpc(complex(xu[i])) for i in idx], pth)
                if xt is not None: d_ = 1; break
        if xt is None: raise RuntimeError(f'tracking failed on edge y {ya:.6f} -> {yc:.6f} for every detour')
        return [complex(v) for v in xt], d_
    idx = set(bad.tolist())
    for i in bad: idx |= set(np.where(np.abs(xu - xu[i]) < 5*o[i, 1])[0].tolist())
    idx = sorted(idx); xt, det = track(idx)
    tgt = xu.copy(); tgt[idx] = xt
    nn = np.argmin(np.abs(Rw[None, :] - tgt[:, None]), 1)
    if len(set(nn)) != n:
        xt, det = track(list(range(n))); tgt = np.array(xt); nn = np.argmin(np.abs(Rw[None, :] - tgt[:, None]), 1)
        assert len(set(nn)) == n, 'assignment not a bijection after full tracking'
    for l in range(n): out[m[l]] = Rw[nn[l]]
    return out, True, det

t0 = time.time()
nb = (len(normal) + BLK - 1)//BLK
if __name__ == '__main__':
    from multiprocessing import Pool
    with Pool(3) as pool:
        fbs = pool.map(rootblock, range(nb))
    R = np.zeros((nS, n), complex)
    for b in range(nb):
        R[normal[b*BLK:(b + 1)*BLK]] = np.load(f'vx_roots/b{b:04d}.npy')
    print(f'{d}: roots at {len(normal)} normal sphere vertices ({nb} blocks), trackmp fallbacks this run {sum(fbs)}, {time.time()-t0:.0f} s', flush=True)
    # ---- labels
    adj = {}
    for e in EP['edges']:
        if e.get('toGroup'): continue
        a, c, m = e['a'], e['c'], e['m']
        if a not in outer or c not in outer: continue
        inv = [0]*n
        for i, j in enumerate(m): inv[j] = i
        adj.setdefault(a, []).append((c, m)); adj.setdefault(c, []).append((a, inv))
    v0 = min(normal, key=lambda v: abs(Y[v] - Y0))
    assert not CS.crossings(Y0, Y[v0]), 'segment y0 -> v0 crosses a cut'
    xe = C.track_path(X0, [Y0, complex(Y[v0])]); hv = CS.h_at(complex(Y[v0]))
    X = np.full((nS, n), np.nan + 0j)
    for j in range(n): X[v0, hv[j]] = complex(xe[j])
    # snap the gauge onto the vertex roots
    for l in range(n): X[v0, l] = R[v0][np.argmin(np.abs(R[v0] - X[v0, l]))]
    assert len(set(np.round(X[v0], 12))) == n
    DET = [0]
    # E95 run 6: level-synchronous BFS, transports in 3 workers, checkpoint per level (vx_labels.npz)
    if os.path.exists('vx_labels.npz'):
        z = np.load('vx_labels.npz', allow_pickle=True); X = z['X']; seen = set(z['seen'].tolist()); front = list(z['front'])
        tree = set(map(tuple, z['tree'].tolist())); ntr = int(z['ntr']); print(f'   resumed: {len(seen)} labelled, front {len(front)}', flush=True)
    else:
        seen = {int(v0)}; front = [int(v0)]; tree = set(); ntr = 0
    from multiprocessing import Pool
    pool = Pool(3)
    while front:
        tasks, tw = [], set()
        for u in front:
            for w, m in adj.get(u, []):
                if w in seen or w in tw: continue
                tw.add(w); tasks.append((u, w, m))
        res = pool.map(transport_task, [(complex(Y[u]), complex(Y[w]), X[u], R[w], m) for u, w, m in tasks], chunksize=8)
        nxt = []
        for (u, w, m), (out, tr, det) in zip(tasks, res):
            X[w] = out; ntr += tr; DET[0] += det; seen.add(w); tree.add((min(u, w), max(u, w))); nxt.append(w)
        front = nxt
        np.savez('vx_labels.tmp.npz', X=X, seen=np.array(sorted(seen)), front=np.array(front, int),
                 tree=np.array(sorted(tree), int).reshape(-1, 2), ntr=ntr); os.replace('vx_labels.tmp.npz', 'vx_labels.npz')
        print(f'   BFS level: +{len(tasks)}, labelled {len(seen)}, tracked {ntr}, {time.time()-t0:.0f} s', flush=True)
    print(f'labels: BFS reached {len(seen)} of {len(outer)} normal vertices with |y| >= {RDISC} ({len(normal) - len(outer)} inside, x NaN); tree edges tracked (ratio test failed) {ntr}, of which detoured {DET[0]}; {time.time()-t0:.0f} s', flush=True)
    bad = ntest = ntr2 = 0
    cand = [e for e in EP['edges'] if not e.get('toGroup') and e['a'] in outer and e['c'] in outer and (min(e['a'], e['c']), max(e['a'], e['c'])) not in tree]
    sub = [cand[i] for i in np.random.default_rng(23).choice(len(cand), min(20000, len(cand)), replace=False)]
    res = pool.map(transport_task, [(complex(Y[e['a']]), complex(Y[e['c']]), X[e['a']], R[e['c']], e['m']) for e in sub], chunksize=8)
    for e, (pr, tr, _) in zip(sub, res):
        ntest += 1; ntr2 += tr
        if np.abs(pr - X[e['c']]).max() > 1e-9: bad += 1
    print(f'   (check (a) on a random subsample of {len(sub)} of {len(cand)} non-tree normal edges)', flush=True)
    print(f'check (a): non-tree normal edges {ntest}, disagreements {bad}, tracked {ntr2}; {time.time()-t0:.0f} s', flush=True)
    Lc = Locator(d); smp = np.random.default_rng(17).choice(np.array(sorted(outer)), 200, replace=False); bad_b = fail_b = 0
    for v in smp:
        xs = Lc.roots_X0(complex(Y[v]))
        if xs is None: fail_b += 1; continue
        h = CS.h_at(complex(Y[v]))
        for j in range(n):
            if abs(X[v, h[j]] - xs[j]) > 1e-8*(1 + abs(xs[j])): bad_b += 1; break
    print(f'check (b): 200 random vertices by continuation from y0: label disagreements {bad_b}, continuation failures {fail_b}; {time.time()-t0:.0f} s', flush=True)
    S = json.load(open('surface_mesh.json'))['V']
    xs = np.empty(len(S), complex); ys = np.empty(len(S), complex)
    grp = {}
    for e in EP['edges']:
        if e.get('toGroup') and e['a'] in outer: grp.setdefault(e['c'], (e['a'], e['m']))
    for i, (v, l) in enumerate(S):
        ys[i] = Y[v]
        if l >= 0: xs[i] = X[v, l]; continue
        if v not in grp: xs[i] = np.nan; continue
        g = -1 - l; nbv, mg = grp[v]
        mean = np.mean([X[nbv, k] for k in range(n) if mg[k] == g])
        rb = np.array([complex(r) for r in C.roots(mpc(complex(Y[v])))])
        xs[i] = rb[np.argmin(np.abs(rb - mean))]
    print(f'surface vertices {len(S)}: x finite {int(np.isfinite(xs).sum())}, NaN (inside the y = 0 disc) {int((~np.isfinite(xs)).sum())}; {time.time()-t0:.0f} s', flush=True)
    np.savez('vertex_x.npz', x=xs, y=ys)
    sys.exit(0 if bad == 0 and bad_b == 0 and len(seen) == len(outer) else 1)
