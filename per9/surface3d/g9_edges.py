#!/usr/bin/env python3
"""E95: overlay edges of G_9 (E91 construction, runs/2026-11-09-E95-per9-probe/G9.json) in E93's overlay_data format: for each
distinct pair of centres joined in G_9, the shortest geodesic on the surface between them = from the centre a (in D) to the
nearest lift of the centre b, found by the tiling BFS of u5_bfs.py (exact when the frontier closes: every tile meeting the ball
B(a, R) is visited, R = the largest current target distance) with the minimizing lift recorded. Centres = the viewer's markers
(u5_centres_ins.npy, fine mesh). Writes g9_edges.json (list of edge dicts: kind, a, b, mult, rank, P, Q, length, certified).
Rows cached in g9_rows/. usage: g9_edges.py [k/N worker]"""
import json, os, pickle, sys, time
import numpy as np
import gmpy2
from gmpy2 import mpc
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here); os.chdir(here)
import fuchs as fu, u5_bfs as ub
gmpy2.get_context().precision = 200
G = json.load(open(os.path.join(here, '..', '..', 'inputs', 'E95-probe', 'G9.json')))
data = json.load(open(os.path.join(here, '..', 'explorer', 'data.json')))
lab2c = {int(c['label']): k for k, c in enumerate(data['centers'])}
gi2c = [lab2c[int(l)] for l in G['labels']]
mult = {}; loops = 0
for e in G['edges']:
    a, b = gi2c[e['a']], gi2c[e['b']]
    if a == b: loops += 1; continue
    k = tuple(sorted((a, b))); mult[k] = mult.get(k, 0) + 1
byp = {}
for a, b in mult: byp.setdefault(a, []).append(b)
pos = np.load('u5_centres_ins.npy')
os.makedirs('g9_rows', exist_ok=True)
T = None
def row(a):
    global T
    fn = f'g9_rows/r{a}.pkl'
    if os.path.exists(fn): return pickle.load(open(fn, 'rb'))
    if T is None: T = ub.Tiling('.')
    t0 = time.time(); bs = byp[a]; Q = pos[bs]; p = pos[a]; pm = mpc(complex(p)); Qh = ub.hyp_from_poinc(Q)
    I = (mpc(1), mpc(0)); r = np.arccosh(np.maximum(1.0, -ub.mink(Qh, ub.hyp_from_poinc(p)[None, :]))); best = [I]*len(bs)
    seen = {fu.mp_key(I)}; front = [I]; nt = 0
    while front:
        nxt = []
        for g in front:
            nt += 1; gi = fu.mp_inv(g)
            pp = ub.mp_hyp((gi[0]*pm + gi[1])/(gi[1].conjugate()*pm + gi[0].conjugate()))
            dq = np.arccosh(np.maximum(1.0, -ub.mink(Qh, pp[None, :])))
            for i in np.where(dq < r)[0]: r[i] = dq[i]; best[i] = g
            R = r.max()
            for j in np.where(T.seg_dist(pp) < R + ub.MARGIN)[0]:
                h = fu.mp_mul(g, T.F[T.E[j]]); k = fu.mp_key(h)
                if k not in seen: seen.add(k); nxt.append(h)
        front = nxt
    lifts = []
    for g, q in zip(best, Q):
        qm = mpc(complex(q)); w = (g[0]*qm + g[1])/(g[1].conjugate()*qm + g[0].conjugate()); lifts.append(complex(w))
    out = dict(a=a, bs=bs, d=[float(x) for x in r], lifts=lifts, tiles=nt, seconds=round(time.time() - t0, 1))
    pickle.dump(out, open(fn + f'.tmp{os.getpid()}', 'wb')); os.replace(fn + f'.tmp{os.getpid()}', fn); return out
if __name__ == '__main__':
    keys = sorted(byp)
    if len(sys.argv) > 1:
        k, N = map(int, sys.argv[1].split('/'))
        for a in keys[k::N]: row(a)
        sys.exit(0)
    edges = []; chk = 0.0
    for a in keys:
        o = row(a)
        for b, d_, q in zip(o['bs'], o['d'], o['lifts']):
            dd = float(fu.hdist(np.array([pos[a]]), np.array([q]))[0]); chk = max(chk, abs(dd - d_))
            edges.append(dict(kind='capture', a=int(a), b=int(b), mult=mult[(a, b)], rank=0, P=[float(pos[a].real), float(pos[a].imag)],
                              Q=[q.real, q.imag], length=float(d_), certified=True))
    L = [e['length'] for e in edges]
    print(f'G_9 capture edges: {sum(mult.values())} with multiplicity (loops skipped {loops}), {len(mult)} distinct pairs, rows {len(keys)}; lengths {min(L):.3f} .. {max(L):.3f}, '
          f'median {np.median(L):.3f}; max |length - distance(P, Q) recomputed in double| {chk:.2e}; every row BFS exhausted (certified shortest)')
    json.dump(edges, open('g9_edges.json', 'w'))
