#!/usr/bin/env python3
"""E94 copy (per-node detours from ../nodes.json, per-vertex branch offsets, 512-bit ../trackmp.py). E93 mesh step 2 (replaces E92's edge_perms.js; same output format): roots per vertex, sheet maps along every edge, groups at
branch vertices. All continuation in 256-bit arithmetic (../trackmp.py), node detours round y = 0 (radius 0.025; mesh vertices
keep >= 0.03 away). usage: edge_perms.py DIR (DIR holds sphere_mesh.json; writes DIR/edge_perms.json).
- Normal vertex v: the 10 roots R[v] at y_v (mpmath.polyroots 60 digits, Newton-polished in 256 bits), in an arbitrary fixed
  order (NOT atlas labels: the gluing only needs consistent per-vertex orders).
- Edge a -> c (both normal): continue R[a] along the edge (straight, or the great-circle arc if max |y| > 1), match to R[c].
  Every edge is done twice (step factor 0.15 and 0.05); the two maps must agree.
- Branch vertex b: roots at the offset point o = b + 1e-4 e^{0.3i}, carried round |y - b| = 1e-4 (128 steps): the cycles of
  that permutation are the points over b (groups). Edge a -> b: continue R[a] to s = b + 1e-4 (y_a - b)/|y_a - b|, carry the
  offset roots from o to s along the short arc, match: root i at a -> group at b."""
import json
import os
import sys
import time
from multiprocessing import Pool
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
import trackmp as tm
from gmpy2 import mpc
DIR = sys.argv[1] if len(sys.argv) > 1 else HERE
B = json.load(open(os.path.join(HERE, '..', 'branch.json')))
BPS = [complex(*b) for b in B['bps']]
NODES = [(complex(*q['c']), q['detour']) for q in json.load(open(os.path.join(HERE, '..', 'nodes.json')))]   # E94: per-node detours
EPS = None                         # E94: per branch vertex, sphere_mesh.json 'eps' (0.3 x innermost ring radius)
C = None; M = None; R = None; OFF = None
def init(r=None, off=None):
    global C, M, R, OFF
    C = tm.Curve(specials=BPS)
    M = json.load(open(os.path.join(DIR, 'sphere_mesh.json')))
    R = r; OFF = off
def vroots(v):
    y = complex(*M['y'][v])
    return [complex(x) for x in C.roots(y)]
def polish(v):
    y = mpc(complex(*M['y'][v])); c = C.cy(y)[0]
    return [C.newton(c, mpc(complex(*x)))[0] for x in R[v]]
def branch_vertex(v):
    EPS = M['eps'][str(v)]
    b = complex(*M['y'][v]); o = b + EPS*np.exp(0.3j)
    xs = C.roots(o)
    circ = [b + EPS*np.exp(1j*(0.3 + 2*np.pi*i/128)) for i in range(129)]
    cur = C.track_path(xs, circ)
    if cur is None: return v, None
    p, w = tm.match(cur, xs)
    g = [-1]*C.n; ng = 0
    for i in range(C.n):
        if g[i] >= 0: continue
        j = i
        while g[j] < 0: g[j] = ng; j = p[j]
        ng += 1
    return v, dict(o=[o.real, o.imag], xs=[[complex(x).real, complex(x).imag] for x in xs], perm=p, g=g, ng=ng, ratio=w)
def edge(k):
    p, q = M['E'][k]
    isB = lambda v: M['bp'][v] >= 0
    if isB(p) and isB(q): return k, None, 'bb'
    a, c = (q, p) if isB(p) else (p, q)
    ya, yc = complex(*M['y'][a]), complex(*M['y'][c])
    xa = polish(a)
    if not isB(c):
        res = []
        for hf in (0.15, 0.05):
            st = {}
            tr = tm.track_edge_nodes(C, xa, ya, yc, NODES, hf, st)
            if tr is None: return k, None, 'track'
            m, w = tm.match(tr, [mpc(complex(*x)) for x in R[c]])
            res.append((m, w, st.get('steps', 0)))
        if sorted(res[0][0]) != list(range(C.n)): return k, None, 'perm'
        return k, dict(a=a, c=c, m=res[0][0], ratio=max(res[0][1], res[1][1]), steps=res[0][2], agree=res[0][0] == res[1][0]), 'ok'
    off = OFF[c]; b = yc; EPS = M['eps'][str(c)]
    s = b + EPS*(ya - b)/abs(ya - b)
    res = []
    for hf in (0.15, 0.05):
        tr = tm.track_edge_nodes(C, xa, ya, s, NODES, hf)
        if tr is None: return k, None, 'track'
        o = complex(*off['o']); a0, a1 = np.angle(o - b), np.angle(s - b)
        da = (a1 - a0 + np.pi) % (2*np.pi) - np.pi
        arc = [b + EPS*np.exp(1j*(a0 + da*i/32)) for i in range(33)]
        cur = C.track_path([mpc(complex(*x)) for x in off['xs']], arc, hf)
        if cur is None: return k, None, 'track-arc'
        m, w = tm.match(tr, cur)
        res.append(([off['g'][i] for i in m], w))
    return k, dict(a=a, c=c, m=res[0][0], toGroup=True, ratio=max(res[0][1], res[1][1]), agree=res[0][0] == res[1][0]), 'ok'
def _cache(name):
    d = os.path.join(DIR, 'edge_cache'); os.makedirs(d, exist_ok=True); return os.path.join(d, name)
if __name__ == '__main__':
    # E94: checkpointed (container restarts kill long runs): vertex roots, branch vertices and blocks of 4000 edges are cached in
    # DIR/edge_cache/ and reused on restart; results are identical to a single pass (each item is computed independently).
    t0 = time.time()
    init()
    nV = len(M['y']); isB = [b >= 0 for b in M['bp']]
    normal = [v for v in range(nV) if not isB[v]]
    if os.path.exists(_cache('roots.json')):
        Rall = json.load(open(_cache('roots.json'))); print('vertex roots: loaded from cache', flush=True)
    else:
        with Pool(4, initializer=init) as pool:
            Rl = pool.map(vroots, normal, chunksize=64)
        Rall = [None]*nV
        for v, r in zip(normal, Rl): Rall[v] = [[x.real, x.imag] for x in r]
        json.dump(Rall, open(_cache('roots.json') + '.tmp', 'w')); os.replace(_cache('roots.json') + '.tmp', _cache('roots.json'))
    sepmin = min(min(abs(complex(*x) - complex(*z)) for i, x in enumerate(r) for z in r[i+1:]) for r in Rall if r is not None)
    print(f'vertex roots: {len(normal)} normal vertices, {time.time()-t0:.1f} s; min root separation at a vertex {sepmin:.3e}', flush=True)
    if os.path.exists(_cache('branch.json')):
        bl = [(int(v), g) for v, g in json.load(open(_cache('branch.json')))]
    else:
        with Pool(4, initializer=init) as pool:
            bl = pool.map(branch_vertex, [v for v in range(nV) if isB[v]])
        json.dump(bl, open(_cache('branch.json'), 'w'))
    group = {v: g for v, g in bl}
    fails = [v for v, g in bl if g is None]
    print('branch vertices', len(bl), '| failed', fails, '| points over each:', ' '.join(str(g['ng']) for v, g in bl if g),
          '| worst loop match ratio', f"{max(g['ratio'] for v, g in bl if g):.2e}", flush=True)
    t1 = time.time(); nE = len(M['E']); BS = 4000; out = []
    for b0 in range(0, nE, BS):
        f = _cache(f'edges_{b0}.json')
        if os.path.exists(f):
            out += [tuple(r) for r in json.load(open(f))]; continue
        with Pool(4, initializer=init, initargs=(Rall, group)) as pool:
            blk = pool.map(edge, range(b0, min(nE, b0 + BS)), chunksize=16)
        json.dump(blk, open(f + '.tmp', 'w')); os.replace(f + '.tmp', f); out += blk
        print(f'   edges {min(nE, b0 + BS)} of {nE} | {time.time()-t1:.0f} s', flush=True)
    status = {}
    for k, e, s_ in out: status[s_] = status.get(s_, 0) + 1
    edges = [None]*nE
    for k, e, s_ in out: edges[k] = e
    good = [e for e in edges if e]
    print(f'edges {nE} in {time.time()-t1:.1f} s (this pass): status {status}; step factor 0.15 vs 0.05 disagreements '
          f'{sum(1 for e in good if not e["agree"])}; worst match ratio {max(e["ratio"] for e in good):.2e}; '
          f'mean steps per normal edge {np.mean([e["steps"] for e in good if "steps" in e]):.1f}', flush=True)
    json.dump(dict(R=Rall, group={str(v): dict(g=g['g'], ng=g['ng'], perm=g['perm']) for v, g in group.items() if g},
                   edges=[dict(a=e['a'], c=e['c'], m=e['m'], **({'toGroup': True} if e.get('toGroup') else {})) if e else None for e in edges]),
              open(os.path.join(DIR, 'edge_perms.json'), 'w'))
    print(f'wrote {DIR}/edge_perms.json; total {time.time()-t0:.1f} s')
