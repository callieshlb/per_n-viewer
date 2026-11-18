#!/usr/bin/env python3
"""E92-U overlays for per6_live.html: the real ovals and the graph G_6 in the conformal model.
Real ovals (U6: they are closed geodesics): the fine-mesh lifts (oval_lifts.npz) are folded into the Dirichlet polygon D;
each maximal run with one fold element is a piece, drawn as the geodesic through its end points (the full geodesic meets
the convex D in exactly that piece). Check: max distance of the piece's points to its geodesic.
G_6 (E91 G6_sources.json): for each pair of centres, the shortest geodesic on the surface between them (min over all
explored group elements g of d(p, g q), certified when d(0,p) + d(0,q) + best <= explored bound); a second-shortest one
for capture pairs of multiplicity 2. This is a stand-in for the capture path, not the path through the shared capture
component. Shared punctures: geodesics centre -> puncture -> centre."""
import json
import os
import sys
import numpy as np
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
from develop import Mq, mob_ap, hd0, hdist, to_klein
D = np.load('develop.npz'); FG = list(D['face_g']); E = D['els_all']; B = float(D['els_bound'])
U_ = np.array([np.exp(1j*np.angle(mob_ap(g, 0))) for g in FG]); Tt = np.array([np.tanh(hd0(mob_ap(g, 0))/2) for g in FG])
Gi = [np.linalg.inv(g) for g in FG]
def fold_track(w):
    """fold points into D; return folded points and an integer key of the accumulated fold element per point."""
    w = w.copy(); key = np.zeros((len(w), 0), int).tolist(); F = np.tile(np.eye(2, dtype=complex), (len(w), 1, 1))
    for it in range(300):
        k = to_klein(w); v = (k[:, None]*np.conj(U_[None, :])).real - Tt[None, :]
        a = np.argmax(v, 1); out = v[np.arange(len(w)), a] > 1e-13
        if not out.any(): break
        for f in np.unique(a[out]):
            m = out & (a == f); G = Gi[f]
            w[m] = (G[0, 0]*w[m] + G[0, 1])/(G[1, 0]*w[m] + G[1, 1]); F[m] = np.einsum('ij,njk->nik', G, F[m])
    o = F[:, 0, 1]/F[:, 1, 1]
    return w, np.round(np.stack([o.real, o.imag], 1)*1e6).astype(np.int64)
def ideal_ends(P, Q):
    Mi = np.linalg.inv(Mq(P)); q = mob_ap(Mi, Q); e = np.exp(1j*np.angle(q)); M = Mq(P)
    return mob_ap(M, -e), mob_ap(M, e)
def geod_dist(z, xi1, xi2):
    al, be = np.angle(xi1), np.angle(xi2); half = (be - al)/2
    if abs(np.cos(half)) < 1e-9: return np.arcsinh(2*np.abs((z*np.exp(-1j*al)).imag)/(1 - np.abs(z)**2))
    C = np.exp(1j*(al + be)/2)/np.cos(half); rho = abs(np.tan(half))
    return np.arcsinh(np.abs(np.abs(z - C)**2 - rho**2)/(rho*(1 - np.abs(z)**2)))
L = np.load('oval_lifts.npz'); ovals = []
for i in range(2):
    pts = L[f'lifts{i}']; w, key = fold_track(pts)
    brk = np.where(np.any(key[1:] != key[:-1], axis=1))[0] + 1
    runs = np.split(np.arange(len(w)), brk); pieces = []; worst = 0.0
    for r in runs:
        if len(r) < 2: continue
        P, Q = w[r[0]], w[r[-1]]
        if hdist(P, Q) < 1e-3: continue
        x1, x2 = ideal_ends(P, Q); dmax = float(geod_dist(w[r], x1, x2).max()); worst = max(worst, dmax)
        pieces.append(dict(xi=[[x1.real, x1.imag], [x2.real, x2.imag]], n=int(len(r)), dmax=dmax))
    # merge pieces that are the same geodesic (the oval may revisit a piece across the starting point)
    uniq = []
    for pc in pieces:
        a = complex(*pc['xi'][0]); b = complex(*pc['xi'][1])
        if any(min(abs(a - complex(*u['xi'][0])) + abs(b - complex(*u['xi'][1])), abs(a - complex(*u['xi'][1])) + abs(b - complex(*u['xi'][0]))) < 1e-3 for u in uniq): continue
        uniq.append(pc)
    print(f'oval {i}: {len(pts)} lifted points, {len(runs)} runs in D, {len(uniq)} geodesic pieces; max distance of oval points to their piece geodesic {worst:.2e}')
    ovals.append(dict(pieces=uniq, maxdev=worst))
# G_6 pairs
CF = open('conformal_data.js').read(); CF = json.loads(CF[CF.index('{'):CF.rindex('}') + 1])
cent = [m for m in CF['markers'] if m['kind'] == 'centre']; punc = [m for m in CF['markers'] if m['kind'] == 'puncture']
G6S = json.load(open('../../inputs/E91-Gn-graph/G6_sources.json')); G6 = json.load(open('../../inputs/E91-Gn-graph/G6.json'))
assert [int(m['label']) for m in cent] == G6S['labels'], 'centre order differs from E91 labels'
mult = {}
for e in G6['edges']: k = tuple(sorted((e['a'], e['b']))); mult[k] = mult.get(k, 0) + 1
a_, b_, c_, d_ = E[:, 0, 0], E[:, 0, 1], E[:, 1, 0], E[:, 1, 1]
def shortest(p, q, nbest=1):
    W = (a_*q + b_)/(c_*q + d_); ds = hdist(p, W); o = np.argsort(ds); out = []
    for j in o:
        if any(abs(W[j] - W[k_]) < 1e-6 for k_ in [x[2] for x in out]): continue
        out.append((float(ds[j]), W[j], j))
        if len(out) == nbest: break
    return [(dd, w, hd0(p) + hd0(q) + dd <= B + 1e-12) for dd, w, _ in out]
def pos(m): return complex(*m['p'])
edges, ncert, ntot = [], 0, 0
for kind in ('capture', 'mating'):
    for a, b in G6S[kind]:
        k = tuple(sorted((a, b))); m = mult.get(k, 1) if kind == 'capture' else 1
        for r, (dd, q, cert) in enumerate(shortest(pos(cent[a]), pos(cent[b]), m)):
            ntot += 1; ncert += cert
            edges.append(dict(kind=kind, a=int(a), b=int(b), rank=r, P=[pos(cent[a]).real, pos(cent[a]).imag], Q=[q.real, q.imag], length=dd, certified=bool(cert)))
plab = {m['label']: m for m in punc}
for a, b in G6S['puncture']:
    la, lb = G6S['labels'][a], G6S['labels'][b]
    nm = f'{min(la, lb)}/{max(la, lb)}'
    if nm not in plab: print('no puncture marker named', nm); continue
    pp = pos(plab[nm])
    for c in (a, b):
        dd, q, cert = shortest(pp, pos(cent[c]))[0]; ntot += 1; ncert += cert
        edges.append(dict(kind='puncture', a=int(c), b=-1, rank=0, P=[pp.real, pp.imag], Q=[q.real, q.imag], length=dd, certified=bool(cert)))
for kind in ('capture', 'mating', 'puncture'):
    ek = [e for e in edges if e['kind'] == kind]
    print(f'{kind}: {len(ek)} geodesic segments (lengths {min(e["length"] for e in ek):.3f} .. {max(e["length"] for e in ek):.3f}), '
          f'certified shortest {sum(e["certified"] for e in ek)}/{len(ek)}')
json.dump(dict(ovals=ovals, edges=edges), open('overlay_data.json', 'w'))
print('wrote overlay_data.json')
