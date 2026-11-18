#!/usr/bin/env python3
"""E93-U overlays for per7_live.html (E92 overlay_export.py for genus 24): the 4 real ovals and the graph G_7 (E91 G7.json).
Real ovals (U6: closed geodesics): oval_lifts.npz (oval_check.py) holds the oval points folded into D, the axis end points of the
holonomy conjugated into the local frame (per point) and the run id (a run = consecutive points in one tile). Each run is a piece
of the geodesic through its axis end points (the full geodesic meets the convex D in exactly that piece). Check: max distance of
the run's points to its geodesic.
G_7 (E91, capture edges with multiplicity): for each distinct pair of centres, the shortest geodesic on the surface between them:
the nearest lift of q to p is fold_p(q) in the Dirichlet polygon D_p about p (200 bits, fuchs.dirichlet_at_mp; certified
by area as in U4/U5). A stand-in for the capture path, NOT the path through the shared capture component; multiplicity-2 pairs
are drawn once (E92 also drew a second-shortest geodesic; not done here)."""
import json
import os
import sys
import time
from multiprocessing import Pool
import numpy as np
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import fuchs as fu
from fuchs import Mq, mob_ap, hd0, hdist, to_klein
from meshloc import load_D
GEN = json.load(open(os.path.join(here, '..', 'branch.json')))['genus']; AREA = 4*np.pi*(GEN - 1)
def geod_dist(z, xi1, xi2):
    al, be = np.angle(xi1), np.angle(xi2); half = (be - al)/2
    if abs(np.cos(half)) < 1e-9: return np.arcsinh(2*np.abs((z*np.exp(-1j*al)).imag)/(1 - np.abs(z)**2))
    C = np.exp(1j*(al + be)/2)/np.cos(half); rho = abs(np.tan(half))
    return np.arcsinh(np.abs(np.abs(z - C)**2 - rho**2)/(rho*(1 - np.abs(z)**2)))
_W = {}
def _init():
    res, D = load_D(here)
    _W.update(F0=fu.load_face_mp(D), rD=float(hd0(fu.to_poinc(res['verts'])).max()))
def _lifts(args):
    i, p, qs = args
    Di = fu.dirichlet_at_mp(p, _W['F0'], AREA, _W['rD'], step=2.5)
    q = mob_ap(np.linalg.inv(Di['Mp']), np.array(qs))
    qf, _, _ = fu.folder(Di)(q)
    return i, Di['certified'], mob_ap(Di['Mp'], qf), hd0(qf)
if __name__ == '__main__':
    t0 = time.time()
    L = np.load('oval_lifts.npz'); ovals = []
    n_ov = len([k for k in L.files if k.startswith('run')])
    for i in range(n_ov):
        f, x1, x2, run = L[f'f{i}'], L[f'x1_{i}'], L[f'x2_{i}'], L[f'run{i}']
        pieces = []; worst = 0.0
        for r in np.unique(run):
            idx = np.where(run == r)[0]
            if len(idx) < 2: continue
            a, b = x1[idx[0]], x2[idx[0]]
            dmax = float(geod_dist(f[idx], a, b).max()); worst = max(worst, dmax)
            pieces.append(dict(xi=[[a.real, a.imag], [b.real, b.imag]], n=int(len(idx)), dmax=dmax))
        uniq = []
        for pc in pieces:
            a = complex(*pc['xi'][0]); b = complex(*pc['xi'][1])
            if any(min(abs(a - complex(*u['xi'][0])) + abs(b - complex(*u['xi'][1])), abs(a - complex(*u['xi'][1])) + abs(b - complex(*u['xi'][0]))) < 1e-6 for u in uniq): continue
            uniq.append(pc)
        print(f'oval {i}: {len(f)} folded points, {len(np.unique(run))} runs in D, {len(uniq)} geodesic pieces; max distance of oval points '
              f'to their piece geodesic {worst:.2e}', flush=True)
        ovals.append(dict(pieces=uniq, maxdev=worst))
    CF = open('conformal_data.js').read(); CF = json.loads(CF[CF.index('{'):CF.rindex('}') + 1])
    cent = [m for m in CF['markers'] if m['kind'] == 'centre']
    data = json.load(open(os.path.join(here, '..', 'explorer', 'data.json')))
    G7 = json.load(open(os.path.join(here, '..', '..', 'inputs', 'E91-Gn-graph', 'G7.json')))
    # G7 index i -> census centre: the centre whose label equals G7['labels'][i]
    lab2c = {int(c['label']): k for k, c in enumerate(data['centers'])}
    gi2c = [lab2c[int(l)] for l in G7['labels']]
    mult = {}
    for e in G7['edges']:
        a, b = gi2c[e['a']], gi2c[e['b']]
        if a == b: continue
        k = tuple(sorted((a, b))); mult[k] = mult.get(k, 0) + 1
    pos = np.array([complex(*m['p']) for m in cent])
    byp = {}
    for a, b in mult: byp.setdefault(a, []).append(b)
    jobs = [(a, pos[a], [pos[b] for b in bs]) for a, bs in byp.items()]
    edges = []; ncert = 0
    with Pool(4, initializer=_init) as pool:
        for a, cert, Q, dd in pool.imap_unordered(_lifts, jobs):
            for b, q, d_ in zip(byp[a], Q, dd):
                edges.append(dict(kind='capture', a=int(a), b=int(b), mult=mult[(a, b)], rank=0, P=[pos[a].real, pos[a].imag],
                                  Q=[q.real, q.imag], length=float(d_), certified=bool(cert)))
                ncert += cert
    print(f'G_7 capture: {sum(mult.values())} edges with multiplicity, {len(mult)} distinct pairs (loops skipped), lengths '
          f'{min(e["length"] for e in edges):.3f} .. {max(e["length"] for e in edges):.3f}, certified shortest {ncert}/{len(edges)}; '
          f'{time.time()-t0:.0f} s', flush=True)
    json.dump(dict(ovals=ovals, edges=edges), open('overlay_data.json', 'w'))
    print('wrote overlay_data.json')
