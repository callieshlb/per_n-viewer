#!/usr/bin/env python3
"""E93-U6 (E92 oval_check.py, lifting reorganised for genus 24): are the 4 real ovals of Per_7(0) closed geodesics of the
computed hyperbolic metric? usage: oval_check.py DIR (mesh directory with develop.npz from dirichlet_only.py).
Ovals: oval_points.json (trace_ovals.py, E90 arcs re-traced). Each point: exact location in the glued mesh (meshloc, 256-bit
continuation; ties = lifts within 10x the best x-mismatch), Klein-barycentric placement in the developed disk, folding into the
Dirichlet polygon D. Lifting: consecutive points f_k, f_{k+1} in D are joined by g in {I} u {face elements of D} (every tile
meeting a thin neighbourhood of D is a face neighbour: all vertex cycles have length 3) minimising d(f_k, g f_{k+1}); the lift
is W_k = H_k(f_k), H_{k+1} = H_k g, accumulated in 200 bits (E92 kept lifts in double; here the ovals are long and lifts would
leave double range). Holonomy h = H_N; geodesic length 2 arccosh(|tr h|/2); polyline length = sum of steps (local); distance
of W_k to the axis of h = distance of f_k to the axis of H_k^-1 h H_k (local frame). Writes DIR/oval_check.json and
DIR/oval_lifts.npz (folded points, axis end points per point, run ids) for overlay_export.py."""
import json
import os
import sys
import time
import numpy as np
import gmpy2
from gmpy2 import mpc
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import fuchs as fu
from meshloc import Locator, place, load_D
DIR = sys.argv[1] if len(sys.argv) > 1 else '.'


def axis_dist(z, xi1, xi2):
    al, be = np.angle(xi1), np.angle(xi2); half = (be - al)/2
    if abs(np.cos(half)) < 1e-9:
        return float(np.arcsinh(2*abs((z*np.exp(-1j*al)).imag)/(1 - abs(z)**2)))
    C = np.exp(1j*(al + be)/2)/np.cos(half); rho = abs(np.tan(half))
    return float(np.arcsinh(abs(abs(z - C)**2 - rho**2)/(rho*(1 - abs(z)**2))))


def fixed_points_mp(e):
    """fixed points on the unit circle of a hyperbolic SU(1,1) element (a, b): conj(b) z^2 + (conj(a) - a) z - b = 0."""
    a, b = e
    A, Bq, Cq = b.conjugate(), a.conjugate() - a, -b
    disc = gmpy2.sqrt(Bq*Bq - 4*A*Cq)
    return complex((-Bq + disc)/(2*A)), complex((-Bq - disc)/(2*A))


if __name__ == '__main__':
    t0 = time.time()
    ovals = json.load(open(os.path.join(here, 'oval_points.json')))
    print('ovals:', [(len(o['arcs']), len(o['pts'])) for o in ovals], '(arcs, points)', flush=True)
    res, D = load_D(DIR); z, Tn = D['z'], D['Tn']
    FG = res['face_g']; Fmp = fu.load_face_mp(D)
    fold = fu.folder(res)
    L = Locator(DIR)
    I2 = np.eye(2, dtype=complex)
    cands_g = [I2] + list(FG)
    def candidates(y, x):
        r = L.locate(y, x)
        if r is None: return []
        best = r[0]
        return [c for c in L.last_cands if c[0] <= max(10*best, 1e-7)]
    results, SAVE = [], {}
    def thin(p, d=1e-3):
        """E94 (registered before U6 ran): keep a traced point when it is >= 1e-3 from the last kept one in chordal y AND/OR x
        (max of the two); the traced arcs are very dense near punctures (oval 0: 395 953 points)."""
        ch = lambda a, b: abs(a - b)/np.sqrt((1 + abs(a)**2)*(1 + abs(b)**2))
        keep = [0]; last = p[0]
        for i in range(1, len(p)):
            if max(ch(p[i][0], last[0]), ch(p[i][1], last[1])) >= d: keep.append(i); last = p[i]
        return [p[i] for i in keep]
    for oi, ov in enumerate(ovals):
        pts = np.array(thin(ov['pts']))
        print(f'   oval {oi}: traced points {len(ov["pts"])}, kept after thinning {len(pts)}', flush=True)
        k0 = next(k for k in range(len(pts)) if abs(pts[k, 0]) < 2 and abs(pts[k, 1]) < 2 and len(candidates(complex(pts[k, 0]), complex(pts[k, 1]))) == 1)
        pts = np.concatenate([pts[k0:], pts[:k0]])
        Pl = list(pts) + [pts[0]]
        H = (mpc(1), mpc(0)); prev = None; steps = []; xmis = []; amb = 0; F_list = []; H_list = []; nfail = 0; runs = []; run = 0
        for k, (yr, xr) in enumerate(Pl):
            cs = candidates(complex(yr), complex(xr))
            if not cs:
                nfail += 1
                if nfail <= 40: print(f'   oval {oi} point {k}: location failed at y = {yr:.6g}, x = {xr:.6g}', flush=True)
                continue
            if len(cs) > 1: amb += 1
            best = None
            for dx, f, wf in cs:
                w = fold(np.array([place(z, Tn, f, wf)]))[0][0]
                if prev is None:
                    cand = (0.0, w, -1, dx)
                else:
                    imgs = np.array([fu.mob_ap(g, w) for g in cands_g])
                    ds = fu.hdist(prev, imgs); j = int(np.argmin(ds))
                    cand = (float(ds[j]), w, j - 1, dx)
                if best is None or cand[0] < best[0]: best = cand
            dist, w, j, dx = best
            if prev is not None:
                steps.append(dist)
                if j >= 0:
                    H = fu.mp_mul(H, Fmp[j]); run += 1
                if dist > 0.2:
                    print(f'   oval {oi} step {k}: y = {yr:.6g}, x = {xr:.6g}, candidates {len(cs)}, step {dist:.3e}', flush=True)
            xmis.append(dx); F_list.append(w); H_list.append(H); runs.append(run); prev = w
        h = H                                                       # W_N = H_N(f_N), f_N = f_0: holonomy
        tr = complex(2*h[0].real)
        trace = float(2*h[0].real)
        Lg = float(2*gmpy2.acosh(abs(h[0].real))) if abs(h[0].real) > 1 else 0.0
        Lp = float(np.sum(steps))
        dists, X1, X2 = [], [], []
        for f, Hk in zip(F_list, H_list):
            c = fu.mp_mul(fu.mp_inv(Hk), fu.mp_mul(h, Hk))
            xi1, xi2 = fixed_points_mp(c)
            dists.append(axis_dist(f, xi1, xi2)); X1.append(xi1); X2.append(xi2)
        dists = np.array(dists)
        r = dict(oval=oi, arcs=len(ov['arcs']), points=len(F_list), location_failures=nfail, ambiguous_points=amb,
                 max_lift_step=float(np.max(steps)), max_x_mismatch=float(np.max(xmis)), trace=trace,
                 hyperbolic=bool(abs(h[0].real) > 1), geodesic_length=Lg, polyline_length=Lp,
                 length_ratio_minus_1=Lp/Lg - 1 if Lg else None, face_crossings=run,
                 max_dist_to_axis=float(dists.max()), median_dist_to_axis=float(np.median(dists)))
        results.append(r)
        print('U6 oval', oi, json.dumps(r), f'| {time.time()-t0:.0f} s', flush=True)
        SAVE[f'f{oi}'] = np.array(F_list); SAVE[f'x1_{oi}'] = np.array(X1); SAVE[f'x2_{oi}'] = np.array(X2); SAVE[f'run{oi}'] = np.array(runs)
    json.dump(results, open(os.path.join(DIR, 'oval_check.json'), 'w'), indent=1)
    np.savez(os.path.join(DIR, 'oval_lifts.npz'), **SAVE)
