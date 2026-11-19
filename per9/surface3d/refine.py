#!/usr/bin/env python3
"""E93-U: iterative refinement of the uniformization in 160-bit arithmetic (gmpy2), run in a mesh directory after uniformize.py.
Why: genus 24 side pairings displace the base point by up to ~24; composing group elements multiplies their errors by ~e^d, so the
double-precision metric (angle sums to ~1e-13, development consistent to ~1e-11) gives elements wrong by ~1e-4 (E93 log).
Same equations as uniformize.py (BPS hyperbolic, the double l0 = chordal lengths taken as exact data): residual F(u) =
angle sums - 2 pi evaluated in 160 bits, correction from the double Jacobian (factorized once), u kept in 160 bits.
Writes uniform_u_hi.npy, uniform_u_lo.npy (u = hi + lo, ~32 digits) and refine_log.json."""
import json
import os
import sys
import time
import numpy as np
import gmpy2
from gmpy2 import mpfr
from scipy.sparse.linalg import splu, cg
import pyamg
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from uniformize import load_mesh, setup, residual, jacobian
gmpy2.get_context().precision = 160
TWO_PI = 2*gmpy2.const_pi()


def mp_angles_sums(u, T, E, l0, nV):
    """angle sums per vertex (list of mpfr) and per-face (a, b, c, A, B, C) mpfr lists, from u (list of mpfr)."""
    half = mpfr('0.5')
    L = [2*gmpy2.asinh(half*mpfr(float(l0[k]))*gmpy2.exp(half*(u[int(i)] + u[int(j)]))) for k, (i, j) in enumerate(E)]
    fe = setup(T)
    S = [mpfr(0)]*nV
    area = mpfr(0)          # E95: area accumulated here; the per-face list (8M x 6 mpfr) exhausted memory on the fine Per_9 mesh
    for f in range(len(T)):
        a, b, c = L[fe[f, 0]], L[fe[f, 1]], L[fe[f, 2]]
        s = (a + b + c)/2
        ss, sa, sb, sc = gmpy2.sinh(s), gmpy2.sinh(s - a), gmpy2.sinh(s - b), gmpy2.sinh(s - c)
        A = 2*gmpy2.atan(gmpy2.sqrt(sb*sc/(ss*sa)))
        B = 2*gmpy2.atan(gmpy2.sqrt(sc*sa/(ss*sb)))
        C = 2*gmpy2.atan(gmpy2.sqrt(sa*sb/(ss*sc)))
        i, j, k = T[f]
        S[i] += A; S[j] += B; S[k] += C
        area += gmpy2.const_pi() - A - B - C
    return S, area


if __name__ == '__main__':
    t0 = time.time()
    T, SV, Y, P, E, l0 = load_mesh(); nV = len(SV)
    u0 = np.load('uniform_u.npy')
    fe = setup(T)
    F, l, tri = residual(u0, T, fe, l0, E, nV)
    J = jacobian(T, fe, l, tri, E, nV)
    # E95: CG + pyamg SA instead of SuperLU (memory; canaried on E94: 1e-13); rtol 1e-14 per correction
    A_ = (-J).tocsr(); A_ = 0.5*(A_ + A_.T); ml = pyamg.smoothed_aggregation_solver(A_, symmetry='symmetric'); Mp = ml.aspreconditioner()
    print(f'double residual at the start: max {np.abs(F).max():.3e} | AMG set up in {time.time()-t0:.1f} s', flush=True)
    u = [mpfr(float(x)) for x in u0]
    hist = []; it0 = 0
    if os.path.exists('refine_ckpt.npz'):                   # E95: resume (restarts); u = hi + lo in 160 bits
        z = np.load('refine_ckpt.npz'); u = [mpfr(float(a)) + mpfr(float(b)) for a, b in zip(z['hi'], z['lo'])]
        it0 = int(z['it']); hist = list(z['hist']); print('resumed from refine_ckpt.npz at refinement', it0, flush=True)
    for it in range(it0, 4):
        S, _ = mp_angles_sums(u, T, E, l0, nV)
        Fm = [s - TWO_PI for s in S]
        Fd = np.array([float(x) for x in Fm])
        hist.append(float(np.abs(Fd).max()))
        print(f'refinement {it}: max |angle sum - 2 pi| (160-bit evaluation) {hist[-1]:.3e} | {time.time()-t0:.1f} s', flush=True)
        if hist[-1] < 1e-30:
            break
        # E95: CG checkpointed every 25 iterations (key = refinement index + residual), warm start after restarts
        key = float(it) + float(np.abs(Fd).max()); x0 = None; nit = [0]
        if os.path.exists('refine_cg_ckpt.npz'):
            z_ = np.load('refine_cg_ckpt.npz')
            if float(z_['key']) == key: x0 = z_['x']; print(f'   CG warm start ({int(z_["it"])} iterations)', flush=True)
        def cb(xk):
            nit[0] += 1
            if nit[0] % 25 == 0:
                np.savez('refine_cg_ckpt.tmp.npz', x=xk, key=key, it=nit[0]); os.replace('refine_cg_ckpt.tmp.npz', 'refine_cg_ckpt.npz')
        du, info = cg(A_, Fd, x0=x0, rtol=1e-14, maxiter=5000, M=Mp, callback=cb)
        print(f'   CG info {info}, relative residual {np.linalg.norm(A_ @ du - Fd)/np.linalg.norm(Fd):.2e}', flush=True)
        u = [ui + mpfr(float(d)) for ui, d in zip(u, du)]
        _hi = np.array([float(x) for x in u]); _lo = np.array([float(x - mpfr(float(x))) for x in u])
        np.savez('refine_ckpt.tmp.npz', hi=_hi, lo=_lo, it=it + 1, hist=np.array(hist)); os.replace('refine_ckpt.tmp.npz', 'refine_ckpt.npz')
    hi = np.array([float(x) for x in u]); lo = np.array([float(x - mpfr(float(x))) for x in u])
    np.save('uniform_u_hi.npy', hi); np.save('uniform_u_lo.npy', lo)
    S, area = mp_angles_sums(u, T, E, l0, nV)
    print(f'final: max |angle sum - 2 pi| {max(abs(float(s - TWO_PI)) for s in S):.3e} | area {float(area):.15f} '
          f'| max |u - u_double| {np.abs(hi - u0).max():.3e}', flush=True)
    json.dump(dict(hist=hist, final=max(abs(float(s - TWO_PI)) for s in S), area=float(area), du_max=float(np.abs(hi - u0).max())),
              open('refine_log.json', 'w'), indent=1)
    if not (max(abs(float(s - TWO_PI)) for s in S) < 1e-25): sys.exit(3)     # E94: fail loudly
