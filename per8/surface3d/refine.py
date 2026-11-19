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
from scipy.sparse.linalg import splu
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
    faces = []
    for f in range(len(T)):
        a, b, c = L[fe[f, 0]], L[fe[f, 1]], L[fe[f, 2]]
        s = (a + b + c)/2
        ss, sa, sb, sc = gmpy2.sinh(s), gmpy2.sinh(s - a), gmpy2.sinh(s - b), gmpy2.sinh(s - c)
        A = 2*gmpy2.atan(gmpy2.sqrt(sb*sc/(ss*sa)))
        B = 2*gmpy2.atan(gmpy2.sqrt(sc*sa/(ss*sb)))
        C = 2*gmpy2.atan(gmpy2.sqrt(sa*sb/(ss*sc)))
        i, j, k = T[f]
        S[i] += A; S[j] += B; S[k] += C
        faces.append((a, b, c, A, B, C))
    return S, faces


if __name__ == '__main__':
    t0 = time.time()
    T, SV, Y, P, E, l0 = load_mesh(); nV = len(SV)
    u0 = np.load('uniform_u.npy')
    fe = setup(T)
    F, l, tri = residual(u0, T, fe, l0, E, nV)
    J = jacobian(T, fe, l, tri, E, nV)
    lu = splu((-J).tocsc(), permc_spec='COLAMD')
    print(f'double residual at the start: max {np.abs(F).max():.3e} | factorized in {time.time()-t0:.1f} s', flush=True)
    u = [mpfr(float(x)) for x in u0]
    hist = []
    for it in range(4):
        S, _ = mp_angles_sums(u, T, E, l0, nV)
        Fm = [s - TWO_PI for s in S]
        Fd = np.array([float(x) for x in Fm])
        hist.append(float(np.abs(Fd).max()))
        print(f'refinement {it}: max |angle sum - 2 pi| (160-bit evaluation) {hist[-1]:.3e} | {time.time()-t0:.1f} s', flush=True)
        if hist[-1] < 1e-30:
            break
        du = lu.solve(Fd)
        u = [ui + mpfr(float(d)) for ui, d in zip(u, du)]
    hi = np.array([float(x) for x in u]); lo = np.array([float(x - mpfr(float(x))) for x in u])
    np.save('uniform_u_hi.npy', hi); np.save('uniform_u_lo.npy', lo)
    S, faces = mp_angles_sums(u, T, E, l0, nV)
    area = sum(gmpy2.const_pi() - A - B - C for a, b, c, A, B, C in faces)
    print(f'final: max |angle sum - 2 pi| {max(abs(float(s - TWO_PI)) for s in S):.3e} | area {float(area):.15f} '
          f'| max |u - u_double| {np.abs(hi - u0).max():.3e}', flush=True)
    json.dump(dict(hist=hist, final=max(abs(float(s - TWO_PI)) for s in S), area=float(area), du_max=float(np.abs(hi - u0).max())),
              open('refine_log.json', 'w'), indent=1)
    if not (max(abs(float(s - TWO_PI)) for s in S) < 1e-25): sys.exit(3)     # E94: fail loudly
