"""E92 shared pieces: P7 coefficients in x as functions of y, batched roots, sheet labels by
continuation from the base point y0 along straight segments (K1 inner steps up to radius 12)."""
import os
import itertools
import json

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
B = json.load(open(HERE + 'branch.json'))
y0 = complex(*B['y0'])
x0 = np.array([complex(*v) for v in B['x0']])
bps = [complex(*b) for b in B['bps']]
xs, ys_ = sp.symbols('x y')
P = sp.sympify(open(os.path.dirname(os.path.abspath(__file__)) + '/../inputs/E90-real-circles/per6_xy_poly.txt').read().replace('^', '**'))
coef = [sp.lambdify(ys_, c, 'numpy') for c in sp.Poly(P, xs).all_coeffs()]
d = len(coef) - 1


assert sp.expand(sp.Poly(P, xs).all_coeffs()[0] - (ys_-1)**2) == 0   # leading coefficient is (y-1)^2


def coeffs(t):
    t = np.asarray(t, dtype=complex)
    out = np.stack([np.broadcast_to(np.asarray(c(t), dtype=complex), t.shape) for c in coef], -1)
    out[..., 0] = (t-1)**2          # avoid cancellation in y^2-2y+1 near the puncture locus y = 1
    return out


def roots_batch(t):
    c = coeffs(t)
    c = c / c[..., :1]
    N = c.shape[0]
    M = np.zeros((N, d, d), dtype=complex)
    M[:, 0, :] = -c[:, 1:]
    M[:, np.arange(1, d), np.arange(0, d-1)] = 1
    return np.linalg.eigvals(M)


PERMS = np.array(list(itertools.permutations(range(d))))


def match(prev, new):
    dist = np.abs(prev[:, :, None] - new[:, None, :])
    idx = dist.argmin(2)
    srt = np.sort(dist, 2)
    ok = (np.sort(idx, 1) == np.arange(d)).all(1) & (4*srt[:, :, 0] < srt[:, :, 1]).all(1)
    out = np.take_along_axis(new, idx, 1)
    bad = np.where(~ok)[0]
    if len(bad):
        cost = np.abs(prev[bad][:, None, :] - new[bad][:, PERMS]).sum(2)
        best = PERMS[cost.argmin(1)]
        out[bad] = np.take_along_axis(new[bad], best, 1)
    return out, len(bad)


def _continue_paths(paths):
    """paths: (N, M) complex, all starting at y0; returns labelled roots at the path ends."""
    cur = np.tile(x0, (paths.shape[0], 1))
    nb = 0
    for k in range(1, paths.shape[1]):
        cur, b = match(cur, roots_batch(paths[:, k]))
        nb += b
    return cur, nb


R0 = 0.12      # y = 0 is a node (trivial monodromy): detour around it on the segment's own side


def label(targets, K1=700, K2=120, detour=True):
    targets = np.asarray(targets, dtype=complex)
    if not detour:
        return _label_straight(targets, K1, K2)
    v = targets - y0
    tstar = np.clip(-np.real(np.conj(v)*y0)/np.maximum(np.abs(v)**2, 1e-300), 0, 1)
    close = (np.abs(y0 + tstar*v) < R0) & (tstar > 0) & (tstar < 1) & (np.abs(targets) > R0)
    out = np.empty((len(targets), d), dtype=complex)
    nb = 0
    if (~close).any():
        out[~close], nb = _label_straight(targets[~close], K1, K2)
    idx = np.where(close)[0]
    if len(idx):
        paths = []
        for i in idx:
            yt = targets[i]
            vv = yt - y0
            a, b, c = abs(vv)**2, 2*np.real(np.conj(vv)*y0), abs(y0)**2 - R0**2
            disc = np.sqrt(max(b*b - 4*a*c, 0.0))
            t1, t2 = (-b - disc)/(2*a), (-b + disc)/(2*a)
            p1, p2 = y0 + t1*vv, y0 + t2*vv
            f1, f2 = np.angle(p1), np.angle(p2)
            df = (f2 - f1 + np.pi) % (2*np.pi) - np.pi          # minor arc: the segment's own side
            seg1 = y0 + (p1 - y0)*np.linspace(0, 1, K1)
            arc = R0*np.exp(1j*(f1 + df*np.linspace(0, 1, 2000)))
            seg2 = p2 + (yt - p2)*np.linspace(0, 1, K1)
            paths.append(np.concatenate([seg1, arc[1:], seg2[1:]]))
        X2, b2 = _continue_paths(np.array(paths))
        far = np.abs(targets[idx]) > 12.0
        out[idx] = X2
        nb += b2
    return out, nb


def _label_straight(targets, K1=700, K2=120):
    """continue x0 from y0 to each target y along the segment; returns (N, d) labelled roots."""
    targets = np.asarray(targets, dtype=complex)
    L = np.abs(targets - y0)
    inner = np.minimum(L, 12.0)
    cur = np.tile(x0, (len(targets), 1))
    nb = 0
    for k in range(1, K1+1):
        rho = inner*k/K1
        pts = y0 + (targets-y0)/np.where(L > 0, L, 1)*rho
        cur, b = match(cur, roots_batch(pts))
        nb += b
    far = L > 12.0
    if far.any():
        for k in range(1, K2+1):
            rho = 12.0*(L[far]/12.0)**(k/K2)
            pts = y0 + (targets[far]-y0)/L[far]*rho
            cur[far], b = match(cur[far], roots_batch(pts))
            nb += b
    return cur, nb


