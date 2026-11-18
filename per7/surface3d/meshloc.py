"""E93-U (E92 meshloc.py generalized to d sheets; continuation in 256 bits): locate a point (y, x) of Per_7(0) in the glued mesh
and in the developed disk. Sphere face by gnomonic barycentrics (convex-hull triangulation of the y-sphere); sheet by exact
continuation: the roots at a normal corner of the triangle (R, polished in 256 bits) are continued to y (node detours) and the
lift whose corner root arrives at x is chosen (E92: Locator.locate). Position in the developed disk by Klein-model barycentric
interpolation (E92: place)."""
import json
import os
import sys
import numpy as np
from scipy.spatial import cKDTree
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here); sys.path.insert(0, os.path.join(here, '..'))
from fuchs import Mq, mob_ap, to_klein, to_poinc
import trackmp as tm
from gmpy2 import mpc
B = json.load(open(os.path.join(here, '..', 'branch.json')))
NODES = [(complex(*b), 0.025) for b, k in zip(B['bps'], B['kinds']) if k == 'node']
CURVE = tm.Curve(specials=[complex(*b) for b in B['bps']])


def chord(a, b): return np.abs(a - b)/np.sqrt((1 + np.abs(a)**2)*(1 + np.abs(b)**2))


class Locator:
    def __init__(self, d='.'):
        S = json.load(open(f'{d}/surface_mesh.json')); M = json.load(open(f'{d}/sphere_mesh.json')); EP = json.load(open(f'{d}/edge_perms.json'))
        self.T = np.array(S['T']); self.SV = np.array(S['V'])
        self.R = EP['R']; self.n = len(next(r for r in self.R if r is not None))
        Y = np.array([complex(*y) for y in M['y']]); r2 = np.abs(Y)**2
        self.Ys = Y; self.isB = np.array(M['bp']) >= 0
        self.U = np.stack([2*Y.real, 2*Y.imag, r2 - 1], -1)/(r2 + 1)[:, None]; self.TS = np.array(M['T'])
        yt = Y[self.SV[self.T, 0]]; sel = np.abs(yt).max(1) < 1.5
        sg = np.sign(((yt[sel, 1] - yt[sel, 0]).conjugate()*(yt[sel, 2] - yt[sel, 0])).imag)
        if (sg < 0).sum() > (sg > 0).sum(): self.T = self.T[:, [0, 2, 1]]
        xv = np.zeros(len(self.SV), complex); br = []
        for k, (v, i) in enumerate(self.SV):
            if i >= 0: xv[k] = complex(*self.R[v][i])
            else: br.append(k)
        for k in br:
            f = np.where((self.T == k).any(1))[0][:6]; nb = [c for ff in f for c in self.T[ff] if c != k and self.SV[c, 1] >= 0]
            xv[k] = np.mean([xv[c] for c in nb])
        self.xv = xv; self.yv = Y[self.SV[:, 0]]
        cen = self.U[self.TS].mean(1); cen /= np.linalg.norm(cen, axis=1)[:, None]; self.kd = cKDTree(cen)

    def face_of(self, y):
        q = abs(y)**2; p = np.array([2*y.real, 2*y.imag, q - 1])/(q + 1)
        for t in self.kd.query(p, 16)[1]:
            w = np.linalg.solve(self.U[self.TS[t]].T, p)
            if (w >= -1e-12).all(): return t, w/w.sum()
        return None, None

    def track_from(self, v, y):
        """the n roots at normal sphere vertex v (R order) continued to y; list of complex (None on failure)."""
        yv = complex(self.Ys[v]); c = CURVE.cy(mpc(yv))[0]
        xs = [CURVE.newton(c, mpc(complex(*x)))[0] for x in self.R[v]]
        r = tm.track_line(CURVE, xs, yv, complex(y), NODES)
        return None if r is None else [complex(z) for z in r]

    def locate(self, y, x):
        """(chordal x-mismatch, surface face, barycentric weights in its corner order) of the best lift; ties in last_cands."""
        y = complex(y); x = complex(x)
        t, w = self.face_of(y)
        if t is None: return None
        vs = [v for v in self.TS[t] if not self.isB[v]]
        v = min(vs, key=lambda vv: abs(self.Ys[vv] - y))
        roots = self.track_from(v, y)
        if roots is None: return None
        cand = []
        for i in range(self.n):
            f = self.n*t + i; c = self.T[f]
            wf = np.array([w[list(self.TS[t]).index(self.SV[cc, 0])] for cc in c])
            k = next(cc for cc in c if self.SV[cc, 0] == v)
            cand.append((float(chord(roots[self.SV[k, 1]], x)), f, wf))
        cand.sort(key=lambda r: r[0])
        self.last_ratio = cand[1][0]/max(cand[0][0], 1e-300); self.last_cands = cand
        return cand[0]


def place(z, Tn, f, wf):
    a, b, c = z[Tn[f]]
    Mi = np.linalg.inv(Mq(a)); Kb, Kc = to_klein(mob_ap(Mi, b)), to_klein(mob_ap(Mi, c))
    return mob_ap(Mq(a), to_poinc(wf[1]*Kb + wf[2]*Kc))


def load_D(d='.'):
    """the certified Dirichlet polygon of develop.npz as a fuchs result dict."""
    D = np.load(f'{d}/develop.npz')
    return dict(verts=D['verts'], face_g=D['face_g'], bounded=True), D
