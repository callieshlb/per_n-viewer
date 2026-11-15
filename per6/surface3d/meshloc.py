"""E92-U: locate a point (y, x) of Per_6(0) in the glued mesh and in the developed disk (shared by the U-step scripts).
Sphere face by gnomonic barycentrics (convex-hull triangulation of the y-sphere), sheet by the closest x over the 5 lifts
(surface face 5t + i lies over sphere face t), position in the developed disk by Klein-model barycentric interpolation."""
import json
import numpy as np
from scipy.spatial import cKDTree
from develop import Mq, mob_ap, to_klein, to_poinc
def chord(a, b): return np.abs(a - b)/np.sqrt((1 + np.abs(a)**2)*(1 + np.abs(b)**2))
class Locator:
    def __init__(self, d='.'):
        S = json.load(open(f'{d}/surface_mesh.json')); M = json.load(open(f'{d}/sphere_mesh.json')); EP = json.load(open(f'{d}/edge_perms.json'))
        self.T = np.array(S['T']); self.SV = np.array(S['V'])
        Y = np.array([complex(*y) for y in M['y']]); r2 = np.abs(Y)**2
        self.U = np.stack([2*Y.real, 2*Y.imag, r2 - 1], -1)/(r2 + 1)[:, None]; self.TS = np.array(M['T'])
        yt = Y[self.SV[self.T, 0]]; sel = np.abs(yt).max(1) < 1.5
        sg = np.sign(((yt[sel, 1] - yt[sel, 0]).conjugate()*(yt[sel, 2] - yt[sel, 0])).imag)
        if (sg < 0).sum() > (sg > 0).sum(): self.T = self.T[:, [0, 2, 1]]
        R = EP['R']; xv = np.zeros(len(self.SV), complex); br = []
        for k, (v, i) in enumerate(self.SV):
            if i >= 0: xv[k] = complex(*R[v][i])
            else: br.append(k)
        for k in br:
            f = np.where((self.T == k).any(1))[0][:6]; nb = [c for ff in f for c in self.T[ff] if c != k and self.SV[c, 1] >= 0]
            xv[k] = np.mean([xv[c] for c in nb])
        self.xv = xv; self.yv = Y[self.SV[:, 0]]
        cen = self.U[self.TS].mean(1); cen /= np.linalg.norm(cen, axis=1)[:, None]; self.kd = cKDTree(cen)
    COEF = None
    def roots(self, y):
        if Locator.COEF is None:
            import os
            Locator.COEF = np.array(json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'explorer', 'data.json')))['coef'], float)
        c = [np.polyval(Locator.COEF[k][::-1], y) for k in range(len(Locator.COEF))]
        return np.roots(c[::-1])
    def track(self, y0, x0, y1, n=120):
        """follow the root x(y) from (y0, x0) to y1 (straight in 1/y when both are large, else in y), nearest-root steps."""
        big = abs(y0) > 2 and abs(y1) > 2
        a, b = (1/y0, 1/y1) if big else (y0, y1)
        x = x0
        for s in np.linspace(0, 1, n + 1)[1:]:
            w = a + s*(b - a); yy = 1/w if big else w
            r = self.roots(yy); x = r[np.argmin(chord(r, x))]
        return x
    def locate(self, y, x):
        """sheet chosen by continuation: the lift whose corner root, followed to y inside the triangle, arrives at x."""
        q = abs(y)**2; p = np.array([2*y.real, 2*y.imag, q - 1])/(q + 1)
        for t in self.kd.query(p, 12)[1]:
            w = np.linalg.solve(self.U[self.TS[t]].T, p)
            if (w >= -1e-12).all(): break
        else: return None
        w = w/w.sum(); cand = []
        for i in range(5):
            f = 5*t + i; c = self.T[f]
            wf = np.array([w[list(self.TS[t]).index(self.SV[cc, 0])] for cc in c])
            cc = next(k for k in c if self.SV[k, 1] >= 0)
            xe = self.track(self.yv[cc], self.xv[cc], y)
            cand.append((float(chord(xe, x)), f, wf))
        cand.sort(key=lambda r: r[0])
        self.last_ratio = cand[1][0]/max(cand[0][0], 1e-300); self.last_cands = cand
        return cand[0]
    def locate_sheet(self, y, x, sheet):
        """as locate, but when two lifts arrive at x equally well (nodes, x ~ y ~ infinity), prefer the lift whose corners carry
        the atlas sheet label `sheet` (root index = atlas label in edge_perms R)."""
        best = self.locate(y, x)
        if self.last_ratio >= 10: return best, 'continuation'
        tied = [c for c in self.last_cands if c[0] <= 10*max(best[0], 1e-300) + 1e-12]
        score = lambda c: sum(int(self.SV[k, 1] == sheet) for k in self.T[c[1]])
        pick = max(tied, key=lambda c: (score(c), -c[0]))
        return pick, f'atlas sheet {sheet} (tied lifts {len(tied)}, corners labelled {score(pick)}/3)'
    def locate_barycentric(self, y, x):
        q = abs(y)**2; p = np.array([2*y.real, 2*y.imag, q - 1])/(q + 1)
        for t in self.kd.query(p, 12)[1]:
            w = np.linalg.solve(self.U[self.TS[t]].T, p)
            if (w >= -1e-12).all(): break
        else: return None
        w = w/w.sum(); best = None
        for i in range(5):
            f = 5*t + i; c = self.T[f]
            wf = np.array([w[list(self.TS[t]).index(self.SV[cc, 0])] for cc in c])
            dx = float(np.sum(wf*chord(self.xv[c], x)))
            if best is None or dx < best[0]: best = (dx, f, wf)
        return best
def place(z, Tn, f, wf):
    a, b, c = z[Tn[f]]
    Mi = np.linalg.inv(Mq(a)); Kb, Kc = to_klein(mob_ap(Mi, b)), to_klein(mob_ap(Mi, c))
    return mob_ap(Mq(a), to_poinc(wf[1]*Kb + wf[2]*Kc))
