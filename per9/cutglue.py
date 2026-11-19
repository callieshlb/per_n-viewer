#!/usr/bin/env python3
"""E95: cut-system gluing (sheet maps of every mesh edge from the P1 lasso permutations; no root tracking).
Cuts = the P1 lasso tails: polylines y0 -> st (trackmp.waypoints with the P1 node detours) -> b. The star of tails plus a ray
y0 -> inf cuts the sphere into a simply connected region holding no candidate; sheets are labelled there by continuation from
y0 (label i = root X0[i]). The loop at infinity is the identity (P1), so the ray is invisible. A lasso (tail, counterclockwise
circle, back) maps label i to sigma_b[i]; the circle leaves the tail to its right and returns from its left, so crossing the
tail of b from its LEFT to its RIGHT (left = i x direction y0 -> b) applies sigma_b, right to left sigma_b^-1. An edge's sheet map
is the composite over its crossings in path order: label l -> p_k[...p_1[l]].
Functions are vectorized over tail segments (numpy)."""
import cmath
import math
import numpy as np
import trackmp as tm


def tails_from_branch(B, nodedet=None, m=96):
    """cuts = STRAIGHT segments y0 -> b. The P1 tails detour round nodes on minor arcs; a chord and its minor arc are homotopic in
    the complement of all candidates (the swept region never contains the disk's centre, and holds no other candidate), so each
    straight segment is homotopic to its P1 tail and the lasso permutations apply to it unchanged. (First canary run used the
    detoured polylines: tails passing 0 share one arc, so the cuts overlapped and crossing order was ambiguous: 3404 / 29 446
    disagreements.) The straight segments form a star from y0 that is non-crossing unless two candidates are collinear with y0."""
    y0 = complex(*B['y0']); out = []
    for b in B['bps']:
        out.append(np.array([y0, complex(*b)], complex))
    return y0, out


class CutSystem:
    def __init__(self, B, nodedet=None):
        self.y0, self.tails = tails_from_branch(B, nodedet)
        self.sigma = [list(p) for p in B['perms']]
        self.d = d = len(self.sigma[0])
        inv = lambda p: [p.index(i) for i in range(d)]
        # E95 fix (canary run 2: the triangle round y0 did not close; 3404 / 29 446 disagreements): y0 lies on every cut, so the
        # slit labels differ from X0 by sector. Sectors in counterclockwise order of phase(b - y0) mod 2 pi (as P1's product
        # check); S0 = the sector containing phase 0+, where slit label = X0 index. h = X0 index -> slit label in the current
        # sector; crossing tail k counterclockwise round y0 (right -> left) maps slit labels by tau_RL[k]; the lasso gives
        # sigma_k = h_k^-1 tau_LR[k] h_k (h_k = h in the sector just before tail k), so tau_LR[k] = h_k sigma_k h_k^-1.
        self.order = sorted(range(len(self.tails)), key=lambda k: cmath.phase(self.tails[k][-1] - self.y0) % (2*math.pi))
        self.perms = [None]*len(self.sigma); self.inv = [None]*len(self.sigma); self.hsec = []
        h = list(range(d))
        for k in self.order:
            self.hsec.append((cmath.phase(self.tails[k][-1] - self.y0) % (2*math.pi), list(h)))
            hi = inv(h); sg = self.sigma[k]
            tLR = [h[sg[hi[l]]] for l in range(d)]
            self.perms[k] = tLR; self.inv[k] = inv(tLR)
            h = [self.inv[k][h[i]] for i in range(d)]
        assert h == list(range(d)), 'tau product round y0 is not the identity'
        A, Bs, K = [], [], []
        for k, T in enumerate(self.tails):
            A.append(T[:-1]); Bs.append(T[1:]); K.append(np.full(len(T) - 1, k))
        self.A = np.concatenate(A); self.B = np.concatenate(Bs); self.K = np.concatenate(K)
        self.lo = np.minimum(self.A.real, self.B.real); self.hi = np.maximum(self.A.real, self.B.real)
        self.loi = np.minimum(self.A.imag, self.B.imag); self.hii = np.maximum(self.A.imag, self.B.imag)
        self.last = np.zeros(len(self.A), bool)
        o = 0
        for T in self.tails:
            o += len(T) - 1; self.last[o - 1] = True        # final segment st -> b of each tail

    def crossings(self, p, q, end_b=None):
        """crossings of the segment p -> q with the tails: list of (t along p->q, tail k, +1 left->right / -1 right->left).
        end_b: index of a tail whose endpoint b is q itself (edge into a branch vertex): its touching at q is ignored."""
        p, q = complex(p), complex(q)
        m = ((np.maximum(p.real, q.real) >= self.lo) & (np.minimum(p.real, q.real) <= self.hi) &
             (np.maximum(p.imag, q.imag) >= self.loi) & (np.minimum(p.imag, q.imag) <= self.hii))
        idx = np.where(m)[0]
        out = []
        if len(idx) == 0: return out
        a, b, k = self.A[idx], self.B[idx], self.K[idx]
        r, s = q - p, b - a
        den = (r.conjugate()*s).imag
        ok = den != 0
        w = a - p
        t = np.where(ok, (w.conjugate()*s).imag/np.where(ok, den, 1), -1)     # along p -> q
        u = np.where(ok, (w.conjugate()*r).imag/np.where(ok, den, 1), -1)     # along the tail segment
        hit = ok & (t > 0) & (t < 1) & (u >= 0) & (u < 1)
        for j in np.where(hit)[0]:
            if end_b is not None and k[j] == end_b and t[j] > 1 - 1e-9: continue
            # side: tail direction s; crossing direction r; left->right iff (s x r) < 0  (r points to the right of s)
            sgn = 1 if (s[j].conjugate()*r).imag < 0 else -1
            out.append((float(t[j]), int(k[j]), sgn))
        return out

    def edge_perm(self, path, end_b=None):
        """sheet map along a polyline (list of complex points); returns list l -> label."""
        lab = list(range(self.d)); n = len(path) - 1
        for i, (p, q) in enumerate(zip(path[:-1], path[1:])):
            cr = sorted(self.crossings(p, q, end_b if i == n - 1 else None))
            for t, k, sg in cr:
                P = self.perms[k] if sg > 0 else self.inv[k]
                lab = [P[l] for l in lab]
        return lab

    def h_at(self, z):
        """X0 index -> slit label in the sector of y0 containing direction z - y0 (for a point near y0)."""
        a = cmath.phase(complex(z) - self.y0) % (2*math.pi); h = list(range(self.d))
        for ph, hk in self.hsec:
            if ph < a: continue
            return hk
        return h

    def groups(self, k):
        g = [-1]*self.d; ng = 0; p = self.perms[k]
        for i in range(self.d):
            if g[i] >= 0: continue
            j = i
            while g[j] < 0: g[j] = ng; j = p[j]
            ng += 1
        return g, ng


def edge_path(ya, yc, nodes):
    """the polyline E94's track_edge_nodes follows: great-circle arc (16 pieces) when max |y| > 1, else straight; node detours."""
    ya, yc = complex(ya), complex(yc)
    pts = tm.sphere_arc(ya, yc, 16) if max(abs(ya), abs(yc)) > 1 else [ya, yc]
    full = [pts[0]]
    for p, q in zip(pts[:-1], pts[1:]):
        full += tm.waypoints(p, q, nodes)[1:]
    return tm.dedup(full)
