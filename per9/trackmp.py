"""E94 copy of E93 trackmp.py (default curve per9_xy_poly.txt here; E95 copy of E94). E93: multiprecision root continuation on the Per_7 curve P(x, y) = 0 (gmpy2 mpc, 256 bits ~ 77 digits).

Why: over y = 0 ten branches meet in two clusters (x -> 0 and x -> 1) mutually tangent to order 4 (root separation ~ |y|^4,
chordal 2e-9 at |y| = 1e-2); double precision cannot separate them for |y| < ~0.3 (E92's np.roots matching failed there).
All values are exact-integer-coefficient evaluations in 256-bit arithmetic.

Tracker (straight segment in y): second-order predictor x + x' h + x'' h^2/2 (implicit derivatives), Newton corrector.
A step is accepted only if every root's Newton converges, every correction is <= 0.1 x (distance from that root to the
nearest other root), and the roots stay distinct; otherwise the step is halved. Step cap 0.15 x distance to the nearest
special point (all discriminant / leading-coefficient roots), grown by 1.5 after success.
Edges with max |y| > 1 follow the great-circle arc of the y-sphere (16 pieces), as E92's edge_perms.js."""
import json
import os
import numpy as np
import gmpy2
from gmpy2 import mpc, mpfr
gmpy2.get_context().precision = 512      # E94: 512 bits (A3 run 1 at 256 bits could not continue to 19 loop points next to y = 0)
HERE = os.path.dirname(os.path.abspath(__file__))
TOL = mpfr('1e-36')
TOL_STAG = mpfr('1e-28')


def load_coef(path=None):
    """coef[k][j]: integer coefficient of x^k y^j of P (default per7_xy_poly.txt; E94: path or env HURWITZ_CURVE)."""
    import sympy as sp
    x, y = sp.symbols('x y')
    path = path or os.environ.get('HURWITZ_CURVE') or os.path.join(HERE, 'per9_xy_poly.txt')   # E95: the Per_9 curve
    P = sp.sympify(open(path).read().replace('^', '**'))
    Pp = sp.Poly(P, x, y)
    dx, dy = sp.Poly(P, x).degree(), sp.Poly(P, y).degree()
    coef = [[0]*(dy + 1) for _ in range(dx + 1)]
    for (kx, jy), c in Pp.terms():
        coef[kx][jy] = int(c)
    return coef


class Curve:
    def __init__(self, coef=None, specials=None):
        self.coef_int = coef if coef is not None else load_coef()
        self.n = len(self.coef_int) - 1                      # degree in x (number of sheets)
        self.C = [[mpc(c) for c in row] for row in self.coef_int]
        self.specials = np.array(specials if specials is not None else [], complex)

    def cy(self, y):
        """c_k(y), c_k'(y), c_k''(y) for k = 0..n."""
        c, c1, c2 = [], [], []
        for row in self.C:
            p = mpc(0); d1 = mpc(0); d2 = mpc(0)
            for a in reversed(row):
                d2 = d2*y + d1; d1 = d1*y + p; p = p*y + a
            c.append(p); c1.append(d1); c2.append(2*d2)
        return c, c1, c2

    @staticmethod
    def horner2(c, x):
        """sum c_k x^k and its first two x-derivatives."""
        p = mpc(0); d1 = mpc(0); d2 = mpc(0)
        for a in reversed(c):
            d2 = d2*x + d1; d1 = d1*x + p; p = p*x + a
        return p, d1, 2*d2

    @staticmethod
    def horner1(c, x):
        p = mpc(0); d1 = mpc(0)
        for a in reversed(c):
            d1 = d1*x + p; p = p*x + a
        return p, d1

    def derivs(self, cc, x):
        """x'(y), x''(y) on the branch through x (implicit differentiation)."""
        c, c1, c2 = cc
        P, Px, Pxx = self.horner2(c, x)
        Py, Pxy = self.horner1(c1, x)
        Pyy = self.horner1(c2, x)[0]
        xp = -Py/Px
        xpp = -(Pyy + 2*Pxy*xp + Pxx*xp*xp)/Px
        return xp, xpp

    def newton(self, c, x, maxit=40):
        """converged when |dx| <= 1e-36 (1 + |x|), or at round-off stagnation below 1e-28 (1 + |x|) (clustered roots)."""
        prev = None
        for it in range(maxit):
            P, Px = self.horner1(c, x)
            if Px == 0:
                return x, -1
            dx = P/Px
            x = x - dx
            adx = abs(dx); sc = 1 + abs(x)
            if adx <= TOL*sc or (prev is not None and adx <= TOL_STAG*sc and adx > 0.25*prev):
                return x, it + 1
            prev = adx
        return x, -1

    def roots(self, y):
        """all n roots at y, 256-bit. Initial values from mpmath.polyroots (60 digits), polished by Newton."""
        import mpmath as mp
        mp.mp.dps = 60
        yy = mp.mpc(complex(y)) if not isinstance(y, mpc) else mp.mpc(mp.mpf(str(y.real)), mp.mpf(str(y.imag)))
        cs = [mp.polyval(row[::-1], yy) for row in self.coef_int]
        r0 = mp.polyroots(cs[::-1], maxsteps=800, extraprec=800)
        yg = y if isinstance(y, mpc) else mpc(complex(y))
        c = self.cy(yg)[0]
        out = []
        for r in r0:
            x, it = self.newton(c, mpc(mpfr(str(r.real)), mpfr(str(r.imag))))
            out.append(x)
        return out

    def sep(self, xs):
        n = len(xs)
        s = [float('inf')]*n
        for i in range(n):
            for j in range(i + 1, n):
                d = float(abs(xs[i] - xs[j]))
                if d < s[i]: s[i] = d
                if d < s[j]: s[j] = d
        return s

    def dist_special(self, y):
        if not len(self.specials):
            return 1.0
        return float(np.abs(self.specials - complex(y)).min())

    def track_segment(self, xs, a, b, hfac=0.15, stats=None):
        """continue the roots xs (list of mpc at y = a) along the straight segment a -> b; returns roots at b or None."""
        a = a if isinstance(a, mpc) else mpc(complex(a)); b = b if isinstance(b, mpc) else mpc(complex(b))
        v = b - a; L = float(abs(v))
        if L == 0:
            return list(xs)
        t = mpfr(0); cur = list(xs); cc = self.cy(a)
        der = [self.derivs(cc, x) for x in cur]
        dt = None; nsteps = 0; nfail = 0
        while t < 1:
            yc = a + t*v
            ac = abs(complex(yc))
            hmax = min(hfac*self.dist_special(yc), max(0.5, 0.1*ac))
            if dt is None:
                dt = hmax/L
            dt = min(dt, hmax/L, float(1 - t))
            ok = False
            for tries in range(60):
                t2 = t + dt if dt < float(1 - t) else mpfr(1)
                h = (t2 - t)*v
                y2 = a + t2*v
                cc2 = self.cy(y2)
                nx, good = [], True
                for x, (xp, xpp) in zip(cur, der):
                    pred = x + xp*h + xpp*h*h/2
                    xn, it = self.newton(cc2[0], pred, 30)
                    if it < 0:
                        good = False; break
                    nx.append((xn, pred))
                if good:
                    new = [q[0] for q in nx]
                    s = self.sep(new)
                    for (xn, pred), si in zip(nx, s):
                        if not (si > 0) or float(abs(xn - pred)) > 0.1*si:
                            good = False; break
                if good:
                    ok = True; t = t2; cur = new; cc = cc2
                    der = [self.derivs(cc, x) for x in cur]
                    nsteps += 1; dt *= 1.5
                    break
                dt /= 2; nfail += 1
                if dt*L < 1e-14:
                    break
            if not ok:
                if stats is not None: stats['fail'] = stats.get('fail', 0) + 1
                return None
        if stats is not None:
            stats['steps'] = stats.get('steps', 0) + nsteps; stats['halvings'] = stats.get('halvings', 0) + nfail
        return cur

    def track_path(self, xs, pts, hfac=0.15, stats=None):
        cur = xs
        for p, q in zip(pts[:-1], pts[1:]):
            cur = self.track_segment(cur, p, q, hfac, stats)
            if cur is None:
                return None
        return cur

    def track_edge(self, xs, a, c, hfac=0.15, stats=None):
        """mesh edge: straight in y, or along the great-circle arc of the y-sphere when max |y| > 1 (16 pieces)."""
        a, c = complex(a), complex(c)
        if max(abs(a), abs(c)) <= 1:
            return self.track_segment(xs, a, c, hfac, stats)
        return self.track_path(xs, sphere_arc(a, c, 16), hfac, stats)


def to_u(y):
    r2 = abs(y)**2
    return np.array([2*y.real, 2*y.imag, r2 - 1])/(r2 + 1)


def to_y(u):
    return complex(u[0], u[1])/(1 - u[2])


def sphere_arc(a, c, K):
    ua, uc = to_u(a), to_u(c)
    om = np.arccos(np.clip(ua @ uc, -1, 1))
    pts = [a]
    for k in range(1, K):
        t = k/K
        u = (np.sin((1 - t)*om)*ua + np.sin(t*om)*uc)/np.sin(om)
        pts.append(to_y(u))
    pts.append(c)
    return pts


def match(tracked, target):
    """index map i -> j (tracked root i is target root j) by nearest, with the worst ratio nearest / runner-up."""
    m, worst = [], 0.0
    for x in tracked:
        d = sorted((float(abs(x - z)), j) for j, z in enumerate(target))
        m.append(d[0][1]); worst = max(worst, d[0][0]/max(d[1][0], 1e-300))
    return m, worst


def tofloat(xs):
    return [complex(x) for x in xs]


# ---------------- node detours (generalizes E92's single node at y = 0) ----------------
def arc_pts(c, r, f1, f2, n=24):
    """points on |y - c| = r from angle f1 to f2 the shorter way (excluding f1)."""
    df = (f2 - f1 + np.pi) % (2*np.pi) - np.pi
    return [c + r*np.exp(1j*(f1 + df*i/n)) for i in range(1, n + 1)]


def waypoints(a, b, nodes, depth=0):
    """polyline from a to b homotopic to the straight segment in the complement of the branch points, keeping away from each
    node (c, R): a chord through the disk |y - c| < R is replaced by its minor arc; an endpoint inside the disk is joined along
    the circle through the farther endpoint radius (capped at R), as E92's explorer core. Disks must be disjoint and free of
    branch points (then every replacement is a homotopy across a region holding only the node, of trivial monodromy)."""
    a, b = complex(a), complex(b)
    if depth > 8:
        return [a, b]
    v = b - a; L2 = abs(v)**2
    if L2 == 0:
        return [a, b]
    for c, R in nodes:
        ra, rb = abs(a - c), abs(b - c)
        t = min(1.0, max(0.0, -((a - c).conjugate()*v).real/L2))
        cl = abs(a + t*v - c)
        if ra < R or rb < R:
            if cl >= 0.5*min(ra, rb):
                continue
            r = max(min(max(ra, rb), R), 1e-12)
            A = c + r*(a - c)/ra; Bp = c + r*(b - c)/rb
            path = waypoints(a, A, nodes, depth + 1) if ra > r else [a, A]
            path = path + arc_pts(c, r, np.angle(a - c), np.angle(b - c))
            if rb > r:
                path = path[:-1] + waypoints(Bp, b, nodes, depth + 1)
            elif rb < r:
                path = path + [b]
            return path
        if 0 < t < 1 and cl < R:
            A_, B_, C_ = L2, 2*((a - c).conjugate()*v).real, ra**2 - R**2
            disc = np.sqrt(max(B_*B_ - 4*A_*C_, 0.0))
            t1, t2 = (-B_ - disc)/(2*A_), (-B_ + disc)/(2*A_)
            p1, p2 = a + t1*v, a + t2*v
            return (waypoints(a, p1, nodes, depth + 1) + arc_pts(c, R, np.angle(p1 - c), np.angle(p2 - c))[:-1]
                    + waypoints(p2, b, nodes, depth + 1))
    return [a, b]


def dedup(path):
    out = [path[0]]
    for p in path[1:]:
        if abs(p - out[-1]) > 0:
            out.append(p)
    return out


def track_line(C, xs, a, b, nodes, hfac=0.15, stats=None):
    return C.track_path(xs, dedup(waypoints(a, b, nodes)), hfac, stats)


def track_edge_nodes(C, xs, a, c, nodes, hfac=0.15, stats=None):
    """mesh edge: great-circle arc when max |y| > 1 (16 pieces), else straight; node detours on every piece."""
    a, c = complex(a), complex(c)
    pts = sphere_arc(a, c, 16) if max(abs(a), abs(c)) > 1 else [a, c]
    full = [pts[0]]
    for p, q in zip(pts[:-1], pts[1:]):
        full += waypoints(p, q, nodes)[1:]
    return C.track_path(xs, dedup(full), hfac, stats)
