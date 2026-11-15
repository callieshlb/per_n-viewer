#!/usr/bin/env python3
"""E90 core library: Per_n(0) in the chart F(z)=1+b/z+c/z^2 (orbit 0->inf->1->...).

Sections: A exact angle combinatorics, B centers (spider + Newton), C continuation
primitives (mpmath), D polynomial -> puncture loop, E real arc tracer.
No conclusions are printed by this module; drivers print raw values only.
"""
import cmath
import math
import mpmath as mp
import numpy as np


# ---------------------------------------------------------------- A. combinatorics
def exact_period_angles(n):
    D = 2**n - 1
    out = []
    for k in range(1, D):
        x = k
        per = 0
        while True:
            x = (2*x) % D
            per += 1
            if x == k:
                break
        if per == n:
            out.append(k)
    return out


def kneading(k, n):
    """nu_1..nu_{n-1} for theta=k/D, as a list indexed from 1 (nu[0] unused), nu[n]=None."""
    D = 2**n - 1
    nu = [None]*(n+1)
    x = k
    for j in range(1, n):
        # 2^{j-1} theta in (theta/2,(theta+1)/2) ?
        nu[j] = 1 if (k < 2*x < k + D) else 0
        x = (2*x) % D
    return nu


def internal_address(nu, n):
    S = [1]
    while S[-1] < n:
        s0 = S[-1]
        s = s0 + 1
        while s < n and nu[s] == nu[((s-1) % s0) + 1]:
            s += 1
        S.append(s)
    return S


def rho(r, nu, n):
    s = r + 1
    while s < n and nu[s] == nu[s-r]:
        s += 1
    return s


def denominators(S, nu, n):
    """Bruin-Schleicher: denominators q_k of the arrows S_k -> S_{k+1}."""
    qs = []
    for i in range(len(S)-1):
        Sk, Sn = S[i], S[i+1]
        r = ((Sn-1) % Sk) + 1
        orb = [r]
        while orb[-1] < Sk:
            orb.append(rho(orb[-1], nu, n))
        qs.append((Sn-r)//Sk + (1 if Sk in orb else 2))
    return qs


def rotation_T(p, q):
    """Characteristic angle pair (as fractions of 2^q-1) of the p/q-limb (p/q != 0)."""
    bits = [((j+1)*p)//q - (j*p)//q for j in range(q)]
    N = 0
    for b in bits:
        N = 2*N + b
    D = 2**q - 1
    orb = []
    x = N
    for j in range(q):
        orb.append(x)
        x = (2*x) % D
    orb = sorted(orb)
    best = None
    for i in range(q):
        a, b2 = orb[i], orb[(i+1) % q]
        gap = (b2 - a) % D
        if best is None or gap < best[0]:
            best = (gap, a, b2)
    return best[1], best[2], D


def limb_of(k, n):
    """(p,q) with k/(2^n-1) in the closed p/q-limb."""
    D = 2**n - 1
    for q in range(2, n+1):
        for p in range(1, q):
            if math.gcd(p, q) != 1:
                continue
            a, b, Dq = rotation_T(p, q)
            # compare k/D with a/Dq, b/Dq exactly
            if a*D <= k*Dq <= b*D:
                return p, q
    return None


def label_angle(k, n):
    nu = kneading(k, n)
    S = internal_address(nu, n)
    q = denominators(S, nu, n)
    return dict(addr=S, qs=q, tail12=all(x == 2 for x in q[1:]), limb=limb_of(k, n),
                nu=''.join(str(x) for x in nu[1:n]))


# ---------------------------------------------------------------- B. centers
def spider(k, n, iters=600, R=1e4, J=12, tol=1e-11):
    D = 2**n - 1
    th = [(k*2**j % D)/D for j in range(n)]
    rad = [R**(0.5**i) for i in range(J+1)]
    legs = []
    for j in range(n):
        u = cmath.exp(2j*math.pi*th[j])
        pts = [r*u for r in rad]
        foot = u if j != n-1 else 0j
        legs.append((pts, foot))
    c_old = legs[0][1]
    for it in range(iters):
        c = legs[0][1]
        new = []
        for j in range(n):
            src = (j+1) % n
            pts, foot = legs[src]
            u = cmath.exp(2j*math.pi*th[j])
            seq = [p - c for p in pts]
            # subdivide last stretch to the foot for safe branch continuation
            last = pts[-1]
            seq += [last + (foot-last)*(t/8.0) - c for t in range(1, 9)]
            z = cmath.sqrt(seq[0])
            target = math.sqrt(R)*u
            if abs(z-target) > abs(-z-target):
                z = -z
            lifted = [z]
            for w in seq[1:]:
                y = cmath.sqrt(w)
                if abs(y-z) > abs(-y-z):
                    y = -y
                lifted.append(y)
                z = y
            npts = [R*u] + lifted[:J]
            nfoot = lifted[-1]
            if j == n-1:
                nfoot = 0j
            new.append((npts, nfoot))
        legs = new
        c = legs[0][1]
        if it > 8 and abs(c-c_old) < tol:
            return c, it
        c_old = c
    return legs[0][1], iters


def newton_center(c, n, dps=40, steps=60):
    mp.mp.dps = dps
    c = mp.mpc(c)
    for _ in range(steps):
        z = mp.mpc(0)
        dz = mp.mpc(0)
        for i in range(n):
            dz = 2*z*dz + 1
            z = z*z + c
        d = z/dz
        c = c - d
        if abs(d) < mp.mpf(10)**(-dps+6):
            break
    return c


def orbit_period(c, n):
    z = mp.mpc(0)
    orb = []
    for i in range(n):
        orb.append(z)
        z = z*z + c
    if abs(z) > mp.mpf(10)**(-mp.mp.dps//2):
        return None
    m = min(abs(orb[i]) for i in range(1, n))
    return float(m)


def center_to_chart(c):
    kap = c/(c+1)
    return 2*kap*(kap-1), kap*kap*(1-kap), kap


# ---------------------------------------------------------------- C. primitives
def R_all(b, c, n, one=1):
    """Residual: a_{n-1}^2 + b a_{n-1} + c, with derivatives; generic number type."""
    z = one
    zb = 0*one
    zc = 0*one
    for _ in range(n-3):
        iz = one/z
        iz2 = iz*iz
        dF = -b*iz2 - 2*c*iz2*iz
        zb, zc = iz + dF*zb, iz2 + dF*zc
        z = one + b*iz + c*iz2
    t = 2*z + b
    return z*z + b*z + c, t*zb + z, t*zc + one


def orbit(b, c, n):
    """a_0=0, a_1=None (inf), a_2=1, a_3, ..., a_{n-1}."""
    a = [0*b, None, 1 + 0*b]
    z = a[2]
    for _ in range(n-3):
        z = 1 + b/z + c/(z*z)
        a.append(z)
    return a


def fixed_points(b, c, z=None):
    if z is None:
        r = mp.polyroots([1, -1, -b, -c], maxsteps=200, extraprec=200)
        return list(r)
    # deflate Z^3-Z^2-bZ-c by (Z-z)
    p = z - 1
    q = z*z - z - b
    s = mp.sqrt(p*p - 4*q)
    return [z, (-p+s)/2, (-p-s)/2]


def mult(z, b, c):
    return -b/(z*z) - 2*c/(z*z*z)


def solve3(A, y):
    # Cramer, A given by rows
    (a, b_, c_), (d, e, f), (g, h, i) = A
    ei_fh = e*i - f*h
    fg_di = f*g - d*i
    dh_eg = d*h - e*g
    det = a*ei_fh + b_*fg_di + c_*dh_eg
    y1, y2, y3 = y
    x1 = (y1*ei_fh + b_*(f*y3 - y2*i) + c_*(y2*h - e*y3))/det
    x2 = (a*(y2*i - f*y3) + y1*fg_di + c_*(d*y3 - y2*g))/det
    x3 = (a*(e*y3 - y2*h) + b_*(y2*g - d*y3) + y1*dh_eg)/det
    return x1, x2, x3


def newton_ray(b, c, z, lam, n, tol, maxit=12):
    for it in range(maxit):
        R, Rb, Rc = R_all(b, c, n, mp.mpc(1))
        E2 = z*z*z - z*z - b*z - c
        E3 = lam*z*z*z + b*z + 2*c
        A = ((Rb, Rc, 0), (-z, -1, 3*z*z - 2*z - b), (z, 2, 3*lam*z*z + b))
        db, dc, dz = solve3(A, (-R, -E2, -E3))
        b += db
        c += dc
        z += dz
        if abs(db) <= tol*abs(b) and abs(dc) <= tol*abs(c) and abs(dz) <= tol*abs(z):
            return b, c, z, it+1
    return None


def newton_g(b, c, z, g, n, tol, maxit=12):
    for it in range(maxit):
        R, Rb, Rc = R_all(b, c, n, mp.mpc(1))
        E2 = z*z*z - z*z - b*z - c
        E3 = z*z*z + g*(b*z + 2*c)
        A = ((Rb, Rc, 0), (-z, -1, 3*z*z - 2*z - b), (g*z, 2*g, 3*z*z + g*b))
        db, dc, dz = solve3(A, (-R, -E2, -E3))
        b += db
        c += dc
        z += dz
        if abs(db) <= tol*abs(b) and abs(dc) <= tol*abs(c) and abs(dz) <= tol*abs(z):
            return b, c, z, it+1
    return None


def guard(b0, c0, b1, c1, n, frac=0.5):
    """Branch-jump guard: every pairwise gap among the marked points (a_1=inf included, via
    1/a_i) may change by at most relative frac. Rigid motion of a cluster is not penalised."""
    a0 = orbit(b0, c0, n)
    a1 = orbit(b1, c1, n)
    idx = [0] + list(range(2, n))
    for x, i in enumerate(idx):
        if i != 0:
            if abs(a0[i]/a1[i]-1) > frac:
                return False
        for j in idx[x+1:]:
            d0 = a0[i]-a0[j]
            if d0 == 0:
                return False
            if abs((a1[i]-a1[j])/d0-1) > frac:
                return False
    return True


def small_g(b, c, z):
    """(g, z3): reciprocal of the largest multiplier and its fixed point, given one fixed point z."""
    fp = fixed_points(b, c, z)
    ms = [mult(w, b, c) for w in fp]
    i = max(range(3), key=lambda j: abs(ms[j]))
    return 1/ms[i], fp[i], ms


# ---------------------------------------------------------------- D. poly -> puncture
class Fail(Exception):
    pass


def ray_to_puncture(c0, n, p, q, eps, tol):
    """From center c0 in the p/q-limb follow lambda=r*exp(-2 pi i p/q) until |g|<=eps."""
    b, c, kap = center_to_chart(c0)
    z = kap
    u = mp.expjpi(-2*mp.mpf(p)/q)
    s = mp.mpf(0)
    ds = mp.mpf('0.15')
    steps = 0
    while True:
        s2 = s + ds
        lam = (1-mp.exp(-s2))*u
        res = newton_ray(b, c, z, lam, n, tol)
        ok = res is not None and guard(b, c, res[0], res[1], n)
        if not ok:
            ds /= 2
            if ds < mp.mpf(10)**-14:
                raise Fail(f'ray stalled s={s} ')
            continue
        b, c, z = res[0], res[1], res[2]
        s = s2
        steps += 1
        if res[3] <= 5:
            ds = min(ds*1.3, mp.mpf('0.4'))
        g, z3, ms = small_g(b, c, z)
        if abs(g) <= eps:
            return b, c, z3, g, steps
        if s > 60:
            raise Fail(f'ray reached s={s}: |g|={abs(g)}, not a puncture')


def g_path(b, c, z, g_from, g_to, n, tol, dmax=0.12):
    """Follow log-linear path in g from g_from to g_to (complex logs given)."""
    lf = mp.log(g_from)
    lt = mp.log(g_to)
    # keep the winding that the caller intends: caller may pass (log values) directly
    return g_path_log(b, c, z, lf, lt, n, tol, dmax)


def g_path_log(b, c, z, lf, lt, n, tol, dmax=0.12):
    L = abs(lt-lf)
    tau = mp.mpf(0)
    dt = min(mp.mpf(1), mp.mpf(dmax)/L) if L > 0 else mp.mpf(1)
    steps = 0
    while tau < 1:
        t2 = min(mp.mpf(1), tau+dt)
        g = mp.exp(lf + (lt-lf)*t2)
        res = newton_g(b, c, z, g, n, tol)
        ok = res is not None and guard(b, c, res[0], res[1], n)
        if not ok:
            dt /= 2
            if dt*L < mp.mpf(10)**-12:
                raise Fail('g path stalled')
            continue
        b, c, z = res[0], res[1], res[2]
        tau = t2
        steps += 1
        if res[3] <= 5:
            dt = min(dt*1.25, mp.mpf(dmax)/L)
    return b, c, z, steps


def snap_to_circle(b, c, z, g, eps, n, tol):
    """Move to |g|=eps at the nearest multiple of pi; return (b,c,z,h,offset)."""
    ang = mp.arg(g)
    h = int(mp.nint(ang/mp.pi))
    b, c, z, st = g_path_log(b, c, z, mp.log(g), mp.log(eps)+1j*h*mp.pi, n, tol)
    return b, c, z, h, ang - h*mp.pi


def close(x, y, rtol):
    return abs(x[0]-y[0]) <= rtol*abs(x[0]) and abs(x[1]-y[1]) <= rtol*abs(x[1])


def loop_points(b, c, z, h0, eps, n, tol, rtol, maxhalf=400):
    """Half-turns around the puncture; L[j] is the point with g=eps*exp(i pi (h0+j))."""
    L = [(b, c, z)]
    h = h0
    for j in range(1, maxhalf+1):
        b, c, z, st = g_path_log(b, c, z, mp.log(eps)+1j*h*mp.pi,
                                 mp.log(eps)+1j*(h+1)*mp.pi, n, tol)
        h += 1
        if j % 2 == 0 and close((b, c), L[0][:2], rtol):
            return L
        L.append((b, c, z))
    raise Fail('loop did not close')


def is_real_pt(pt, rtol):
    return abs(mp.im(pt[0])) <= rtol*abs(pt[0]) and abs(mp.im(pt[1])) <= rtol*abs(pt[1])


# ---------------------------------------------------------------- E. real arc tracer
def pattern(a, n):
    fin = sorted((a[i], i) for i in range(n) if i != 1)
    cyc = [i for _, i in fin] + [1]
    k = cyc.index(0)
    cyc = cyc[k:] + cyc[:k]
    rev = [cyc[0]] + cyc[:0:-1]
    return tuple(min(cyc, rev))


def real_g(b, c):
    """|g| and sign for a real map, double precision, chart values possibly huge."""
    r = np.roots([1.0, -1.0, -float(b), -float(c)])
    ms = [(-float(b)/z**2 - 2*float(c)/z**3) for z in r]
    m = max(ms, key=abs)
    return 1/abs(m), m


class Arc:
    pass


def omega_gap(a, b, c, n, pat):
    """Index i such that omega lies between pat[i] and pat[i+1] on RP^1."""
    w = -2*c/b if b != 0 else None
    pos = []
    for i in pat:
        pos.append(math.inf if a[i] is None else float(a[i]))
    # angular coordinate on RP^1
    th = [math.pi/2 if math.isinf(p) else math.atan(p) for p in pos]
    wt = math.pi/2 if w is None or math.isinf(float(w)) else math.atan(float(w))
    m = len(th)
    # orientation of pat on the circle (mod pi)
    s = sum(((th[(i+1) % m]-th[i]) % math.pi) for i in range(m))
    orient = 1 if abs(s-math.pi) < 1e-6 else -1
    for i in range(m):
        x = th[i]
        y = th[(i+1) % m]
        if orient == 1:
            if ((wt-x) % math.pi) <= ((y-x) % math.pi):
                return i
        else:
            if ((x-wt) % math.pi) <= ((x-y) % math.pi):
                return i
    return None


def trace_real(b0, c0, n, eps, dps, direction_out=True, ds0=0.01, maxsteps=400000, log=None, record=None):
    """Real arc from (b0,c0) with |g|=eps outward until |g|<=eps again.

    Runs in doubles while they suffice, otherwise mpmath. Returns dict with raw path summary.
    """
    sc = 1 if c0 > 0 else -1
    mode = 'mp'
    mp.mp.dps = dps
    u = mp.asinh(b0)
    v = mp.log(abs(c0))

    def bc(u_, v_):
        if mode == 'mp':
            return mp.sinh(u_), sc*mp.exp(v_)
        return math.sinh(u_), sc*math.exp(v_)

    def res(u_, v_):
        b, c = bc(u_, v_)
        one = mp.mpf(1) if mode == 'mp' else 1.0
        R, Rb, Rc = R_all(b, c, n, one)
        Ru = Rb*(mp.cosh(u_) if mode == 'mp' else math.cosh(u_))
        Rv = Rc*c
        return R, Ru, Rv

    # initial tangent: choose sign with |g| increasing
    R, Ru, Rv = res(u, v)
    nn = mp.sqrt(Ru*Ru+Rv*Rv)
    t = (-Rv/nn, Ru/nn)
    trial = mp.mpf(10)**-6
    bt, ct = bc(u+trial*t[0], v+trial*t[1])
    g1 = real_g(bt, ct)[0]
    bt, ct = bc(u-trial*t[0], v-trial*t[1])
    g2 = real_g(bt, ct)[0]
    if (g1 < g2) == bool(direction_out):
        t = (-t[0], -t[1])
    prev_t = t
    ds = ds0
    pat0 = pattern(orbit(*bc(u, v), n), n)
    steps = 0
    left = False
    gaps = []
    polys = []
    nreal = []
    prev_pv = None
    mp_steps = 0
    gmax = 0.0
    while steps < maxsteps:
        steps += 1
        # double mode when conditioning allows: decide each step from the point itself
        b, c = bc(u, v)
        fb, fc = float(b), float(c)
        a_now = orbit(fb, fc, n) if abs(fc) > 1e-300 else None
        use_double = False
        if a_now is not None:
            try:
                fin = sorted(float(a_now[i]) for i in range(n) if i != 1)
                mingap = min((fin[i+1]-fin[i])/(1+abs(fin[i])) for i in range(n-2))
                use_double = mingap > 1e-5 and abs(fc) > 1e-4 and max(abs(x) for x in fin) < 1e6
            except (ZeroDivisionError, OverflowError):
                use_double = False
        mode_new = 'fl' if use_double else 'mp'
        if mode_new != mode:
            mode = mode_new
            if mode == 'fl':
                u, v = float(u), float(v)
            else:
                u, v = mp.mpf(u), mp.mpf(v)
        if mode == 'mp':
            mp_steps += 1
        # tangent
        R, Ru, Rv = res(u, v)
        nn2 = Ru*Ru+Rv*Rv
        nn = mp.sqrt(nn2) if mode == 'mp' else math.sqrt(nn2)
        t = (-Rv/nn, Ru/nn)
        if t[0]*prev_t[0]+t[1]*prev_t[1] < 0:
            t = (-t[0], -t[1])
        tol = (mp.mpf(10)**(-dps//2)) if mode == 'mp' else 1e-12
        while True:
            up = u + ds*t[0]
            vp = v + ds*t[1]
            ok = False
            try:
                for it in range(7):
                    R2, Ru2, Rv2 = res(up, vp)
                    n2 = Ru2*Ru2+Rv2*Rv2
                    dd = R2/n2
                    du = dd*Ru2
                    dv = dd*Rv2
                    up -= du
                    vp -= dv
                    if abs(du)+abs(dv) < tol:
                        ok = True
                        break
            except (ZeroDivisionError, OverflowError, ValueError):
                ok = False
            if ok:
                # tangent turning and marked-point guard
                _, Ru3, Rv3 = res(up, vp)
                nn3 = mp.sqrt(Ru3*Ru3+Rv3*Rv3) if mode == 'mp' else math.sqrt(Ru3*Ru3+Rv3*Rv3)
                tn = (-Rv3/nn3, Ru3/nn3)
                dot = abs(tn[0]*t[0]+tn[1]*t[1])
                b1, c1 = bc(up, vp)
                if dot < 0.985 or not guard(b, c, b1, c1, n, 0.3):
                    ok = False
            if ok:
                break
            ds = ds/2
            if ds < 1e-16:
                return dict(status='stalled', b=b, c=c, steps=steps, mp_steps=mp_steps,
                            pattern=pat0, gaps=gaps, polys=polys, gmax=gmax)
        prev_t = t
        u, v = up, vp
        if it <= 2:
            ds = min(ds*1.35, 0.08)
        b, c = bc(u, v)
        if record is not None:
            record.append((b, c) if mode == 'mp' else (mp.mpf(b), mp.mpf(c)))
        a = orbit(b, c, n)
        pi_ = pattern(a, n)
        if pi_ != pat0:
            return dict(status='pattern-change', b=b, c=c, steps=steps, pattern=pat0,
                        pattern2=pi_, gaps=gaps, polys=polys, gmax=gmax)
        gv, mv = real_g(b, c)
        gmax = max(gmax, gv)
        # omega events
        gi = omega_gap(a, b, c, n, pat0)
        if not gaps or gaps[-1] != gi:
            gaps.append(gi)
        # polynomial event: omega fixed
        if b != 0:
            w = -2*c/b
            pv = float((1 + b/w + c/(w*w) - w)/(1+abs(w)))
            if prev_pv is not None and pv*prev_pv < 0 and abs(pv) < 1 and abs(prev_pv) < 1:
                polys.append((float(b), float(c)))
            prev_pv = pv
        else:
            prev_pv = None
        if gv > 3*eps:
            left = True
        if gv > eps*1.05 and not left and steps > 2000:
            left = True
        if left and gv <= eps:
            return dict(status='landed', b=b, c=c, steps=steps, mp_steps=mp_steps, pattern=pat0,
                        gaps=gaps, polys=polys, gmax=gmax, sign=1 if mv > 0 else -1)
    return dict(status='maxsteps', b=b, c=c, steps=steps, mp_steps=mp_steps, pattern=pat0,
                gaps=gaps, polys=polys, gmax=gmax)


def land_real(b, c, sgn_mult, eps, n, tol):
    """Newton onto g=sgn*eps exactly from a nearby real point; returns mp reals."""
    b = mp.mpf(b)
    c = mp.mpf(c)
    fp = mp.polyroots([1, -1, -b, -c], maxsteps=300, extraprec=400)
    ms = [mult(w, b, c) for w in fp]
    i = max(range(3), key=lambda j: abs(ms[j]))
    z = mp.mpc(fp[i])
    g = mp.mpf(sgn_mult)*eps
    g0 = 1/ms[i]
    bb, cc, zz, st = g_path_log(mp.mpc(b), mp.mpc(c), z, mp.log(mp.mpc(g0)), mp.log(mp.mpc(g)), n, tol)
    return mp.re(bb), mp.re(cc)
