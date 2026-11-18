"""E93-U: the development step in 160-bit arithmetic (gmpy2) from the refined metric (uniform_u_hi + uniform_u_lo, refine.py).
Same BFS order and placement rule as E92's devmp.py / develop.develop(); per-face hyperbolic lengths and angles are computed in
160 bits (E92 took them from the double solve). Side pairings, U3 residuals and vertex-cycle products in 160 bits.
Returns positions (double, base point = centroid of the root face moved to 0), generators (double and 45-digit strings)."""
import numpy as np
import gmpy2
from gmpy2 import mpc, mpfr
from devmp import bfs
gmpy2.get_context().precision = 160
ONE = mpfr(1)


def face_geometry(T, E, l0, u):
    """for faces T (corner order as used by the development): (tanh(l_opp/2), exp(i A)) per corner, 160 bits."""
    from uniformize import setup
    half = mpfr('0.5')
    L = [2*gmpy2.asinh(half*mpfr(float(l0[k]))*gmpy2.exp(half*(u[int(i)] + u[int(j)]))) for k, (i, j) in enumerate(E)]
    fe = setup(T)
    th, eA = [], []
    for f in range(len(T)):
        a, b, c = L[fe[f, 0]], L[fe[f, 1]], L[fe[f, 2]]
        s = (a + b + c)/2
        ss, sa, sb, sc = gmpy2.sinh(s), gmpy2.sinh(s - a), gmpy2.sinh(s - b), gmpy2.sinh(s - c)
        A = 2*gmpy2.atan(gmpy2.sqrt(sb*sc/(ss*sa)))
        B = 2*gmpy2.atan(gmpy2.sqrt(sc*sa/(ss*sb)))
        C = 2*gmpy2.atan(gmpy2.sqrt(sa*sb/(ss*sc)))
        th.append((gmpy2.tanh(a/2), gmpy2.tanh(b/2), gmpy2.tanh(c/2)))
        eA.append(tuple(mpc(gmpy2.cos(x), gmpy2.sin(x)) for x in (A, B, C)))
    return th, eA


def hd(z, w): return 2*gmpy2.asinh(abs(z - w)/gmpy2.sqrt((1 - abs(z)**2)*(1 - abs(w)**2)))
def M_to(p):
    s = gmpy2.sqrt(1 - abs(p)**2); return (ONE/s, p/s)                 # (a, b) of [[1, p], [conj p, 1]]/s: 0 -> p
def M_from(p):
    s = gmpy2.sqrt(1 - abs(p)**2); return (ONE/s, -p/s)
def ap(e, z):
    a, b = e; return (a*z + b)/(b.conjugate()*z + a.conjugate())
def mul(e, f):
    a1, b1 = e; a2, b2 = f
    return (a1*a2 + b1*b2.conjugate(), a1*b2 + b1*a2.conjugate())
def inv(e):
    a, b = e; return (a.conjugate(), -b)
def iso_two(p1, q1, p2, q2):
    a = ap(M_from(p1), q1); b = ap(M_from(p2), q2)
    rho = (b/abs(b))/(a/abs(a)); s = gmpy2.sqrt(rho)
    R = (s, mpc(0))                                                        # rotation by rho: [[s, 0], [0, 1/s]], |s| = 1
    return mul(M_to(p2), mul(R, M_from(p1)))
def ident_err(e):
    a, b = e
    return float(min(max(abs(a - 1), abs(b)), max(abs(a + 1), abs(b))))


def mp_stage2(Tn, T, E, l0, u, cyc, sides, pairs, orig, log=print):
    root, order = bfs(Tn)
    th, eA = face_geometry(T, E, l0, u)
    nN = int(Tn.max()) + 1
    z = [None]*nN; err = 0.0; nre = 0
    i, j, k = Tn[root]
    z[i] = mpc(0); z[j] = mpc(th[root][2]); z[k] = th[root][1]*eA[root][0]
    for n_, f in enumerate(order[1:]):
        c = Tn[f]
        for r in range(3):
            p, q, w = c[r], c[(r + 1) % 3], c[(r + 2) % 3]
            if z[p] is not None and z[q] is not None: break
        zp, zq = z[p], z[q]; cp = zp.conjugate()
        qq = (zq - zp)/(ONE - cp*zq)
        wn = (qq/abs(qq))*th[f][(r + 1) % 3]*eA[f][r]
        zw = (wn + zp)/(ONE + cp*wn)
        if z[w] is None: z[w] = zw
        else:
            nre += 1; e = float(hd(z[w], zw))
            if e > err: err = e
    log(f'develop (160 bits): root face {root} | placed {sum(v is not None for v in z)} of {nN} | re-placements {nre}, '
        f'max re-placement error (hyperbolic) {err:.3e}')
    gens, res, sv_bad = [], 0.0, 0
    for s1, s2 in pairs:
        i1 = sides[s1]['idx']; i2 = sides[s2]['idx'][::-1]
        sv_bad += int(np.sum(orig[cyc[i1]] != orig[cyc[i2]]))
        Z1 = [z[cyc[t]] for t in i1]; Z2 = [z[cyc[t]] for t in i2]
        g = iso_two(Z1[0], Z1[-1], Z2[0], Z2[-1])
        res = max(res, max(float(hd(ap(g, a), b)) for a, b in zip(Z1, Z2)))
        gens.append((s1, s2, g))
    log(f'U3 (160 bits): side pairings {len(gens)} | max hyperbolic residual over all side points {res:.3e} | vertex mismatches {sv_bad}')
    gmap = {}
    for s1, s2, g in gens: gmap[s1] = (s2, g); gmap[s2] = (s1, inv(g))
    nS = len(sides); used = set(); cyc_err = 0.0; ncyc = 0; lens = []
    for s0 in range(nS):
        if s0 in used: continue
        g = (mpc(1), mpc(0)); s = s0; L = 0
        while True:
            used.add(s); partner, gs = gmap[s]; g = mul(gs, g); L += 1
            s = (partner + 1) % nS
            if s == s0 or L > 400: break
        cyc_err = max(cyc_err, ident_err(g)); ncyc += 1; lens.append(L)
    log(f'U3 (160 bits): vertex cycles {ncyc} | lengths {sorted(lens)} | max distance of products from +-I {cyc_err:.3e}')
    p0 = (z[Tn[root][0]] + z[Tn[root][1]] + z[Tn[root][2]])/3
    Mb, Mbi = M_from(p0), M_to(p0)
    zc = np.array([complex(ap(Mb, v)) for v in z])
    Gm = [mul(Mb, mul(g, Mbi)) for _, _, g in gens]
    G = [np.array([[complex(a), complex(b)], [complex(b).conjugate(), complex(a).conjugate()]]) for a, b in Gm]
    Gs = [[(f'{a.real:.45e}', f'{a.imag:.45e}'), (f'{b.real:.45e}', f'{b.imag:.45e}'),
           (f'{b.real:.45e}', f'{(-b.imag):.45e}'), (f'{a.real:.45e}', f'{(-a.imag):.45e}')] for a, b in Gm]
    return zc, G, root, dict(dev_err=err, U3_res=res, U3_cycle_err=cyc_err, cycles=ncyc, gens_str=Gs)
