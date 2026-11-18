"""E93: Fuchsian-group tools for genus 24 (where E92's word enumeration up to 2 r_D + 0.25 + r_D is out of reach).

1. closure_mp: the Dirichlet polygon about 0 from the side-pairing generators, in 200-bit arithmetic (gmpy2): the generators of a
   genus-24 development displace 0 by up to ~40, so their products cannot be formed in double. Rounds of 'keep the K shortest
   distinct elements of S u S.S' (E92's stage A), until the polygon cut out by the bisector half-planes of the kept elements is
   bounded with area 4 pi (g - 1). Certification (E92 U4): every kept element is in the group, so that polygon contains the
   Dirichlet domain; equal area means it IS the Dirichlet domain.
2. dirichlet_at (double): the Dirichlet polygon about any point p, from seed elements (e.g. the faces of D), by rounds of
   'add every product of two face elements whose bisector cuts the current polygon' until the area is 4 pi (g - 1) (same
   certification). Used to compute surface distances by folding (fold q into D_p: d(p, fold(q)) = min_g d(p, g q)), which
   replaces E92's 'min over all explored elements'.
Elements are SU(1,1) matrices [[a, b], [conj b, conj a]]."""
import json
import numpy as np
import gmpy2
from gmpy2 import mpc, mpfr
TAU = 2*np.pi


# ---------------- multiprecision closure ----------------
def mp_ctx(prec=200):
    gmpy2.get_context().precision = prec


def mp_from_np(g):
    """complex 2x2 (any normalisation) -> (a, b) gmpy2, SU(1,1)-projected."""
    a = (mpc(complex(g[0, 0])) + mpc(complex(g[1, 1])).conjugate())/2
    b = (mpc(complex(g[0, 1])) + mpc(complex(g[1, 0])).conjugate())/2
    return mp_norm((a, b))


def mp_norm(e):
    a, b = e
    s = gmpy2.sqrt(abs(a)**2 - abs(b)**2)
    a, b = a/s, b/s
    if a.real < 0 or (a.real == 0 and a.imag < 0):
        a, b = -a, -b
    return (a, b)


def mp_mul(e, f):
    a1, b1 = e; a2, b2 = f
    return mp_norm((a1*a2 + b1*b2.conjugate(), a1*b2 + b1*a2.conjugate()))


def mp_inv(e):
    a, b = e
    return mp_norm((a.conjugate(), -b))


def mp_d(e):
    """displacement of 0 (hyperbolic distance d(0, g0)) and direction of g(0) = b / conj(a)."""
    a, b = e
    return float(2*gmpy2.acosh(abs(a))), float(gmpy2.phase(b) + gmpy2.phase(a))


def mp_key(e):
    """(d, theta sinh d) on a 1e-3 grid (E92's HSet idea): distinct orbit points of the surface group are >= the systole apart;
    the same element reached by different words differs by the relation error of the development (U3 cycles ~1e-7, amplified
    along words), so a duplicate can land in a neighbouring cell; such survivors only add redundant half-planes, and the
    near-identity products they create are removed by the identity filter d > 1e-3 (E92)."""
    a, b = e
    d = 2*gmpy2.acosh(abs(a))
    th = gmpy2.phase(b) + gmpy2.phase(a)
    th = th - 2*gmpy2.const_pi()*gmpy2.rint(th/(2*gmpy2.const_pi()))
    return (int(gmpy2.rint(d*10**3)), int(gmpy2.rint(th*gmpy2.sinh(d)*10**3)))


def mp_to_np(e):
    a, b = e
    a, b = complex(a), complex(b)
    return np.array([[a, b], [b.conjugate(), a.conjugate()]])


IDCUT = None


def closure_mp(gens, area_target, K0=800, rounds=40, log=print):
    """gens: list of (a, b). Stage A (E92): the K0 shortest distinct elements of S u S.S, two rounds. Stage B (E93): face-local
    closure: F = face elements of the current polygon (and inverses); add every product f g (f, g in F) whose bisector cuts the
    current polygon (some vertex strictly outside); repeat until the polygon is bounded with area 4 pi (g - 1). Products and
    displacements in 200 bits; the polygon is clipped in double. Returns (elements, Dirichlet result)."""
    global IDCUT
    IDCUT = gmpy2.cosh(mpfr('5e-4'))
    S = {}
    for g in gens:
        for h in (g, mp_inv(g)):
            if mp_d(h)[0] > 1e-3: S[mp_key(h)] = h
    def prods(A, B, window=None):
        out = {}
        dA = [mp_d(e)[0] for e in A]; dB = np.array([mp_d(e)[0] for e in B]); o = np.argsort(dB); B = [B[i] for i in o]; dB = dB[o]
        for e, de in zip(A, dA):
            lo, hi = (0, len(B)) if window is None else (np.searchsorted(dB, de - window), np.searchsorted(dB, de + window, 'right'))
            for f in B[lo:hi]:
                h = mp_mul(e, f)
                if abs(h[0]) < IDCUT: continue
                out[mp_key(h)] = h
        return out
    for it in range(2):
        items = sorted(S.values(), key=lambda e: mp_d(e)[0])[:K0]
        dK = mp_d(items[-1])[0]
        S = dict((mp_key(e), e) for e in items); S.update(prods(items, items, dK))
        items = sorted(S.values(), key=lambda e: mp_d(e)[0])[:K0]
        S = dict((mp_key(e), e) for e in items)
        log(f'U4 stage A{it}: kept {len(items)}, d range {mp_d(items[0])[0]:.4f} .. {mp_d(items[-1])[0]:.4f}')
    res = None
    for it in range(rounds):
        items = list(S.values())
        res = dirichlet_np([mp_to_np(e) for e in items])
        log(f'U4 stage B{it}: elements {len(items)}, faces {res["nfaces"]}, bounded {res["bounded"]}, area {res["area"]:.12f} '
            f'(diff {res["area"] - area_target:.3e}), circumradius {res["rD"]:.6f}', flush=True)
        if res['bounded'] and abs(res['area'] - area_target) < 1e-8:
            break
        fidx = sorted(set(t for t in res['tag'] if t >= 0))
        F = [items[res['ok_index'][t]] for t in fidx]
        F = F + [mp_inv(e) for e in F]
        cand = prods(F, F)
        V = res['verts']
        keys = list(cand.keys()); C = [cand[k] for k in keys]
        U, T, d = halfplanes(np.array([mp_to_np(e) for e in C]))
        viol = ((V[None, :]*np.conj(U[:, None])).real - T[:, None]).max(1)
        cut = [C[i] for i in np.where(viol > 1e-12)[0] if keys[i] not in S]
        log(f'   face elements {len(F)}, products {len(cand)}, cutting the polygon {len(cut)}', flush=True)
        if not cut:
            break
        S = dict((mp_key(e), e) for e in F)
        for e in cut: S[mp_key(e)] = e
    items = list(S.values())
    return items, res


# ---------------- double-precision Dirichlet polygons ----------------
def mob_ap(g, z): return (g[0, 0]*z + g[0, 1])/(g[1, 0]*z + g[1, 1])
def mob_norm(g): return g/np.sqrt(np.linalg.det(g))
def Mq(q): return mob_norm(np.array([[1, q], [np.conj(q), 1]], complex))
def to_klein(z): return 2*z/(1 + np.abs(z)**2)
def to_poinc(k): return k/(1 + np.sqrt(np.maximum(0, 1 - np.abs(k)**2)))
def hd0(w): return 2*np.arctanh(np.minimum(np.abs(w), 1 - 1e-16))
def hdist(z, w): return 2*np.arcsinh(np.abs(z - w)/np.sqrt((1 - np.abs(z)**2)*(1 - np.abs(w)**2)))


def su11(G):
    """project a stack of 2x2 complex matrices onto SU(1,1)."""
    G = np.asarray(G, complex)
    d = G[..., 0, 0]*G[..., 1, 1] - G[..., 0, 1]*G[..., 1, 0]; G = G/np.sqrt(d)[..., None, None]
    a = (G[..., 0, 0] + np.conj(G[..., 1, 1]))/2; b = (G[..., 0, 1] + np.conj(G[..., 1, 0]))/2
    s = np.sqrt(np.abs(a)**2 - np.abs(b)**2); a, b = a/s, b/s
    Q = np.empty_like(G); Q[..., 0, 0] = a; Q[..., 0, 1] = b; Q[..., 1, 0] = np.conj(b); Q[..., 1, 1] = np.conj(a)
    return Q


def halfplanes(G):
    """Klein half-plane {k : Re(k conj u) <= t} of the bisector of 0 and g(0): u = direction, t = tanh(d/2)."""
    G = np.asarray(G)
    a, b = G[:, 0, 0], G[:, 0, 1]
    z = b/np.conj(a)
    d = 2*np.arccosh(np.maximum(np.abs(a), 1.0))
    return np.exp(1j*np.angle(z)), np.tanh(d/2), d


def clip(U, T, R0=1 - 1e-15, nstart=16384):
    """convex polygon (Klein) = regular nstart-gon at Klein radius R0 cut by the half-planes Re(k conj u) <= t; returns vertices
    and the half-plane index of the edge starting at each vertex (-1 = artificial outer boundary). E93: the start polygon's
    inradius is R0 cos(pi/nstart) = 1 - 1.8e-8 (hyperbolic ~9.2); E92's 512-gon at 0.999999 had inradius ~5.8 (E92's
    circumradius 3.6 < 5.8, but genus 24's ~6 is not). Each cut is vectorized (convexity: the inside run is contiguous)."""
    V = np.exp(1j*np.linspace(0, TAU, nstart + 1)[:-1])*R0
    tag = -np.ones(nstart, int)
    for hi in np.argsort(T):
        u, t = U[hi], T[hi]
        f = (V*np.conj(u)).real - t
        ins = f <= 0
        if ins.all():
            continue
        if not ins.any():
            return np.zeros(0, complex), []
        n = len(V)
        nxt = np.roll(ins, -1)
        ex = np.where(ins & ~nxt)[0]; en = np.where(~ins & nxt)[0]
        if len(ex) != 1 or len(en) != 1:
            raise RuntimeError('clip: polygon not convex or degenerate cut')
        ex, en = int(ex[0]), int(en[0])
        e1 = (ex + 1) % n; e2 = (en + 1) % n
        X1 = V[ex] + (V[e1] - V[ex])*(f[ex]/(f[ex] - f[e1]))
        X2 = V[en] + (V[e2] - V[en])*(f[en]/(f[en] - f[e2]))
        run = (np.arange(e2, e2 + ((ex - e2) % n) + 1)) % n           # inside vertices e2 .. ex
        V = np.concatenate([V[run], [X1, X2]])
        tag = np.concatenate([tag[run], [int(hi), tag[en]]])
    n = len(V)
    keep = [a for a in range(n) if tag[a] != tag[a - 1]]
    return V[keep], [int(tag[a]) for a in keep]


def poly_area(verts_k):
    zp = to_poinc(verts_k); n = len(zp); s = 0.0
    for a in range(n):
        M = np.linalg.inv(Mq(zp[a])); u, w = mob_ap(M, zp[a - 1]), mob_ap(M, zp[(a + 1) % n])
        s += abs(np.angle(w/u))
    return (n - 2)*np.pi - s


def dirichlet_np(G):
    """Dirichlet-type polygon about 0 cut out by the bisectors of the elements G (double)."""
    G = su11(np.array(G))
    U, T, d = halfplanes(G)
    ok = d > 1e-3
    okidx = np.where(ok)[0]
    G, U, T, d = G[ok], U[ok], T[ok], d[ok]
    V, tag = clip(U, T)
    bounded = len(tag) > 2 and all(t >= 0 for t in tag)
    area = poly_area(V) if bounded else np.inf
    rD = float(hd0(to_poinc(V)).max()) if bounded else np.inf
    return dict(verts=V, tag=tag, face_g=G[[t for t in tag]] if bounded else None, G=G, bounded=bounded, area=area, rD=rD,
                nfaces=len(tag), ok_index=okidx)


def face_pairing(res):
    """U4: face a with element g is carried by g^-1 onto its partner face; partner element must be g^-1."""
    V = res['verts']; FG = res['face_g']; nf = len(V)
    zp = to_poinc(V)
    P, Q = zp, np.roll(zp, -1)
    perr, unmatched, pair_of = 0.0, 0, []
    for a in range(nf):
        gi = np.linalg.inv(FG[a]); P2, Q2 = mob_ap(gi, P[a]), mob_ap(gi, Q[a])
        dd = np.maximum(hdist(P2, Q), hdist(Q2, P))
        b = int(np.argmin(dd)); pair_of.append(b); perr = max(perr, float(dd[b]))
        e = mob_norm(FG[a] @ FG[b])
        if min(np.abs(e - np.eye(2)).max(), np.abs(e + np.eye(2)).max()) > 1e-6:
            unmatched += 1
    return perr, unmatched, np.array(pair_of)


def vertex_cycles(res, pair_of):
    V = res['verts']; nf = len(V); zp = to_poinc(V)
    par = list(range(nf))
    def fnd(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    for a in range(nf):
        b = pair_of[a]
        for p_, q_ in ((a, (b + 1) % nf), ((a + 1) % nf, b)):
            rp, rq = fnd(p_), fnd(q_)
            if rp != rq: par[rp] = rq
    cls = {}
    for a in range(nf): cls.setdefault(fnd(a), []).append(a)
    ang = []
    for a in range(nf):
        Mi = np.linalg.inv(Mq(zp[a])); u, w = mob_ap(Mi, zp[a - 1]), mob_ap(Mi, zp[(a + 1) % nf]); ang.append(abs(np.angle(w/u)))
    return cls, [sum(ang[a] for a in m) for m in cls.values()]


def products_cut(F, V, chunk=2000):
    """all products f g (f, g in F) whose bisector half-plane cuts the polygon with Klein vertices V (some vertex strictly outside)."""
    F = np.asarray(F); out = []
    for s0 in range(0, len(F), chunk):
        P = su11(np.einsum('nij,mjk->nmik', F[s0:s0 + chunk], F).reshape(-1, 2, 2))
        U, T, d = halfplanes(P)
        ok = d > 1e-3
        P, U, T = P[ok], U[ok], T[ok]
        viol = ((V[None, :]*np.conj(U[:, None])).real - T[:, None]).max(1)
        out.append(P[viol > 1e-12])
    return np.concatenate(out) if out else np.zeros((0, 2, 2), complex)


def dedupe(G, tol=1e-4):
    """distinct elements by orbit point, robust to grid straddling: points (d, sinh d cos th, sinh d sin th) (locally ~ hyperbolic
    distance; distinct orbit points are >= the systole apart) within tol are merged (KD-tree pairs + connected components)."""
    from scipy.spatial import cKDTree
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    G = np.asarray(G)
    if len(G) < 2: return G
    a, b = G[:, 0, 0], G[:, 0, 1]
    d = 2*np.arccosh(np.maximum(np.abs(a), 1.0)); th = np.angle(b/np.conj(a)); sh = np.sinh(d)
    X = np.stack([d, sh*np.cos(th), sh*np.sin(th)], 1)
    pr = cKDTree(X).query_pairs(tol, output_type='ndarray')
    if not len(pr): return G
    A = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(len(G), len(G)))
    _, lab = connected_components(A, directed=False)
    _, idx = np.unique(lab, return_index=True)
    return G[np.sort(idx)]


def products_cut_capped(F, V, dcap, chunk=2000):
    """products f g (f, g in F) whose bisector cuts the polygon V and whose displacement is <= dcap (no far products in double)."""
    F = np.asarray(F); out = []
    for s0 in range(0, len(F), chunk):
        P = su11(np.einsum('nij,mjk->nmik', F[s0:s0 + chunk], F).reshape(-1, 2, 2))
        U, T, d = halfplanes(P)
        ok = (d > 1e-3) & (d <= dcap)
        P, U, T = P[ok], U[ok], T[ok]
        viol = ((V[None, :]*np.conj(U[:, None])).real - T[:, None]).max(1)
        out.append(P[viol > 1e-12])
    return np.concatenate(out) if out else np.zeros((0, 2, 2), complex)


def _closure_at(q, F0, area_target, dcap, rounds=12, log=None):
    """Dirichlet polygon about q from group elements F0 (0-frame): start from F0 u F0^-1 u F0.F0 conjugated to q = 0 and capped
    at displacement dcap; then rounds of capped cutting products of face elements."""
    Mp = Mq(q); Mpi = np.linalg.inv(Mp)
    F = su11(np.einsum('ij,njk,kl->nil', Mpi, np.asarray(F0), Mp))
    F = dedupe(np.concatenate([F, su11(np.linalg.inv(F))]))
    FF = su11(np.einsum('nij,mjk->nmik', F, F).reshape(-1, 2, 2))
    d = halfplanes(FF)[2]
    S = dedupe(np.concatenate([F, FF[(d > 1e-3) & (d <= dcap)]]))
    S = S[halfplanes(S)[2] > 1e-3]
    res = dirichlet_np(S)
    for it in range(rounds):
        if res['bounded'] and abs(res['area'] - area_target) < CERT_TOL:
            break
        Fc = S[[res['ok_index'][t] for t in sorted(set(t for t in res['tag'] if t >= 0))]]
        Fc = dedupe(np.concatenate([Fc, su11(np.linalg.inv(Fc))]))
        add = products_cut_capped(Fc, res['verts'], dcap)
        if log: log(f'   closure at |q| = {float(hd0(np.array([q]))[0]):.3f} round {it}: faces {res["nfaces"]}, bounded {res["bounded"]}, '
                    f'area {res["area"]:.10f}, cutting products {len(add)}')
        if not len(add):
            break
        S = dedupe(np.concatenate([Fc, add]))
        res = dirichlet_np(S)
    res['Mp'] = Mp
    res['certified'] = bool(res['bounded'] and abs(res['area'] - area_target) < CERT_TOL)
    return res


CERT_TOL = 1e-6      # D_p (double precision, conjugated elements): area round-off reaches ~1e-8 (E93 log); U4 itself uses 1e-8


def dirichlet_at(p, F0, area_target, rD, step=0.75, log=None):
    """Dirichlet polygon about p (Poincare disk) from the faces F0 of the Dirichlet polygon about 0 (circumradius rD): the base
    point walks from 0 to p along the geodesic in steps <= step, each polygon seeded by the faces of the previous one; every
    product formed in double has displacement <= 2 rD + 2 (no far products). Certified when bounded with area 4 pi (g - 1) to
    1e-8. Result in the p = 0 frame plus Mp = Mq(p)."""
    dp = float(hd0(np.array([p]))[0])
    k = max(1, int(np.ceil(dp/step)))
    u = p/abs(p) if abs(p) > 0 else 1
    F = np.asarray(F0); dcap = 2*rD + 2
    for j in range(1, k + 1):
        q = u*np.tanh(dp*j/k/2)
        res = _closure_at(q, F, area_target, dcap, log=log)
        if not res['certified']:
            return res
        F = su11(np.einsum('ij,njk,kl->nil', res['Mp'], res['face_g'], np.linalg.inv(res['Mp'])))
    return res


def folder(res):
    """fold function for the polygon of res (in its own coordinates): returns (folded points, number of steps)."""
    FG = res['face_g']; U, T, _ = halfplanes(FG); Gi = su11(np.linalg.inv(FG))
    def fold(w, maxit=2000):
        w = np.array(w, complex, copy=True); Fm = np.tile(np.eye(2, dtype=complex), (len(w), 1, 1))
        for it in range(maxit):
            k = to_klein(w); v = (k[:, None]*np.conj(U[None, :])).real - T[None, :]
            a = np.argmax(v, 1); out = v[np.arange(len(w)), a] > 1e-13
            if not out.any():
                return w, Fm, it
            for f in np.unique(a[out]):
                m = out & (a == f); G = Gi[f]
                w[m] = (G[0, 0]*w[m] + G[0, 1])/(G[1, 0]*w[m] + G[1, 1]); Fm[m] = np.einsum('ij,njk->nik', G, Fm[m])
        raise RuntimeError('fold did not converge')
    return fold


# ---------------- Dirichlet polygon about p in 200 bits (E93 U5 / viewer) ----------------
def mp_Mq(q):
    q = mpc(complex(q)); s = gmpy2.sqrt(1 - abs(q)**2)
    return (mpc(1)/s, q/s)


def mp_conj(Minv, g, M):
    return mp_mul(Minv, mp_mul(g, M))


def mp_halfplanes(L):
    a = [e[0] for e in L]; b = [e[1] for e in L]
    d = np.array([float(2*gmpy2.acosh(abs(x))) for x in a])
    th = np.array([float(gmpy2.phase(y) + gmpy2.phase(x)) for x, y in zip(a, b)])
    return np.exp(1j*th), np.tanh(d/2), d


def mp_dedupe(L):
    out = {}
    for e in L:
        out.setdefault(mp_key(e), e)
    return list(out.values())


def load_face_mp(npz):
    fu_ctx = mp_ctx(200)
    return [mp_norm((mpc(mpfr(a), mpfr(b)), mpc(mpfr(c), mpfr(d)))) for a, b, c, d in json.loads(str(npz['face_g_str']))]


def _closure_at_mp(q, F0, area_target, dcap, rounds=12, log=None):
    """as _closure_at, all products and conjugations in 200 bits (F0: list of (a, b) in the 0-frame)."""
    global IDCUT
    IDCUT = gmpy2.cosh(mpfr('5e-4'))
    M = mp_Mq(q); Mi = mp_inv(M)
    F = mp_dedupe([mp_conj(Mi, g, M) for g in F0] + [mp_conj(Mi, mp_inv(g), M) for g in F0])
    dF = mp_halfplanes(F)[2]
    S = dict((mp_key(e), e) for e in F)
    o = np.argsort(dF); Fs = [F[i] for i in o]; dFs = dF[o]
    for e, de in zip(Fs, dFs):                                    # products with displacement <= dcap (d(fg) >= |d_f - d_g|)
        hi = np.searchsorted(dFs, de + dcap, 'right'); lo = np.searchsorted(dFs, de - dcap)
        for f in Fs[lo:hi]:
            h = mp_mul(e, f)
            if abs(h[0]) < IDCUT: continue
            if float(2*gmpy2.acosh(abs(h[0]))) <= dcap: S.setdefault(mp_key(h), h)
    L = list(S.values())
    res = dirichlet_np([mp_to_np(e) for e in L])
    for it in range(rounds):
        if res['bounded'] and abs(res['area'] - area_target) < 1e-8:
            break
        Fc = [L[res['ok_index'][t]] for t in sorted(set(t for t in res['tag'] if t >= 0))]
        Fc = mp_dedupe(Fc + [mp_inv(e) for e in Fc])
        cand = {}
        for e in Fc:
            for f in Fc:
                h = mp_mul(e, f)
                if abs(h[0]) < IDCUT: continue
                cand.setdefault(mp_key(h), h)
        C = list(cand.values())
        U, T, d = mp_halfplanes(C)
        V = res['verts']
        viol = ((V[None, :]*np.conj(U[:, None])).real - T[:, None]).max(1)
        add = [C[i] for i in np.where((viol > 1e-12) & (d <= 2*dcap))[0]]
        if log: log(f'   mp closure at |q| = {float(hd0(np.array([q]))[0]):.3f} round {it}: faces {res["nfaces"]}, bounded {res["bounded"]}, '
                    f'area {res["area"]:.10f}, cutting products {len(add)}')
        if not add:
            break
        L = mp_dedupe(Fc + add)
        res = dirichlet_np([mp_to_np(e) for e in L])
    res['Mp'] = Mq(q); res['Mp_mp'] = M
    res['face_mp'] = [L[res['ok_index'][t]] for t in res['tag']] if res['bounded'] else None
    res['certified'] = bool(res['bounded'] and abs(res['area'] - area_target) < 1e-8)
    return res


def dirichlet_at_mp(p, F0, area_target, rD, step=0.75, log=None):
    """Dirichlet polygon about p (Poincare disk) from the 200-bit faces F0 of the Dirichlet polygon about 0 (circumradius rD):
    the base point walks from 0 to p in steps <= step, each polygon seeded by the previous faces (all in 200 bits; products
    capped at displacement 2 rD + 2 for seeding). Certified when bounded with area 4 pi (g - 1) to 1e-8 (as U4).
    Result in the p = 0 frame plus Mp = Mq(p)."""
    dp = float(hd0(np.array([p]))[0])
    k = max(1, int(np.ceil(dp/step)))
    u = p/abs(p) if abs(p) > 0 else 1
    F = F0; dcap = 2*rD + 2
    for j in range(1, k + 1):
        q = u*np.tanh(dp*j/k/2)
        res = _closure_at_mp(q, F, area_target, dcap, log=log)
        if not res['certified']:
            return res
        M = res['Mp_mp']; Mi = mp_inv(M)
        F = [mp_conj(M, g, Mi) for g in res['face_mp']]
    return res
