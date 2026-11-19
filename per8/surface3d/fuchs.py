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


def save_els(path, L):
    json.dump([[f'{e[0].real:.60e}', f'{e[0].imag:.60e}', f'{e[1].real:.60e}', f'{e[1].imag:.60e}'] for e in L], open(path + '.tmp', 'w'))
    import os; os.replace(path + '.tmp', path)


def load_els(path):
    return [mp_norm((mpc(mpfr(a), mpfr(b)), mpc(mpfr(c), mpfr(d)))) for a, b, c, d in json.load(open(path))]


def closure_mp(gens, area_target, K0=2000, rounds=200, log=print, ckpt=None):   # E94: K0 2000 (822 faces expected)
    """gens: list of (a, b). Stage A (E92): the K0 shortest distinct elements of S u S.S, two rounds (E94: products screened by
    displacement in double, formed in 200 bits). Stage B (E93): face-local closure: F = face elements of the current polygon
    (and inverses); add every product f g (f, g in F) whose bisector cuts the current polygon; E94 escalation when no pair
    cuts: triple products of faces (triples_cut_screened), then products of faces with the stage-A set. Repeat until the
    polygon is bounded with area 4 pi (g - 1) (the certificate; the escalations only choose candidates). All decisions and
    the polygon in 200 bits; stage-B element sets checkpointed. Returns (elements, Dirichlet result)."""
    global IDCUT
    import os
    IDCUT = gmpy2.cosh(mpfr('5e-4'))
    S = {}
    for g in gens:
        for h in (g, mp_inv(g)):
            if mp_d(h)[0] > 1e-3: S[mp_key(h)] = h
    for it in range(2):
        items = sorted(S.values(), key=lambda e: mp_d(e)[0])[:K0]
        Id = np.array([mp_to_np(e) for e in items]); dK = mp_d(items[-1])[0]
        S = dict((mp_key(e), e) for e in items)
        for s0 in range(0, len(items), 50):
            d = halfplanes(su11(np.einsum('nij,mjk->nmik', Id[s0:s0 + 50], Id).reshape(-1, 2, 2)))[2]
            for k in np.where((d > 1e-3) & (d <= dK + 0.01))[0]:
                a_, b_ = divmod(int(k), len(items)); e = mp_mul(items[s0 + a_], items[b_])
                if abs(e[0]) >= IDCUT: S.setdefault(mp_key(e), e)
        items = sorted(S.values(), key=lambda e: mp_d(e)[0])[:K0]
        S = dict((mp_key(e), e) for e in items)
        log(f'U4 stage A{it}: kept {len(items)}, d range {mp_d(items[0])[0]:.4f} .. {mp_d(items[-1])[0]:.4f}')
    SA = list(S.values())
    if ckpt is not None and os.path.exists(ckpt):                     # E94: stage-B checkpoint (element set of the last round)
        S = dict((mp_key(e), e) for e in load_els(ckpt))
        log(f'U4: resumed from checkpoint {ckpt}: {len(S)} elements')
    res = None
    for it in range(rounds):
        items = list(S.values())
        res = dirichlet_mp(items)                                    # E94: polygon in 200 bits
        V = res['verts']; r = np.abs(V); ins = r < 1
        log(f'U4 stage B{it}: elements {len(items)}, faces {res["nfaces"]}, bounded {res["bounded"]}, area {res["area"]:.12f} '
            f'(diff {res["area"] - area_target:.3e}), circumradius {res["rD"]:.6f} | vertices outside the disk {int((~ins).sum())}, '
            f'max distance of the others {float(np.arctanh(r[ins]).max()) if ins.any() else np.nan:.4f}', flush=True)
        if res['bounded'] and abs(res['area'] - area_target) < 1e-8:
            break
        fidx = sorted(set(t for t in res['tag'] if t >= 0))
        F = [items[res['ok_index'][t]] for t in fidx]
        F = mp_dedupe(F + [mp_inv(e) for e in F])
        cut = [e for e in products_cut_screened(F, res) if mp_key(e) not in S]   # E94: double screen, 200-bit decision
        log(f'   face elements {len(F)}, pair products cutting the polygon {len(cut)}', flush=True)
        if not cut:
            cut = [e for e in triples_cut_screened(F, res, log=log) if mp_key(e) not in S]
        if not cut:
            cut = [e for e in mp_dedupe(products_AB_screened(F, SA, res) + products_AB_screened(SA, F, res)) if mp_key(e) not in S]
            log(f'   faces x stage-A products cutting {len(cut)}', flush=True)
        if not cut:
            break
        S = dict((mp_key(e), e) for e in F)
        for e in cut: S[mp_key(e)] = e
        if ckpt is not None: save_els(ckpt, list(S.values()))
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
    def fold(w, maxit=2000, chunk=40000):
        w = np.array(w, complex, copy=True)
        if len(w) > chunk:                                       # E94: chunked (857k points x 822 faces would need ~16 GB)
            outs = [fold(w[s0:s0 + chunk], maxit, chunk) for s0 in range(0, len(w), chunk)]
            return np.concatenate([o[0] for o in outs]), np.concatenate([o[1] for o in outs]), max(o[2] for o in outs)
        Fm = np.tile(np.eye(2, dtype=complex), (len(w), 1, 1))
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


def support(Vk):
    """support function of the convex polygon with Klein vertices Vk (CCW): returns h(theta) = max_j Re(V_j e^{-i theta}), by
    binary search on the outward edge-normal angles (the extreme vertex and its two neighbours are evaluated)."""
    m = len(Vk)
    ed = np.roll(Vk, -1) - Vk                                        # edge j: V_j -> V_{j+1}
    phi = np.unwrap(np.angle(-1j*ed))                                # outward normal (CCW polygon)
    def h(theta):
        th = np.mod(theta - phi[0], TAU) + phi[0]
        idx = np.searchsorted(phi, th, 'right')
        u = np.exp(-1j*theta)
        return np.max([(Vk[(idx + k) % m]*u).real for k in (-1, 0, 1)], 0)
    return h


def products_cut_screened(Fc, res, dmax=np.inf, margin=1e-6, chunk=100):
    """E94 (U5): products f g (f, g in Fc, 200 bits) whose bisector cuts the polygon of res: all products are screened in double
    (support function, margin 1e-6 >> the double error of products of elements with displacement <= ~20), the survivors are
    recomputed in 200 bits and decided by cutting_mp. Replaces E93's all-pairs 200-bit products (2.7e6 per round at genus 69)."""
    Fd = np.array([mp_to_np(e) for e in Fc])
    dF = np.array([mp_d(e)[0] for e in Fc])
    h = support(res['verts'])
    cand = {}
    for s0 in range(0, len(Fd), chunk):
        P = su11(np.einsum('nij,mjk->nmik', Fd[s0:s0 + chunk], Fd).reshape(-1, 2, 2))
        U, T, d = halfplanes(P)
        mg = np.maximum(margin, 1e-12*np.exp(((dF[s0:s0 + chunk, None] + dF[None, :]).ravel() - d)/2))   # >= 1e4 x the double error
        ok = (d > 1e-3) & (d <= dmax + 1e-6) & ~(h(np.angle(U)) - T <= -mg)
        nb = len(Fd)
        for k in np.where(ok)[0]:
            a_, b_ = divmod(int(k), nb)
            e = mp_mul(Fc[s0 + a_], Fc[b_])
            if abs(e[0]) < IDCUT: continue
            cand.setdefault(mp_key(e), e)
    C = list(cand.values())
    return [C[i] for i in cutting_mp(res, C)]


def products_AB_screened(A, B, res, dmax=np.inf, margin=1e-6, chunk=50):
    """products a b (a in A, b in B; 200 bits) whose bisector cuts the polygon of res: double screen (per-pair margin
    1e-12 e^((d_a + d_b - d)/2)), 200-bit decision (cutting_mp)."""
    Ad = np.array([mp_to_np(e) for e in A]); Bd = np.array([mp_to_np(e) for e in B])
    dA = np.array([mp_d(e)[0] for e in A]); dB = np.array([mp_d(e)[0] for e in B])
    h = support(res['verts']); cand = {}
    for s0 in range(0, len(A), chunk):
        P = su11(np.einsum('nij,mjk->nmik', Ad[s0:s0 + chunk], Bd).reshape(-1, 2, 2))
        U, T, d = halfplanes(P)
        mg = np.maximum(margin, 1e-12*np.exp(((dA[s0:s0 + chunk, None] + dB[None, :]).ravel() - d)/2))
        ok = (d > 1e-3) & (d <= dmax + 1e-6) & ~(h(np.angle(U)) - T <= -mg)
        for k in np.where(ok)[0]:
            a_, b_ = divmod(int(k), len(B))
            e = mp_mul(A[s0 + a_], B[b_])
            if abs(e[0]) < IDCUT: continue
            cand.setdefault(mp_key(e), e)
    C = list(cand.values())
    return [C[i] for i in cutting_mp(res, C)]


def triples_cut_screened(Fc, res, dmax=np.inf, binw=1e-4, per_bin=3, chunk=4, log=None):
    """E94 escalation when pair closure stalls (missing faces that only arise from each other under the length-3 cycle
    relations): triple products f g h (Fc, 200 bits) screened in double against the polygon of res (margin 1e-12 e^((d1+d2+d3-d)/2),
    >= 1e4 x the double error); of the survivors only the per_bin deepest (largest support violation) per direction bin of width
    binw are formed in 200 bits and decided by cutting_mp. A heuristic selection: completeness is certified by the area only."""
    Fd = np.array([mp_to_np(e) for e in Fc]); dF = np.array([mp_d(e)[0] for e in Fc]); n = len(Fc)
    P2 = su11(np.einsum('nij,mjk->nmik', Fd, Fd).reshape(-1, 2, 2)); ds2 = (dF[:, None] + dF[None, :]).ravel()
    h = support(res['verts'])
    TH, VI, II, KK = [], [], [], []
    for s0 in range(0, n, chunk):
        P3 = su11(np.einsum('nij,mjk->nmik', Fd[s0:s0 + chunk], P2).reshape(-1, 2, 2))
        U, T, d = halfplanes(P3)
        ds = (dF[s0:s0 + chunk, None] + ds2[None, :]).ravel()
        th = np.angle(U); vi = h(th) - T
        ok = (d > 1e-3) & (d <= dmax + 1e-6) & ~(vi <= -np.maximum(1e-9, 1e-12*np.exp((ds - d)/2)))
        k = np.where(ok)[0]
        TH.append(th[k]); VI.append(vi[k]); II.append(s0 + k // len(P2)); KK.append(k % len(P2))
    TH, VI, II, KK = (np.concatenate(a) for a in (TH, VI, II, KK))
    if not len(TH): return []
    b = np.floor((TH + np.pi)/binw).astype(np.int64)
    o = np.lexsort((-VI, b)); b = b[o]
    first = np.r_[0, np.where(np.diff(b) != 0)[0] + 1]
    rank = np.arange(len(b)) - np.repeat(first, np.diff(np.r_[first, len(b)]))
    sel = o[rank < per_bin]
    cand = {}
    for i, k in zip(II[sel], KK[sel]):
        a_, c_ = divmod(int(k), n)
        e = mp_mul(Fc[int(i)], mp_mul(Fc[a_], Fc[c_]))
        if abs(e[0]) < IDCUT: continue
        cand.setdefault(mp_key(e), e)
    C = list(cand.values())
    out = [C[i] for i in cutting_mp(res, C)]
    if log: log(f'   triple search: screened survivors {len(TH)}, selected {len(sel)}, distinct {len(C)}, cutting {len(out)}')
    return out


def _closure_at_mp(q, F0, area_target, dcap, rounds=20, log=None):
    """Dirichlet polygon about q from the 200-bit elements F0 (0-frame): seeded (as E93) by F u {f g : d(f g) <= dcap}, F = F0 u F0^-1
    conjugated to q = 0 (E94: products screened in double, formed in 200 bits; seeding by F alone stalled with ideal vertices
    left in the E93 test), then rounds of cutting products of face elements (screened in double, decided in 200 bits); certified
    when bounded with area 4 pi (g - 1) to 1e-8 (as U4)."""
    global IDCUT
    IDCUT = gmpy2.cosh(mpfr('5e-4'))
    M = mp_Mq(q); Mi = mp_inv(M)
    F = mp_dedupe([mp_conj(Mi, g, M) for g in F0] + [mp_conj(Mi, mp_inv(g), M) for g in F0])
    S = dict((mp_key(e), e) for e in F)
    Fd = np.array([mp_to_np(e) for e in F])                        # seed as E93: F u {f g : d(f g) <= dcap}; screened in double
    for s0 in range(0, len(F), 200):
        d = halfplanes(su11(np.einsum('nij,mjk->nmik', Fd[s0:s0 + 200], Fd).reshape(-1, 2, 2)))[2]
        for k in np.where((d > 1e-3) & (d <= dcap + 0.01))[0]:
            a_, b_ = divmod(int(k), len(F))
            h = mp_mul(F[s0 + a_], F[b_])
            if abs(h[0]) < IDCUT: continue
            if float(2*gmpy2.acosh(abs(h[0]))) <= dcap: S.setdefault(mp_key(h), h)
    L = list(S.values())
    res = dirichlet_mp(L)
    for it in range(rounds):
        if res['bounded'] and abs(res['area'] - area_target) < 1e-8:
            break
        Fc = [L[res['ok_index'][t]] for t in sorted(set(t for t in res['tag'] if t >= 0))]
        Fc = mp_dedupe(Fc + [mp_inv(e) for e in Fc])
        known = set(mp_key(e) for e in L)
        add = [e for e in products_cut_screened(Fc, res, dmax=2*dcap) if mp_key(e) not in known]
        if not add:
            add = [e for e in triples_cut_screened(Fc, res, dmax=2*dcap, log=log) if mp_key(e) not in known]
        if log: log(f'   mp closure at |q| = {float(hd0(np.array([q]))[0]):.3f} round {it}: faces {res["nfaces"]}, bounded {res["bounded"]}, '
                    f'area {res["area"]:.10f}, cutting products {len(add)}')
        if not add:
            break
        L = mp_dedupe(Fc + add)
        res = dirichlet_mp(L)
    res['Mp'] = Mq(q); res['Mp_mp'] = M
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


# ---------------- E94: Dirichlet polygons in 200 bits (projective Klein clip) ----------------
# E93's clip ran in double Klein coordinates, starting from a 16384-gon at Klein radius 1 - 1e-15 (inradius ~9.27). The genus-69
# polygon reaches further than that (E94 coarse U4 run 1 stalled with artificial edges left and no cutting product), and near
# the circle double Klein coordinates lose ~1e-8 (relative to 1 - |k| ~ 1e-9), too coarse for the 1e-8 area certificate.
# Here: the polygon is a cyclic list of faces (half-planes); vertex j = face j-1 meet face j, computed from the two face normals
# in 200 bits and stored as Klein-affine coordinates (t = 1). Start: the Euclidean square |Re k|, |Im k| <= 2 in the Klein chart
# (4 artificial faces, tag -1). Bounded = no artificial face left and every vertex strictly inside the unit disk. Area by
# Gauss-Bonnet from the face normals: (m - 2) pi - sum of interior angles, cos(angle) = -<n_{j-1}, n_j> (Minkowski, unit normals).
# A double-precision prefilter (vertex Klein coordinates rounded to double, margin 1e-9) skips half-planes that clearly miss.
HTOL = 1e-12          # a half-plane cuts a vertex if the vertex is outside by sinh(hyperbolic distance) > HTOL


def mp_normal(e):
    """unit outward normal (t, x, y) (Minkowski <X, Y> = -t t' + x x' + y y') of the bisector of 0 and e(0):
    D = {X : <X, n> <= 0}; e(0) = (|a|^2 + |b|^2, 2 a b) on the hyperboloid, n = (e(0) - 0)/|e(0) - 0| = (|b|, a b/|b|)."""
    a, b = e
    mb = abs(b); w = a*b/mb
    return (mb, w.real, w.imag)


_ART = [(mpfr(2), mpfr(1), mpfr(0)), (mpfr(2), mpfr(0), mpfr(1)), (mpfr(2), mpfr(-1), mpfr(0)), (mpfr(2), mpfr(0), mpfr(-1))]


def _meet(n1, n2):
    """Klein-affine point (kx, ky) of the intersection of the lines <X, n1> = 0, <X, n2> = 0 (X = J (n1 x n2), scaled to t = 1)."""
    c0 = n1[1]*n2[2] - n1[2]*n2[1]; c1 = n1[2]*n2[0] - n1[0]*n2[2]; c2 = n1[0]*n2[1] - n1[1]*n2[0]
    t = -c0
    if t == 0:
        raise RuntimeError('meet: lines parallel in the Klein chart')
    return (c1/t, c2/t)


def _sval(k, n):
    """<(1, kx, ky), n>."""
    return -n[0] + k[0]*n[1] + k[1]*n[2]


def _outside(k, n):
    r2 = 1 - k[0]**2 - k[1]**2
    s = _sval(k, n)
    if r2 > 0:
        return s > HTOL*gmpy2.sqrt(r2)            # sinh(signed distance) = s / sqrt(1 - |k|^2) for unit n
    return s > 1e-30*gmpy2.sqrt(n[1]**2 + n[2]**2)   # vertex outside the disk (unbounded polygon): Klein distance > 1e-30


def mp_clip(N, order, U=None, T=None):
    """N: list of unit normals; process half-planes in 'order'. Returns (faces, verts): faces = list of indices into N (-1-q for
    artificial face q), verts[j] = Klein-affine mp point of face j-1 meet face j."""
    faces = [-1, -2, -3, -4]
    nrm = lambda f: N[f] if f >= 0 else _ART[-1 - f]
    verts = [_meet(nrm(faces[j - 1]), nrm(faces[j])) for j in range(4)]
    kd = np.array([complex(float(v[0]), float(v[1])) for v in verts])
    for hi in order:
        if U is not None and ((kd*np.conj(U[hi])).real - T[hi]).max() < -1e-9:
            continue
        n = N[hi]
        if not any(_outside(v, n) for v in verts):                    # cut only if some vertex is outside beyond the tolerance
            continue
        out = [_sval(v, n) > 0 for v in verts]                       # E94 run 3: the clip itself uses exact 200-bit signs
        m = len(verts)
        if all(out):
            raise RuntimeError('mp_clip: half-plane excludes the whole polygon')
        st = [j for j in range(m) if out[j] and not out[j - 1]]
        if len(st) != 1:
            raise RuntimeError(f'mp_clip: outside run not contiguous ({len(st)} runs)')
        p = st[0]; q = p
        while out[(q + 1) % m]:
            q += 1
        q %= m
        fprev = faces[(p - 1) % m]; fq = faces[q]
        X1 = _meet(nrm(fprev), n); X2 = _meet(n, nrm(fq))
        # new cyclic lists starting at face q: faces q .. p-1 (kept), then hi; vertices: X2 (= hi meet q), V_{q+1} .. V_{p-1}, X1
        keepf = [faces[(q + i) % m] for i in range(((p - 1 - q) % m) + 1)]
        keepv = [verts[(q + i) % m] for i in range(1, ((p - 1 - q) % m) + 1)]
        faces = keepf + [hi]
        verts = [X2] + keepv + [X1]
        # degenerate edges (a kept face whose remaining edge is shorter than 1e-25 in the Klein chart: the new line passes through
        # its far vertex): drop the face and recompute the vertex
        while len(faces) > 3:
            m = len(faces); j = m - 1                                   # new face hi is last; its neighbours are m-2 and 0
            e_prev = abs(verts[j][0] - verts[j - 1][0]) + abs(verts[j][1] - verts[j - 1][1])     # edge of face m-2
            e_next = abs(verts[1][0] - verts[0][0]) + abs(verts[1][1] - verts[0][1])             # edge of face 0
            if e_prev < 1e-25:
                del faces[m - 2]; del verts[m - 2]
                verts[-1] = _meet(nrm(faces[-2]), nrm(faces[-1]))
            elif e_next < 1e-25:
                del faces[0]; del verts[0]
                verts[0] = _meet(nrm(faces[-1]), nrm(faces[0]))
            else:
                break
        kd = np.array([complex(float(v[0]), float(v[1])) for v in verts])
    return faces, verts


def _mink(n1, n2): return -n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2]


def dirichlet_mp(L):
    """Dirichlet-type polygon about 0 cut out by the bisectors of the 200-bit elements L (list of (a, b)); all geometry in 200 bits.
    Same fields as dirichlet_np (verts = Klein double, tag = index into the d > 1e-3 elements, face_g double), plus face_mp,
    normals, angles and area in 200 bits, and the minimum edge length."""
    cut = gmpy2.cosh(mpfr('5e-4'))
    okidx = [i for i, e in enumerate(L) if abs(e[0]) > cut]
    Lok = [L[i] for i in okidx]
    N = [mp_normal(e) for e in Lok]
    U, T, d = mp_halfplanes(Lok)
    faces, verts = mp_clip(N, np.argsort(T), U, T)
    m = len(faces)
    inside = [1 - v[0]**2 - v[1]**2 > 0 for v in verts]
    bounded = m > 2 and all(f >= 0 for f in faces) and all(inside)
    Vk = np.array([complex(float(v[0]), float(v[1])) for v in verts])
    out = dict(verts=Vk, tag=faces, bounded=bounded, nfaces=m, ok_index=np.array(okidx), area=np.inf, rD=np.inf, face_g=None,
               G=su11(np.array([mp_to_np(e) for e in Lok])) if Lok else None, face_mp=None, verts_mp=verts)
    if bounded:
        nf = [N[f] for f in faces]
        ang = [gmpy2.acos(-_mink(nf[j - 1], nf[j])) for j in range(m)]
        area = (m - 2)*gmpy2.const_pi() - sum(ang)
        rr = [gmpy2.sqrt(v[0]**2 + v[1]**2) for v in verts]
        X = [(1/gmpy2.sqrt(1 - r**2), v[0]/gmpy2.sqrt(1 - r**2), v[1]/gmpy2.sqrt(1 - r**2)) for v, r in zip(verts, rr)]
        el = [gmpy2.acosh(-_mink(X[j], X[(j + 1) % m])) for j in range(m)]
        out.update(area=float(area), area_mp=area, angles_mp=ang, normals_mp=nf, X_mp=X, rD=float(max(gmpy2.atanh(r) for r in rr)),
                   min_edge=float(min(el)), face_mp=[Lok[f] for f in faces], face_g=np.array([mp_to_np(Lok[f]) for f in faces]))
    return out


def cutting_mp(res, C, margin=1e-9, chunk=20000):
    """indices i of elements C[i] (200 bits) whose bisector half-plane cuts the polygon of res by more than HTOL (sinh of the
    hyperbolic distance) at some vertex; double prefilter with margin, decision in 200 bits."""
    if not len(C): return []
    U, T, d = mp_halfplanes(C)
    Vk = res['verts']; Vmp = res['verts_mp']
    cand = []
    for s0 in range(0, len(C), chunk):
        f = (Vk[None, :]*np.conj(U[s0:s0 + chunk, None])).real - T[s0:s0 + chunk, None]
        for i in np.where(f.max(1) > -margin)[0]:
            cand.append((s0 + int(i), np.where(f[i] > -margin)[0]))
    out = []
    for i, js in cand:
        n = mp_normal(C[i])
        if any(_outside(Vmp[j], n) for j in js):
            out.append(i)
    return out


def face_pairing_mp(res):
    """U4 in 200 bits: face a (element g) is carried by g^-1 onto its partner face b (endpoint mismatch = max hyperbolic distance
    of g^-1(V_a), g^-1(V_{a+1}) from V_{b+1}, V_b); the partner element must be g^-1 (relation error |g_a g_b - id| reported)."""
    X = res['X_mp']; F = res['face_mp']; m = len(F)
    W = [mpc(x[1], x[2])/(1 + x[0]) for x in X]                     # Poincare points
    Wd = np.array([complex(w) for w in W])
    def hd(z, w): return 2*gmpy2.asinh(abs(z - w)/gmpy2.sqrt((1 - abs(z)**2)*(1 - abs(w)**2)))
    def ap(g, w): return (g[0]*w + g[1])/(g[1].conjugate()*w + g[0].conjugate())
    perr = mpfr(0); relerr = mpfr(0); unmatched = 0; pair_of = []
    for a in range(m):
        gi = mp_inv(F[a]); P2, Q2 = ap(gi, W[a]), ap(gi, W[(a + 1) % m])
        dd = np.abs(Wd - complex(Q2)) + np.abs(np.roll(Wd, -1) - complex(P2))
        best = None
        for b in np.argsort(dd)[:3]:
            b = int(b); e = max(hd(P2, W[(b + 1) % m]), hd(Q2, W[b]))
            if best is None or e < best[0]: best = (e, b)
        e, b = best; pair_of.append(b); perr = max(perr, e)
        h = mp_mul(F[a], F[b]); r = max(abs(h[1]), abs(abs(h[0]) - 1))
        relerr = max(relerr, r)
        if r > 1e-6: unmatched += 1
    return float(perr), unmatched, np.array(pair_of), float(relerr)


def vertex_cycles_mp(res, pair_of):
    """vertex cycles (as vertex_cycles) with the interior angles in 200 bits (angle at vertex a = between faces a-1 and a)."""
    m = res['nfaces']; par = list(range(m))
    def fnd(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    for a in range(m):
        b = int(pair_of[a])
        for p_, q_ in ((a, (b + 1) % m), ((a + 1) % m, b)):
            rp, rq = fnd(p_), fnd(q_)
            if rp != rq: par[rp] = rq
    cls = {}
    for a in range(m): cls.setdefault(fnd(a), []).append(a)
    ang = res['angles_mp']
    return cls, [sum(ang[a] for a in c) for c in cls.values()]
