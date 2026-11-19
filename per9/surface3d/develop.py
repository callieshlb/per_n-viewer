#!/usr/bin/env python3
"""E93-U step 2 (E92 develop.py; Dirichlet stage replaced for genus 24, see fuchs.py): develop the uniformized surface into the Poincare disk and extract its Fuchsian group.
Run in a mesh directory holding surface_mesh.json, sphere_mesh.json, uniform_u.npy (from uniformize.py).
- Cut open: tree-cotree (primal BFS from vertex 0, dual BFS), loops of the 2g leftover edges, dangling branches pruned;
  corners joined across uncut edges by connected components (vectorized).
- Develop: BFS over the disk faces from a central face, each new vertex placed from two hyperbolic lengths.
  Re-placements of an already placed vertex give the development error.
- Side pairings: boundary runs between cut-graph branch vertices, paired by their surface edges; each pairing is the
  isometry fixed by the run's two end points, checked on every point of the run (U3); vertex cycles multiply to +-I (U3).
- Dirichlet polygon about a generic base point moved to 0: half-planes in the Klein model from group elements found by
  BFS over words; area must equal 20 pi (U4); faces must close up in pairs (U4).
Writes develop.npz (positions, generators, Dirichlet data) and prints raw values."""
import json
import os
import sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import breadth_first_order, connected_components
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from uniformize import load_mesh, setup, residual, angles
TAU = 2*np.pi
# ---------------- Mobius helpers (2x2 complex, SU(1,1) normalized) ----------------
def mob_ap(g, z): return (g[0, 0]*z + g[0, 1])/(g[1, 0]*z + g[1, 1])
def mob_norm(g): return g/np.sqrt(np.linalg.det(g))
def adjinv(g):
    """E95: exact inverse of a determinant-1 2x2 matrix (adjugate). np.linalg.inv fails for Per_9 generators: with displacement d,
    |a| ~ cosh(d/2) and |a|^2 - |b|^2 = 1 cancels to 0 in double (coarse run 1: 'Singular matrix')."""
    return np.array([[g[1, 1], -g[0, 1]], [-g[1, 0], g[0, 0]]])
def Mq(q): return mob_norm(np.array([[1, q], [np.conj(q), 1]], complex))               # 0 -> q
def rotm(ph): return np.array([[np.exp(.5j*ph), 0], [0, np.exp(-.5j*ph)]])
def hdist(z, w): return 2*np.arcsinh(np.abs(z - w)/np.sqrt((1 - np.abs(z)**2)*(1 - np.abs(w)**2)))
def iso_two(p1, q1, p2, q2):
    """orientation-preserving isometry with p1 -> p2 and q1 -> (direction of) q2."""
    a = mob_ap(np.linalg.inv(Mq(p1)), q1); b = mob_ap(np.linalg.inv(Mq(p2)), q2)
    return mob_norm(Mq(p2) @ rotm(np.angle(b) - np.angle(a)) @ np.linalg.inv(Mq(p1)))
def ident_err(g):
    g = mob_norm(g); return min(np.abs(g - np.eye(2)).max(), np.abs(g + np.eye(2)).max())
# ---------------- cut open ----------------
def cut_open(T, nV):
    nF = len(T)
    Es = np.sort(np.concatenate([T[:, [1, 2]], T[:, [2, 0]], T[:, [0, 1]]]), axis=1)
    Eu, inv = np.unique(Es, axis=0, return_inverse=True); nE = len(Eu)
    fe = inv.reshape(3, nF).T
    ef = np.zeros((nE, 2), int); cnt = np.zeros(nE, int)
    for j in range(3):
        for f, k in zip(range(nF), fe[:, j]):
            ef[k, cnt[k]] = f; cnt[k] += 1
    A = sp.coo_matrix((np.ones(nE), (Eu[:, 0], Eu[:, 1])), shape=(nV, nV)); A = (A + A.T).tocsr()
    _, pred = breadth_first_order(A, 0, directed=False, return_predecessors=True)
    ekey = {(int(a), int(b)): k for k, (a, b) in enumerate(Eu)}
    tree = np.zeros(nE, bool)
    for v in range(nV):
        p = pred[v]
        if p >= 0: tree[ekey[(min(v, p), max(v, p))]] = True
    nt = np.where(~tree)[0]
    D = sp.coo_matrix((np.arange(1, len(nt) + 1), (ef[nt, 0], ef[nt, 1])), shape=(nF, nF)).tocsr(); D = D + D.T
    _, dpred = breadth_first_order(D, 0, directed=False, return_predecessors=True)
    cot = np.zeros(nE, bool)
    pair_edge = {}
    for k in nt: pair_edge[(ef[k, 0], ef[k, 1])] = k; pair_edge[(ef[k, 1], ef[k, 0])] = k
    for f in range(nF):
        p = dpred[f]
        if p >= 0: cot[pair_edge[(f, p)]] = True
    left = np.where(~tree & ~cot)[0]
    cut = set(left.tolist())
    def to_root(v):
        out = []
        while pred[v] >= 0:
            p = pred[v]; out.append(ekey[(min(v, p), max(v, p))]); v = p
        return out
    for k in left:
        a, b = Eu[k]; cut |= set(to_root(a)) | set(to_root(b))
    deg = np.zeros(nV, int)
    for k in cut: deg[Eu[k]] += 1
    changed = True
    while changed:
        changed = False
        for k in list(cut):
            a, b = Eu[k]
            if deg[a] == 1 or deg[b] == 1: cut.discard(k); deg[a] -= 1; deg[b] -= 1; changed = True
    iscut = np.zeros(nE, bool); iscut[list(cut)] = True
    cdeg = np.zeros(nV, int)
    for k in cut: cdeg[Eu[k]] += 1
    # corners (f, j) -> id 3f + j; join across uncut edges
    pos_in = lambda f, v: 3*f + np.argmax(T[f] == v[:, None], axis=1)
    un = np.where(~iscut)[0]
    f1, f2 = ef[un, 0], ef[un, 1]
    r, c = [], []
    for col in (0, 1):
        v = Eu[un, col]
        r.append(3*f1 + np.argmax(T[f1] == v[:, None], axis=1)); c.append(3*f2 + np.argmax(T[f2] == v[:, None], axis=1))
    G = sp.coo_matrix((np.ones(2*len(un)), (np.concatenate(r), np.concatenate(c))), shape=(3*nF, 3*nF))
    ncomp, lab = connected_components(G, directed=False)
    Tn = lab.reshape(nF, 3)
    orig = np.zeros(ncomp, int); orig[lab] = T.reshape(-1)
    info = dict(nE=nE, left=len(left), cut_edges=int(iscut.sum()), branch_nodes=int((cdeg >= 3).sum()), disk_vertices=int(ncomp))
    return Tn, orig, cdeg, info, fe
def boundary(Tn, orig, cdeg):
    nF = len(Tn)
    En = np.sort(np.concatenate([Tn[:, [0, 1]], Tn[:, [1, 2]], Tn[:, [2, 0]]]), axis=1)
    eu, cnt = np.unique(En, axis=0, return_counts=True)
    bset = {(int(a), int(b)) for a, b in eu[cnt == 1]}
    succ = {}
    for j in range(3):
        a, b = Tn[:, j], Tn[:, (j + 1) % 3]
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        for f in range(nF):
            if (int(lo[f]), int(hi[f])) in bset: succ[int(a[f])] = int(b[f])
    start = next(iter(succ)); cyc = [start]
    while True:
        w = succ[cyc[-1]]
        if w == start: break
        cyc.append(w)
    cyc = np.array(cyc); cv = orig[cyc]; nb = len(cyc)
    mark = [i for i, v in enumerate(cv) if cdeg[v] >= 3 or v == 0]
    sides = []
    for s in range(len(mark)):
        i0, i1 = mark[s], mark[(s + 1) % len(mark)]
        idx = list(range(i0, i1 + 1)) if i1 > i0 else list(range(i0, nb)) + list(range(0, i1 + 1))
        es = frozenset((min(int(cv[a]), int(cv[b])), max(int(cv[a]), int(cv[b]))) for a, b in zip(idx[:-1], idx[1:]))
        sides.append(dict(idx=idx, es=es))
    pairs = {}
    for si, s in enumerate(sides): pairs.setdefault(s['es'], []).append(si)
    return cyc, sides, list(pairs.values()), len(bset), chi_disk(Tn, eu)
def chi_disk(Tn, eu): return int(Tn.max() + 1 - len(eu) + len(Tn))
# ---------------- develop ----------------
def develop(Tn, lface, Aface):
    """lface[f, j] = hyperbolic length of the side opposite corner j; Aface[f, j] = angle at corner j."""
    nF = len(Tn); nN = int(Tn.max()) + 1
    En = np.sort(np.concatenate([Tn[:, [1, 2]], Tn[:, [2, 0]], Tn[:, [0, 1]]]), axis=1)
    eu, inv, cnt = np.unique(En, axis=0, return_inverse=True, return_counts=True)
    fe = inv.reshape(3, nF).T
    ef = -np.ones((len(eu), 2), int); k_ = np.zeros(len(eu), int)
    for j in range(3):
        for f, k in zip(range(nF), fe[:, j]): ef[k, k_[k]] = f; k_[k] += 1
    intr = cnt == 2
    Dg = sp.coo_matrix((np.ones(intr.sum()), (ef[intr, 0], ef[intr, 1])), shape=(nF, nF)); Dg = (Dg + Dg.T).tocsr()
    bf = np.unique(ef[~intr, 0])                                    # boundary faces
    # central face = farthest (dual BFS) from the boundary
    src = nF
    Db = sp.bmat([[Dg, sp.coo_matrix((np.ones(len(bf)), (bf, np.zeros(len(bf), int))), shape=(nF, 1))],
                  [sp.coo_matrix((np.ones(len(bf)), (np.zeros(len(bf), int), bf)), shape=(1, nF)), None]]).tocsr()
    order_b, pred_b = breadth_first_order(Db, src, directed=False, return_predecessors=True)
    depth = np.zeros(nF + 1, int)
    for v in order_b[1:]: depth[v] = depth[pred_b[v]] + 1
    root = int(np.argmax(depth[:nF]))
    order, pred = breadth_first_order(Dg, root, directed=False, return_predecessors=True)
    z = np.full(nN, np.nan + 0j); err = 0.0
    i, j, k = Tn[root]
    z[i] = 0; z[j] = np.tanh(lface[root, 2]/2); z[k] = np.tanh(lface[root, 1]/2)*np.exp(1j*Aface[root, 0])
    for f in order[1:]:
        c = Tn[f]
        known = ~np.isnan(z[c])
        # rotate so that corners 0, 1 are placed (the parent shares an edge)
        for r in range(3):
            cc = np.roll(c, -r)
            if not np.isnan(z[cc[0]]) and not np.isnan(z[cc[1]]): break
        p, q, w = cc; jr = (2 + r) % 3                              # corner index (in f) of w; side opposite corner of p is ...
        jp = r % 3
        lpw = lface[f, (jp + 1) % 3]                                 # side (p, w) is opposite corner q = (jp+1)%3
        ap_ = Aface[f, jp]
        M = np.linalg.inv(Mq(z[p])); qq = mob_ap(M, z[q]); th = np.angle(qq)
        wn = np.tanh(lpw/2)*np.exp(1j*(th + ap_))
        zw = mob_ap(Mq(z[p]), wn)
        if np.isnan(z[w]): z[w] = zw
        else: err = max(err, float(hdist(z[w], zw)))
    return z, root, err
# ---------------- Dirichlet polygon ----------------
def to_klein(z): return 2*z/(1 + np.abs(z)**2)
def to_poinc(k): return k/(1 + np.sqrt(np.maximum(0, 1 - np.abs(k)**2)))
def gkey(g):
    g = mob_norm(g); f = g.reshape(-1)
    k = next(i for i in range(4) if abs(f[i]) > 1e-9)
    if f[k].real < 0 or (f[k].real == 0 and f[k].imag < 0): f = -f
    return tuple(np.round(np.concatenate([f.real, f.imag])*1e6).astype(np.int64))
def hd0(w): return 2*np.arctanh(np.minimum(np.abs(w), 1 - 1e-16))
def group_elements(gens, bpts, R, maxel=40000):
    """all g with g(Delta) meeting B(0, R) (Delta = developed disk containing 0; tested on its sampled boundary),
    found by BFS over adjacent tiles g -> g s."""
    I2 = np.eye(2, dtype=complex)
    els = [I2]; seen = {gkey(I2)}; i = 0
    while i < len(els) and len(els) < maxel:
        g = els[i]; i += 1
        for s in gens:
            h = mob_norm(g @ s); k = gkey(h)
            if k in seen: continue
            seen.add(k)
            if hd0(mob_ap(h, bpts)).min() <= R: els.append(h)
    return els
def words_bfs(S, bound, keep, maxel=400000):
    """all products of elements of S whose orbit point stays within `bound` (BFS); returns those within `keep`."""
    I2 = np.eye(2, dtype=complex)
    els = [I2]; seen = {gkey(I2)}; i = 0
    while i < len(els) and len(els) < maxel:
        g = els[i]; i += 1
        for s in S:
            h = mob_norm(g @ s); k = gkey(h)
            if k in seen: continue
            seen.add(k)
            if hd0(mob_ap(h, 0)) <= bound: els.append(h)
    return [g for g in els if hd0(mob_ap(g, 0)) <= keep], len(els), len(els) >= maxel
def vnorm(P):
    """normalize to det 1 and project onto SU(1,1) exactly (a, b; conj b, conj a), killing round-off drift."""
    d = P[:, 0, 0]*P[:, 1, 1] - P[:, 0, 1]*P[:, 1, 0]; P = P/np.sqrt(d)[:, None, None]
    a = (P[:, 0, 0] + np.conj(P[:, 1, 1]))/2; b = (P[:, 0, 1] + np.conj(P[:, 1, 0]))/2
    s = np.sqrt(np.abs(a)**2 - np.abs(b)**2); a, b = a/s, b/s
    Q = np.empty_like(P); Q[:, 0, 0] = a; Q[:, 0, 1] = b; Q[:, 1, 0] = np.conj(b); Q[:, 1, 1] = np.conj(a); return Q
def vprod(A, B):
    return vnorm(np.einsum('nij,mjk->nmik', A, B).reshape(-1, 2, 2))
def vorb(P): return P[:, 0, 1]/P[:, 1, 1]
def vkeys(z): return np.stack([np.round(z.real*1e8), np.round(z.imag*1e8)], 1).astype(np.int64)
class HSet:
    """set of orbit points up to hyperbolic distance ~EPS (cells in (d, theta*sinh d) of size EPS, 3x3 neighbourhood
    checked, theta wrap handled). Distinct orbit points of a surface group are >= 2*injectivity radius apart (>> EPS)."""
    EPS = 1e-4
    def __init__(self): self.c = set()
    @staticmethod
    def _cell(z):
        d = float(hd0(z)); th = float(np.angle(z)); s = np.sinh(d)
        return int(np.floor(d/HSet.EPS)), th, s
    def _keys(self, z):
        kd, th, s = self._cell(z); out = []
        ths = [th] + ([th - 2*np.pi*np.sign(th)] if abs(th) > np.pi - 1e-3 else [])
        for t in ths:
            kt = int(np.floor(t*s/HSet.EPS))
            out += [(kd + a, kt + b) for a in (-1, 0, 1) for b in (-1, 0, 1)]
        return out, (kd, int(np.floor(th*s/HSet.EPS)))
    def add_if_new(self, z):
        ks, own = self._keys(z)
        if any(k in self.c for k in ks): return False
        self.c.add(own); return True
def vuniq(P):
    hs = HSet(); keep = [i for i, z in enumerate(vorb(P)) if hs.add_if_new(z)]; return P[keep]
def shorten(S0, K=1500, iters=12, log=print):
    """keep the K shortest distinct elements of S u S.S, repeatedly (all are group elements)."""
    S = vuniq(np.array(S0))
    for it in range(iters):
        P = np.concatenate([S, vprod(S, S)]); P = vuniq(P)
        d = hd0(vorb(P)); o = np.argsort(d); P, d = P[o], d[o]
        P, d = P[d > 1e-3], d[d > 1e-3]
        S = P[:K]
        verts, ftag, H = dirichlet(list(S[:400]), 99.0, pre=True)
        bounded = all(t >= 0 for t in ftag)
        area = poly_area(verts)[0] if bounded else np.inf
        rN = float(hd0(to_poinc(verts)).max()) if bounded else np.inf
        log(f'U4 stage A{it}: kept {len(S)}, shortest {d[0]:.4f}, 400th {d[min(399, len(d)-1)]:.4f}, faces {len(ftag)}, bounded {bounded}, '
            f'area {area:.10f}, circumradius {rN:.4f}')
        if bounded and abs(area - 20*np.pi) < 1e-6: return S, rN
    return S, rN
def words_bfs_vec(S, bound, keep, maxel=2000000):
    S = np.array(S); I2 = np.eye(2, dtype=complex)[None]
    hs = HSet(); hs.add_if_new(0j); allE = [I2]; front = I2; n = 1
    while len(front) and n < maxel:
        P = vprod(front, S); z = vorb(P); m = hd0(z) <= bound; P, z = P[m], z[m]
        new = np.array([hs.add_if_new(zz) for zz in z], bool)
        P = P[new]
        allE.append(P); front = P; n += len(P)
    E_ = np.concatenate(allE); d = hd0(vorb(E_))
    global LAST_EXPLORED; LAST_EXPLORED = (E_, bound)
    return list(E_[d <= keep]), n, n >= maxel
def dirichlet(els, R, pre=False):
    H = []
    for g in (els if pre else els[1:]):
        z0 = mob_ap(g, 0); d = hd0(z0)
        if d < 1e-3: continue
        if d <= R: H.append((np.exp(1j*np.angle(z0)), np.tanh(d/2), g))
    H.sort(key=lambda h: h[1])
    poly = [np.exp(1j*t)*0.999999 for t in np.linspace(0, TAU, 513)[:-1]]
    tag = [-1]*len(poly)
    for hi, (u, t, g) in enumerate(H):
        new, ntag = [], []
        n = len(poly)
        for a in range(n):
            P, Q = poly[a], poly[(a + 1) % n]; tP = tag[a]
            fP, fQ = (P*np.conj(u)).real - t, (Q*np.conj(u)).real - t
            if fP <= 0: new.append(P); ntag.append(tP)
            if (fP <= 0) != (fQ <= 0):
                X = P + (Q - P)*(fP/(fP - fQ)); new.append(X); ntag.append(tP if fP > 0 else hi)
        poly, tag = new, ntag
    V = np.array(poly); n = len(V)
    keep = [a for a in range(n) if tag[a] != tag[a - 1]]
    return V[keep], [tag[a] for a in keep], H
def poly_area(verts_k):
    """hyperbolic area of a convex geodesic polygon given by Klein vertices (Gauss-Bonnet)."""
    zp = to_poinc(verts_k); n = len(zp); s = 0.0
    for a in range(n):
        M = np.linalg.inv(Mq(zp[a])); u, w = mob_ap(M, zp[a - 1]), mob_ap(M, zp[(a + 1) % n])
        s += abs(np.angle(w/u))
    return (n - 2)*np.pi - s, zp
if __name__ == '__main__':
    T, SV, Y, P, E, l0 = load_mesh(); nV = len(SV)
    yt = Y[SV[T, 0]]; sel = np.abs(yt).max(1) < 1.5
    sg = np.sign(((yt[sel, 1] - yt[sel, 0]).conjugate()*(yt[sel, 2] - yt[sel, 0])).imag)
    print('orientation in the y-plane (|y| < 1.5): counterclockwise', int((sg > 0).sum()), 'clockwise', int((sg < 0).sum()))
    if (sg < 0).sum() > (sg > 0).sum(): T = T[:, [0, 2, 1]]; print('flipped mesh orientation to counterclockwise in y')
    u = np.load('uniform_u.npy')
    fe0 = setup(T)
    F, l, tri = residual(u, T, fe0, l0, E, nV)
    a, b, c, A, B, C, ok = tri
    print('U1 recheck: max |angle sum - 2pi|', f'{np.abs(F).max():.3e}', '| degenerate', int((~ok).sum()))
    lface = np.stack([a, b, c], 1); Aface = np.stack([A, B, C], 1)
    Tn, orig, cdeg, info, _ = cut_open(T, nV)
    print('cut open:', info)
    cyc, sides, pairs, nbe, chid = boundary(Tn, orig, cdeg)
    print('disk chi', chid, '| boundary edges', nbe, '| sides', len(sides), '| pairs', len(pairs), '| unpaired', sum(len(p) != 2 for p in pairs))
    if 'mp2' in sys.argv:                 # E93: 160-bit development from the refined metric (refine.py)
        import gmpy2
        from devmp2 import mp_stage2
        uh, ul = np.load('uniform_u_hi.npy'), np.load('uniform_u_lo.npy')
        u_mp = [gmpy2.mpfr(float(a)) + gmpy2.mpfr(float(b)) for a, b in zip(uh, ul)]
        zc, G, root, mpinfo = mp_stage2(Tn, T, E, l0, u_mp, cyc, sides, pairs, orig)
        np.savez('develop_mp2_raw.npz', zc=zc, G=np.array(G), root=root, gens_str=np.array(json.dumps(mpinfo.get('gens_str'))))  # E95: save before anything else can fail
        _ad = [abs(complex(g[0, 0])) for g in G]; print(f'generators: {len(G)} | max |a| {max(_ad):.3e} (displacement {2*np.arccosh(max(_ad)):.3f})', flush=True)
        G2 = G + [adjinv(g) for g in G]
        derr, res, ncyc, cyc_err = mpinfo['dev_err'], mpinfo['U3_res'], mpinfo['cycles'], mpinfo['U3_cycle_err']
    elif 'mp' in sys.argv:
        from devmp import mp_stage
        zc, G, root, mpinfo = mp_stage(Tn, lface, Aface, cyc, sides, pairs, orig)
        G2 = G + [adjinv(g) for g in G]
        derr, res, ncyc, cyc_err = mpinfo['dev_err'], mpinfo['U3_res'], mpinfo['cycles'], mpinfo['U3_cycle_err']
    else:
        z, root, derr = develop(Tn, lface, Aface)
        print('develop: root face', root, '| placed', int(np.isfinite(z).sum()), 'of', len(z), '| max re-placement error (hyperbolic)', f'{derr:.3e}',
              '| max distance from root', f'{float(np.nanmax(hdist(0, z))):.4f}')
        gens, res, sv_bad = [], 0.0, 0
        for s1, s2 in pairs:
            i1 = sides[s1]['idx']; i2 = sides[s2]['idx'][::-1]
            sv_bad += int(np.sum(orig[cyc[i1]] != orig[cyc[i2]]))
            Z1, Z2 = z[cyc[i1]], z[cyc[i2]]
            g = iso_two(Z1[0], Z1[-1], Z2[0], Z2[-1])
            res = max(res, float(np.max(hdist(mob_ap(g, Z1), Z2))))
            gens.append((s1, s2, g))
        print('U3: side pairings', len(gens), '| max hyperbolic residual over all side points', f'{res:.3e}', '| vertex mismatches', sv_bad)
        gmap = {}
        for s1, s2, g in gens: gmap[s1] = (s2, g); gmap[s2] = (s1, np.linalg.inv(g))
        nS = len(sides); used = set(); cyc_err = 0.0; ncyc = 0; lens = []
        for s0 in range(nS):
            if s0 in used: continue
            g = np.eye(2, dtype=complex); s = s0; L = 0
            while True:
                used.add(s); partner, gs = gmap[s]; g = gs @ g; L += 1
                s = (partner + 1) % nS
                if s == s0 or L > 400: break
            cyc_err = max(cyc_err, ident_err(g)); ncyc += 1; lens.append(L)
        print('U3: vertex cycles', ncyc, '| lengths', sorted(lens), '| max distance of products from +-I', f'{cyc_err:.3e}')
        # base point: centroid of the root face, moved to 0
        p0 = z[Tn[root]].mean()
        Mb = np.linalg.inv(Mq(p0))
        zc = mob_ap(Mb, z)
        G = [mob_norm(Mb @ g @ np.linalg.inv(Mb)) for _, _, g in gens]
        G2 = G + [np.linalg.inv(g) for g in G]
    GEN = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'branch.json')))['genus']
    AREA = 4*np.pi*(GEN - 1)
    gens_str = mpinfo.get('gens_str') if ('mp' in sys.argv or 'mp2' in sys.argv) else None
    np.savez('develop_pre.npz', zc=zc, G=np.array(G), cyc=cyc, Tn=Tn, orig=orig, root=root,
             gens_str=np.array(json.dumps(gens_str)))
    print('base point moved to 0; developed disk: max distance from 0', f'{float(hd0(zc[np.isfinite(zc)]).max()):.4f}',
          '| max generator displacement', f'{max(hd0(mob_ap(g, 0)) for g in G):.4f}', flush=True)
    import dirichlet_only
    dirichlet_only.run(AREA, json.load(open('develop_log.json')) if os.path.exists('develop_log.json') else {},
                       dict(U1=float(np.abs(F).max()), dev_err=derr, U3_res=res, U3_cycles=ncyc, U3_cycle_err=cyc_err))
