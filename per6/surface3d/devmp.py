"""E92-U: the development step of develop.py in 30-digit arithmetic (mpmath), for `python3 develop.py mp`.
Same BFS order and placement rule as develop.develop(); per-face lengths and angles come from the double-precision
solve (relative error ~1e-16, applied in each face's local frame), while every position, Mobius map, side-pairing fit
and vertex-cycle product is computed with 30 digits. Returns positions and generators in double after moving the base
point (centroid of the root face) to 0."""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import breadth_first_order
import mpmath as mp
mp.mp.dps = 30
ONE = mp.mpf(1)
def hd_mp(z, w): return 2*mp.asinh(abs(z - w)/mp.sqrt((1 - abs(z)**2)*(1 - abs(w)**2)))
def M_to(p): return mp.matrix([[1, p], [mp.conj(p), 1]])/mp.sqrt(1 - abs(p)**2)          # 0 -> p
def M_from(p): return mp.matrix([[1, -p], [-mp.conj(p), 1]])/mp.sqrt(1 - abs(p)**2)       # p -> 0
def apm(M, z): return (M[0, 0]*z + M[0, 1])/(M[1, 0]*z + M[1, 1])
def iso_two_mp(p1, q1, p2, q2):
    a = apm(M_from(p1), q1); b = apm(M_from(p2), q2)
    rho = (b/abs(b))/(a/abs(a)); s = mp.sqrt(rho)
    return M_to(p2)*mp.matrix([[s, 0], [0, 1/s]])*M_from(p1)
def mnorm(g):
    d = g[0, 0]*g[1, 1] - g[0, 1]*g[1, 0]; return g/mp.sqrt(d)
def ident_err_mp(g):
    g = mnorm(g); I = mp.eye(2)
    e1 = max(abs(g[i, j] - I[i, j]) for i in range(2) for j in range(2)); e2 = max(abs(g[i, j] + I[i, j]) for i in range(2) for j in range(2))
    return float(min(e1, e2))
def to_np(g): return np.array([[complex(g[0, 0]), complex(g[0, 1])], [complex(g[1, 0]), complex(g[1, 1])]])
def bfs(Tn):
    nF = len(Tn)
    En = np.sort(np.concatenate([Tn[:, [1, 2]], Tn[:, [2, 0]], Tn[:, [0, 1]]]), axis=1)
    eu, inv, cnt = np.unique(En, axis=0, return_inverse=True, return_counts=True)
    fe = inv.reshape(3, nF).T
    ef = -np.ones((len(eu), 2), int); k_ = np.zeros(len(eu), int)
    for j in range(3):
        for f, k in zip(range(nF), fe[:, j]): ef[k, k_[k]] = f; k_[k] += 1
    intr = cnt == 2
    Dg = sp.coo_matrix((np.ones(intr.sum()), (ef[intr, 0], ef[intr, 1])), shape=(nF, nF)); Dg = (Dg + Dg.T).tocsr()
    bf = np.unique(ef[~intr, 0])
    Db = sp.bmat([[Dg, sp.coo_matrix((np.ones(len(bf)), (bf, np.zeros(len(bf), int))), shape=(nF, 1))],
                  [sp.coo_matrix((np.ones(len(bf)), (np.zeros(len(bf), int), bf)), shape=(1, nF)), None]]).tocsr()
    order_b, pred_b = breadth_first_order(Db, nF, directed=False, return_predecessors=True)
    depth = np.zeros(nF + 1, int)
    for v in order_b[1:]: depth[v] = depth[pred_b[v]] + 1
    root = int(np.argmax(depth[:nF]))
    order, _ = breadth_first_order(Dg, root, directed=False, return_predecessors=True)
    return root, order
def mp_stage(Tn, lface, Aface, cyc, sides, pairs, orig, log=print):
    root, order = bfs(Tn)
    nN = int(Tn.max()) + 1
    th = np.tanh(lface/2); eA = np.exp(1j*Aface)
    z = [None]*nN; err = 0.0; nre = 0
    i, j, k = Tn[root]
    z[i] = mp.mpc(0); z[j] = mp.mpc(mp.tanh(mp.mpf(lface[root, 2])/2))
    z[k] = mp.tanh(mp.mpf(lface[root, 1])/2)*mp.expj(mp.mpf(Aface[root, 0]))
    for n_, f in enumerate(order[1:]):
        c = Tn[f]
        for r in range(3):
            p, q, w = c[r], c[(r + 1) % 3], c[(r + 2) % 3]
            if z[p] is not None and z[q] is not None: break
        zp, zq = z[p], z[q]; cp = mp.conj(zp)
        qq = (zq - zp)/(ONE - cp*zq)
        wn = (qq/abs(qq))*mp.mpc(complex(th[f, (r + 1) % 3]*eA[f, r]))
        zw = (wn + zp)/(ONE + cp*wn)
        if z[w] is None: z[w] = zw
        else:
            nre += 1; e = float(hd_mp(z[w], zw))
            if e > err: err = e
        if n_ % 200000 == 0: log(f'   mp development: {n_} of {len(order) - 1} faces')
    log(f'develop (30 digits): root face {root} | placed {sum(v is not None for v in z)} of {nN} | re-placements {nre}, '
        f'max re-placement error (hyperbolic) {err:.3e}')
    gens, res, sv_bad = [], 0.0, 0
    for s1, s2 in pairs:
        i1 = sides[s1]['idx']; i2 = sides[s2]['idx'][::-1]
        sv_bad += int(np.sum(orig[cyc[i1]] != orig[cyc[i2]]))
        Z1 = [z[cyc[t]] for t in i1]; Z2 = [z[cyc[t]] for t in i2]
        g = iso_two_mp(Z1[0], Z1[-1], Z2[0], Z2[-1])
        res = max(res, max(float(hd_mp(apm(g, a), b)) for a, b in zip(Z1, Z2)))
        gens.append((s1, s2, g))
    log(f'U3 (30 digits): side pairings {len(gens)} | max hyperbolic residual over all side points {res:.3e} | vertex mismatches {sv_bad}')
    gmap = {}
    for s1, s2, g in gens: gmap[s1] = (s2, g); gmap[s2] = (s1, g**-1)
    nS = len(sides); used = set(); cyc_err = 0.0; ncyc = 0; lens = []
    for s0 in range(nS):
        if s0 in used: continue
        g = mp.eye(2); s = s0; L = 0
        while True:
            used.add(s); partner, gs = gmap[s]; g = gs*g; L += 1
            s = (partner + 1) % nS
            if s == s0 or L > 400: break
        cyc_err = max(cyc_err, ident_err_mp(g)); ncyc += 1; lens.append(L)
    log(f'U3 (30 digits): vertex cycles {ncyc} | lengths {sorted(lens)} | max distance of products from +-I {cyc_err:.3e}')
    p0 = (z[Tn[root][0]] + z[Tn[root][1]] + z[Tn[root][2]])/3
    Mb = M_from(p0); Mbi = M_to(p0)
    zc = np.array([complex(apm(Mb, v)) for v in z])
    G = [to_np(mnorm(Mb*g*Mbi)) for _, _, g in gens]
    return zc, G, root, dict(dev_err=err, U3_res=res, U3_cycle_err=cyc_err, cycles=ncyc)
