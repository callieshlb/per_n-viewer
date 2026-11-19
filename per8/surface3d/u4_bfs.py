#!/usr/bin/env python3
"""E94-U4 seed (replaces the word closure of fuchs.closure_mp, which at genus 69 needed long words in long generators: their
products amplified the development's relation error into garbage, E94 U4 run 5). Nearest-lift walk over the surface mesh:
every surface triangle (a triangle of the cut-open disk, positions zc in the base-point frame of develop_pre.npz) gets one lift
g_t (T_t) near the base point 0, by a Dijkstra-type search over triangle lifts ordered by the hyperbolic distance of the lift's
centroid from 0, crossing cut edges with their side-pairing generators (identified in double, products formed in 200 bits). Each
g_t maps a disk triangle (distance <= 15.6) to near 0, so every partial product stays short (no amplification by long words).
For adjacent triangles t, u whose chosen lifts are not adjacent in H^2 the element gamma = g_t s g_u^-1 (s = side pairing or
id) is a Voronoi neighbour of 0: its bisector passes between the two lifts. These candidates (and inverses) seed the Dirichlet
polygon; certification is unchanged (bounded, area 4 pi (g - 1) to 1e-8)."""
import heapq
import json
import numpy as np
import gmpy2
from gmpy2 import mpc, mpfr
import fuchs as fu


def bfs_candidates(npz, log=print):
    zc, Tn, orig, root = npz['zc'], npz['Tn'], npz['orig'], int(npz['root'])
    gs = json.loads(str(npz['gens_str']))
    gens = []
    for g in gs:
        m = [mpc(mpfr(re), mpfr(im)) for re, im in g]
        gens.append(fu.mp_norm(((m[0] + m[3].conjugate())/2, (m[1] + m[2].conjugate())/2)))
    G2 = gens + [fu.mp_inv(g) for g in gens]
    Gd = np.array([fu.mp_to_np(g) for g in G2])
    nT = len(Tn)
    # directed edges -> triangle
    de = {}
    for t in range(nT):
        a, b, c = (int(x) for x in Tn[t])
        de[(a, b)] = t; de[(b, c)] = t; de[(c, a)] = t
    bnd = {}
    for (a, b), t in de.items():
        if (b, a) not in de:
            bnd[(int(orig[a]), int(orig[b]))] = (t, a, b)
    nbr = [[] for _ in range(nT)]          # (u, generator index or -1)
    worst = 0.0; nb = 0
    for t in range(nT):
        a, b, c = (int(x) for x in Tn[t])
        for (p, q) in ((a, b), (b, c), (c, a)):
            if (q, p) in de:
                nbr[t].append((de[(q, p)], -1)); continue
            t2, p2, q2 = bnd[(int(orig[q]), int(orig[p]))]           # partner edge q' -> p' (orig q, orig p)
            # side pairing s with s(z[p2... ]) : s maps z[q2] -> z[p] and z[p2] -> z[q]
            w1 = (Gd[:, 0, 0]*zc[q2] + Gd[:, 0, 1])/(Gd[:, 1, 0]*zc[q2] + Gd[:, 1, 1])
            w2 = (Gd[:, 0, 0]*zc[p2] + Gd[:, 0, 1])/(Gd[:, 1, 0]*zc[p2] + Gd[:, 1, 1])
            err = np.maximum(fu.hdist(w1, zc[p]), fu.hdist(w2, zc[q]))
            k = int(np.argmin(err)); worst = max(worst, float(err[k])); nb += 1
            nbr[t].append((t2, k))
    log(f'U4 seed: {nT} triangles, {nb} cut-edge crossings, side pairings identified with max endpoint error {worst:.2e}')
    cen = zc[Tn].mean(1)
    E_mp = [(mpc(1), mpc(0))]; E_d = [np.eye(2, dtype=complex)]; cache = {}
    def child(e, k):
        key = (e, k)
        if key not in cache:
            m = fu.mp_mul(E_mp[e], G2[k]); E_mp.append(m); E_d.append(fu.mp_to_np(m)); cache[key] = len(E_mp) - 1
        return cache[key]
    def dist(e, t):
        g = E_d[e]; w = (g[0, 0]*cen[t] + g[0, 1])/(g[1, 0]*cen[t] + g[1, 1])
        return 2*np.arctanh(min(abs(w), 1 - 1e-16))
    el = -np.ones(nT, int); H = [(dist(0, root), root, 0)]
    while H:
        d, t, e = heapq.heappop(H)
        if el[t] >= 0: continue
        el[t] = e
        for u, k in nbr[t]:
            if el[u] < 0:
                e2 = e if k < 0 else child(e, k)
                heapq.heappush(H, (dist(e2, u), u, e2))
    dmax = max(dist(int(el[t]), t) for t in range(0, nT, 7))
    log(f'U4 seed: lifts chosen for {int((el >= 0).sum())} of {nT} triangles; tile elements used {len(set(el.tolist()))}; '
        f'max centroid distance (sampled) {dmax:.4f}')
    cand = {}; seen = set()
    for t in range(nT):
        for u, k in nbr[t]:
            e_t, e_u = int(el[t]), int(el[u])
            if k < 0 and e_t == e_u: continue
            if k >= 0 and cache.get((e_t, k)) == e_u: continue
            key = (e_t, k, e_u)
            if key in seen: continue
            seen.add(key)
            gt = E_mp[e_t] if k < 0 else fu.mp_mul(E_mp[e_t], G2[k])
            gam = fu.mp_mul(gt, fu.mp_inv(E_mp[e_u]))
            if abs(gam[0]) < fu.IDCUT: continue
            for h in (gam, fu.mp_inv(gam)):
                cand.setdefault(fu.mp_key(h), h)
    C = list(cand.values())
    dd = np.array([fu.mp_d(h)[0] for h in C])
    log(f'U4 seed: Voronoi-neighbour candidates {len(C)} (with inverses), displacement {dd.min():.4f} .. {dd.max():.4f}')
    return C


def dirichlet_from_seed(C, area_target, rounds=20, log=print):
    """Dirichlet polygon from the seed C (200 bits), then rounds of face-local pair products whose bisector cuts the polygon,
    with displacement capped at 2 rD + 1 (rD = circumradius of the current bounded polygon; a face element of D displaces 0 by at
    most 2 rD, so the cap loses no face once the polygon is bounded); certified when bounded with area 4 pi (g - 1) to 1e-8."""
    L = list(C); res = fu.dirichlet_mp(L)
    for it in range(rounds):
        log(f'U4 round {it}: elements {len(L)}, faces {res["nfaces"]}, bounded {res["bounded"]}, area diff '
            f'{res["area"] - area_target:.3e}, circumradius {res["rD"]:.6f}', flush=True)
        if res['bounded'] and abs(res['area'] - area_target) < 1e-8:
            break
        Lok = [L[i] for i in res['ok_index']]
        F = fu.mp_dedupe([Lok[t] for t in res['tag'] if t >= 0] + [fu.mp_inv(Lok[t]) for t in res['tag'] if t >= 0])
        known = set(fu.mp_key(e) for e in L)
        cap = 2*res['rD'] + 1 if res['bounded'] else np.inf
        add = [e for e in fu.products_cut_screened(F, res, dmax=cap) if fu.mp_key(e) not in known]
        log(f'   face elements {len(F)}, cap {cap:.4f}, cutting products {len(add)}', flush=True)
        if not add:
            break
        L = fu.mp_dedupe(F + add); res = fu.dirichlet_mp(L)
    return L, res
