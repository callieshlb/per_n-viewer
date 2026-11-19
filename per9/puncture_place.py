#!/usr/bin/env python3
"""E95 A3 repair (2026-11-12): every puncture is the boundary point of multiplier e^{-2 pi i p/q} of the hyperbolic
component of a polynomial centre in the p/q-limb. For each census puncture P (E90 data_n9 json), its first centre (polys[0], which
defines loop point 0) is joined to the loop point by the census's own path, re-walked here with every step recorded:
   pern.ray_to_puncture (lambda = (1 - e^-s) e^{-2 pi i p/q}, Newton on (b, c, z)) until |g| <= eps, then pern.snap_to_circle.
Along the path (b, c) gives y = a_3 and x = a_8 by the dynamics (pern.orbit): no root-finding on the curve. The SHEET is carried
combinatorially: slit label at the centre (512-bit continuation y0 -> y_centre, which worked for all 252 centres in U5) and then
the cut crossings of the y-polyline (cutglue.CutSystem.edge_perm, as the gluing). Steps are refined until |dy| <= 0.25 x the
distance from y to the nearest branch candidate, so each chord is homotopic to the true path in the complement of the candidates.
Place of the puncture: if the loop point is within RB of a branch candidate b: (b, cycle of the label under tau_b); otherwise
(sheet = X0 index of the label, y rounded 1e-2) as E94's A3.
Canary (exit status): at the punctures whose loop point WAS continued in atlas_A3 (sheet not None), the sheet found here must equal
the continuation's sheet. Writes puncture_place.json; --locate adds the mesh location (face, weights) on the fine mesh.
usage: puncture_place.py [worker k/N]   (results cached per puncture in pp_cache/)."""
import json, os, pickle, sys, time
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, 'surface3d'))
E90 = os.path.join(HERE, '..', 'inputs', 'E90-real-circles'); sys.path.insert(0, E90)
import pern
from meshloc import Locator, CS, B as BR
n = 9; RB = 0.01
census = json.load(open(os.path.join(E90, 'data_n9_eps1e-3_dps60.json')))
dps = int(census['dps']); eps = mp.mpf(census['eps']); mp.mp.dps = dps; tol = mp.mpf(10)**(-dps//2)
BPS = np.array([complex(*b) for b in BR['bps']])
A3 = {r['i']: r for r in (json.loads(l) for l in open(os.path.join(HERE, 'atlas_A3_ckpt.jsonl')))}
os.makedirs(os.path.join(HERE, 'pp_cache'), exist_ok=True)


MOVES = np.array([[CS.perms[k][l] != l for l in range(CS.d)] for k in range(len(BPS))])     # candidate k ramifies the sheet with slit label l
def dist_c(y, labs=None):
    """distance from y to the nearest candidate that moves one of the given slit labels (all candidates if labs is None). A real ray
    passes exactly through real candidates on sheets that are unramified there; those do not constrain the chord."""
    m = np.ones(len(BPS), bool) if labs is None else MOVES[:, list(labs)].any(1)
    return float(np.abs(BPS[m] - complex(y)).min()) if m.any() else np.inf
def xy(b, c):
    a = pern.orbit(b, c, n); return complex(a[n - 1]), complex(a[3])


def walk(cen, lab):
    """recorded path centre -> loop point; returns list of (x, y) and the end (b, c)."""
    p, q = cen['limb']; c0 = mp.mpc(mp.mpf(cen['c'][0]), mp.mpf(cen['c'][1]))
    b, c, kap = pern.center_to_chart(c0); z = kap
    u = mp.expjpi(-2*mp.mpf(p)/q); s = mp.mpf(0); ds = mp.mpf('0.05'); pts = [xy(b, c)]; steps = rej = 0; mind = np.inf
    while True:                                               # pern.ray_to_puncture with recording and the |dy| rule
        s2 = s + ds; lam = (1 - mp.exp(-s2))*u
        res = pern.newton_ray(b, c, z, lam, n, tol)
        ok = res is not None and pern.guard(b, c, res[0], res[1], n)
        if ok:
            xn, yn = xy(res[0], res[1]); yo = pts[-1][1]
            lab2 = CS.edge_perm([yo, yn])[lab]
            dd_ = min(dist_c(yn, (lab, lab2)), dist_c(yo, (lab, lab2)))
            if abs(yn - yo) > 0.25*dd_ and dd_ > 1e-9: ok = False; rej += 1     # within 1e-9 of a ramifying candidate the path has arrived over it (place = label cycle; double y cannot resolve further)
        if not ok:
            ds /= 2
            if ds < mp.mpf(10)**-14: raise pern.Fail(f'ray stalled s={mp.nstr(s, 8)} y={pts[-1][1]}')
            continue
        b, c, z = res[0], res[1], res[2]; s = s2; steps += 1; pts.append((xn, yn)); lab = lab2; mind = min(mind, dd_)
        if res[3] <= 5: ds = min(ds*1.3, mp.mpf('0.2'))
        g, z3, ms = pern.small_g(b, c, z)
        if abs(g) <= eps: break
        if s > 60: raise pern.Fail(f'ray reached s={s}')
    ang = mp.arg(g); h = int(mp.nint(ang/mp.pi))              # pern.snap_to_circle = g_path_log(log g -> log eps + i h pi), recorded
    lf, lt = mp.log(g), mp.log(eps) + 1j*h*mp.pi; L = abs(lt - lf); tau = mp.mpf(0)
    dt = min(mp.mpf(1), mp.mpf('0.12')/L) if L > 0 else mp.mpf(1); z = z3
    while tau < 1:
        t2 = min(mp.mpf(1), tau + dt); gg = mp.exp(lf + (lt - lf)*t2)
        res = pern.newton_g(b, c, z, gg, n, tol)
        ok = res is not None and pern.guard(b, c, res[0], res[1], n)
        if ok:
            xn, yn = xy(res[0], res[1]); yo = pts[-1][1]
            lab2 = CS.edge_perm([yo, yn])[lab]
            dd_ = min(dist_c(yn, (lab, lab2)), dist_c(yo, (lab, lab2)))
            if abs(yn - yo) > 0.25*dd_ and dd_ > 1e-9: ok = False; rej += 1     # within 1e-9 of a ramifying candidate the path has arrived over it (place = label cycle; double y cannot resolve further)
        if not ok:
            dt /= 2
            if dt*L < mp.mpf(10)**-12: raise pern.Fail('snap stalled')
            continue
        b, c, z = res[0], res[1], res[2]; tau = t2; steps += 1; pts.append((xn, yn)); lab = lab2; mind = min(mind, dd_)
        if res[3] <= 5: dt = min(dt*1.25, mp.mpf('0.12')/L)
    return pts, b, c, steps, rej, lab, mind


LOC = None
def one(zi):
    global LOC
    fn = os.path.join(HERE, 'pp_cache', f'p{zi}.pkl')
    if os.path.exists(fn): return pickle.load(open(fn, 'rb'))
    if LOC is None: LOC = Locator(os.path.join(HERE, 'surface3d'))
    mp.mp.dps = dps; t0 = time.time()
    P = census['punctures'][zi]; polys = P['polys'] if not isinstance(P['polys'], str) else eval(P['polys'])
    ci = int(polys[0]['center']); cen = census['centers'][ci]
    out = dict(i=zi, centre=ci, ks=cen['ks'], limb=cen['limb'], status='ok')
    try:
        c0 = mp.mpc(mp.mpf(cen['c'][0]), mp.mpf(cen['c'][1])); b0, c0_, _ = pern.center_to_chart(c0); xc, yc = xy(b0, c0_)
        roots = LOC.roots_X0(yc)
        if roots is None: raise pern.Fail('centre continuation failed')
        dd = sorted((abs(r - xc), j) for j, r in enumerate(roots)); jc = dd[0][1]
        lab_c = CS.h_at(yc)[jc]
        pts, b, c, steps, rej, lab_e, mind = walk(cen, lab_c)
        xe, ye = pts[-1]
        assert CS.edge_perm([p_[1] for p_ in pts])[lab_c] == lab_e
        he = CS.h_at(ye); sheet = he.index(lab_e)
        kb = int(np.argmin(np.abs(BPS - ye))); db = float(abs(BPS[kb] - ye))
        grp = CS.groups(kb)[0][lab_e] if db < RB else None
        out.update(x=[xe.real, xe.imag], y=[ye.real, ye.imag], b=mp.nstr(b, 25), c=mp.nstr(c, 25), steps=steps, rejected=rej, npts=len(pts),
                   centre_ratio=float(dd[1][0]/max(dd[0][0], 1e-300)), label_centre=int(lab_c), label=int(lab_e), sheet=int(sheet),
                   near_b=kb if db < RB else None, dist_b=db, group=None if grp is None else int(grp),
                   ymax=float(max(abs(p_[1]) for p_ in pts)), min_dist_c=float(mind), seconds=round(time.time() - t0, 1))
    except pern.Fail as ex:
        out.update(status='fail', why=str(ex), seconds=round(time.time() - t0, 1))
    pickle.dump(out, open(fn + f'.tmp{os.getpid()}', 'wb')); os.replace(fn + f'.tmp{os.getpid()}', fn)
    return out


if __name__ == '__main__':
    NP = len(census['punctures'])
    if len(sys.argv) > 1 and '/' in sys.argv[1]:
        k, N = map(int, sys.argv[1].split('/'))
        for zi in list(range(NP))[k::N]:
            r = one(zi); print(f'   puncture {zi}: {r["status"]} {r.get("why", "")} sheet {r.get("sheet")} steps {r.get("steps")} ({r["seconds"]} s)', flush=True)
        sys.exit(0)
    R = [one(zi) if os.path.exists(os.path.join(HERE, 'pp_cache', f'p{zi}.pkl')) or not os.environ.get('PP_SKIP_PENDING') else dict(i=zi, status='pending', why='still running') for zi in range(NP)]
    ok = [r for r in R if r['status'] == 'ok']
    print(f'punctures {NP}: placed {len(ok)}, failed {NP - len(ok)} {[(r["i"], r["why"][:40]) for r in R if r["status"] != "ok"]}')
    print(f'path steps median {int(np.median([r["steps"] for r in ok]))}, max {max(r["steps"] for r in ok)}; steps rejected by the |dy| rule: total {sum(r["rejected"] for r in ok)}; '
          f'centre label: min runner-up/nearest ratio {min(r["centre_ratio"] for r in ok):.2e}; min distance of a path point to a candidate that ramifies its sheet {min(r["min_dist_c"] for r in ok):.2e}')
    cmpd = agree = 0; dis = []; nref_bad = []
    for r in ok:
        a = A3.get(r['i'])
        if a is None or a['sheet'] is None: continue
        if a['droot'] > 1e-6*(1 + abs(complex(*a['x']))): nref_bad.append(r['i']); continue   # the A3 point itself is not on a root there (puncture 21: nearest root 309.6 away): no reference
        if abs(complex(*a['y']) - complex(*r['y'])) > 1e-12*(1 + abs(complex(*r['y']))): nref_bad.append(r['i']); continue   # same point only: atlas_A3 moved loop points within 1e-12 of y = 1 by 1e-9 and took the nearest root there (not a continuation)
        cmpd += 1
        if a['sheet'] == r['sheet']: agree += 1
        else: dis.append((r['i'], a['sheet'], r['sheet'], r['y']))
    print(f'CANARY vs the A3 continuation (same loop point, sheet continued there): compared {cmpd}, agree {agree}, disagree {len(dis)} {dis[:8]}; excluded (A3 reference point not on a root, or moved off the loop point): {nref_bad}')
    places = {}
    for r in ok:
        key = ('b', r['near_b'], r['group']) if r['near_b'] is not None else ('s', r['sheet'], round(r['y'][0], 2), round(r['y'][1], 2))
        places.setdefault(key, []).append(r['i'])
    dup = {str(k): v for k, v in places.items() if len(v) > 1}
    nb = sum(1 for r in ok if r['near_b'] is not None)
    print(f'A3 (label method): {len(ok)} punctures placed, distinct places {len(places)}; over a branch candidate (within {RB}): {nb}, elsewhere {len(ok) - nb}; shared places: {dup}')
    from collections import Counter
    print('punctures per branch candidate (index: y: count):', [(k, str(np.round(BPS[k], 6)), v) for k, v in Counter(r['near_b'] for r in ok if r['near_b'] is not None).most_common()])
    json.dump(R, open(os.path.join(HERE, 'puncture_place.json'), 'w'), indent=1)
    sys.exit(0 if len(dis) == 0 and cmpd > 0 else 1)
