#!/usr/bin/env python3
"""E99: 40-digit repair of edge_perms.json in the zone Z = {|y| < 0.05} u {0 < |y - 1| < 1e-6} (REGISTRATION.md log).
A port of the explorer tracker (explorer/per6_core_mesh.js: coeffs, Newton, waypoints with the node detour, trackSegment,
track, globalRoots) and of edge_perms.js (sphereTrack, trackEdge with the +-20% midpoint retry, branch-vertex offset loop and
toGroup edges) to mpmath. Differences from the double code, all registered: step factor 0.05 (not 0.15), Newton tolerance
1e-24 (not 1e-11), distinctness chordal 1e-25 (not 1e-10), mp.dps 40.
Recomputes R[v] for every normal vertex in Z and every edge with an endpoint in Z (or a null double result anywhere).
Second route for the labels (independent path): node: double roots at q = 0.08 y/|y| (explorer core, via node), then a
40-digit straight track q -> y; y = 1: from the nearest normal vertex outside Z (double R) by a 40-digit straight track.
Usage: python3 edge_perms_mp.py [NPROC]. Reads sphere_mesh.json, edge_perms.json; writes edge_perms.json (backup kept as
edge_perms_double.json) and edge_perms_mp.log.json."""
import json, os, sys, time, subprocess
import mpmath as mp
from multiprocessing import Pool
mp.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, '..', 'explorer', 'data.json')))
COEF = [[mp.mpf(int(c)) for c in row] for row in D['coef']]
assert all(float(c) == int(c) for row in D['coef'] for c in row)
SPECIALS = [mp.mpc(*b) for b in D['bps']]
NODE_R = mp.mpf('0.12'); FLOOR = mp.mpf('1e-14'); SFAC = mp.mpf('0.05'); NTOL = mp.mpf('1e-24'); DIST = mp.mpf('1e-25')
def chord(a, b):
    return abs(a - b)/mp.sqrt((1 + abs(a)**2)*(1 + abs(b)**2))
def coeffs(y):
    c, dc = [], []
    for row in COEF:
        p = mp.mpc(0); q = mp.mpc(0)
        for a in reversed(row):
            q = q*y + p; p = p*y + a
        c.append(p); dc.append(q)
    e = y - 1; c[-1] = e*e; dc[-1] = 2*e
    return c, dc
def peval(x, cc):
    c, dc = cc; p = px = py = mp.mpc(0)
    for k in range(len(c) - 1, -1, -1):
        px = px*x + p; p = p*x + c[k]; py = py*x + dc[k]
    return p, px, py
def newton(x, cc, maxit):
    for it in range(maxit):
        p, px, _ = peval(x, cc); dx = p/px; x = x - dx
        if abs(dx) <= NTOL*(1 + abs(x)): return x, it + 1
    return x, -1
def dist_special(y): return min(abs(y - s) for s in SPECIALS)
def arc_to(r, f1, f2, n):
    df = ((f2 - f1 + mp.pi) % (2*mp.pi) + 2*mp.pi) % (2*mp.pi) - mp.pi
    return [mp.mpc(r*mp.cos(f1 + df*i/n), r*mp.sin(f1 + df*i/n)) for i in range(1, n + 1)]
def waypoints(a, b):
    v = b - a; L2 = abs(v)**2
    if L2 == 0: return [a, b]
    ra, rb = abs(a), abs(b)
    if rb < NODE_R or ra < NODE_R:
        tt = min(1, max(0, -(a.real*v.real + a.imag*v.imag)/L2))
        cl = abs(a + tt*v)
        if cl >= 0.5*min(ra, rb): return [a, b]
        r = max(min(max(ra, rb), NODE_R), mp.mpf('1e-12'))
        A = r*a/ra; Bp = r*b/rb
        path = [a]
        if ra > r: path = waypoints(a, A)
        elif ra < r: path.append(A)
        path = path + arc_to(r, mp.atan2(a.imag, a.real), mp.atan2(b.imag, b.real), 24)
        if rb != r: path.append(b)
        return path
    t = -(a.real*v.real + a.imag*v.imag)/L2; cl = a + t*v
    if not (0 < t < 1 and abs(cl) < NODE_R and abs(a) > NODE_R and abs(b) > NODE_R): return [a, b]
    A = L2; Bq = 2*(a.real*v.real + a.imag*v.imag); C = abs(a)**2 - NODE_R**2
    disc = mp.sqrt(max(Bq*Bq - 4*A*C, 0)); t1 = (-Bq - disc)/(2*A); t2 = (-Bq + disc)/(2*A)
    p1 = a + t1*v; p2 = a + t2*v
    f1 = mp.atan2(p1.imag, p1.real); f2 = mp.atan2(p2.imag, p2.real)
    df = ((f2 - f1 + mp.pi) % (2*mp.pi) + 2*mp.pi) % (2*mp.pi) - mp.pi
    out = [a, p1] + [mp.mpc(NODE_R*mp.cos(f1 + df*i/24), NODE_R*mp.sin(f1 + df*i/24)) for i in range(1, 25)] + [b]
    return out
def dxdy_of(cur, cc):
    out = []
    for x in cur:
        p, px, py = peval(x, cc); out.append(-py/px)
    return out
DEADLINE = [None]
class Budget(Exception): pass
def track_segment(xs, a, b):
    v = b - a; L = abs(v)
    if L == 0: return list(xs)
    t = mp.mpf(0); cur = list(xs); cc = coeffs(a); dxdy = dxdy_of(cur, cc); guard = 0
    while t < 1 and guard < 2000000:
        guard += 1
        if DEADLINE[0] is not None and (guard & 63) == 0 and time.time() > DEADLINE[0]: raise Budget()
        ycur = a + t*v
        h = min(max(mp.mpf('0.5'), mp.mpf('0.02')*abs(ycur)), max(FLOOR, SFAC*dist_special(ycur)))   # run 5: cap 0.5 -> max(0.5, 0.02|y|)
        dt = min(1 - t, h/L); ok = False
        for tries in range(40):
            t2 = t + dt; y2 = a + t2*v; dy = v*dt; cc2 = coeffs(y2); nx = []; good = True
            for r in range(len(cur)):
                pred = cur[r] + dxdy[r]*dy
                xr, it = newton(pred, cc2, 12)
                if it < 0 or abs(xr - pred) > mp.mpf('0.3')*(abs(dxdy[r]*dy) + mp.mpf('1e-9')*(1 + abs(cur[r]))) + mp.mpf('1e-12'):
                    good = False; break
                nx.append(xr)
            if good:
                for i in range(len(nx)):
                    for j in range(i + 1, len(nx)):
                        if chord(nx[i], nx[j]) < DIST: good = False
            if good:
                ok = True; t = t2; cur = nx; cc = cc2; dxdy = dxdy_of(cur, cc); break
            dt /= 2
        if not ok: return None
    return cur
def track(xs, a, b):
    w = waypoints(a, b); cur = xs
    for i in range(1, len(w)):
        if cur is None: return None
        cur = track_segment(cur, w[i - 1], w[i])
    return cur
X0 = None
def global_roots(y):
    global X0
    if X0 is None:
        y0 = mp.mpc(*D['y0']); cc = coeffs(y0); X0 = [newton(mp.mpc(*x), cc, 50)[0] for x in D['x0']]
    return track(X0, mp.mpc(*D['y0']), y)
def to_u(y):
    r2 = abs(y)**2; return [2*y.real/(r2 + 1), 2*y.imag/(r2 + 1), (r2 - 1)/(r2 + 1)]
def to_y(u): return mp.mpc(u[0]/(1 - u[2]), u[1]/(1 - u[2]))
def track_straight(xs, a, c, stats):
    r = track(xs, a, c)
    if r is not None: return r
    stats['retried'] += 1
    v = c - a; mid = (a + c)/2; res = []
    for s in (mp.mpf('0.2'), mp.mpf('-0.2')):
        m = mid + mp.mpc(-s*v.imag, s*v.real)
        r1 = track(xs, a, m); r2 = track(r1, m, c) if r1 is not None else None
        if r2 is not None: res.append(r2)
    if len(res) == 2 and any(chord(x, y) > mp.mpf('1e-8') for x, y in zip(res[0], res[1])): stats['retry_disagree'] += 1
    return res[0] if res else None
def sphere_track(xs, a, c, stats):
    ua, uc = to_u(a), to_u(c)
    dot = max(-1, min(1, sum(p*q for p, q in zip(ua, uc)))); om = mp.acos(dot); K = 16
    cur, prev = xs, a
    for k in range(1, K + 1):
        if cur is None: break
        t = mp.mpf(k)/K; s1 = mp.sin((1 - t)*om)/mp.sin(om); s2 = mp.sin(t*om)/mp.sin(om)
        u = [s1*p + s2*q for p, q in zip(ua, uc)]
        q = c if k == K else to_y(u)
        cur = track_straight(cur, prev, q, stats); prev = q
    return cur
def track_edge(xs, a, c, stats):
    if max(abs(a), abs(c)) > 1: return sphere_track(xs, a, c, stats)
    return track_straight(xs, a, c, stats)
def nearest(x, xs):
    ds = [chord(x, y) for y in xs]; i = min(range(len(xs)), key=lambda k: ds[k]); return i, ds[i]
def second(ds):
    s = sorted(ds); return s[1] if len(s) > 1 else mp.inf

# ---------------- worker tasks ----------------
def task_vertex(args):
    v, y = args; y = mp.mpc(*y); t0 = time.time()
    r = global_roots(y)
    if r is None: return v, None, None, time.time() - t0
    # route B
    if abs(y) < mp.mpf('0.05'):
        q = mp.mpf('0.08')*y/abs(y)
        rb = BQ.get(v)
        rbq = [mp.mpc(*z) for z in rb] if rb is not None else None
        rB = track(rbq, q, y) if rbq is not None else None
    else:
        a_y, a_R = NB.get(v, (None, None))
        rB = track([mp.mpc(*z) for z in a_R], mp.mpc(*a_y), y) if a_R is not None else None
    lab = None
    if rB is not None:
        lab = [nearest(x, r)[0] for x in rB]
    sep = min(chord(r[i], r[j]) for i in range(5) for j in range(i + 1, 5))
    return v, [[float(z.real), float(z.imag)] for z in r], dict(labB=lab, sep=float(sep), rmp=[[mp.nstr(z.real, 40), mp.nstr(z.imag, 40)] for z in r]), time.time() - t0
OFF = {}
def task_edge(args):
    k, a, c, ya, yc, Ra, Rc, grp, off = args
    stats = dict(retried=0, retry_disagree=0)
    ya, yc = mp.mpc(*ya), mp.mpc(*yc); Ra = [mp.mpc(*z) for z in Ra]
    if grp is None:
        Rc = [mp.mpc(*z) for z in Rc]
        sa, sc = 1/ya, 1/yc; vv = sc - sa
        tt = min(1, max(0, -(sa.real*vv.real + sa.imag*vv.imag)/abs(vv)**2)) if abs(vv) > 0 else 0
        near_pole = abs(ya) > 1e3 and abs(yc) > 1e3 and abs(sa + tt*vv) < 0.25*min(abs(sa), abs(sc))
        tr = None if near_pole else track_edge(Ra, ya, yc, stats)
        if tr is None and abs(ya) > 1e3 and abs(yc) > 1e3:      # near y = inf (pinf = id): detour via m = i ya, same |y|
            stats['inf_detour'] = 1; mm = 1j*ya
            t1 = track_edge(Ra, ya, mm, stats); tr = track_edge(t1, mm, yc, stats) if t1 is not None else None
        if tr is None: return k, None, stats, None
        m, worst, marg = [], 0, mp.inf
        for x in tr:
            i, d = nearest(x, Rc); m.append(i); worst = max(worst, d)
            marg = min(marg, second([chord(x, z) for z in Rc])/max(d, mp.mpf('1e-300')))
        return k, dict(a=a, c=c, m=m), stats, (float(worst), float(marg))
    b = yc; dvec = ya - b; L = abs(dvec)
    E2 = mp.mpf(off['eps']); s = b + E2*dvec/L
    tr = track_edge(Ra, ya, s, stats)
    if tr is None: return k, None, stats, None
    o = mp.mpc(*off['o']); a0 = mp.atan2((o - b).imag, (o - b).real); a1 = mp.atan2((s - b).imag, (s - b).real)
    da = ((a1 - a0 + mp.pi) % (2*mp.pi) + 2*mp.pi) % (2*mp.pi) - mp.pi
    cur = [mp.mpc(*z) for z in off['xs']]; prev = o
    for i in range(1, 33):
        if cur is None: break
        qq = b + E2*mp.expj(a0 + da*i/32); cur = track(cur, prev, qq); prev = qq
    if cur is None: return k, None, stats, None
    m, worst = [], 0
    for x in tr:
        i, d = nearest(x, cur); m.append(grp['g'][i]); worst = max(worst, d)
    return k, dict(a=a, c=c, m=m, toGroup=True), stats, (float(worst), None)
def task_branch(args):
    v, b = args; b = mp.mpc(*b)
    for eps in ('1e-4', '3e-4', '1e-3', '3e-3'):
        eps = mp.mpf(eps)
        for th in ('0.3', '1.9', '3.4', '5.0'):
            th = mp.mpf(th); o = b + eps*mp.expj(th)
            xs = global_roots(o)
            if xs is None: continue
            cur, prev = xs, o
            for i in range(1, 65):
                if cur is None: break
                q = b + eps*mp.expj(th + 2*mp.pi*i/64); cur = track(cur, prev, q); prev = q
            if cur is None: continue
            p = [nearest(x, xs)[0] for x in cur]
            g = [-1]*5; ng = 0
            for i in range(5):
                if g[i] >= 0: continue
                j = i
                while g[j] < 0: g[j] = ng; j = p[j]
                ng += 1
            return v, dict(o=[float(o.real), float(o.imag)], xs=[[float(z.real), float(z.imag)] for z in xs], eps=float(eps)), \
                dict(g=g, ng=ng, perm=p), dict(o_mp=(mp.nstr(o.real, 40), mp.nstr(o.imag, 40)))
    return v, None, None, None
def task_edge_router(j):
    if j[3] is None or (j[7] is not None and j[8] is None): return j[0], None, dict(retried=0, retry_disagree=0, timeout=0), None
    DEADLINE[0] = time.time() + 120
    try:
        r = task_edge(j); r[2]['timeout'] = 0; return r
    except Budget:
        return j[0], None, dict(retried=0, retry_disagree=0, timeout=1), None
    finally:
        DEADLINE[0] = None
BQ, NB = {}, {}
def init_worker(bq, nb):
    global BQ, NB; BQ, NB = bq, nb
if __name__ == '__main__':
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    M = json.load(open('sphere_mesh.json')); EP = json.load(open('edge_perms.json'))
    if not os.path.exists('edge_perms_double.json'): os.replace('edge_perms.json', 'edge_perms_double.json'); json.dump(EP, open('edge_perms.json', 'w'))
    Y = [complex(*y) for y in M['y']]; isB = [b >= 0 for b in M['bp']]; nV = len(Y)
    inZ = [(not isB[v]) and (abs(Y[v]) < 0.05 or 0 < abs(Y[v] - 1) < 1e-6 or 0 < abs(Y[v] - 0.5) < 1e-6) for v in range(nV)]   # run 3: + y = 1/2 (cluster conditioning)
    Zv = [v for v in range(nV) if inZ[v]]
    print(f'zone Z: {len(Zv)} normal vertices (node {sum(abs(Y[v]) < 0.05 for v in Zv)}, y=1 {sum(abs(Y[v] - 1) < 1e-6 for v in Zv)}, y=1/2 {sum(abs(Y[v] - 0.5) < 1e-6 for v in Zv)}); '
          f'double failures inside Z {sum(EP["R"][v] is None for v in Zv)}, outside Z {sum(EP["R"][v] is None for v in range(nV) if not isB[v] and not inZ[v])}', flush=True)
    # route-B starting data: node -> double roots at q = 0.08 y/|y| via the explorer core (node); y=1 -> nearest outside-Z neighbour
    nodeZ = [v for v in Zv if abs(Y[v]) < 0.05]
    qs = [[0.08*Y[v].real/abs(Y[v]), 0.08*Y[v].imag/abs(Y[v])] for v in nodeZ]
    json.dump(qs, open('_q.json', 'w'))
    js = ("const C=require('../explorer/per6_core.js');C.init(require('../explorer/data.json'));"
          "const q=require('./_q.json');require('fs').writeFileSync('_qR.json',JSON.stringify(q.map(y=>C.globalRoots(y))));")
    subprocess.run(['node', '-e', js], check=True)
    qR = json.load(open('_qR.json')); bq = {v: r for v, r in zip(nodeZ, qR)}
    print(f'route B node starts: {len(nodeZ)}, double failures at q {sum(r is None for r in qR)}', flush=True)
    adj = {v: set() for v in Zv}
    for p, q in M['E']:
        if p in adj: adj[p].add(q)
        if q in adj: adj[q].add(p)
    nb = {}
    for v in Zv:
        if abs(Y[v]) < 0.05: continue
        # BFS outward through Z to the nearest normal vertex outside Z with double roots
        seen, fr, found = {v}, [v], None
        while fr and found is None:
            nf = []
            for u in fr:
                for w in sorted(adj[u]):
                    if w in seen: continue
                    seen.add(w)
                    if not isB[w] and not inZ[w] and EP['R'][w] is not None: found = w; break
                    if inZ[w]: nf.append(w)
                if found is not None: break
            fr = nf
        if found is not None: nb[v] = (M['y'][found], EP['R'][found])
    print(f'route B y=1 starts: {len(nb)} of {len(Zv) - len(nodeZ)}', flush=True)
    t0 = time.time(); logd = dict(vertices={}, edges={})
    CK = 'edge_perms_mp_vertices.json'
    if os.path.exists(CK): VR = {int(k): tuple(x) for k, x in json.load(open(CK)).items()}; print(f'vertex results loaded from {CK}: {len(VR)}', flush=True)
    else: VR = None
    with Pool(nproc, initializer=init_worker, initargs=(bq, nb)) as pool:
        nfail = 0; nB_mis = 0; nB_none = 0; dbl_dis = 0; dbl_lab = 0; seps = []
        if VR is None: VR = {}
        todo = [v for v in Zv if v not in VR]
        print(f'vertices to compute: {len(todo)}', flush=True)
        for v, R, info, dtv in pool.imap_unordered(task_vertex, [(v, M['y'][v]) for v in todo], chunksize=8): VR[v] = (R, info)
        json.dump({str(k): list(x) for k, x in VR.items()}, open(CK, 'w'))
        for v in Zv:
            R, info = VR[v]
            old = EP['R'][v]
            if R is None: nfail += 1; EP['R'][v] = None; continue
            if info['labB'] is None: nB_none += 1
            elif info['labB'] != list(range(5)): nB_mis += 1; logd['vertices'][v] = dict(labB=info['labB'])
            if old is not None:
                d = max(abs(complex(*a) - complex(*b)) for a, b in zip(old, R))
                if d > 1e-6: dbl_dis += 1
                lab = [min(range(5), key=lambda i: abs(complex(*o) - complex(*R[i]))) for o in old]
                if lab != list(range(5)): dbl_lab += 1
            EP['R'][v] = R; seps.append(info['sep'])
        print(f'vertices (40 digits): {len(Zv)}, failures {nfail}; route B: label mismatches {nB_mis}, not evaluable {nB_none}; '
              f'min chordal root separation {min(seps):.3e}; double results inside Z differing by > 1e-6: {dbl_dis}, with permuted labels: {dbl_lab}; '
              f'{time.time() - t0:.0f} s', flush=True)
        # branch vertex at y = 1 (offset data + groups) in 40 digits, compared with the double groups
        bvs = [v for v in range(nV) if isB[v] and (abs(Y[v] - 1) < 1e-12 or abs(Y[v] - 0.5) < 1e-12)]   # run 4: both cone branch vertices
        off = {}
        for v, o, grp, extra in pool.imap_unordered(task_branch, [(v, M['y'][v]) for v in bvs]):
            gd = EP['group'][str(v)]
            print(f'branch vertex {v} (y = {Y[v]}): 40-digit groups {grp}; double groups {gd}; offset eps {o["eps"] if o else None}', flush=True)
            off[v] = o
            if grp is not None: EP['group'][str(v)] = grp
        # edges
        jobs = []
        for k, (p, q) in enumerate(M['E']):
            if isB[p] and isB[q]: continue
            a, c = (q, p) if isB[p] else (p, q)
            ek = EP['edges'][k]
            nonperm = ek is not None and not ek.get('toGroup') and len(set(ek['m'])) != 5
            if nonperm: print('double non-permutation edge', k, 'y', M['y'][a], M['y'][c], 'in Z', inZ[a], inZ[c], flush=True)
            if not (inZ[a] or inZ[c] or ek is None or nonperm): continue
            if EP['R'][a] is None or (not isB[c] and EP['R'][c] is None): jobs.append((k, a, c, None, None, None, None, None, None)); continue
            if isB[c]:
                o = off.get(c)
                if o is None:
                    print('edge to a branch vertex without 40-digit offset data (outside y = 1):', k, a, c, flush=True); o = None
                jobs.append((k, a, c, M['y'][a], M['y'][c], EP['R'][a], None, EP['group'][str(c)], o))
            else:
                jobs.append((k, a, c, M['y'][a], M['y'][c], EP['R'][a], EP['R'][c], None, None))
        print(f'edges to recompute in 40 digits: {len(jobs)} (double nulls among them {sum(EP["edges"][j[0]] is None for j in jobs)})', flush=True)
        t1 = time.time(); efail = 0; changed = 0; worst = 0.0; marg = float('inf'); retried = 0; rdis = 0
        ntime = 0; done = 0; failed_ids = []
        for k, e, st, wm in pool.imap_unordered(task_edge_router, jobs, chunksize=4):
            done += 1
            if done % 2000 == 0: print(f'   edges done {done} of {len(jobs)} ({time.time() - t1:.0f} s)', flush=True)
            retried += st['retried']; rdis += st['retry_disagree']; ntime += st.get('timeout', 0); ninf = globals().setdefault('_ninf', [0]); ninf[0] += st.get('inf_detour', 0)
            if e is None: efail += 1; EP['edges'][k] = None; failed_ids.append(k); continue
            if EP['edges'][k] is None or EP['edges'][k].get('m') != e['m']: changed += 1
            EP['edges'][k] = e; worst = max(worst, wm[0])
            if wm[1] is not None: marg = min(marg, wm[1])
        bad_perm = sum(1 for j in jobs if EP['edges'][j[0]] is not None and not EP['edges'][j[0]].get('toGroup') and len(set(EP['edges'][j[0]]['m'])) != 5)
        logd['edge_failures'] = failed_ids
        print(f'edges (40 digits): {len(jobs)}, failures {efail} (of them 120 s timeouts {ntime}; ids in edge_perms_mp.log.json), infinity detours {_ninf[0]}, maps changed vs double {changed}, non-permutation maps {bad_perm}; '
              f'retried {retried} (both sides disagree {rdis}); worst chordal match {worst:.2e}; min runner-up/best match ratio {marg:.2e}; '
              f'{time.time() - t1:.0f} s', flush=True)
    json.dump(EP, open('edge_perms.json', 'w'))
    json.dump(logd, open('edge_perms_mp.log.json', 'w'))
    for f in ('_q.json', '_qR.json'): os.remove(f)
    print('wrote edge_perms.json', flush=True)
