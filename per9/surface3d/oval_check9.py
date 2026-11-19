#!/usr/bin/env python3
"""E95-U6 (E94 oval_check.py with the location step replaced; lifting, holonomy and axis test unchanged): are the traced real
ovals of Per_9(0) closed geodesics of the computed hyperbolic metric?
Ovals: oval_points.json (trace_ovals.py: the census's CLOSED real circles, a lower bound of the real locus: 8 real arcs are
unlanded; points per arc recorded in arc_end). Thinning as E94 (>= 1e-3 chordal in y or x between kept points), per arc.
Location (E94 used a 256-bit continuation per point, which fails near y = 0 here):
  - one ANCHOR per arc: a point with |y|, |x| < 2 whose sheet is found by the 512-bit continuation y0 -> y (unique: runner-up /
    nearest root distance >= 1e3); its slit label = h_at(y)[sheet]; mesh-independent, cached in oval_anchors.json;
  - the other points of the arc by CONTINUITY on the surface, outward from the anchor: of the 42 points of the surface over y
    (one per sheet, placed in the developed disk and folded into D), the one nearest to the previous point (distance modulo the
    face pairings). No root-finding. The margin second-nearest / nearest is recorded.
Lifting (E94): consecutive folded points are joined by g in {I} u {face elements} minimising the distance; H accumulated in 200
bits; holonomy h = H_N; geodesic length 2 arccosh(|tr h| / 2); distance of every point to the axis of h (in its local frame).
U6 thresholds (E93/E94, per oval, fine mesh): max axis distance < 2e-2, |polyline length / geodesic length - 1| < 1e-2, and both
smaller on the fine mesh than on the coarse one. usage: oval_check9.py DIR. Writes DIR/oval_check.json, DIR/oval_lifts.npz."""
import json, os, signal, sys, time
import numpy as np
import gmpy2
from gmpy2 import mpc
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import fuchs as fu
from meshloc import Locator, place, load_D, CS
from oval_check import axis_dist, fixed_points_mp
DIR = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else here)
gmpy2.get_context().precision = 200
t0 = time.time()
ovals = json.load(open(os.path.join(here, 'oval_points.json')))
ch = lambda a, b: abs(a - b)/np.sqrt((1 + abs(a)**2)*(1 + abs(b)**2))
def thin_arc(p, d=1e-3):
    keep = [0]; last = p[0]
    for i in range(1, len(p)):
        if max(ch(p[i][0], last[0]), ch(p[i][1], last[1])) >= d: keep.append(i); last = p[i]
    return [p[i] for i in keep]
ARCS = []; CATS = []                                   # per oval: list of arcs, each an array of (y, x); catalog id of each arc's start
for ov in ovals:
    a0 = 0; arcs = []; cats = []
    for (c0, c1), e in zip(ov['arcs'], ov['arc_end']):
        if e > a0: arcs.append(np.array(thin_arc(ov['pts'][a0:e]))); cats.append(c0)
        a0 = e
    ARCS.append(arcs); CATS.append(cats)
# run 2: anchors at the arc starts (oval_arc_labels.json: slit label at the census loop point the arc starts from, carried round
# the puncture's g-circle from the internal-ray placement), then oval_anchors_wide/, then the continuation anchors below. An arc
# with no anchor is NOT located (run 1 located it by continuity across the puncture, which picked wrong sheets).
AL = json.load(open(os.path.join(here, '..', 'oval_arc_labels.json')))
def wide(oi, ai):
    f = os.path.join(here, 'oval_anchors_wide', f'{oi}_{ai}.json')
    return json.load(open(f)) if os.path.exists(f) else None
print('ovals:', [(len(o['arcs']), len(o['pts']), sum(len(a) for a in A)) for o, A in zip(ovals, ARCS)], '(arcs, traced points, kept)', flush=True)

# ---------------------------------------------------------------- anchors (mesh-independent)
class _TO(Exception): pass
def _al(*a): raise _TO()
signal.signal(signal.SIGALRM, _al)
AP = os.path.join(here, 'oval_anchors.json')
anch = json.load(open(AP)) if os.path.exists(AP) else {}
Lf = None
def anchor(oi, ai, arc):
    global Lf
    key = f'{oi}_{ai}'
    if key in anch: return anch[key]
    if Lf is None: Lf = Locator(here)
    order = sorted(range(len(arc)), key=lambda k: abs(k - len(arc)//2))
    order = [k for k in order if abs(arc[k][0]) < 2 and abs(arc[k][1]) < 2][:12]
    out = None
    for k in order:
        y, x = complex(arc[k][0]), complex(arc[k][1])
        signal.alarm(25)
        try: roots = Lf.roots_X0(y)
        except _TO: roots = None
        signal.alarm(0)
        if roots is None: continue
        dd = sorted((abs(r - x), j) for j, r in enumerate(roots))
        if dd[0][0] < 1e-8*(1 + abs(x)) and dd[1][0] >= 1e3*max(dd[0][0], 1e-300):
            out = dict(k=int(k), label=int(CS.h_at(y)[dd[0][1]]), mismatch=float(dd[0][0]), ratio=float(dd[1][0]/max(dd[0][0], 1e-300))); break
    anch[key] = out; json.dump(anch, open(AP, 'w')); return out
na = 0; nmiss = 0
for oi, arcs in enumerate(ARCS):
    for ai, arc in enumerate(arcs):
        a = anchor(oi, ai, arc); na += a is not None; nmiss += a is None
print(f'anchors: {na} arcs anchored by continuation, {nmiss} without an anchor (located by continuity from the previous arc)', f'| {time.time()-t0:.0f} s', flush=True)

# ---------------------------------------------------------------- location on this mesh
res, D = load_D(DIR); Z, Tn = D['z'], D['Tn']; FG = res['face_g']; Fmp = fu.load_face_mp(D); fold = fu.folder(res)
L = Locator(DIR); n = L.n
GA = np.concatenate([np.eye(2, dtype=complex)[None], FG])                # candidates g: identity, then the face elements
# run 4: faces with a corner at a ramified branch vertex b (ramification e = surface faces at the corner / sphere faces at b > 1)
# are interpolated in the local coordinate u = (y - b)^(1/e) (the developing map is conformal in u, not in y): barycentric
# weights of u(y) in the triangle (0, 1, u_Q), u = ((y - b)/(y_P - b))^(1/e), principal branch. U6_SQRT=0 restores run 3.
SQRT = os.environ.get('U6_SQRT', '1') == '1'
cntS = np.bincount(L.T.ravel(), minlength=len(L.SV)); cntY = np.bincount(L.TS.ravel(), minlength=len(L.Ys))
NSQ = [0]
def WF(t, w, f, y):
    cs = L.T[f]; wf = np.array([w[list(L.TS[t]).index(L.SV[cc, 0])] for cc in cs])
    if not SQRT: return wf
    for i in range(3):
        b = L.SV[cs[i], 0]
        if not L.isB[b]: continue
        e = int(round(cntS[cs[i]]/cntY[b])); yb = complex(L.Ys[b])
        if e < 2 or not np.isfinite(yb) or abs(yb) > 1e6: continue
        iP, iQ = (i + 1) % 3, (i + 2) % 3; yP, yQ = complex(L.Ys[L.SV[cs[iP], 0]]), complex(L.Ys[L.SV[cs[iQ], 0]])
        if L.isB[L.SV[cs[iP], 0]] or L.isB[L.SV[cs[iQ], 0]] or y == yb: continue
        uQ = ((yQ - yb)/(yP - yb))**(1.0/e); u = ((y - yb)/(yP - yb))**(1.0/e)
        wP, wQ = np.linalg.solve(np.array([[1.0, uQ.real], [0.0, uQ.imag]]), np.array([u.real, u.imag]))
        out = np.empty(3); out[i] = 1 - wP - wQ; out[iP] = wP; out[iQ] = wQ; NSQ[0] += 1; return out
    return wf
def sheets(y):
    """folded disk positions of the 42 points of the surface over y, and the slit label of each at the nearest normal corner."""
    t, w = L.face_of(y)
    if t is None: return None, None
    ws = []
    for i in range(n):
        f = n*t + i; wf = WF(t, w, f, y); ws.append(place(Z, Tn, f, wf))
    return fold(np.array(ws))[0], t
# run 3: the sheet of each point is chosen by x on the FINE mesh (vertex_x.npz interpolated in each of the 42 faces over y;
# run 2's nearest-point rule retraced the same sheet at the real fold points of y). Accepted when the best chordal x distance is
# < 1/2 of the runner-up; otherwise (x NaN in the mask / at branch corners, or ambiguous) nearest to the geodesic extrapolation
# of the two previous points among the candidates. The slit label (at y, cut-system frame; mesh-independent) of every chosen
# point is saved to oval_point_labels.json; any other mesh places the points by these labels (no choice made there).
PL_PATH = os.path.join(here, 'oval_point_labels.json'); FINE = DIR == here
VX = np.load(os.path.join(here, 'vertex_x.npz'))['x'] if FINE else None
PLAB = None if FINE else json.load(open(PL_PATH))
def sheets_x(y):
    t, w = L.face_of(y)
    vs = [v for v in L.TS[t] if not L.isB[v]]; v = min(vs, key=lambda vv: abs(L.Ys[vv] - y))
    m = CS.edge_perm([y, complex(L.Ys[v])]); inv = {m[l]: l for l in range(n)}
    ws, xs, labs = [], [], []
    for i in range(n):
        f = n*t + i; wf = WF(t, w, f, y); ws.append(place(Z, Tn, f, wf))
        xs.append(np.dot(wf, VX[L.T[f]])); k = next(cc for cc in L.T[f] if L.SV[cc, 0] == v); labs.append(inv[int(L.SV[k, 1])])
    return fold(np.array(ws))[0], np.array(xs), labs
def gimg(p): return (GA[:, 0, 0]*p + GA[:, 0, 1])/(GA[:, 1, 0]*p + GA[:, 1, 1])
def extrap(p2, p1):
    """geodesic reflection of p2 through p1 (p2 taken in p1's frame: nearest image)."""
    im = gimg(p2); q = im[np.argmin(fu.hdist(p1, im))]
    u = (q - p1)/(1 - np.conj(p1)*q); return (-u + p1)/(1 - np.conj(p1)*u)
# run 5: the sheet is carried by its slit label along the real segment y_prev -> y (CS.edge_perm; exact away from the branch
# points of that sheet). It is given up only for a sheet within disk distance 0.5 of the previous point (modulo faces) whose
# interpolated x is clearly better: dx(carried) > max(3 dx(best near), 1e-4). That is the switch at a real fold of y; sheets
# with the same x to 1e-5 elsewhere on the surface (run 4's flips, 4 to 12 away in the disk) are not candidates.
STAT = dict(carried=0, switched=0, carried_no_x=0)
def choose(prev, prev2, y, x, yprev, labprev, yprev2=None):
    cand, xs, labs = sheets_x(y)
    dx = np.abs(xs - x)/np.sqrt((1 + np.abs(xs)**2)*(1 + abs(x)**2)); dx[~np.isfinite(dx)] = np.inf
    it = labs.index(int(CS.edge_perm([yprev, y])[labprev]))
    img = (GA[:, 0, 0][:, None]*cand[None, :] + GA[:, 0, 1][:, None])/(GA[:, 1, 0][:, None]*cand[None, :] + GA[:, 1, 1][:, None])
    dd = fu.hdist(prev, img).min(0); near = np.where(dd < 0.5)[0]
    # run 7 (run 6's rule lagged at the folds: 155 backward moves on oval 3 = its whole length excess): a real fold of y is where
    # y turns round along the arc. At a turning step (sign(y - yprev) != sign(yprev - yprev2)) the sheet switches to the partner:
    # the near sheet (disk distance < 0.5 of the previous point) other than the carried one that is closest in the disk.
    # Elsewhere the label is carried, and given up only on strong x evidence (run 5: dx(carried) > max(3 dx(best near), 1e-4)).
    turn = yprev2 is not None and (y.real - yprev.real)*(yprev.real - yprev2.real) < 0
    others = near[near != it]
    if turn and len(others):
        pick = int(others[np.argmin(dd[others])]); STAT['turn_switched'] = STAT.get('turn_switched', 0) + 1
        if np.isfinite(dx[pick]) and np.isfinite(dx[it]): STAT['turn_x_agrees'] = STAT.get('turn_x_agrees', 0) + int(dx[pick] <= dx[it])
    else:
        if turn: STAT['turn_no_partner'] = STAT.get('turn_no_partner', 0) + 1
        best = int(near[np.argmin(dx[near])]) if len(near) else it
        if not np.isfinite(dx[it]) or not np.isfinite(dx[best]): STAT['carried_no_x'] += 1; pick = it
        elif best != it and dx[it] > max(3*dx[best], 1e-4): STAT['switched'] += 1; pick = best
        else: STAT['carried'] += 1; pick = it
    return cand[pick], labs[pick], float(dd[pick])
def by_label(y, lab):
    t, w = L.face_of(y)
    vs = [v for v in L.TS[t] if not L.isB[v]]; v = min(vs, key=lambda vv: abs(L.Ys[vv] - y))
    lv = CS.edge_perm([y, complex(L.Ys[v])])[lab]
    for i in range(n):
        f = n*t + i; k = next(cc for cc in L.T[f] if L.SV[cc, 0] == v)
        if L.SV[k, 1] == lv:
            wf = WF(t, w, f, y); return fold(np.array([place(Z, Tn, f, wf)]))[0][0]
    return None
def nearest(prev, cand):
    """index of the candidate nearest to prev modulo {I} u faces, its distance, the runner-up candidate's distance."""
    img = (GA[:, 0, 0][:, None]*cand[None, :] + GA[:, 0, 1][:, None])/(GA[:, 1, 0][:, None]*cand[None, :] + GA[:, 1, 1][:, None])
    d = fu.hdist(prev, img).min(0); o = np.argsort(d)
    return int(o[0]), float(d[o[0]]), float(d[o[1]])
results, SAVE = [], {}; XCHK, XRES, NOANCH = [], [], []; PL_OUT = {}
ONLY = [int(v) for v in os.environ.get('U6_OVALS', '').split(',') if v != '']
for oi, arcs in enumerate(ARCS):
    if ONLY and oi not in ONLY: continue
    F_all = []; margins = []; nfail = 0; nanch = 0; SRC = {}
    prev_end = None
    for ai, arc in enumerate(arcs):
        pos = [None]*len(arc); al = AL.get(str(CATS[oi][ai])); src = None; lab_k = [None]*len(arc)
        if not FINE:
            for k in range(len(arc)):
                lb = PLAB[f'{oi}_{ai}'][k]
                if lb is not None: pos[k] = by_label(complex(arc[k][0]), lb)
                if pos[k] is None: nfail += 1
            F_all += [p for p in pos if p is not None]; nanch += 1; continue
        if al is not None:
            yl = complex(*al['y']); y1 = complex(arc[0][0])
            p0 = by_label(y1, CS.edge_perm([yl, y1])[al['label']])
            if p0 is not None: pos[0] = p0; src = 'start'; lab_k[0] = int(CS.edge_perm([yl, y1])[al['label']])
        for a, nm in ((wide(oi, ai), 'wide'), (anch.get(f'{oi}_{ai}'), 'cont')):
            if a is None: continue
            p0 = by_label(complex(arc[a['k']][0]), a['label'])
            if p0 is None: continue
            if src is None: pos[a['k']] = p0; src = nm; lab_k[a['k']] = int(a['label'])
            else: XCHK.append((oi, ai, nm, a['k'], p0))                    # second anchor: cross-check after continuity
        if src is None: nfail += len(arc); NOANCH.append((oi, ai, len(arc))); continue
        nanch += 1; SRC[src] = SRC.get(src, 0) + 1
        k0 = next(k for k in range(len(arc)) if pos[k] is not None)
        for rng in (range(k0 + 1, len(arc)), range(k0 - 1, -1, -1)):
            prev = pos[k0]; prev2 = None; kp = k0; kp2 = None
            for k in rng:
                p_, lb, mg_ = choose(prev, prev2, complex(arc[k][0]), complex(arc[k][1]), complex(arc[kp][0]), lab_k[kp],
                                     None if kp2 is None else complex(arc[kp2][0]))
                pos[k] = p_; lab_k[k] = lb; prev2 = prev; prev = p_; margins.append(mg_); kp2 = kp; kp = k
        PL_OUT[f'{oi}_{ai}'] = lab_k
        got = [p for p in pos if p is not None]; F_all += got
        for (o_, a_, nm, k, p0) in [x for x in XCHK if x[0] == oi and x[1] == ai]:
            if pos[k] is not None:
                dx = float(fu.hdist(pos[k], (GA[:, 0, 0]*p0 + GA[:, 0, 1])/(GA[:, 1, 0]*p0 + GA[:, 1, 1])).min()); XRES.append((oi, ai, nm, dx))
    F_all = np.array(F_all)
    xr = [x[3] for x in XRES if x[0] == oi]
    print(f'   oval {oi}: anchored arcs {nanch}/{len(arcs)} by source {SRC}; arcs without anchor (oval, arc, points) {[x for x in NOANCH if x[0] == oi]}; '
          f'second-anchor cross-checks {len(xr)}, max distance {max(xr, default=0):.3e}, over 1e-2: {sum(1 for v in xr if v > 1e-2)}', flush=True)
    if len(F_all) < 3: results.append(dict(oval=oi, arcs=len(arcs), anchored_arcs=nanch, points=len(F_all), location_failures=nfail)); continue
    # ---- E94 lifting and axis test
    Pl = list(F_all) + [F_all[0]]; H = (mpc(1), mpc(0)); steps = []; H_list = [H]; runs = [0]; run = 0
    for k in range(1, len(Pl)):
        w = Pl[k]; img = (GA[:, 0, 0]*w + GA[:, 0, 1])/(GA[:, 1, 0]*w + GA[:, 1, 1])
        ds = fu.hdist(Pl[k - 1], img); j = int(np.argmin(ds)); steps.append(float(ds[j]))
        if j > 0: H = fu.mp_mul(H, Fmp[j - 1]); run += 1
        if k < len(Pl) - 1: H_list.append(H); runs.append(run)
    h = H
    trace = float(2*h[0].real); Lg = float(2*gmpy2.acosh(abs(h[0].real))) if abs(h[0].real) > 1 else 0.0; Lp = float(np.sum(steps))
    dists, X1, X2 = [], [], []
    for f, Hk in zip(F_all, H_list):
        c = fu.mp_mul(fu.mp_inv(Hk), fu.mp_mul(h, Hk)); xi1, xi2 = fixed_points_mp(c)
        dists.append(axis_dist(f, xi1, xi2)); X1.append(xi1); X2.append(xi2)
    dists = np.array(dists); mg = np.array(margins if margins else [np.inf])
    r = dict(oval=oi, arcs=len(arcs), anchored_arcs=nanch, xcheck_n=len(xr), xcheck_max=max(xr, default=0.0), xcheck_over_1e2=sum(1 for v in xr if v > 1e-2), points=len(F_all), location_failures=nfail, max_lift_step=float(np.max(steps)),
             steps_over_0p2=int((np.array(steps) > 0.2).sum()), margin_min=float(mg.min()), margin_below_2=int((mg < 2).sum()), trace=trace,
             hyperbolic=bool(abs(h[0].real) > 1), geodesic_length=Lg, polyline_length=Lp, length_ratio_minus_1=Lp/Lg - 1 if Lg else None,
             face_crossings=run, max_dist_to_axis=float(dists.max()), median_dist_to_axis=float(np.median(dists)))
    results.append(r)
    print('U6 oval', oi, json.dumps(r), f'| {time.time()-t0:.0f} s', flush=True)
    SAVE[f'f{oi}'] = F_all; SAVE[f'x1_{oi}'] = np.array(X1); SAVE[f'x2_{oi}'] = np.array(X2); SAVE[f'run{oi}'] = np.array(runs)
print('sheet choices:', STAT, '| face placements in the root coordinate:', NSQ[0], '(U6_SQRT', SQRT, ')', flush=True)
if FINE and not ONLY: json.dump(PL_OUT, open(PL_PATH, 'w'))
if ONLY: DIR = os.path.join(DIR, 'only_' + '_'.join(map(str, ONLY))); os.makedirs(DIR, exist_ok=True); json.dump(PL_OUT, open(os.path.join(DIR, 'oval_point_labels.json'), 'w'))
json.dump(results, open(os.path.join(DIR, 'oval_check.json'), 'w'), indent=1)
np.savez(os.path.join(DIR, 'oval_lifts.npz'), **SAVE)
