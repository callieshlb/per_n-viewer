#!/usr/bin/env python3
"""E95: puncture markers for the viewer from puncture_place.json (internal-ray placement). Each placed puncture's loop point
(y, slit label) is located on the FINE mesh without root-finding: sphere face of y, the label carried to a normal corner by cut
crossings (as meshloc.Locator.locate, with the label known), Klein-barycentric placement in the developed face, fold into D.
Appends 'puncture' markers to surface3d/conformal_data.js, items to surface3d/marker_info.json and 'punctures' to
explorer/data.json (idempotent: existing puncture entries are replaced). Marker name = labels of the centres sharing the puncture
(E94 convention). Raw checks printed: located count, distinct disk positions, distance from each marker to the folded centre it
was reached from (must be < the ray's hyperbolic length; printed range only)."""
import json, os, sys
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__)); S3 = os.path.join(HERE, 'surface3d'); sys.path.insert(0, HERE); sys.path.insert(0, S3)
E90 = os.path.join(HERE, '..', 'inputs', 'E90-real-circles'); sys.path.insert(0, E90)
import fuchs as fu
from meshloc import Locator, place, load_D, CS
census = json.load(open(os.path.join(E90, 'data_n9_eps1e-3_dps60.json')))
PP = json.load(open(os.path.join(HERE, 'puncture_place.json')))
def ev(v): return eval(v) if isinstance(v, str) else v
cen = census['centers']
lab = lambda i: next((v for v in ev(cen[i]['ks']) if v % 2 == 0), min(ev(cen[i]['ks'])))
L = Locator(S3); res, D = load_D(S3); Z, Tn = D['z'], D['Tn']; fold = fu.folder(res); n = L.n
def cs(z, k=10):
    z = complex(z); return f'{z.real:.{k}g} {"-" if z.imag < 0 else "+"} {abs(z.imag):.{k}g}i'
mk, items, dj, fail = [], [], [], []
for r in PP:
    P = census['punctures'][r['i']]; polys = ev(P['polys'])
    if r['status'] != 'ok': fail.append((r['i'], r.get('why', '')[:30])); continue
    y = complex(*r['y']); x = complex(*r['x'])
    t, w = L.face_of(y)
    vs = [v for v in L.TS[t] if not L.isB[v]]; v = min(vs, key=lambda vv: abs(L.Ys[vv] - y))
    lv = CS.edge_perm([y, complex(L.Ys[v])])[r['label']]
    f = None
    for i in range(n):
        ff = n*t + i; c = L.T[ff]; k = next(cc for cc in c if L.SV[cc, 0] == v)
        if L.SV[k, 1] == lv: f = ff; break
    if f is None: fail.append((r['i'], 'no face')); continue
    wf = np.array([w[list(L.TS[t]).index(L.SV[cc, 0])] for cc in L.T[f]])
    p = fold(np.array([place(Z, Tn, f, wf)]))[0][0]
    names = sorted(lab(int(q['center'])) for q in polys); name = '/'.join(str(v) for v in names)
    q_, k_, m_ = int(ev(P['q'])), int(ev(P['k'])), float(ev(P['m'])); real = bool(ev(P['real_pts']))
    c0 = cen[r['centre']]; p_, qq = ev(c0['limb'])
    shar = [f'{lab(int(q["center"]))} (limb {ev(cen[int(q["center"])]["limb"])[0]}/{ev(cen[int(q["center"])]["limb"])[1]}, internal address {"-".join(str(t_) for t_ in ev(cen[int(q["center"])]["addr"]))}, j {q["j"]}, off {float(q["off"]):+.4f})' for q in polys]
    lines = [f'puncture {name}: k = {k_}, q = {q_}, multiplicity m = {m_:g}, real: {real}',
             f'shared by {len(polys)} centre(s): ' + '; '.join(shar),
             f'boundary point of multiplier e^(-2 pi i {p_}/{qq}) of the hyperbolic component of centre {lab(r["centre"])} (placed along its internal ray, {r["steps"]} steps)',
             f'census loop point (|g| = eps = {census["eps"]}): b = {r["b"][:22]}, c = {r["c"][:22]}',
             f'at the loop point: a3 = y = {cs(y)}, a8 = x = {cs(x)}',
             (f'the loop point lies {r["dist_b"]:.2e} from the branch value #{r["near_b"]} of y; sheet label cycle {r["group"]} there' if r['near_b'] is not None
              else f'sheet (continuation index) {r["sheet"]}; nearest branch value of y at distance {r["dist_b"]:.3f}'),
             'the marker sits at the loop point, next to the puncture']
    mk.append(dict(kind='puncture', label=name, p=[float(p.real), float(p.imag)], y=r['y'], x=r['x']))
    items.append(dict(kind='puncture', title=f'puncture {name}', lines=lines, y=r['y'], x=r['x'], label=name, m=m_, q=q_, k=k_, real=real,
                      centres=names, census_index=r['i']))
    dj.append(dict(x=r['x'], y=r['y'], sheet=r['sheet'], name=name, real=real, q=q_, m=m_, centres=names))
P_ = np.array([complex(*m['p']) for m in mk]); dmin = min(fu.hdist(P_[i:i+1], np.delete(P_, i)).min() for i in range(len(P_)))
names_all = [m['label'] for m in mk]
print(f'puncture markers: {len(mk)} of {len(PP)} (not placed / not located: {fail}); all inside the unit disk {bool((np.abs(P_) < 1).all())}; '
      f'min hyperbolic distance between two markers {dmin:.4f}; distinct names {len(set(names_all))}; real {sum(i["real"] for i in items)}; '
      f'multiplicity histogram {sorted({m: sum(1 for i in items if i["m"] == m) for m in set(i["m"] for i in items)}.items())}')
# ---- write (idempotent)
cf = open(os.path.join(S3, 'conformal_data.js')).read(); CF = json.loads(cf[cf.index('{'):cf.rindex('}') + 1])
CF['markers'] = [m for m in CF['markers'] if m['kind'] != 'puncture'] + mk
open(os.path.join(S3, 'conformal_data.js'), 'w').write('window.CF = ' + json.dumps(CF) + ';\n')
MI = json.load(open(os.path.join(S3, 'marker_info.json')))
MI['items'] = [i for i in MI['items'] if i['kind'] != 'puncture'] + items
json.dump(MI, open(os.path.join(S3, 'marker_info.json'), 'w'))
assert len(CF['markers']) == len(MI['items']) and all(a['kind'] == b['kind'] and a['label'] == b['label'] for a, b in zip(CF['markers'], MI['items']))
DJ = json.load(open(os.path.join(HERE, 'explorer', 'data.json'))); DJ['punctures'] = dj
json.dump(DJ, open(os.path.join(HERE, 'explorer', 'data.json'), 'w'))
print(f'wrote conformal_data.js ({len(CF["markers"])} markers), marker_info.json ({len(MI["items"])} items), explorer/data.json ({len(dj)} punctures)')
