#!/usr/bin/env python3
"""E94 copy (n = 8). E93 viewer: details shown when a marker (polynomial centre or puncture) is clicked in per7_live.html.
One entry per marker, in the order of conformal_data.js's markers (63 centres in census order, then 40 punctures).
Centre (the polynomial z^2 + c with a period-7 critical point): label (even numerator k of angle k/127), characteristic angles,
c (25 digits, E90 census), internal address, limb, kneading sequence, real, rabbit, its point (a_3, a_6) and marked orbit, atlas
sheet, G_7 capture degree. Puncture: name (labels of the centres sharing it), k, q, multiplicity m, real, the sharing centres
(label, j, off), the census loop point (b, c) and its (a_3 .. a_6) (the marker sits at the loop point, |g| = eps), atlas sheet.
Writes marker_info.json."""
import json
import os
import re
import sys
import mpmath as mp
here = os.path.dirname(os.path.abspath(__file__))
E90 = os.path.join(here, '..', '..', 'inputs', 'E90-real-circles'); sys.path.insert(0, E90)
import pern  # noqa: E402
census = json.load(open(os.path.join(E90, 'data_n8_eps1e-3_dps60.json')))
data = json.load(open(os.path.join(here, '..', 'explorer', 'data.json')))
AP = json.load(open(os.path.join(here, '..', 'atlas_points.json')))
CF = open(os.path.join(here, 'conformal_data.js')).read(); CF = json.loads(CF[CF.index('{'):CF.rindex('}') + 1])
lab_of = {min(e['ks']): next(v for v in e['ks'] if v % 2 == 0) for e in census['centers']}   # E94: even numerator (no G_8 file)
deg = None
def cs(z, n=10):
    z = complex(z); return f'{z.real:.{n}g} {"-" if z.imag < 0 else "+"} {abs(z.imag):.{n}g}i'
def mpc_str(s):
    m = re.fullmatch(r'\(?\s*([-+]?[0-9.eE+-]+?)\s*([-+])\s*([0-9.eE+-]+)j\s*\)?', s.strip())
    return mp.mpc(mp.mpf(m.group(1)), mp.mpf(m.group(2) + m.group(3)))
mp.mp.dps = 60
def polish(b, c):                                   # as atlas_points.py (E94): Newton on c onto the curve, 60 digits
    for _ in range(60):
        r, rb, rc = pern.R_all(b, c, 8); dc = r/rc; c = c - dc
        if abs(dc) < mp.mpf(10)**-50*(1 + abs(c)): break
    return b, c
items = []
for e, c, a in zip(census['centers'], data['centers'], AP['centers']):
    b_, c_, _ = pern.center_to_chart(mp.mpc(mp.mpf(e['c'][0]), mp.mpf(e['c'][1])))
    orb = pern.orbit(b_, c_, 8)
    k1, k2 = e['ks']
    lines = [f'polynomial z^2 + c, critical point 0 of period 8',
             f'c = {e["c"][0]} {"-" if e["c"][1].startswith("-") else "+"} {e["c"][1].lstrip("-")}i',
             f'characteristic angles {k1}/255, {k2}/255   (label {lab_of[min(e["ks"])]})',
             f'internal address {"-".join(str(t) for t in e["addr"])}, limb {e["limb"][0]}/{e["limb"][1]}, kneading {e["nu"]}',
             f'real: {e["real"]}, rabbit: {e["rabbit"]}',
             f'on Per_8(0): a3 = y = {cs(orb[3])}, a4 = {cs(orb[4])}, a5 = {cs(orb[5])}, a6 = {cs(orb[6])}, a7 = x = {cs(orb[7])}',
             f'atlas sheet {a["sheet"]}']
    items.append(dict(kind='centre', title=f'centre {c["label"]}', lines=lines, y=c['y'], x=c['x'], ks=[int(k1), int(k2)], label=str(c['label'])))
for P, p, a in zip(census['punctures'], data['punctures'], AP['punctures']):
    b_, c_ = (mpc_str(s) for s in P['loop0'])
    b_, c_ = polish(b_, c_)
    orb = pern.orbit(b_, c_, 8)
    shar = [f'{lab_of[min(census["centers"][q["center"]]["ks"])]} (j {q["j"]}, off {q["off"]:+.4f})' for q in P['polys']]
    lines = [f'puncture {p["name"]}: k = {P["k"]}, q = {P["q"]}, multiplicity m = {P["m"]:g}, real: {bool(P["real_pts"])}',
             f'shared by centres: {", ".join(shar)}',
             f'census loop point (|g| = eps): b = {cs(b_, 12)}, c = {cs(c_, 12)}',
             f'at the loop point: a3 = y = {cs(orb[3])}, a4 = {cs(orb[4])}, a5 = {cs(orb[5])}, a6 = {cs(orb[6])}, a7 = x = {cs(orb[7])}',
             f'atlas sheet {a["sheet"]}   (the marker sits at the loop point, next to the puncture)']
    items.append(dict(kind='puncture', title=f'puncture {p["name"]}', lines=lines, y=p['y'], x=p['x'], label=p['name']))
mk = CF['markers']
assert len(mk) == len(items) and all(m['kind'] == it['kind'] and it['title'].endswith(' ' + m['label']) for m, it in zip(mk, items)), 'marker order'
json.dump(dict(items=items), open(os.path.join(here, 'marker_info.json'), 'w'))
print('marker_info.json:', len(items), 'entries;', items[0]['lines'][1], '|', items[-1]['lines'][0])
