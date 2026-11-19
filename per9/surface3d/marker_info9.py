#!/usr/bin/env python3
"""E95 viewer (E94 marker_info.py, n = 9, centres only): details shown when a marker is clicked in per9_live.html. One entry per
marker of conformal_data.js (252 centres, order of explorer/data.json); census entries matched by their angles ks (E90
data_n9_eps1e-3_dps60.json). No puncture entries (A3 not met). Writes marker_info.json and overlay_data.json (no ovals: U6 not run)."""
import json, os, sys
import mpmath as mp
here = os.path.dirname(os.path.abspath(__file__))
E90 = os.path.join(here, '..', '..', 'inputs', 'E90-real-circles'); sys.path.insert(0, E90)
import pern  # noqa: E402
census = json.load(open(os.path.join(E90, 'data_n9_eps1e-3_dps60.json')))
data = json.load(open(os.path.join(here, '..', 'explorer', 'data.json')))
CF = open(os.path.join(here, 'conformal_data.js')).read(); CF = json.loads(CF[CF.index('{'):CF.rindex('}') + 1])
def L(v): return eval(v) if isinstance(v, str) and v[:1] in '[(' else v
byks = {tuple(sorted(int(k) for k in L(e['ks']))): e for e in census['centers']}
def cs(z, n=10):
    z = complex(z); return f'{z.real:.{n}g} {"-" if z.imag < 0 else "+"} {abs(z.imag):.{n}g}i'
mp.mp.dps = 60
items = []; worst = 0.0
for c in data['centers']:
    e = byks[tuple(sorted(c['ks']))]
    cc = L(e['c']); addr = L(e['addr']); limb = L(e['limb'])
    b_, c_, _ = pern.center_to_chart(mp.mpc(mp.mpf(cc[0]), mp.mpf(cc[1])))
    orb = pern.orbit(b_, c_, 9)
    worst = max(worst, abs(complex(orb[3]) - complex(*c['y'])), abs(complex(orb[8]) - complex(*c['x'])))
    k1, k2 = sorted(c['ks'])
    lines = ['polynomial z^2 + c, critical point 0 of period 9',
             f'c = {cc[0]} {"-" if str(cc[1]).startswith("-") else "+"} {str(cc[1]).lstrip("-")}i',
             f'characteristic angles {k1}/511, {k2}/511   (label {c["label"]})',
             f'internal address {"-".join(str(t) for t in addr)}, limb {limb[0]}/{limb[1]}, kneading {e["nu"]}',
             f'real: {e["real"]}, rabbit: {e["rabbit"]}',
             f'on Per_9(0): a3 = y = {cs(orb[3])}, a4 = {cs(orb[4])}, a5 = {cs(orb[5])}, a6 = {cs(orb[6])}, a7 = {cs(orb[7])}, a8 = x = {cs(orb[8])}']
    items.append(dict(kind='centre', title=f'centre {c["label"]}', lines=lines, y=c['y'], x=c['x'], ks=[int(k1), int(k2)], label=str(c['label'])))
mk = CF['markers']
assert len(mk) == len(items) and all(m['kind'] == it['kind'] and it['title'].endswith(' ' + m['label']) for m, it in zip(mk, items)), 'marker order'
json.dump(dict(items=items), open(os.path.join(here, 'marker_info.json'), 'w'))
json.dump(dict(ovals=[], edges=[]), open(os.path.join(here, 'overlay_data.json'), 'w'))
print('marker_info.json:', len(items), 'entries; max |census orbit point - data.json (y, x)|', f'{worst:.2e}', '|', items[0]['lines'][1], '|', items[-1]['lines'][2])
print('overlay_data.json: ovals 0, edges 0')
