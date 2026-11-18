#!/usr/bin/env python3
"""E99 fu1 (REGISTRATION.md, 6/58 implementation (i)): build fu1/site/ = the round-0 Per_6 site files plus re-normalised chart
layers for the 6/58 cone. Appends to chart_grids.bin one copy of every chart served (filled chart map) within hyperbolic
R_SERVE of the cone point, stored as Y' = (y - y_c)/(y - b_y), X' = (x - x_c)/(x - b_x) from the float64 grids; existing
layers byte-identical. chart_meta.json: NC += m, charts += m entries (ay = y_c, ax = x_c, same b), and
cones = [{label, p, R, yc, xc, layer (chart -> new layer, -1 if none), Q (6 x 8, coefficient of dx^k dy^j)}]."""
import base64, json, os, shutil, sys
import numpy as np, sympy as sp
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__)); R0 = os.path.join(HERE, '..', 'surface3d'); sys.path.insert(0, R0)
from develop import mob_ap, to_klein, hdist
OUT = os.path.join(HERE, 'site'); os.makedirs(OUT, exist_ok=True)
CONES = [dict(label='6/58', yc=0.5, xc=0.5, R=0.15)]
R_SERVE = 0.17
meta = json.load(open(f'{R0}/chart_meta.json')); S, NC, MS = meta['S'], meta['NC'], meta['mapSize']
G32 = np.fromfile(f'{R0}/chart_grids.bin', '<f4').reshape(NC, S, S, 4); G64 = np.load(f'{HERE}/chart_grids64.npy')
assert (G64.astype(np.float32) == G32).all()
cmap = np.array(Image.open(f'{R0}/chart_map.png').convert('RGBA'))
D = np.load(f'{R0}/develop.npz'); z = D['z']; orig = D['orig']; FG = list(D['face_g'])
U_ = np.array([np.exp(1j*np.angle(mob_ap(g, 0))) for g in FG]); Tt = np.array([abs(mob_ap(g, 0)) for g in FG]); Gi = [np.linalg.inv(g) for g in FG]
def fold(w):
    for _ in range(1000):
        v = (to_klein(np.array([w]))[0]*np.conj(U_)).real - Tt; a = int(np.argmax(v))
        if v[a] <= 1e-13: return w
        w = mob_ap(Gi[a], w)
s = open(f'{R0}/surface_data.js').read(); SD = json.loads(s[s.index('{'):s.rindex('}') + 1])
yv = np.frombuffer(base64.b64decode(SD['y']), np.float64).reshape(-1, 2); yv = yv[:, 0] + 1j*yv[:, 1]
xv = np.frombuffer(base64.b64decode(SD['x']), np.float64).reshape(-1, 2); xv = xv[:, 0] + 1j*xv[:, 1]
coef = meta['coef']; X_, Y_ = sp.symbols('x y')
P7 = sum(sp.Integer(int(c))*X_**k*Y_**j for k, row in enumerate(coef) for j, c in enumerate(row))
yy = np.mgrid[0:MS, 0:MS]; pix = ((yy[1] + .5)/MS*2 - 1) + 1j*(1 - (yy[0] + .5)/MS*2)
newG, newC, cones = [], [], []
for c in CONES:
    yc, xc = c['yc'], c['xc']
    cand = np.where((yv == yc) & (np.abs(xv - xc) < 1e-2))[0]          # the cone (branch) vertex: y exactly y_c, x ~ x_c
    assert len(cand) == 1, cand
    dv = np.where(orig == cand[0])[0]; ps = sorted({complex(np.round(fold(z[n]), 12)) for n in dv}, key=abs)
    p = ps[0]
    sel = (np.abs(pix) < 1) & (hdist(pix, p) < R_SERVE) & (cmap[..., 3] == 255)
    pl = np.unique(cmap[..., 0][sel].astype(int) + 256*cmap[..., 1][sel].astype(int))
    charts = sorted({meta['placements'][i]['chart'] for i in pl})
    layer = [-1]*NC
    for k in charts:
        ch = meta['charts'][k]; by, bx = complex(*ch['by']), complex(*ch['bx'])
        Y = G64[k, :, :, 0] + 1j*G64[k, :, :, 1]; X = G64[k, :, :, 2] + 1j*G64[k, :, :, 3]
        y = (complex(*ch['ay']) - by*Y)/(1 - Y); x = (complex(*ch['ax']) - bx*X)/(1 - X)
        YA = (y - yc)/(y - by); XA = (x - xc)/(x - bx)
        newG.append(np.stack([YA.real, YA.imag, XA.real, XA.imag], -1)); layer[k] = NC + len(newC)
        newC.append(dict(ch, ay=[yc, 0.0], ax=[xc, 0.0], cone=c['label'], base_chart=k))
    Q = sp.Poly(sp.expand(P7.subs({X_: X_ + sp.nsimplify(xc), Y_: Y_ + sp.nsimplify(yc)})), X_, Y_); Qa = np.zeros((6, 8))
    for (i, j), v in Q.terms(): Qa[i, j] = float(v)
    print(f"cone {c['label']}: cone vertex surface id {int(cand[0])}, disk copies {len(dv)}, point in D {p:.12f} (copies folded: {len(ps)}); "
          f"placements served within {R_SERVE}: {len(pl)}, charts {charts}; Q nonzero terms {int((Qa != 0).sum())}, Q[0][0] {Qa[0][0]}")
    cones.append(dict(label=c['label'], p=[p.real, p.imag], R=c['R'], yc=yc, xc=xc, layer=layer, Q=Qa.reshape(-1).tolist()))
G = np.concatenate([G32, np.array(newG).astype(np.float32)])
assert (G[:NC] == G32).all()
G.astype('<f4').tofile(f'{OUT}/chart_grids.bin')
meta2 = dict(meta, NC=NC + len(newC), charts=meta['charts'] + newC, cones=cones)
json.dump(meta2, open(f'{OUT}/chart_meta.json', 'w'))
for f in ('conformal_data.js', 'overlay_data.json', 'chart_map.png'): shutil.copy(f'{R0}/{f}', f'{OUT}/{f}')
print(f'wrote site/: layers {NC} -> {NC + len(newC)}, chart_grids.bin {os.path.getsize(f"{OUT}/chart_grids.bin")/1e6:.1f} MB; first {NC} layers byte-identical')
