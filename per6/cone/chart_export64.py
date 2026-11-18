#!/usr/bin/env python3
"""E92-U viewer data for per6_live.html: an atlas of conformal charts on compact Per_6(0) storing (y, x), so the dynamics
can be computed per screen pixel.
- Chart centres: farthest-point sampling on the surface (distance = min over group elements of the hyperbolic distance
  between orbit points) until every surface point is within R_COV of a centre.
- Chart k: Poincare coordinate u with u = 0 at the centre (u = Mq(c_k)^-1 (p)); grid S x S on [-R_U, R_U]^2,
  R_U = tanh(1.45/2), so every point within 1.45 of the centre is on the grid.
- Values: y = a_3 and x = a_5, stored through a per-chart Mobius change Y = (y - a)/(y - b) (b chosen far, chordally,
  from every value y takes on the chart, a its antipode), so the stored functions are bounded holomorphic functions of u.
  Grid values come from a weighted holomorphic quadratic least-squares fit to the mesh vertices near each grid point.
- Placements: (chart k, group element g) whose centre g(c_k) is nearest to each point of the Dirichlet polygon D; a
  1024^2 chart map over D stores the placement index.
Checks printed: coverage radius, fit residuals at the mesh vertices, and a hold-out test of grid interpolation."""
import base64
import io
import json
import os
import sys
import numpy as np
from PIL import Image
from scipy.spatial import cKDTree
here = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'round0_surface3d'); sys.path.insert(0, here); os.chdir(here)   # E99 fu1: run on round 0
from develop import Mq, mob_ap, hd0, hdist, to_klein, to_poinc
R_COV, R_U_H, S, K = 1.2, 1.45, 160, 16
R_U = np.tanh(R_U_H/2)
D = np.load('develop.npz'); z, Tn, orig = D['z'], D['Tn'], D['orig']; FG = list(D['face_g']); E = D['els_all']; verts = D['verts']
I2 = np.eye(2, dtype=complex)
U_ = np.array([np.exp(1j*np.angle(mob_ap(g, 0))) for g in FG]); Tt = np.array([np.tanh(hd0(mob_ap(g, 0))/2) for g in FG])
Ginv = [np.linalg.inv(g) for g in FG]
rD = float(hd0(to_poinc(verts)).max())
def vfold(w, maxit=200):
    w = w.copy()
    for it in range(maxit):
        k = to_klein(w); v = (k[:, None]*np.conj(U_[None, :])).real - Tt[None, :]
        a = np.argmax(v, 1); out = v[np.arange(len(w)), a] > 1e-13
        if not out.any(): return w
        for f in np.unique(a[out]):
            m = out & (a == f); Gf = Ginv[f]; w[m] = (Gf[0, 0]*w[m] + Gf[0, 1])/(Gf[1, 0]*w[m] + Gf[1, 1])
    raise RuntimeError('vfold')
def inD(w):
    k = to_klein(w); return ((k[:, None]*np.conj(U_[None, :])).real - Tt[None, :] <= 1e-12).all(1)
src = open('surface_data.js').read(); Dd = json.loads(src[src.index('{'):src.rindex('}') + 1])
yv = np.frombuffer(base64.b64decode(Dd['y']), np.float64).reshape(-1, 2); yv = yv[:, 0] + 1j*yv[:, 1]
xv = np.frombuffer(base64.b64decode(Dd['x']), np.float64).reshape(-1, 2); xv = xv[:, 0] + 1j*xv[:, 1]
nV = len(yv)
first = np.full(nV, -1)
for n_, o in enumerate(orig):
    if first[o] < 0: first[o] = n_
P = vfold(z[first])                                  # one point of D per surface vertex
Yv, Xv = yv, xv
print(f'surface vertices {nV}, all folded into D: {bool(inD(P).all())}; Dirichlet circumradius {rD:.4f}', flush=True)
# group elements near the identity (for surface distances)
dE = hd0(E[:, 0, 1]/E[:, 1, 1])
H = E[dE <= 2*rD + 0.6]
print(f'group elements with orbit point within {2*rD + 0.6:.3f}: {len(H)}', flush=True)
def apE(G, w):
    return (G[:, 0, 0][:, None]*w[None, :] + G[:, 0, 1][:, None])/(G[:, 1, 0][:, None]*w[None, :] + G[:, 1, 1][:, None])
# farthest-point sampling of chart centres
rng = np.random.default_rng(3)
samp = P[rng.choice(nV, 30000, replace=False)]
centres = [0j]; dmin = np.full(len(samp), np.inf)
while True:
    c = centres[-1]
    imgs = (H[:, 0, 0]*c + H[:, 0, 1])/(H[:, 1, 0]*c + H[:, 1, 1])
    for s0 in range(0, len(samp), 5000):
        ss = samp[s0:s0 + 5000]
        dmin[s0:s0 + 5000] = np.minimum(dmin[s0:s0 + 5000], hdist(ss[None, :], imgs[:, None]).min(0))
    cov = float(dmin.max())
    if cov <= R_COV or len(centres) >= 48: break
    centres.append(samp[int(np.argmax(dmin))])
centres = np.array(centres)
print(f'chart centres: {len(centres)}; coverage radius (max distance from a sample point to the nearest centre) {cov:.4f}', flush=True)
NC = len(centres)
# placements and chart map over D
size = 1024
yy, xx = np.mgrid[0:size, 0:size]; wpix = (((xx + 0.5)/size*2 - 1) + 1j*(1 - (yy + 0.5)/size*2)).reshape(-1)
m = np.where(np.abs(wpix) < 1)[0]; m = m[inD(wpix[m])]
cand = []
for k, c in enumerate(centres):
    imgs = (E[:, 0, 0]*c + E[:, 0, 1])/(E[:, 1, 0]*c + E[:, 1, 1])
    for j in np.where(hd0(imgs) <= rD + R_COV + 0.3)[0]: cand.append((k, j, imgs[j]))
cpos = np.array([c_[2] for c_ in cand])
best = np.full(len(m), -1); bd = np.full(len(m), np.inf)
for s0 in range(0, len(cpos), 50):
    dd = hdist(wpix[m][None, :], cpos[s0:s0 + 50, None])
    a = np.argmin(dd, 0); v = dd[a, np.arange(len(m))]
    upd = v < bd; best[upd] = s0 + a[upd]; bd[upd] = v[upd]
used = np.unique(best); remap = {int(u): i for i, u in enumerate(used)}
cmap = np.zeros((size*size, 4), np.uint8)
idx = np.array([remap[int(b)] for b in best])
cmap[m, 0] = idx & 255; cmap[m, 1] = idx >> 8; cmap[m, 3] = 255
placements = []
for u in used:
    k, j, _ = cand[int(u)]
    Mp = np.linalg.inv(Mq(centres[k])) @ np.linalg.inv(E[j])           # disk point p -> chart coordinate u
    Mp = Mp/np.sqrt(np.linalg.det(Mp))
    placements.append(dict(chart=int(k), M=[Mp[0, 0].real, Mp[0, 0].imag, Mp[0, 1].real, Mp[0, 1].imag,
                                             Mp[1, 0].real, Mp[1, 0].imag, Mp[1, 1].real, Mp[1, 1].imag]))
print(f'placements used {len(placements)}; texels of D {len(m)}; max distance to the chosen centre {bd.max():.4f} (grid reaches {R_U_H})', flush=True)
# per-chart grids
def chord(a, b): return np.abs(a - b)/np.sqrt((1 + np.abs(a)**2)*(1 + np.abs(b)**2))
fib = np.arange(4000) + 0.5; th = np.arccos(1 - 2*fib/4000); ph = np.pi*(1 + 5**0.5)*fib
sph = np.stack([np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)], 1)
cand_pts = (sph[:, 0] + 1j*sph[:, 1])/(1 - sph[:, 2] + 1e-15)
def far_point(vals):
    vs = vals[np.random.default_rng(5).choice(len(vals), min(4000, len(vals)), replace=False)]
    d = np.array([chord(vs, b).min() for b in cand_pts]); b = cand_pts[int(np.argmax(d))]
    a = -1/np.conj(b) if abs(b) > 1e-12 else 1e12
    return a, b, float(d.max())
ug = np.linspace(-R_U, R_U, S); UG = (ug[None, :] + 1j*ug[::-1, None]).reshape(-1)     # row 0 = top (+Im u)
grids = np.zeros((NC, S, S, 4), np.float64); meta = []; worst_fit = 0.0; ho = []   # fu1: keep float64
for k, c in enumerate(centres):
    imgsE = (E[:, 0, 0]*c + E[:, 0, 1])/(E[:, 1, 0]*c + E[:, 1, 1])
    # tiles g(D) meeting the data disk |u| < 0.8 (radius 2.2 about c) have d(c, g(0)) <= rD + 2.2; complete since
    # d(0, g(0)) <= 2 rD + 2.2 < explored bound
    sel = np.where(hdist(E[:, 0, 1]/E[:, 1, 1], c) <= rD + 2.25)[0]
    assert 2*rD + 2.25 <= float(D['els_bound']), 'explored elements do not cover the chart disks'
    pts, vy, vx = [], [], []
    for j in sel:
        G = np.linalg.inv(Mq(c)) @ E[j]
        w = (G[0, 0]*P + G[0, 1])/(G[1, 0]*P + G[1, 1]); keep = np.abs(w) < 0.8
        pts.append(w[keep]); vy.append(Yv[keep]); vx.append(Xv[keep])
    pts, vy, vx = np.concatenate(pts), np.concatenate(vy), np.concatenate(vx)
    ay, by, sy = far_point(vy); ax, bx, sx = far_point(vx)
    with np.errstate(all='ignore'):
        Yf = (vy - ay)/(vy - by); Xf = (vx - ax)/(vx - bx)
    tree = cKDTree(np.stack([pts.real, pts.imag], 1))
    dd, nb = tree.query(np.stack([UG.real, UG.imag], 1), K)
    sig = 1.5*dd[:, 5:6] + 1e-12
    wgt = np.exp(-(dd/sig)**2)
    dz = pts[nb] - UG[:, None]
    A = np.stack([np.ones_like(dz), dz, dz*dz], -1)                          # holomorphic quadratic
    AhW = np.conj(A).transpose(0, 2, 1)*wgt[:, None, :]
    N_ = AhW @ A
    for col, vals in ((0, Yf), (2, Xf)):
        rhs = (AhW @ vals[nb][:, :, None])[:, :, 0]
        coef = np.linalg.solve(N_, rhs[..., None])[..., 0]   # E99: numpy-2 solve broadcasting (as E93)
        g = coef[:, 0].reshape(S, S)
        grids[k, :, :, col] = g.real; grids[k, :, :, col + 1] = g.imag
        res = np.abs((A @ coef[:, :, None])[:, :, 0] - vals[nb])[:, :6]
        worst_fit = max(worst_fit, float(np.nanmax(res)))
    meta.append(dict(c=[c.real, c.imag], ay=[complex(ay).real, complex(ay).imag], by=[by.real, by.imag],
                     ax=[complex(ax).real, complex(ax).imag], bx=[bx.real, bx.imag], ysep=sy, xsep=sx, npts=int(len(pts))))
    # hold-out: bilinear on the grid at mesh points inside |u| < R_U*0.95, compare y (chordal)
    inside = np.where(np.abs(pts) < R_U*0.95)[0][::97]
    fx = (pts[inside].real + R_U)/(2*R_U)*(S - 1); fy = (R_U - pts[inside].imag)/(2*R_U)*(S - 1)
    i0, j0 = np.floor(fy).astype(int), np.floor(fx).astype(int); ty, tx = fy - i0, fx - j0
    Gy = grids[k, :, :, 0] + 1j*grids[k, :, :, 1]
    val = (Gy[i0, j0]*(1 - tx)*(1 - ty) + Gy[i0, j0 + 1]*tx*(1 - ty) + Gy[i0 + 1, j0]*(1 - tx)*ty + Gy[i0 + 1, j0 + 1]*tx*ty)
    yrec = (ay - by*val)/(1 - val)
    ho.append(np.median(chord(yrec, vy[inside]))); ho_max = float(chord(yrec, vy[inside]).max())
    print(f'chart {k}: centre {c:.4f}, mesh points {len(pts)}, Mobius separation (chordal) y {sy:.3f} x {sx:.3f}, '
          f'hold-out bilinear y error median {ho[-1]:.2e} max {ho_max:.2e}', flush=True)
print(f'fit residual at the 6 nearest mesh points, max over all charts and grid points (stored variables) {worst_fit:.2e}', flush=True)
def evalc(k, u):
    fx = (u.real + R_U)/(2*R_U)*(S - 1); fy = (R_U - u.imag)/(2*R_U)*(S - 1); i0, j0 = int(fy), int(fx); ty, tx = fy - i0, fx - j0
    g = grids[k, i0, j0]*(1-tx)*(1-ty) + grids[k, i0, j0+1]*tx*(1-ty) + grids[k, i0+1, j0]*(1-tx)*ty + grids[k, i0+1, j0+1]*tx*ty
    Y = complex(g[0], g[1]); a, b = complex(*meta[k]['ay']), complex(*meta[k]['by']); return (a - b*Y)/(1 - Y)
PM = [(pl['chart'], np.array([[pl['M'][0] + 1j*pl['M'][1], pl['M'][2] + 1j*pl['M'][3]], [pl['M'][4] + 1j*pl['M'][5], pl['M'][6] + 1j*pl['M'][7]]])) for pl in placements]
errs = []
for v in np.random.default_rng(11).choice(nV, 1500, replace=False):
    for k, Mp in PM:
        u = mob_ap(Mp, P[v])
        if abs(u) < R_U*0.97: errs.append(float(chord(evalc(k, u), yv[v])))
errs = np.array(errs)
print(f'CONSISTENCY: {len(errs)} (mesh vertex, covering placement) evaluations; chordal y error median {np.median(errs):.2e}, '
      f'99th pct {np.percentile(errs, 99):.2e}, max {errs.max():.2e}, fraction > 1e-3 {(errs > 1e-3).mean():.4f}', flush=True)
np.save(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'chart_grids64.npy'), grids)   # fu1: float64 grids; shipped files NOT rewritten
bi = io.BytesIO(); Image.fromarray(cmap.reshape(size, size, 4), 'RGBA').save(bi, 'PNG', optimize=True)
pass   # fu1: chart_map.png not rewritten
coef = json.load(open('../explorer/data.json'))['coef']
json.dump(dict(NC=NC, S=S, RU=R_U, RUH=R_U_H, RCOV=R_COV, coverage=cov, charts=meta, placements=placements, mapSize=size, coef=coef),
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'chart_meta_fu1.json'), 'w'))
print(f'wrote chart_grids.bin ({grids.nbytes/1e6:.1f} MB), chart_map.png ({len(bi.getvalue())/1e6:.2f} MB), chart_meta.json', flush=True)
