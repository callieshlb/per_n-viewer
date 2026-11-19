#!/usr/bin/env python3
"""E93-U viewer data for per7_live.html (E92 chart_export.py, element sets rebuilt for genus 24): an atlas of conformal charts on
compact Per_7(0) storing (y, x) = (a_3, a_6), so the dynamics can be computed per screen pixel.
- Chart centres: farthest-point sampling over 20 000 mesh vertices folded into D; distance = min over the elements
  H = {I} u F u F.F (F = faces of D; products in 200 bits) of d(s, g c). This is an UPPER bound of the surface distance, so the
  printed coverage radius is an upper bound of the true one (safe for the stop rule coverage <= 1.2).
- Chart k: Poincare coordinate u with u = 0 at the centre; grid S x S on [-R_U, R_U]^2, R_U = tanh(1.45/2).
- Data: every lift of every mesh vertex within |u| < 0.8 (hyperbolic 2.2) of the centre: the tiles g D meeting that disk are
  found by BFS over face neighbours, keeping tiles whose distance lower bound (max over faces of the signed distance of
  g^-1(c) to the face line) is <= 2.3 (tiles meeting the disk are edge-connected; the kept set contains them).
- Values: y and x stored through a per-chart Mobius change Y = (y - a)/(y - b) (bounded holomorphic functions of u), grid values
  by weighted holomorphic quadratic least squares on the 16 nearest mesh lifts (E92).
- Placements: (chart k, g in H) whose centre image g(c_k) is nearest to each texel of D; 1024^2 chart map over D.
Checks printed: coverage radius, hold-out interpolation error, and the CONSISTENCY check C (registered): chordal y error of the
chart value vs the mesh vertex value over (mesh vertex, covering placement) pairs; median < 1e-4 and fraction > 1e-3 < 1e-3."""
import base64
import io
import json
import os
import sys
import numpy as np
from PIL import Image
from scipy.spatial import cKDTree
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import fuchs as fu
from fuchs import Mq, mob_ap, hd0, hdist, to_klein, to_poinc
R_COV, R_U_H, S, K = 1.2, 1.45, 160, 16
R_U = np.tanh(R_U_H/2)
D = np.load('develop.npz'); z, Tn, orig = D['z'], D['Tn'], D['orig']; FG = D['face_g']; verts = D['verts']
res = dict(verts=verts, face_g=FG, bounded=True)
U_, Tt, _ = fu.halfplanes(FG)
fold = fu.folder(res)
rD = float(hd0(to_poinc(verts)).max())
SD = json.load(open('surface_data.json'))
yv = np.frombuffer(base64.b64decode(SD['y']), np.float64).reshape(-1, 2); yv = yv[:, 0] + 1j*yv[:, 1]
xv = np.frombuffer(base64.b64decode(SD['x']), np.float64).reshape(-1, 2); xv = xv[:, 0] + 1j*xv[:, 1]
nV = len(yv)
first = np.full(nV, -1)
first[orig[::-1]] = np.arange(len(orig))[::-1]
P = fold(z[first])[0]
def inD(w):
    w = np.asarray(w); out = np.zeros(len(w), bool)                 # E94: chunked (845k points x 822 faces)
    for s0 in range(0, len(w), 40000):
        k = to_klein(w[s0:s0 + 40000]); out[s0:s0 + 40000] = ((k[:, None]*np.conj(U_[None, :])).real - Tt[None, :] <= 1e-12).all(1)
    return out
print(f'surface vertices {nV}, all folded into D: {bool(inD(P).all())}; Dirichlet circumradius {rD:.4f}', flush=True)
# H = {I} u F u F.F in 200 bits
Fmp = fu.load_face_mp(D)
Fall = fu.mp_dedupe(Fmp + [fu.mp_inv(e) for e in Fmp])
prods = {}
for e in Fall:
    for f in Fall:
        h = fu.mp_mul(e, f)
        if abs(h[0]) < fu.gmpy2.cosh(fu.mpfr('5e-4')): continue
        prods.setdefault(fu.mp_key(h), h)
for e in Fall: prods.setdefault(fu.mp_key(e), e)
H = np.array([np.eye(2, dtype=complex)] + [fu.mp_to_np(e) for e in prods.values()])
dH = fu.halfplanes(H[1:])[2]
print(f'H = I u F u F.F: {len(H)} elements, displacement up to {dH.max():.3f}', flush=True)
def images(Hs, c):
    return (Hs[:, 0, 0]*c + Hs[:, 0, 1])/(Hs[:, 1, 0]*c + Hs[:, 1, 1])
rng = np.random.default_rng(3)
samp = P[rng.choice(nV, 20000, replace=False)]
centres = [0j]; dmin = np.full(len(samp), np.inf)
while True:
    c = centres[-1]
    im = images(H, c); im = im[hd0(im) <= rD + 2.0]
    for s0 in range(0, len(samp), 4000):
        ss = samp[s0:s0 + 4000]
        dmin[s0:s0 + 4000] = np.minimum(dmin[s0:s0 + 4000], hdist(ss[None, :], im[:, None]).min(0))
    cov = float(dmin.max())
    if cov <= R_COV or len(centres) >= 1000: break
    centres.append(samp[int(np.argmax(dmin))])
centres = np.array(centres); NC = len(centres)
print(f'chart centres: {NC}; coverage radius (upper bound: max over samples of the distance to the nearest centre) {cov:.4f}', flush=True)
# placements and chart map over D
size = 1024
yy, xx = np.mgrid[0:size, 0:size]; wpix = (((xx + 0.5)/size*2 - 1) + 1j*(1 - (yy + 0.5)/size*2)).reshape(-1)
m = np.where(np.abs(wpix) < 1)[0]; m = m[inD(wpix[m])]
cand = []
for k, c in enumerate(centres):
    im = images(H, c)
    for j in np.where(hd0(im) <= rD + R_COV + 0.3)[0]: cand.append((k, j, im[j]))
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
    Mp = np.linalg.inv(Mq(centres[k])) @ np.linalg.inv(H[j])
    Mp = Mp/np.sqrt(np.linalg.det(Mp))
    placements.append(dict(chart=int(k), M=[Mp[0, 0].real, Mp[0, 0].imag, Mp[0, 1].real, Mp[0, 1].imag,
                                             Mp[1, 0].real, Mp[1, 0].imag, Mp[1, 1].real, Mp[1, 1].imag]))
print(f'placements used {len(placements)}; texels of D {len(m)}; max distance to the chosen centre {bd.max():.4f} (grid reaches {R_U_H})', flush=True)
# tiles near a chart centre
def face_lower_bound(w):
    k = to_klein(w); den = np.sqrt((1 - np.abs(k)**2)[:, None]*(1 - Tt**2)[None, :])
    s = np.arcsinh(((k[:, None]*np.conj(U_[None, :])).real - Tt[None, :])/den)
    return np.maximum(s.max(1), 0)
def tiles_near(c, R=2.3, cap=5000):
    I2 = (fu.mpc(1), fu.mpc(0))
    seen = {fu.mp_key(I2): I2}; front = [I2]; keep = [I2]
    while front and len(seen) < cap:
        nxt = []
        for g in front:
            for f in Fmp:
                h = fu.mp_mul(g, f); kk = fu.mp_key(h)
                if kk in seen: continue
                seen[kk] = h
                hn = fu.mp_to_np(h)
                w = mob_ap(np.linalg.inv(hn), c)
                if face_lower_bound(np.array([w]))[0] <= R:
                    keep.append(h); nxt.append(h)
        front = nxt
    return [fu.mp_to_np(h) for h in keep]
def chord(a, b): return np.abs(a - b)/np.sqrt((1 + np.abs(a)**2)*(1 + np.abs(b)**2))
fib = np.arange(4000) + 0.5; th = np.arccos(1 - 2*fib/4000); ph = np.pi*(1 + 5**0.5)*fib
sph = np.stack([np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)], 1)
cand_pts = (sph[:, 0] + 1j*sph[:, 1])/(1 - sph[:, 2] + 1e-15)
def far_point(vals):
    vs = vals[np.random.default_rng(5).choice(len(vals), min(4000, len(vals)), replace=False)]
    d = np.array([chord(vs, b).min() for b in cand_pts]); b = cand_pts[int(np.argmax(d))]
    a = -1/np.conj(b) if abs(b) > 1e-12 else 1e12
    return a, b, float(d.max())
ug = np.linspace(-R_U, R_U, S); UG = (ug[None, :] + 1j*ug[::-1, None]).reshape(-1)
grids = np.zeros((NC, S, S, 4), np.float32); meta = []; ho = []; ntiles = []
for k, c in enumerate(centres):
    tl = tiles_near(c); ntiles.append(len(tl))
    Mc = np.linalg.inv(Mq(c))
    pts, vy, vx = [], [], []
    for g in tl:
        G = Mc @ g
        w = (G[0, 0]*P + G[0, 1])/(G[1, 0]*P + G[1, 1]); keep = np.abs(w) < 0.8
        pts.append(w[keep]); vy.append(yv[keep]); vx.append(xv[keep])
    pts, vy, vx = np.concatenate(pts), np.concatenate(vy), np.concatenate(vx)
    ay, by, sy = far_point(vy); ax, bx, sx = far_point(vx)
    with np.errstate(all='ignore'):
        Yf = (vy - ay)/(vy - by); Xf = (vx - ax)/(vx - bx)
    tree = cKDTree(np.stack([pts.real, pts.imag], 1))
    dd, nb = tree.query(np.stack([UG.real, UG.imag], 1), K)
    sig = 1.5*dd[:, 5:6] + 1e-12
    wgt = np.exp(-(dd/sig)**2)
    dz = pts[nb] - UG[:, None]
    A = np.stack([np.ones_like(dz), dz, dz*dz], -1)
    AhW = np.conj(A).transpose(0, 2, 1)*wgt[:, None, :]
    N_ = AhW @ A
    for col, vals in ((0, Yf), (2, Xf)):
        rhs = (AhW @ vals[nb][:, :, None])[:, :, 0]
        coef = np.linalg.solve(N_, rhs[..., None])[..., 0]          # numpy 2: batched right-hand side needs a trailing axis
        g_ = coef[:, 0].reshape(S, S)
        grids[k, :, :, col] = g_.real; grids[k, :, :, col + 1] = g_.imag
    meta.append(dict(c=[c.real, c.imag], ay=[complex(ay).real, complex(ay).imag], by=[by.real, by.imag],
                     ax=[complex(ax).real, complex(ax).imag], bx=[bx.real, bx.imag], ysep=sy, xsep=sx, npts=int(len(pts)), tiles=len(tl)))
    inside = np.where(np.abs(pts) < R_U*0.95)[0][::97]
    fx = (pts[inside].real + R_U)/(2*R_U)*(S - 1); fy = (R_U - pts[inside].imag)/(2*R_U)*(S - 1)
    i0, j0 = np.floor(fy).astype(int), np.floor(fx).astype(int); ty, tx = fy - i0, fx - j0
    Gy = grids[k, :, :, 0] + 1j*grids[k, :, :, 1]
    val = (Gy[i0, j0]*(1 - tx)*(1 - ty) + Gy[i0, j0 + 1]*tx*(1 - ty) + Gy[i0 + 1, j0]*(1 - tx)*ty + Gy[i0 + 1, j0 + 1]*tx*ty)
    yrec = (ay - by*val)/(1 - val)
    ho.append(float(np.median(chord(yrec, vy[inside]))))
    if k % 10 == 0 or k == NC - 1:
        print(f'chart {k}: centre {c:.4f}, tiles {len(tl)}, mesh lifts {len(pts)}, Mobius separation y {sy:.3f} x {sx:.3f}, '
              f'hold-out bilinear y error median {ho[-1]:.2e}', flush=True)
print(f'tiles per chart {min(ntiles)} .. {max(ntiles)}; hold-out median over charts: max {max(ho):.2e}', flush=True)
def cw(t):
    t2 = t*t; t3 = t2*t; return np.array([-.5*t3 + t2 - .5*t, 1.5*t3 - 2.5*t2 + 1, -1.5*t3 + 2*t2 + .5*t, .5*t3 - .5*t2])
def evalc(k, u):
    """bicubic (Catmull-Rom) read as in the viewer, then decode y."""
    fx = (u.real + R_U)/(2*R_U)*(S - 1); fy = (R_U - u.imag)/(2*R_U)*(S - 1); ix, iy = int(np.floor(fx)), int(np.floor(fy))
    wx, wy = cw(fx - ix), cw(fy - iy); g = np.zeros(4)
    for j in range(4):
        for i in range(4):
            g += wx[i]*wy[j]*grids[k, min(S - 1, max(0, iy + j - 1)), min(S - 1, max(0, ix + i - 1))]
    Y = complex(g[0], g[1]); a, b = complex(*meta[k]['ay']), complex(*meta[k]['by']); return (a - b*Y)/(1 - Y)
PM = [(pl['chart'], np.array([[pl['M'][0] + 1j*pl['M'][1], pl['M'][2] + 1j*pl['M'][3]], [pl['M'][4] + 1j*pl['M'][5], pl['M'][6] + 1j*pl['M'][7]]])) for pl in placements]
def evalc_bl(k, u):
    """bilinear read (E92 chart_export.py's consistency check)."""
    fx = (u.real + R_U)/(2*R_U)*(S - 1); fy = (R_U - u.imag)/(2*R_U)*(S - 1); i0, j0 = int(fy), int(fx); ty, tx = fy - i0, fx - j0
    g = grids[k, i0, j0]*(1-tx)*(1-ty) + grids[k, i0, j0+1]*tx*(1-ty) + grids[k, i0+1, j0]*(1-tx)*ty + grids[k, i0+1, j0+1]*tx*ty
    Y = complex(g[0], g[1]); a, b = complex(*meta[k]['ay']), complex(*meta[k]['by']); return (a - b*Y)/(1 - Y)
errs, errs_bc = [], []
for v in np.random.default_rng(11).choice(nV, 1500, replace=False):
    for k, Mp in PM:
        u = mob_ap(Mp, P[v])
        if abs(u) < R_U*0.97:
            errs.append(float(chord(evalc_bl(k, u), yv[v]))); errs_bc.append(float(chord(evalc(k, u), yv[v])))
errs, errs_bc = np.array(errs), np.array(errs_bc)
print(f'CONSISTENCY (registered, bilinear as E92): {len(errs)} (mesh vertex, covering placement) evaluations; chordal y error median '
      f'{np.median(errs):.2e}, 99th pct {np.percentile(errs, 99):.2e}, max {errs.max():.2e}, fraction > 1e-3 {(errs > 1e-3).mean():.4f}', flush=True)
print(f'consistency with the viewer\'s bicubic read: median {np.median(errs_bc):.2e}, 99th pct {np.percentile(errs_bc, 99):.2e}, '
      f'max {errs_bc.max():.2e}, fraction > 1e-3 {(errs_bc > 1e-3).mean():.4f}', flush=True)
# E94: split into parts < 95 MB (the viewer concatenates chart_grids_0.bin, chart_grids_1.bin, ...)
raw = grids.astype('<f4').tobytes(); per = S*S*4*4; cpp = max(1, (60_000_000)//per)   # 60 MB raw: base64 copies for the offline page stay < 95 MB
parts = []
for i, k0 in enumerate(range(0, NC, cpp)):
    open(f'chart_grids_{i}.bin', 'wb').write(raw[k0*per:min(NC, k0 + cpp)*per]); parts.append(f'chart_grids_{i}.bin')
print('chart grid parts', parts, flush=True)
bi = io.BytesIO(); Image.fromarray(cmap.reshape(size, size, 4), 'RGBA').save(bi, 'PNG', optimize=True)
open('chart_map.png', 'wb').write(bi.getvalue())
coef = json.load(open(os.path.join(here, '..', 'explorer', 'data.json')))['coef']
json.dump(dict(parts=parts, NC=NC, S=S, RU=R_U, RUH=R_U_H, RCOV=R_COV, coverage=cov, charts=meta, placements=placements, mapSize=size, coef=coef,
               leadQ=[0, 0, 0, 0, 1], leadPow=6, n=8,
               consistency=dict(n=len(errs), median=float(np.median(errs)), p99=float(np.percentile(errs, 99)), max=float(errs.max()),
                                frac_gt_1e3=float((errs > 1e-3).mean()))), open('chart_meta.json', 'w'))
print(f'wrote {len(parts)} chart grid parts ({grids.nbytes/1e6:.1f} MB), chart_map.png ({len(bi.getvalue())/1e6:.2f} MB), chart_meta.json', flush=True)
