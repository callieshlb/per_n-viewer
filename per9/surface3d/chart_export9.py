#!/usr/bin/env python3
"""E95 viewer data for per9_live.html: E94's chart_export.py rebuilt for genus 188 (same chart definition, grid layout, Mobius
value encoding, least-squares fit, checks and file formats; the element sets come from the tiling BFS of u5_bfs.py instead of
H = I u F u F.F, which would be ~2e7 elements here).
Inputs: develop.npz (D, 200-bit faces), P_folded.npy (fold_vertices9.py: every surface vertex folded into D), vertex_x.npz
(vertex_x.py: x, y per surface vertex; x NaN inside the masked disc |y| < 0.1).
- tiles_near(c, R): the tiles g D meeting the ball B(c, R) (BFS over face neighbours, crossing a face only when its hyperbolic
  distance from g^-1 c is < R; 200-bit words, exact D geometry of u5_bfs.Tiling).
- Chart centres: farthest-point sampling over NS mesh vertices (plus the 252 centres); distance sample -> centre = min over
  tiles_near(c, 2.7) of d(g s, c) (the surface distance whenever it is < 2.7; an upper bound otherwise). Stop at coverage 1.2.
- Chart k: Poincare coordinate u = Mq(c_k)^-1 (.), grid S x S on [-R_U, R_U]^2, R_U = tanh(1.45/2) (E94). Data: every lift of
  every mesh vertex with |u| < 0.8, from tiles_near(c_k, 2.3). Values as E94: per-chart Mobius change of y and of x, weighted
  holomorphic quadratic least squares on the 16 nearest lifts (x: lifts with finite x only).
- Neighbour table (new): for chart k every centre lift within 2.6 of c_k: (j, T) with u_j = T u_k; the viewer walks from the
  chart-map guess to the nearest centre (the 1024^2 chart map cannot resolve D beyond hyperbolic radius ~5.5).
- Chart map: nearest centre lift per texel of D (hyperboloid inner products), other texels filled from the nearest assigned one.
- Check C (registered, E94 definition): chordal y error of the chart value (bilinear) vs the mesh vertex value over (vertex,
  covering placement) pairs, 1500 random vertices; median < 1e-4 and fraction > 1e-3 below 1e-3. Exit 0 iff C passes.
Stages are cached (chart9_centres.pkl, chart9_grids.npy + chart9_done.json) so a restart resumes. Raw values only."""
import json, os, pickle, sys, time
import numpy as np
import gmpy2
from gmpy2 import mpc
from PIL import Image
from scipy.spatial import cKDTree
from scipy.ndimage import distance_transform_edt
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import fuchs as fu, u5_bfs as ub
from fuchs import Mq, mob_ap, hd0, to_klein, to_poinc
gmpy2.get_context().precision = 200
R_COV, R_U_H, S, K = 1.2, 1.45, int(os.environ.get('CHART_S', '96')), 16
R_U = float(np.tanh(R_U_H/2)); NS = 30000; R_NB = 2.6; MAPSIZE = 1024
NW = int(os.environ.get('CHART_WORKERS', '3'))
t0 = time.time()
def log(*a): print(*a, f'| {time.time()-t0:.0f} s', flush=True)

D = np.load('develop.npz'); FG = D['face_g']; verts = D['verts']
U_, Tt, _ = fu.halfplanes(FG); Ur = np.stack([U_.real, U_.imag])
rD = float(hd0(to_poinc(verts)).max())
TL = ub.Tiling('.')
P = np.load('P_folded.npy'); nV = len(P)
VX = np.load('vertex_x.npz'); yv, xv = VX['y'], VX['x']; assert len(yv) == nV
xok = np.isfinite(xv)
# E95 v3 payload: the charts store the two zeros (w, x) of f(z) = (z - x)(z - w)/z^2 instead of (y, x): w = (x + y - 1)/(x - 1),
# y = (1 - x)(1 - w). The viewer needs the map, and (y, x) -> w divides by x - 1 (ill-conditioned in float32 near the punctures
# with x -> 1); products are well-conditioned. Channel pair 1 = w (meta ay, by), pair 2 = x (ax, bx); meta payload = 'wx'.
with np.errstate(all='ignore'):
    wv = (xv + yv - 1)/(xv - 1)
xok &= np.isfinite(wv)
log(f'surface vertices {nV}; x finite {int(xok.sum())}, NaN {int((~xok).sum())} (max |y| among NaN {np.abs(yv[~xok]).max() if (~xok).any() else 0:.4f}); '
    f'Dirichlet circumradius {rD:.4f}; faces {len(FG)}; edge->element match worst {TL.match_worst:.2e}')


def inD(w):
    w = np.asarray(w); out = np.zeros(len(w), bool)
    for s0 in range(0, len(w), 40000):
        k = to_klein(w[s0:s0 + 40000]); out[s0:s0 + 40000] = ((np.stack([k.real, k.imag], 1) @ Ur) - Tt[None, :] <= 1e-12).all(1)
    return out


def tiles_near(c, R):
    """numpy 2x2 elements g with g(D) meeting B(c, R) (c in D); the identity first."""
    cm = mpc(complex(c)); I2 = (mpc(1), mpc(0)); seen = {fu.mp_key(I2)}; front = [I2]; out = []
    while front:
        nxt = []
        for g in front:
            out.append(fu.mp_to_np(g))
            gi = fu.mp_inv(g)
            pp = ub.mp_hyp((gi[0]*cm + gi[1])/(gi[1].conjugate()*cm + gi[0].conjugate()))
            for j in np.where(TL.seg_dist(pp) < R + 1e-3)[0]:
                h = fu.mp_mul(g, TL.F[TL.E[j]]); k = fu.mp_key(h)
                if k not in seen: seen.add(k); nxt.append(h)
        front = nxt
    return np.array(out)


def mob(G, w): return (G[..., 0, 0]*w + G[..., 0, 1])/(G[..., 1, 0]*w + G[..., 1, 1])
def minv(G):
    out = np.empty_like(G); out[..., 0, 0] = G[..., 1, 1]; out[..., 1, 1] = G[..., 0, 0]; out[..., 0, 1] = -G[..., 0, 1]; out[..., 1, 0] = -G[..., 1, 0]
    return out/np.sqrt(G[..., 0, 0]*G[..., 1, 1] - G[..., 0, 1]*G[..., 1, 0])[..., None, None]
def to_c(c, w): return (w - c)/(1 - np.conj(c)*w)            # Mq(c)^-1 w
def adj(g): return np.array([[g[1, 1], -g[0, 1]], [-g[1, 0], g[0, 0]]])


# ---------------------------------------------------------------- 1. chart centres (cached)
if os.path.exists('chart9_centres.pkl'):
    centres, tiles27, cov, samp = pickle.load(open('chart9_centres.pkl', 'rb'))
    log(f'chart centres from cache: {len(centres)}, coverage {cov:.4f}')
else:
    rng = np.random.default_rng(3)
    samp = np.concatenate([P[rng.choice(nV, NS, replace=False)], np.load('u5_centres_ins.npy')])
    centres = [0j]; tiles27 = []; dmin = np.full(len(samp), np.inf)
    while True:
        c = centres[-1]; tl = tiles_near(c, 2.7); tiles27.append(tl)
        for g in tl:
            dmin = np.minimum(dmin, hd0(to_c(c, mob(g, samp))))
        cov = float(dmin.max())
        if len(centres) % 100 == 0: log(f'   {len(centres)} centres, coverage so far {cov:.3f}, tiles for the last {len(tl)}')
        if cov <= R_COV or len(centres) >= 3000: break
        centres.append(samp[int(np.argmax(dmin))])
    centres = np.array(centres)
    pickle.dump((centres, tiles27, cov, samp), open('chart9_centres.pkl', 'wb'))
NC = len(centres)
nt = [len(t) for t in tiles27]
log(f'chart centres: {NC}; coverage radius over {len(samp)} samples (upper bound) {cov:.4f}; tiles meeting B(c, 2.7): min {min(nt)}, median {int(np.median(nt))}, max {max(nt)}')

# ---------------------------------------------------------------- 2. neighbour table
nbl = []; nmax = 0
for k in range(NC):
    ck = centres[k]; ent = []
    for g in tiles27[k]:
        q = to_c(ck, mob(g, centres))                         # lifts g c_j in chart k's coordinate
        dq = hd0(q)
        for j in np.where(dq < R_NB)[0]:
            if dq[j] < 1e-9: continue                         # the centre itself
            Tm = np.linalg.inv(Mq(centres[j])) @ adj(g) @ Mq(ck)   # adj(g) = g^-1 exactly (det 1; an LU inverse loses ~1e-7 x |a| here)
            Tm = fu.su11(Tm[None])[0]
            ent.append((float(dq[j]), int(j), Tm[0, 0], Tm[0, 1]))
    ent.sort(key=lambda e: e[0]); nbl.append(ent); nmax = max(nmax, len(ent))
NB = min(nmax, 40)
nb = np.zeros((NC, NB, 5)); nb[:, :, 0] = -1
for k, ent in enumerate(nbl):
    for i, (dq, j, a, b) in enumerate(ent[:NB]): nb[k, i] = [j, a.real, a.imag, b.real, b.imag]
cnt = [len(e) for e in nbl]
log(f'neighbour table: centre lifts within {R_NB} per chart: min {min(cnt)}, median {int(np.median(cnt))}, max {max(cnt)}; NB {NB} (truncated charts {sum(c > NB for c in cnt)}); '
    f'nearest neighbour distance: min {min(e[0][0] for e in nbl):.3f}, max {max(e[0][0] for e in nbl):.3f}')

# ---------------------------------------------------------------- 3. placements and chart map
size = MAPSIZE
yy, xx = np.mgrid[0:size, 0:size]; wpix = (((xx + 0.5)/size*2 - 1) + 1j*(1 - (yy + 0.5)/size*2)).reshape(-1)
m = np.where(np.abs(wpix) < 1)[0]; m = m[inD(wpix[m])]
cand_k, cand_g, cand_pos = [], [], []
for k in range(NC):
    for g in tiles27[k]:
        pos = mob(minv(g[None])[0], centres[k])              # lift g^-1 c_k, the one near D when g^-1 D... (all kept; cheap)
        if inD(np.array([pos]))[0] or TL.seg_dist(ub.hyp_from_poinc(np.array([pos]))[0]).min() <= 1.3:   # lifts within 1.3 of D
            cand_k.append(k); cand_g.append(g); cand_pos.append(pos)
cand_pos = np.array(cand_pos); CH = ub.hyp_from_poinc(cand_pos); CHJ = CH*np.array([-1.0, 1, 1])
best = np.full(len(m), -1); bd = np.full(len(m), np.inf)
CHK = max(200, int(2e8//max(1, len(cand_pos))))
for s0 in range(0, len(m), CHK):
    XH = ub.hyp_from_poinc(wpix[m[s0:s0 + CHK]])
    ip = -(XH @ CHJ.T)                                        # cosh of the distance
    a = np.argmin(ip, 1); best[s0:s0 + CHK] = a; bd[s0:s0 + CHK] = np.arccosh(np.maximum(1, ip[np.arange(len(a)), a]))
used = np.unique(best); remap = np.full(len(cand_pos), -1); remap[used] = np.arange(len(used))
assert len(used) < 65536
idxmap = np.full(size*size, -1); idxmap[m] = remap[best]
hole = (idxmap < 0).reshape(size, size)
ii = distance_transform_edt(hole, return_distances=False, return_indices=True)
idxfull = idxmap.reshape(size, size)[ii[0], ii[1]]
cmap = np.zeros((size, size, 4), np.uint8); cmap[..., 0] = idxfull & 255; cmap[..., 1] = idxfull >> 8; cmap[..., 3] = 255
placements = []
for u in used:
    k, g = cand_k[int(u)], cand_g[int(u)]
    Mp = np.linalg.inv(Mq(centres[k])) @ g; Mp = Mp/np.sqrt(np.linalg.det(Mp))
    placements.append(dict(chart=int(k), M=[Mp[0, 0].real, Mp[0, 0].imag, Mp[0, 1].real, Mp[0, 1].imag, Mp[1, 0].real, Mp[1, 0].imag, Mp[1, 1].real, Mp[1, 1].imag]))
log(f'chart map {size}: candidates {len(cand_pos)}, placements used {len(placements)}; texels of D {len(m)} (others filled from the nearest: {int(hole.sum())}); '
    f'distance texel centre -> chosen centre lift: median {np.median(bd):.3f}, 99% {np.percentile(bd, 99):.3f}, max {bd.max():.3f} (the viewer walks to the nearest centre)')

# ---------------------------------------------------------------- 4. chart grids (parallel, cached)
XHv = None
def _init():
    global XHv
    XHv = ub.hyp_from_poinc(P)*np.array([-1.0, 1, 1])        # J X_v
def chord(a, b): return np.abs(a - b)/np.sqrt((1 + np.abs(a)**2)*(1 + np.abs(b)**2))
fib = np.arange(4000) + 0.5; th_ = np.arccos(1 - 2*fib/4000); ph_ = np.pi*(1 + 5**0.5)*fib
sph = np.stack([np.sin(th_)*np.cos(ph_), np.sin(th_)*np.sin(ph_), np.cos(th_)], 1)
cand_pts = (sph[:, 0] + 1j*sph[:, 1])/(1 - sph[:, 2] + 1e-15)
def far_point(vals):
    vs = vals[np.random.default_rng(5).choice(len(vals), min(4000, len(vals)), replace=False)]
    dd = np.array([chord(vs, b).min() for b in cand_pts]); b = cand_pts[int(np.argmax(dd))]
    a = -1/np.conj(b) if abs(b) > 1e-12 else 1e12
    return a, b, float(dd.max())
ug = np.linspace(-R_U, R_U, S); UG = (ug[None, :] + 1j*ug[::-1, None]).reshape(-1)
CUT = float(np.cosh(2.25))
def fit(pts, vals):
    tree = cKDTree(np.stack([pts.real, pts.imag], 1))
    dd, nbr = tree.query(np.stack([UG.real, UG.imag], 1), K)
    sig = 1.5*dd[:, 5:6] + 1e-12; wgt = np.exp(-(dd/sig)**2)
    dz = pts[nbr] - UG[:, None]
    A = np.stack([np.ones_like(dz), dz, dz*dz], -1)
    AhW = np.conj(A).transpose(0, 2, 1)*wgt[:, None, :]
    N_ = AhW @ A + 1e-30*np.eye(3)
    rhs = (AhW @ vals[nbr][:, :, None])[:, :, 0]
    return np.linalg.solve(N_, rhs[..., None])[..., 0][:, 0].reshape(S, S)
def chart(k):
    c = centres[k]; tl = tiles27[k]
    pts, ids = [], []
    for g in tl:
        q = mob(minv(g[None])[0], c)                          # vertices v with d(g P_v, c) = d(P_v, g^-1 c) < 2.25
        Q = ub.hyp_from_poinc(np.array([q]))[0]
        idx = np.where(-(XHv @ Q) < CUT)[0]
        if not len(idx): continue
        w = to_c(c, mob(g, P[idx])); keep = np.abs(w) < 0.8
        pts.append(w[keep]); ids.append(idx[keep])
    pts = np.concatenate(pts); ids = np.concatenate(ids); vy = yv[ids]; vx = xv[ids]; vw = wv[ids]; fx_ = xok[ids]
    assert len(pts) >= K, ('chart with too few mesh lifts', k, len(pts))
    gr = np.zeros((S, S, 4), np.float32)
    if fx_.sum() >= K:
        ay, by, sy = far_point(vw[fx_]); ax, bx, sx = far_point(vx[fx_])
        with np.errstate(all='ignore'):
            Wf = (vw[fx_] - ay)/(vw[fx_] - by); Xf = (vx[fx_] - ax)/(vx[fx_] - bx)
        gy = fit(pts[fx_], Wf); gx = fit(pts[fx_], Xf)
    else:                                                     # no vertex with x in this chart: placeholder (x = w = 0)
        ay, by, sy, ax, bx, sx = 0j, 1e6 + 0j, 0.0, 0j, 1e6 + 0j, 0.0; gy = np.zeros((S, S), complex); gx = np.zeros((S, S), complex)
    gr[..., 0] = gy.real; gr[..., 1] = gy.imag; gr[..., 2] = gx.real; gr[..., 3] = gx.imag
    inside = np.where((np.abs(pts) < R_U*0.95) & fx_)[0][::97]
    fx = (pts[inside].real + R_U)/(2*R_U)*(S - 1); fy = (R_U - pts[inside].imag)/(2*R_U)*(S - 1)
    i0, j0 = np.floor(fy).astype(int), np.floor(fx).astype(int); ty, tx = fy - i0, fx - j0
    bl = lambda G_: (G_[i0, j0]*(1 - tx)*(1 - ty) + G_[i0, j0 + 1]*tx*(1 - ty) + G_[i0 + 1, j0]*(1 - tx)*ty + G_[i0 + 1, j0 + 1]*tx*ty)
    vW, vX = bl(gy), bl(gx)
    with np.errstate(all='ignore'):
        yrec = (1 - (ax - bx*vX)/(1 - vX))*(1 - (ay - by*vW)/(1 - vW))
    meta = dict(c=[c.real, c.imag], ay=[complex(ay).real, complex(ay).imag], by=[by.real, by.imag], ax=[complex(ax).real, complex(ax).imag],
                bx=[bx.real, bx.imag], ysep=sy, xsep=sx, npts=int(len(pts)), nx=int(fx_.sum()), tiles=int(len(tl)))
    return k, gr, meta, float(np.median(chord(yrec, vy[inside])))

gpath = 'chart9_grids.npy'; dpath = 'chart9_done.json'
if os.path.exists(gpath) and os.path.exists(dpath):
    done = json.load(open(dpath))
    if done.get('NC') != NC or done.get('S') != S: done = dict(NC=NC, S=S, meta={})
else:
    done = dict(NC=NC, S=S, meta={})
grids = np.lib.format.open_memmap(gpath, mode='r+' if (os.path.exists(gpath) and done['meta']) else 'w+', dtype=np.float32, shape=(NC, S, S, 4))
todo = [k for k in range(NC) if str(k) not in done['meta']]
log(f'chart grids: {NC - len(todo)} from cache, {len(todo)} to compute ({NW} workers, S = {S})')
if todo:
    from multiprocessing import Pool
    with Pool(NW, initializer=_init) as pool:
        for n_, (k, gr, meta, ho_) in enumerate(pool.imap_unordered(chart, todo)):
            grids[k] = gr; meta['ho'] = ho_; done['meta'][str(k)] = meta
            if n_ % 25 == 0 or n_ == len(todo) - 1:
                grids.flush(); json.dump(done, open(dpath + '.tmp', 'w')); os.replace(dpath + '.tmp', dpath)
                log(f'   chart {k} ({n_ + 1}/{len(todo)}): tiles {meta["tiles"]}, mesh lifts {meta["npts"]} (x finite {meta["nx"]}), Mobius separation y {meta["ysep"]:.3f} x {meta["xsep"]:.3f}, hold-out y median {ho_:.2e}')
meta = [done['meta'][str(k)] for k in range(NC)]
ho = [mt['ho'] for mt in meta]
log(f'mesh lifts per chart {min(mt["npts"] for mt in meta)} .. {max(mt["npts"] for mt in meta)}; hold-out bilinear y error, median per chart: median {np.median(ho):.2e}, max {max(ho):.2e}; '
    f'grid finite {bool(np.isfinite(grids).all())}')

# ---------------------------------------------------------------- 5. check C (registered) and x consistency
def cw(t):
    t2 = t*t; t3 = t2*t; return np.array([-.5*t3 + t2 - .5*t, 1.5*t3 - 2.5*t2 + 1, -1.5*t3 + 2*t2 + .5*t, .5*t3 - .5*t2])
def read(k, u, cubic):
    fx = (u.real + R_U)/(2*R_U)*(S - 1); fy = (R_U - u.imag)/(2*R_U)*(S - 1)
    if cubic:
        ix, iy = int(np.floor(fx)), int(np.floor(fy)); wx, wy = cw(fx - ix), cw(fy - iy); g = np.zeros(4)
        for j in range(4):
            for i in range(4): g += wx[i]*wy[j]*grids[k, min(S - 1, max(0, iy + j - 1)), min(S - 1, max(0, ix + i - 1))]
    else:
        i0, j0 = int(fy), int(fx); ty, tx = fy - i0, fx - j0
        g = grids[k, i0, j0]*(1-tx)*(1-ty) + grids[k, i0, j0+1]*tx*(1-ty) + grids[k, i0+1, j0]*(1-tx)*ty + grids[k, i0+1, j0+1]*tx*ty
    mt = meta[k]; Wc = complex(g[0], g[1]); X = complex(g[2], g[3])
    w_ = (complex(*mt['ay']) - complex(*mt['by'])*Wc)/(1 - Wc); x_ = (complex(*mt['ax']) - complex(*mt['bx'])*X)/(1 - X)
    return (1 - x_)*(1 - w_), x_          # y = (1 - x)(1 - w)
ey, eyc, ex, exm, near = [], [], [], [], []
for v in np.random.default_rng(11).choice(nV, 1600, replace=False):
    best_u = None
    for g in tiles_near(P[v], R_U_H):
        u = to_c(centres, mob(minv(g[None])[0], P[v]))        # chart coordinates relative to the lifts g c_k
        if not xok[v]: continue                               # vertices over branch values have no stored (x, w)
        for k in np.where(np.abs(u) < R_U*0.97)[0]:
            yb, xb = read(int(k), u[k], False); yc, xc = read(int(k), u[k], True)
            ey.append(float(chord(yb, yv[v]))); eyc.append(float(chord(yc, yv[v])))
            if xok[v]: (ex if abs(yv[v]) >= 0.1 else exm).append(float(chord(xc, xv[v])))
            if best_u is None or abs(u[k]) < best_u[0]: best_u = (abs(u[k]), float(chord(yc, yv[v])))
    if best_u is not None: near.append(best_u[1])
ey, eyc, ex, near = np.array(ey), np.array(eyc), np.array(ex), np.array(near)
log(f'CONSISTENCY C (registered, bilinear as E92/E94): {len(ey)} (mesh vertex, covering placement) evaluations; chordal y error median {np.median(ey):.2e}, '
    f'99th pct {np.percentile(ey, 99):.2e}, max {ey.max():.2e}, fraction > 1e-3 {(ey > 1e-3).mean():.4f}')
log(f'viewer read (bicubic): y median {np.median(eyc):.2e}, 99th pct {np.percentile(eyc, 99):.2e}, max {eyc.max():.2e}, fraction > 1e-3 {(eyc > 1e-3).mean():.4f}; '
    f'nearest-centre chart only ({len(near)} vertices; vertices with no covering chart {1500 - len(near)}): median {np.median(near):.2e}, 99th pct {np.percentile(near, 99):.2e}, max {near.max():.2e}')
log(f'x (bicubic, vertices with |y| >= 0.1, {len(ex)} evaluations): chordal error median {np.median(ex):.2e}, 99th pct {np.percentile(ex, 99):.2e}, max {ex.max():.2e}, fraction > 1e-3 {(ex > 1e-3).mean():.4f}')
C_ok = bool(np.median(ey) < 1e-4 and (ey > 1e-3).mean() < 1e-3)

# ---------------------------------------------------------------- 6. files
raw = np.ascontiguousarray(grids).astype('<f4').tobytes(); per = S*S*4*4; cpp = max(1, 60_000_000//per)
parts = []
for i, s0 in enumerate(range(0, NC, cpp)):
    fn = f'chart_grids_{i}.bin'; open(fn, 'wb').write(raw[s0*per:(s0 + cpp)*per]); parts.append(fn)
Image.fromarray(cmap, 'RGBA').save('chart_map.png', optimize=True)
coef = json.load(open(os.path.join(here, '..', 'explorer', 'data.json'))).get('coef')
json.dump(dict(S=S, NC=NC, RU=R_U, mapSize=size, placements=placements, charts=meta, parts=parts, NB=NB, nb=[round(float(v), 12) for v in nb.reshape(-1)],
               coef=coef, payload='wx', mask_y=(float(os.environ['MASK_Y']) if os.environ.get('MASK_Y') else (0.1 if int((~xok & (np.abs(yv) < 0.1)).sum()) > 2000 else None)), coverage=cov, check_C=dict(n=len(ey), median=float(np.median(ey)), frac_gt_1e3=float((ey > 1e-3).mean()), passed=C_ok)),
          open('chart_meta.json', 'w'))
log(f'wrote {len(parts)} chart grid parts ({len(raw)/1e6:.1f} MB), chart_map.png ({os.path.getsize("chart_map.png")/1e6:.2f} MB), chart_meta.json ({os.path.getsize("chart_meta.json")/1e6:.2f} MB); check C passed: {C_ok}')
sys.exit(0 if C_ok else 1)
