#!/usr/bin/env python3
"""E92-U viewer data (per6_conformal.html): the conformal model of compact Per_6(0) on its Dirichlet polygon.
Inputs: develop.npz (developed disk z, Dirichlet faces), surface_data.js (colours), tiling_data.js (per-vertex y, x,
regions), explorer/data.json (markers). Tiles g(Delta) meeting the Dirichlet polygon D are found by BFS over the side
pairings; their triangles are rendered into a texture in Poincare coordinates (the painting is moved by isometries only,
so it stays conformal). Also: id raster, face data (Klein half-planes, fold maps g^-1, pair ids), vertex cycles
(check: V - E + F = -10 and angle sum 2 pi per cycle), markers placed by exact point location."""
import base64
import io
import json
import os
import sys
import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import distance_transform_edt
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
from develop import Mq, mob_ap, mob_norm, hdist, to_klein, to_poinc, gkey, hd0
from meshloc import Locator, place
D = np.load('develop.npz'); z = D['z']; Tn = D['Tn']; orig = D['orig']
G = list(D['gens']); G2 = G + [np.linalg.inv(g) for g in G]
vk = D['verts']; face_g = list(D['face_g']); pair_of = list(D['pair_of']); nf = len(vk)
I2 = np.eye(2, dtype=complex)
U_ = np.array([np.exp(1j*np.angle(mob_ap(g, 0))) for g in face_g]); Tt = np.array([np.tanh(hd0(mob_ap(g, 0))/2) for g in face_g])
Ginv = [np.linalg.inv(g) for g in face_g]
def viol(w):
    k = to_klein(np.atleast_1d(w)); return (k[:, None]*np.conj(U_[None, :])).real - Tt[None, :]
def fold(w):
    w = complex(w)
    for _ in range(1000):
        v = viol(w)[0]; a = int(np.argmax(v))
        if v[a] <= 1e-13: return w
        w = mob_ap(Ginv[a], w)
    raise RuntimeError('fold did not converge')
zp = to_poinc(vk)
print('Dirichlet polygon: faces', nf, '| vertex radii (Poincare) max', f'{np.abs(zp).max():.6f}', '| circumradius', f'{hd0(zp).max():.4f}')
# vertex cycles
par = list(range(nf))
def fnd(a):
    while par[a] != a: par[a] = par[par[a]]; a = par[a]
    return a
for a in range(nf):
    b = pair_of[a]
    for p_, q_ in ((a, (b + 1) % nf), ((a + 1) % nf, b)):
        rp, rq = fnd(p_), fnd(q_)
        if rp != rq: par[rp] = rq
cls = {}
for a in range(nf): cls.setdefault(fnd(a), []).append(a)
ang = []
for a in range(nf):
    Mi = np.linalg.inv(Mq(zp[a])); u, w = mob_ap(Mi, zp[a - 1]), mob_ap(Mi, zp[(a + 1) % nf]); ang.append(abs(np.angle(w/u)))
cyc_sums = [sum(ang[a] for a in m) for m in cls.values()]
print('U4 vertex cycles:', len(cls), '| V - E + F =', len(cls) - nf//2 + 1, '(genus 6: -10) | angle sums / 2pi:',
      [round(s/(2*np.pi), 9) for s in cyc_sums])
# fold every developed vertex into D (vectorized), remembering the element F_v with F_v(z_v) in D
def vfold(w, maxit=200):
    w = w.copy(); Fm = np.tile(I2, (len(w), 1, 1))
    for it in range(maxit):
        k = to_klein(w); v = (k[:, None]*np.conj(U_[None, :])).real - Tt[None, :]
        a = np.argmax(v, 1); out = v[np.arange(len(w)), a] > 1e-13
        if not out.any(): return w, Fm, it
        for f in np.unique(a[out]):
            m = out & (a == f); Gf = Ginv[f]
            w[m] = (Gf[0, 0]*w[m] + Gf[0, 1])/(Gf[1, 0]*w[m] + Gf[1, 1]); Fm[m] = np.einsum('ij,njk->nik', Gf, Fm[m])
    raise RuntimeError('vfold did not converge')
zf, Fm, nit = vfold(z)
okey = np.round(np.stack([(Fm[:, 0, 1]/Fm[:, 1, 1]).real, (Fm[:, 0, 1]/Fm[:, 1, 1]).imag], 1)*1e6).astype(np.int64)
ukeys, vid = np.unique(okey, axis=0, return_inverse=True); vid = vid.reshape(-1)
rep = {int(i): Fm[int(np.argmax(vid == i))] for i in range(len(ukeys))}
tri_ids = vid[Tn]
print(f'folded {len(z)} developed vertices into D in {nit} steps; distinct group elements used {len(ukeys)}; '
      f'triangles straddling an edge of D {int((tri_ids.min(1) != tri_ids.max(1)).sum())}', flush=True)
def raster(size, fill):
    img = Image.new('RGBA', (size, size), (0, 0, 0, 0)); dr = ImageDraw.Draw(img); nd = 0
    for i, M in rep.items():
        fs = np.where((tri_ids == i).any(1))[0]
        if not len(fs): continue
        w = (M[0, 0]*z[Tn[fs]] + M[0, 1])/(M[1, 0]*z[Tn[fs]] + M[1, 1])
        Q = np.stack([(0.5 + 0.5*w.real)*size, (0.5 - 0.5*w.imag)*size], -1); nd += len(fs)
        for t, f in enumerate(fs):
            dr.polygon([tuple(Q[t, 0]), tuple(Q[t, 1]), tuple(Q[t, 2])], fill=fill(f) + (255,))
    A = np.array(img); hole = A[..., 3] == 0
    yy, xx = np.mgrid[0:size, 0:size]; wpix = ((xx + 0.5)/size*2 - 1) + 1j*(1 - (yy + 0.5)/size*2)
    inD = np.zeros(size*size, bool); wf_ = wpix.reshape(-1); m = np.where(np.abs(wf_) < 1)[0]
    for s0 in range(0, len(m), 200000):
        mm = m[s0:s0 + 200000]; inD[mm] = (viol(wf_[mm]) <= 0).all(1)
    inD = inD.reshape(size, size)
    print(f'raster {size}: triangles drawn {nd}; pixels inside D {int(inD.sum())}, of them uncovered {int((inD & hole).sum())}', flush=True)
    idx = distance_transform_edt(hole, return_distances=False, return_indices=True)
    return Image.fromarray(A[idx[0], idx[1]][..., :3], 'RGB')
src = open('surface_data.js').read(); Dd = json.loads(src[src.index('{'):src.rindex('}') + 1])
col = np.frombuffer(base64.b64decode(Dd['col']), np.uint8).reshape(-1, 3)
tex = raster(4096, lambda f: tuple(int(x) for x in col[orig[Tn[f, 0]]]))
bt = io.BytesIO(); tex.save(bt, 'WEBP', quality=88, method=4); tex.save('conformal_texture.png')
ids = raster(1024, lambda f: ((int(orig[Tn[f, 0]]) >> 16) & 255, (int(orig[Tn[f, 0]]) >> 8) & 255, int(orig[Tn[f, 0]]) & 255))
bi = io.BytesIO(); ids.save(bi, 'PNG', optimize=True)
# markers by exact point location
L = Locator('.')
D2 = json.load(open('../explorer/data.json')); mk = []; worst = 0.0
for kind, items in (('centre', D2['centers']), ('puncture', D2['punctures'])):
    for it in items:
        loc, how = L.locate_sheet(complex(*it['y']), complex(*it['x']), it['sheet'])
        dx, f, wf = loc; worst = max(worst, dx)
        if how != 'continuation': print(f"  marker {kind} {it.get('label', it.get('name'))}: {how}")
        p = fold(place(z, Tn, f, wf))
        mk.append(dict(kind=kind, label=str(it.get('label', it.get('name'))), p=[p.real, p.imag], y=it['y'], x=it['x']))
print('markers', len(mk), '| max x-mismatch after continuation (chordal)', f'{worst:.2e}')
# faces for the viewer
pid = {}
for a in range(nf):
    k = min(a, pair_of[a]); pid.setdefault(k, len(pid))
faces = []
for a in range(nf):
    gi = Ginv[a]; P, Q = zp[a], zp[(a + 1) % nf]
    Mi = np.linalg.inv(Mq(P)); qq = mob_ap(Mi, Q); mid = mob_ap(Mq(P), np.tanh(np.arctanh(abs(qq))/2)*np.exp(1j*np.angle(qq)))
    faces.append(dict(u=[U_[a].real, U_[a].imag], t=float(Tt[a]), gi=[gi[0, 0].real, gi[0, 0].imag, gi[0, 1].real, gi[0, 1].imag,
                      gi[1, 0].real, gi[1, 0].imag, gi[1, 1].real, gi[1, 1].imag], pair=pid[min(a, pair_of[a])], partner=int(pair_of[a]),
                      P=[P.real, P.imag], Q=[Q.real, Q.imag], mid=[mid.real, mid.imag]))
td = open('tiling_data.js').read(); TD = json.loads(td[td.index('{'):td.rindex('}') + 1])
vcls = [0]*nf
for n_, m in enumerate(cls.values()):
    for a in m: vcls[a] = n_
out = dict(faces=faces, nPairs=len(pid), vcls=vcls, nCycles=len(cls),
           tex='data:image/webp;base64,' + base64.b64encode(bt.getvalue()).decode(), texSize=4096,
           ids='data:image/png;base64,' + base64.b64encode(bi.getvalue()).decode(), idSize=1024,
           Y=TD['Y'], X=TD['X'], piece=TD['piece'], sheet=TD['sheet'], regions=TD['regions'], markers=mk)
open('conformal_data.js', 'w').write('window.CF = ' + json.dumps(out) + ';\n')
print(f'texture {len(bt.getvalue())/1e6:.1f} MB, ids {len(bi.getvalue())/1e6:.1f} MB; wrote conformal_data.js')
