#!/usr/bin/env python3
"""E94 copy (n = 8: region name, genus from branch.json). E93-U viewer data (E92 conformal_export.py for genus 24): the conformal model of compact Per_7(0) on its Dirichlet polygon D.
Inputs: develop.npz (developed disk z, D's faces), surface_data.json (colour_export.py), explorer/data.json (markers).
Every developed vertex is folded into D (vectorized), triangles are drawn per fold element (keys by (d, theta sinh d), as
fuchs.dedupe: Poincare-coordinate rounding would merge far elements at genus 24) into a texture in Poincare coordinates.
Also: id raster, face data (Klein half-planes, fold maps g^-1, pair ids), vertex cycles (check V - E + F = -46, angle sums
2 pi), markers placed by exact point location (meshloc). Writes conformal_data.js."""
import base64
import io
import json
import os
import sys
import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import distance_transform_edt
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import fuchs as fu
from fuchs import Mq, mob_ap, to_klein, to_poinc, hd0
from meshloc import Locator, place
D = np.load('develop.npz'); z = D['z']; Tn = D['Tn']; orig = D['orig']
vk = D['verts']; face_g = list(D['face_g']); pair_of = list(D['pair_of']); nf = len(vk)
res = dict(verts=vk, face_g=D['face_g'], bounded=True)
U_, Tt, _ = fu.halfplanes(D['face_g'])
Ginv = [np.linalg.inv(g) for g in face_g]
fold = fu.folder(res)
def viol(w):
    k = to_klein(np.atleast_1d(w)); return (k[:, None]*np.conj(U_[None, :])).real - Tt[None, :]
zp = to_poinc(vk)
print('Dirichlet polygon: faces', nf, '| vertex radii (Poincare) max', f'{np.abs(zp).max():.6f}', '| circumradius', f'{hd0(zp).max():.4f}', flush=True)
cls, sums = fu.vertex_cycles(res, np.array(pair_of))
GEN_ = json.load(open(os.path.join(here, '..', 'branch.json')))['genus']
print('U4 vertex cycles:', len(cls), '| V - E + F =', len(cls) - nf//2 + 1, f'(genus {GEN_:g}: {2 - 2*GEN_:g}) | angle sums / 2pi in',
      f'[{min(sums)/(2*np.pi):.9f}, {max(sums)/(2*np.pi):.9f}]', flush=True)
zf, Fm, nit = fold(z)
a_, b_ = Fm[:, 0, 0], Fm[:, 0, 1]
dF = 2*np.arccosh(np.maximum(np.abs(a_), 1.0)); thF = np.angle(b_/np.conj(a_))
okey = np.stack([np.round(dF*1e5), np.round(thF*np.sinh(dF)*1e5)], 1).astype(np.int64)
ukeys, vid = np.unique(okey, axis=0, return_inverse=True); vid = vid.reshape(-1)
first = np.full(len(ukeys), -1); first[vid[::-1]] = np.arange(len(vid))[::-1]
rep = {int(i): Fm[int(first[i])] for i in range(len(ukeys))}
tri_ids = vid[Tn]
print(f'folded {len(z)} developed vertices into D in {nit} steps; distinct fold elements {len(ukeys)} (max displacement '
      f'{dF.max():.3f}); triangles straddling an edge of D {int((tri_ids.min(1) != tri_ids.max(1)).sum())}', flush=True)
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
    for s0 in range(0, len(m), 40000):                                  # E94: 40k chunks (822 faces)
        mm = m[s0:s0 + 40000]; inD[mm] = (viol(wf_[mm]) <= 0).all(1)
    inD = inD.reshape(size, size)
    print(f'raster {size}: triangles drawn {nd}; pixels inside D {int(inD.sum())}, of them uncovered {int((inD & hole).sum())}', flush=True)
    idx = distance_transform_edt(hole, return_distances=False, return_indices=True)
    return Image.fromarray(A[idx[0], idx[1]][..., :3], 'RGB')
SD = json.load(open('surface_data.json'))
col = np.frombuffer(base64.b64decode(SD['col']), np.uint8).reshape(-1, 3)
yv = np.frombuffer(base64.b64decode(SD['y']), np.float64).reshape(-1, 2)
xv = np.frombuffer(base64.b64decode(SD['x']), np.float64).reshape(-1, 2)
tex = raster(4096, lambda f: tuple(int(x) for x in col[orig[Tn[f, 0]]]))
bt = io.BytesIO(); tex.save(bt, 'WEBP', quality=88, method=4); tex.save('conformal_texture.png')
ids = raster(1024, lambda f: ((int(orig[Tn[f, 0]]) >> 16) & 255, (int(orig[Tn[f, 0]]) >> 8) & 255, int(orig[Tn[f, 0]]) & 255))
bi = io.BytesIO(); ids.save(bi, 'PNG', optimize=True)
L = Locator('.')
D2 = json.load(open(os.path.join(here, '..', 'explorer', 'data.json'))); mk = []; worst = 0.0; ties = []
for kind, items in (('centre', D2['centers']), ('puncture', D2['punctures'])):
    for it in items:
        loc = L.locate(complex(*it['y']), complex(*it['x']))
        dx, f, wf = loc; worst = max(worst, dx)
        if L.last_ratio < 10: ties.append((kind, it.get('label', it.get('name')), L.last_ratio))
        p = fold(np.array([place(z, Tn, f, wf)]))[0][0]
        mk.append(dict(kind=kind, label=str(it.get('label', it.get('name'))), p=[p.real, p.imag], y=it['y'], x=it['x']))
print('markers', len(mk), '| max x-mismatch after continuation (chordal)', f'{worst:.2e}', '| lifts within 10x of the best:', ties, flush=True)
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
vcls = [0]*nf
for n_, m in enumerate(cls.values()):
    for a in m: vcls[a] = n_
b64 = lambda a: base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()
out = dict(faces=faces, nPairs=len(pid), vcls=vcls, nCycles=len(cls),
           tex='data:image/webp;base64,' + base64.b64encode(bt.getvalue()).decode(), texSize=4096,
           ids='data:image/png;base64,' + base64.b64encode(bi.getvalue()).decode(), idSize=1024,
           Y=b64(yv.astype(np.float32)), X=b64(xv.astype(np.float32)), piece=b64(np.zeros(len(yv), np.int8)), sheet=None,
           regions={'0': 'Per_8(0)'}, markers=mk)
open('conformal_data.js', 'w').write('window.CF = ' + json.dumps(out) + ';\n')
print(f'texture {len(bt.getvalue())/1e6:.1f} MB, ids {len(bi.getvalue())/1e6:.1f} MB; wrote conformal_data.js')
