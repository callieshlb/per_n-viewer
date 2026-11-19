#!/usr/bin/env python3
"""E95 viewer data (E94 conformal_export.py without the mesh texture / id raster: 4 M vertices): faces of the certified Dirichlet
polygon D (fine mesh develop.npz), exactly as E94 builds them, vertex cycles, and the 252 centres as markers (p = u5_centres_ins.npy,
the inserted placement adopted for U5, already folded into D). Writes conformal_data.js (window.CF)."""
import base64, io, json, os, sys
import numpy as np
from PIL import Image
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import fuchs as fu
from fuchs import Mq, mob_ap, to_klein, to_poinc, hd0
D = np.load(os.path.join(here, 'develop.npz'))
vk = D['verts']; FG = D['face_g']; pair_of = [int(v) for v in D['pair_of']]; nf = len(vk)
res = dict(verts=vk, face_g=FG, bounded=True)
U_, Tt, _ = fu.halfplanes(FG)
Ginv = [np.linalg.inv(g) for g in FG]
zp = to_poinc(vk)
print('Dirichlet polygon: faces', nf, '| vertex radii (Poincare) max', f'{np.abs(zp).max():.9f}', '| circumradius', f'{hd0(zp).max():.4f}', flush=True)
cls, sums = fu.vertex_cycles(res, np.array(pair_of))
GEN_ = json.load(open(os.path.join(here, '..', 'branch.json')))['genus']
print('U4 vertex cycles:', len(cls), '| V - E + F =', len(cls) - nf//2 + 1, f'(genus {GEN_:g}: {2 - 2*GEN_:g}) | angle sums / 2pi in',
      f'[{min(sums)/(2*np.pi):.9f}, {max(sums)/(2*np.pi):.9f}] (double-precision vertices)', flush=True)
D2 = json.load(open(os.path.join(here, '..', 'explorer', 'data.json')))
P = np.load(os.path.join(here, 'u5_centres_ins.npy')); assert len(P) == len(D2['centers']) and np.isfinite(P).all()
viol = ((to_klein(P)[:, None]*np.conj(U_[None, :])).real - Tt[None, :]).max(1)
assert (viol <= 1e-9).all(), ('centre outside D', float(viol.max()))
mk = [dict(kind='centre', label=str(c['label']), p=[float(p.real), float(p.imag)], y=c['y'], x=c['x']) for c, p in zip(D2['centers'], P)]
print('markers', len(mk), '| max face violation of a centre (Klein)', f'{viol.max():.2e}', '| max hyperbolic radius', f'{hd0(P).max():.4f}', flush=True)
pid = {}
for a in range(nf):
    k = min(a, pair_of[a]); pid.setdefault(k, len(pid))
faces = []
for a in range(nf):
    gi = Ginv[a]; Pa, Q = zp[a], zp[(a + 1) % nf]
    Mi = np.linalg.inv(Mq(Pa)); qq = mob_ap(Mi, Q); mid = mob_ap(Mq(Pa), np.tanh(np.arctanh(abs(qq))/2)*np.exp(1j*np.angle(qq)))
    faces.append(dict(u=[U_[a].real, U_[a].imag], t=float(Tt[a]), gi=[gi[0, 0].real, gi[0, 0].imag, gi[0, 1].real, gi[0, 1].imag,
                      gi[1, 0].real, gi[1, 0].imag, gi[1, 1].real, gi[1, 1].imag], pair=pid[min(a, pair_of[a])], partner=int(pair_of[a]),
                      P=[Pa.real, Pa.imag], Q=[Q.real, Q.imag], mid=[mid.real, mid.imag]))
assert all(np.isfinite(v) for f in faces for key in ('u', 'gi', 'P', 'Q', 'mid') for v in f[key]), 'non-finite face data'
vcls = [0]*nf
for n_, m in enumerate(cls.values()):
    for a in m: vcls[a] = n_
def png1(rgb):
    b = io.BytesIO(); Image.new('RGB', (1, 1), rgb).save(b, 'PNG'); return 'data:image/png;base64,' + base64.b64encode(b.getvalue()).decode()
b64 = lambda a: base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()
out = dict(faces=faces, nPairs=len(pid), vcls=vcls, nCycles=len(cls), tex=png1((244, 243, 239)), texSize=1, ids=png1((0, 0, 0)), idSize=1,
           Y=b64(np.zeros(2, np.float32)), X=b64(np.zeros(2, np.float32)), piece=b64(np.zeros(1, np.int8)), sheet=None,
           regions={'0': 'Per_9(0)'}, markers=mk)
s = 'window.CF = ' + json.dumps(out) + ';\n'
open(os.path.join(here, 'conformal_data.js'), 'w').write(s)
print(f'wrote conformal_data.js ({len(s)/1e6:.2f} MB): faces {nf}, pairs {len(pid)}, cycles {len(cls)}, markers {len(mk)}; no mesh texture / id raster')
