#!/usr/bin/env python3
"""E93/E94 chart map fix (2026-11-11; external report chart_map_defect.md, reproduced). chart_export.py assigned a placement
only to pixels whose CENTRE lies in the Dirichlet tile D; the other pixels are RGBA (0, 0, 0, 0), which the shader decodes as
placement 0. Pixels whose centre is just outside D but whose square overlaps D therefore send a strip of D along its boundary
to placement 0 (|u|/RU up to 1.6, outside the chart grid). Fill: every alpha-0 pixel within 2 pixels of an assigned pixel
(and touching the unit disk) gets the placement that minimises max |u|/RU over a 5 x 5 sample of its square (points with
|z| < 1), among the placements of the assigned pixels within 2 pixels, else among all placements; if no placement reaches
max < 0.98 on the whole square, the samples inside D only are used; if that fails too, the pixel is left and counted.
Assigned pixels are not touched. usage: chart_map_fill.py [IN [OUT]] (default chart_map.png -> chart_map.png)"""
import io, json, os, sys
import numpy as np
from PIL import Image
from scipy.ndimage import binary_dilation
here = os.path.dirname(os.path.abspath(__file__))
fin = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, 'chart_map.png'); fout = sys.argv[2] if len(sys.argv) > 2 else fin
meta = json.load(open(os.path.join(here, 'chart_meta.json'))); RU, MS = meta['RU'], meta['mapSize']
PM = np.array([[complex(p['M'][0], p['M'][1]), complex(p['M'][2], p['M'][3]), complex(p['M'][4], p['M'][5]), complex(p['M'][6], p['M'][7])] for p in meta['placements']])
txt = open(os.path.join(here, 'conformal_data.js')).read(); txt = txt[txt.index('{'):txt.rindex('}') + 1]
CF = json.loads(txt); FU = np.array([f['u'] for f in CF['faces']], float); FT = np.array([f['t'] for f in CF['faces']], float)
img = np.array(Image.open(fin).convert('RGBA')); assert img.shape == (MS, MS, 4)
asg = img[..., 3] == 255; pidx = img[..., 0].astype(int) + 256*img[..., 1].astype(int)
h = 2/MS; s = np.linspace(0, 1, 5)
cand = binary_dilation(asg, np.ones((5, 5), bool)) & ~asg
rr, cc = np.nonzero(cand)
x0 = cc*h - 1; y1 = 1 - rr*h                                   # pixel square [x0, x0 + h] x [y1 - h, y1]
near = np.hypot(np.clip(0, x0, x0 + h), np.clip(0, y1 - h, y1)) < 1
rr, cc = rr[near], cc[near]
def mu_max(P, z):                                               # P: placement indices (m,), z: samples (k,) -> (m,)
    M = PM[P]; u = (M[:, 0, None]*z[None, :] + M[:, 1, None])/(M[:, 2, None]*z[None, :] + M[:, 3, None])
    return np.abs(u).max(1)/RU
def inD(z):
    k = np.stack([2*z.real, 2*z.imag], 1)/(1 + np.abs(z)**2)[:, None]
    return (k @ FU.T - FT[None, :]).max(1) <= 0
out = img.copy(); st = dict(candidates=len(rr), filled_local=0, filled_all=0, filled_Donly=0, left_overlaps_D=0, left_no_D_sample=0, worst_mu=0.0)
for r, c in zip(rr, cc):
    X, Y = np.meshgrid(c*h - 1 + s*h, 1 - r*h - s*h); z = (X + 1j*Y).ravel(); z = z[np.abs(z) < 1]
    if len(z) == 0: continue
    loc = np.unique(pidx[max(0, r - 2):r + 3, max(0, c - 2):c + 3][asg[max(0, r - 2):r + 3, max(0, c - 2):c + 3]])
    best, how = None, None
    for P, tag, zz in ((loc, 'filled_local', z), (np.arange(len(PM)), 'filled_all', z)):
        if len(P) == 0: continue
        m = mu_max(P, zz); j = int(np.argmin(m))
        if m[j] < 0.98: best, how, bm = int(P[j]), tag, float(m[j]); break
    if best is None:
        zD = z[inD(z)]
        if len(zD) == 0: st['left_no_D_sample'] += 1; continue
        m = mu_max(np.arange(len(PM)), zD); j = int(np.argmin(m))
        if m[j] < 0.98: best, how, bm = j, 'filled_Donly', float(m[j])
        else: st['left_overlaps_D'] += 1; continue
    out[r, c] = (best & 255, best >> 8, 0, 255); st[how] += 1; st['worst_mu'] = max(st['worst_mu'], bm)
assert np.array_equal(out[asg], img[asg])
print('chart_map_fill', os.path.basename(os.path.dirname(here)), st, flush=True)
bi = io.BytesIO(); Image.fromarray(out, 'RGBA').save(bi, 'PNG', optimize=True); open(fout, 'wb').write(bi.getvalue())
print('wrote', fout, f'{len(bi.getvalue())/1e6:.2f} MB', flush=True)
