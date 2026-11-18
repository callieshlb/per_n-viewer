#!/usr/bin/env python3
"""E99: per-surface-vertex fields for the viewer's click info (replaces the E92 tiling_data.js fields Y, X, piece, sheet,
regions, which are indexed by E92 vertices). Y, X, sheet from this mesh's surface_data.js; piece (E92 pants-region label,
click text only) transferred from the nearest E92 surface vertex in (y, x) (stereographic embedding of both). Writes
vertex_fields.js. Prints the transfer distances (raw)."""
import base64, json, numpy as np
from scipy.spatial import cKDTree
E92 = '../../2026-11-07-E92-per6-atlas/surface3d'
def load(d):
    s = open(f'{d}/surface_data.js').read(); Dd = json.loads(s[s.index('{'):s.rindex('}') + 1])
    y = np.frombuffer(base64.b64decode(Dd['y']), np.float64).reshape(-1, 2); x = np.frombuffer(base64.b64decode(Dd['x']), np.float64).reshape(-1, 2)
    return y[:, 0] + 1j*y[:, 1], x[:, 0] + 1j*x[:, 1], Dd
def emb(z):
    z = np.where(np.isfinite(z), z, 1e300); r2 = np.abs(z)**2
    with np.errstate(over='ignore', invalid='ignore'):
        u = np.stack([2*z.real/(1 + r2), 2*z.imag/(1 + r2), (r2 - 1)/(r2 + 1)], 1)
    u[~np.isfinite(u).all(1)] = [0, 0, 1]; return u
yn, xn, Dn = load('.'); yo, xo, _ = load(E92)
plab = np.load(f'{E92}/pieces.npy').astype(np.int8); assert len(plab) == len(yo)
t = open(f'{E92}/tiling_data.js').read(); TD = json.loads(t[t.index('{'):t.rindex('}') + 1])
kd = cKDTree(np.concatenate([emb(yo), emb(xo)], 1))
d, j = kd.query(np.concatenate([emb(yn), emb(xn)], 1))
piece = plab[j]
print(f'vertices {len(yn)} (E92 {len(yo)}); transfer distance (6-D chordal) median {np.median(d):.2e}, 99% {np.quantile(d, .99):.2e}, max {d.max():.2e}')
print('piece label counts:', dict(zip(*[a.tolist() for a in np.unique(piece, return_counts=True)])))
b64 = lambda a: base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()
sheet = base64.b64decode(Dn['sheet']); assert len(sheet) == len(yn)
out = dict(Y=b64(np.stack([yn.real, yn.imag], 1).astype(np.float32)), X=b64(np.stack([xn.real, xn.imag], 1).astype(np.float32)),
           piece=b64(piece), sheet=Dn['sheet'], regions=TD['regions'])
open('vertex_fields.js', 'w').write('window.VF = ' + json.dumps(out) + ';\n')
print('wrote vertex_fields.js')
