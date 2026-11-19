#!/usr/bin/env python3
"""E95: x on the sheets over the masked disc |y| < 0.1 (198 normal sphere vertices of the fine mesh), to un-mask the viewer.
Sheets collapse there (chordal min separation in x: 3.9e-9 at |y| = 0.1, 1.6e-13 at 0.03, 1.3e-21 at 0.003), so labels are carried
inward from the labelled outer vertices along mesh edges with the glued label maps, by vertex_x.transport_task (512-bit tracking of
the ambiguous roots); when tracking fails or exceeds INNER_CAP = 15 s the nearest root is taken instead (COUNTED: such a root is the right one
only up to the width of its collapsing cluster). Roots then Newton-polished at 512 bits. Writes vx_labels_full.npz {X, seen,
approx (vertices that needed a nearest-root fallback)}. Raw counts printed. NOT a certified labelling inside the disc."""
import json, os, signal, sys, time
import numpy as np
here = os.path.dirname(os.path.abspath(__file__)); os.chdir(here); sys.path.insert(0, here)
sys.argv = [sys.argv[0], '.']
import vertex_x as vx
from gmpy2 import mpc
Y, bp, EP, n, C = vx.Y, vx.bp, vx.EP, vx.n, vx.C
normal = np.where(bp < 0)[0]; nb_ = (len(normal) + vx.BLK - 1)//vx.BLK
R = np.zeros((len(Y), n), complex)
for b in range(nb_): R[normal[b*vx.BLK:(b + 1)*vx.BLK]] = np.load(f'vx_roots/b{b:04d}.npy')
z = np.load('vx_labels_polished.npz'); X = z['X'].copy(); seen = set(int(v) for v in z['seen'])
inner = [int(v) for v in normal if int(v) not in seen]
adj = {}
for e in EP['edges']:
    if e.get('toGroup'): continue
    a, c, m = e['a'], e['c'], e['m']; inv = [0]*n
    for i, j in enumerate(m): inv[j] = i
    adj.setdefault(a, []).append((c, m)); adj.setdefault(c, []).append((a, inv))
CAP = int(os.environ.get('INNER_CAP', '15'))
class _TO(Exception): pass
def _al(*a): raise _TO()
signal.signal(signal.SIGALRM, _al)
def polish(v, xs):
    cc = C.cy(mpc(complex(Y[v])))[0]; return np.array([complex(C.newton(cc, mpc(complex(x)))[0]) for x in xs])
t0 = time.time(); approx = []; tracked = 0; front = sorted(seen & set(u for v in inner for u, _ in adj.get(v, [])))
todo = set(inner); order = sorted(inner, key=lambda v: -abs(Y[v]))          # from the rim inward
while todo:
    prog = False
    for w in order:
        if w not in todo: continue
        src = [(u, m) for u, m in adj.get(w, []) if u in seen]
        if not src: continue
        u, m_uw = max(src, key=lambda t: abs(Y[t[0]]))                       # the labelled neighbour farthest from y = 0
        minv = [0]*n
        for i, j in enumerate(m_uw): minv[j] = i                             # adj[w] holds maps w -> u; we need u -> w
        Rw = polish(w, R[w])
        signal.alarm(CAP)
        try:
            out, tr, det = vx.transport_task((complex(Y[u]), complex(Y[w]), X[u], Rw, minv)); tracked += tr
        except (_TO, RuntimeError, AssertionError):
            out = np.empty(n, complex); nn = np.argmin(np.abs(Rw[None, :] - X[u][:, None]), 1)
            for l in range(n): out[minv[l]] = Rw[nn[l]]
            approx.append(w)
        signal.alarm(0)
        X[w] = out; seen.add(w); todo.discard(w); prog = True
        if (len(inner) - len(todo)) % 20 == 0: print(f'   {len(inner) - len(todo)} of {len(inner)} done, fallbacks {len(approx)}, |y| now {abs(Y[w]):.4f}, {time.time()-t0:.0f} s', flush=True)
    if not prog: break
dup = sum(1 for w in inner if w in seen and len(set(np.round(X[w], 14))) < n)
print(f'inner normal vertices {len(inner)}: labelled {len(inner) - len(todo)}, unreachable {len(todo)}; by tracking {len(inner) - len(todo) - len(approx)} (edges tracked {tracked}), '
      f'nearest-root fallback {len(approx)} (min |y| among them {min((abs(Y[w]) for w in approx), default=0):.4f}, max {max((abs(Y[w]) for w in approx), default=0):.4f}); '
      f'vertices where two labels received the same root (to 1e-14) {dup}; {time.time()-t0:.0f} s', flush=True)
np.savez('vx_labels_full.npz', X=X, seen=np.array(sorted(seen)), approx=np.array(approx, int))
