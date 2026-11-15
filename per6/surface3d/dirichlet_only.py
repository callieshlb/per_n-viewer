#!/usr/bin/env python3
"""E92-U4: rerun only the Dirichlet stage of develop.py from develop_pre.npz (positions + generators already computed),
writing develop.npz in the same format. Usage: python3 dirichlet_only.py (in a mesh directory)."""
import json
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import develop as dv
from develop import hd0, vorb, shorten, words_bfs_vec, dirichlet, poly_area, to_poinc, mob_ap, hdist, ident_err
D = np.load('develop_pre.npz'); zc, G, cyc, Tn, orig, root = D['zc'], list(D['G']), D['cyc'], D['Tn'], D['orig'], int(D['root'])
G2 = G + [np.linalg.inv(g) for g in G]
print('loaded develop_pre.npz: generators', len(G), '| max generator displacement', f'{max(hd0(mob_ap(g, 0)) for g in G):.4f}')
S, rA = shorten(G2)
for it in range(6):
    Rb = 2*rA + 0.25
    dS = hd0(vorb(S)); SB = list(S[dS <= Rb])
    els, nexp, capped = words_bfs_vec(SB, Rb + rA + 0.25, Rb)
    verts, ftag, H = dirichlet(els, Rb)
    bounded = all(t >= 0 for t in ftag)
    area = poly_area(verts)[0] if bounded else np.inf
    rN = float(hd0(to_poinc(verts)).max()) if bounded else np.inf
    expect = (np.cosh(Rb + rA + 0.25) - 1)/10
    print(f'U4 stage B{it}: generators {len(SB)} (kept elements with d <= {Rb:.3f}), words explored {nexp} (area estimate {expect:.0f}; capped {capped}), '
          f'kept {len(els)}, half-planes {len(H)}, faces {len(ftag)}, area {area:.12f} (diff {area - 20*np.pi:.2e}), circumradius {rN:.6f}')
    if bounded and abs(area - 20*np.pi) < 1e-8 and rN <= rA + 1e-9: break
    rA = rN; S = np.array(els[1:])
nf = len(verts); faces = []
for a_ in range(nf):
    g = H[ftag[a_]][2]; faces.append(dict(g=g, P=to_poinc(verts[a_]), Q=to_poinc(verts[(a_ + 1) % nf])))
perr, unmatched, pair_of = 0.0, 0, [-1]*nf
for a_, fa in enumerate(faces):
    gi = np.linalg.inv(fa['g']); P2, Q2 = mob_ap(gi, fa['P']), mob_ap(gi, fa['Q'])
    best, bi = np.inf, -1
    for b_, fb in enumerate(faces):
        d = max(hdist(P2, fb['Q']), hdist(Q2, fb['P']))
        if d < best: best, bi = d, b_
    pair_of[a_] = bi; perr = max(perr, best)
    if ident_err(fa['g'] @ faces[bi]['g']) > 1e-6: unmatched += 1
print('U4: faces', nf, '| max endpoint mismatch after g^-1', f'{perr:.3e}', '| faces whose partner element is not g^-1:', unmatched)
np.savez('develop.npz', els_all=dv.LAST_EXPLORED[0], els_bound=dv.LAST_EXPLORED[1], z=zc, Tn=Tn, orig=orig, root=root, gens=np.array(G),
         verts=verts, ftag=np.array(ftag), face_g=np.array([f['g'] for f in faces]), pair_of=np.array(pair_of), cyc=cyc)
lg = json.load(open('develop_log.json')) if os.path.exists('develop_log.json') else {}
lg.update(U4_area=float(area), U4_faces=nf, U4_pair_err=float(perr), U4_unmatched=unmatched, U4_explored=nexp)
json.dump(lg, open('develop_log.json', 'w'), indent=1)
print('wrote develop.npz')
