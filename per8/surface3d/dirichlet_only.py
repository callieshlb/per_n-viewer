#!/usr/bin/env python3
"""E93-U4: the Dirichlet stage, from develop_pre.npz (positions + side-pairing generators as 30-digit strings), run in a mesh
directory. Replaces E92's word search (complete to 2 r_D + 0.25 + r_D: ~1e6+ elements at genus 24) by fuchs.closure_mp
(200 bits; certification as E92: area = 4 pi (g - 1) of the polygon cut out by group elements). U4 face checks as E92.
Writes develop.npz in E92's format; els_all = the kept short elements (NOT a complete list: els_bound = -1; downstream E93
scripts use folding into Dirichlet polygons instead of 'all elements up to a bound')."""
import json
import os
import sys
import time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import fuchs as fu
from gmpy2 import mpc, mpfr


def run(AREA, lg, extra=None):
    t0 = time.time()
    D = np.load('develop_pre.npz'); zc, G, cyc, Tn, orig, root = D['zc'], list(D['G']), D['cyc'], D['Tn'], D['orig'], int(D['root'])
    gs = json.loads(str(D['gens_str']))
    fu.mp_ctx(200)
    gens = []
    for g in gs:
        m = [mpc(mpfr(re), mpfr(im)) for re, im in g]
        gens.append(fu.mp_norm(((m[0] + m[3].conjugate())/2, (m[1] + m[2].conjugate())/2)))
    dg = [fu.mp_d(g)[0] for g in gens]
    print(f'loaded develop_pre.npz: generators {len(gens)} | displacement min {min(dg):.4f} max {max(dg):.4f}', flush=True)
    import gmpy2, u4_bfs                                             # E94: nearest-lift seed + capped pair rounds (u4_bfs.py)
    fu.IDCUT = gmpy2.cosh(mpfr('5e-4'))
    lg_ = lambda *a, **k: print(*a, flush=True)
    keep, res = u4_bfs.dirichlet_from_seed(u4_bfs.bfs_candidates(D, log=lg_), AREA, log=lg_)
    print(f'closure done in {time.time()-t0:.1f} s', flush=True)
    nf = res['nfaces']
    if not (res['bounded'] and abs(res['area'] - AREA) < 1e-8):
        print('U4: closure did not reach a bounded polygon of area 4 pi (g - 1):', res['bounded'], res['area'])
        print('U4: not certified - develop.npz NOT written', flush=True); sys.exit(5)     # E94 run 5: never write an uncertified D
    if not res['bounded']:
        print('U4: polygon unbounded - no face checks possible', flush=True); sys.exit(4)
    perr, unmatched, pair_of, relerr = fu.face_pairing_mp(res)          # E94: face checks in 200 bits
    cls, sums = fu.vertex_cycles_mp(res, pair_of)
    sums = [float(x) for x in sums]
    print(f'U4 (200 bits): area - 4 pi (g - 1) = {float(res["area_mp"] - AREA):.3e} | min edge length {res["min_edge"]:.3e} | '
          f'max relation error |g_a g_partner - id| {relerr:.3e} | angle sums / 2 pi - 1 in [{min(sums)/(2*np.pi) - 1:.3e}, '
          f'{max(sums)/(2*np.pi) - 1:.3e}]', flush=True)
    print(f'U4: faces {nf}, area {res["area"]:.12f} (diff {res["area"] - AREA:.3e}), circumradius {res["rD"]:.6f} | max endpoint '
          f'mismatch after g^-1 {perr:.3e} | faces whose partner element is not g^-1: {unmatched} | vertex cycles {len(cls)}, '
          f'V - E + F = {len(cls) - nf//2 + 1}, angle sums / 2 pi in [{min(sums)/(2*np.pi):.10f}, {max(sums)/(2*np.pi):.10f}], '
          f'cycle lengths {sorted(set(len(m) for m in cls.values()))}', flush=True)
    els = np.array([fu.mp_to_np(e) for e in keep])
    fmp = res['face_mp']                 # face elements in 200 bits (for fuchs.dirichlet_at_mp)
    face_str = [[f'{e[0].real:.50e}', f'{e[0].imag:.50e}', f'{e[1].real:.50e}', f'{e[1].imag:.50e}'] for e in fmp]
    assert np.abs(np.array([fu.mp_to_np(e) for e in fmp]) - res['face_g']).max() < 1e-6*np.abs(res['face_g']).max()
    np.savez('develop.npz', els_all=els, els_bound=-1.0, z=zc, Tn=Tn, orig=orig, root=root, gens=np.array(G), verts=res['verts'],
             ftag=np.array(res['tag']), face_g=res['face_g'], pair_of=pair_of, cyc=cyc,
             face_g_str=np.array(json.dumps(face_str)))
    if extra: lg.update(extra)
    lg.update(U4_area=float(res['area']), U4_area_diff=float(res['area'] - AREA), U4_faces=nf, U4_pair_err=float(perr),
              U4_unmatched=unmatched, U4_cycles=len(cls), U4_rD=res['rD'], kept=len(keep), U4_relerr=relerr,
              U4_min_edge=res['min_edge'], U4_anglesum_dev=float(max(abs(x/(2*np.pi) - 1) for x in sums)))
    json.dump(lg, open('develop_log.json', 'w'), indent=1)
    print('wrote develop.npz')


if __name__ == '__main__':
    GEN = json.load(open(os.path.join(HERE, '..', 'branch.json')))['genus']
    run(4*np.pi*(GEN - 1), json.load(open('develop_log.json')) if os.path.exists('develop_log.json') else {})
