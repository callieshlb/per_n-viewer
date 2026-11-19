#!/usr/bin/env python3
"""E94 step P1a (E93 disc.py with the n = 8 curve): discriminant of P (per8_xy_poly.txt here) in x, factored over Q; cached in
disc_factors.txt."""
import time
import sympy as sp
x, y = sp.symbols('x y')
P = sp.sympify(open('per8_xy_poly.txt').read().replace('^', '**'))
print('deg_x', sp.Poly(P, x).degree(), 'deg_y', sp.Poly(P, y).degree(), 'terms', len(sp.Poly(P, x, y).terms()), flush=True)
for k, c in enumerate(sp.Poly(P, x).all_coeffs()[:2] + sp.Poly(P, x).all_coeffs()[-1:]): print('coefficient', ['x^20', 'x^19', 'x^0'][k], sp.factor(c), flush=True)
t0 = time.time()
disc = sp.discriminant(sp.Poly(P, x))
disc = sp.Poly(disc.as_expr(), y)
print('discriminant degree in y', disc.degree(), 'time', round(time.time()-t0, 1), flush=True)
c, facs = sp.factor_list(disc.as_expr())
print('content', c, 'factors', [(sp.Poly(f, y).degree(), e) for f, e in facs], 'time', round(time.time()-t0, 1), flush=True)
with open('disc_factors.txt', 'w') as fh:
    for f, e in facs:
        fh.write(f'{e} {str(sp.expand(f))}\n')
