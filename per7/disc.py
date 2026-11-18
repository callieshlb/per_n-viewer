#!/usr/bin/env python3
"""E93 step A0: discriminant of P (per7_xy_poly.txt) in x, factored over Q; cached in disc_factors.txt."""
import time
import sympy as sp
x, y = sp.symbols('x y')
P = sp.sympify(open('../inputs/E90-real-circles/per7_xy_poly.txt').read().replace('^', '**'))
t0 = time.time()
disc = sp.discriminant(sp.Poly(P, x))
disc = sp.Poly(disc.as_expr(), y)
print('discriminant degree in y', disc.degree(), 'time', round(time.time()-t0, 1), flush=True)
c, facs = sp.factor_list(disc.as_expr())
print('content', c, 'factors', [(sp.Poly(f, y).degree(), e) for f, e in facs], 'time', round(time.time()-t0, 1), flush=True)
with open('disc_factors.txt', 'w') as fh:
    for f, e in facs:
        fh.write(f'{e} {str(sp.expand(f))}\n')
