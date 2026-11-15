#!/usr/bin/env python3
"""E92 step A: branch points of y over P7, lasso monodromy at a base point, monodromy at infinity,
Riemann-Hurwitz genus. Writes branch.json. Raw values only."""
import os
import cmath, json, math
import numpy as np
import sympy as sp
x, y = sp.symbols('x y')
P = sp.sympify(open(os.path.dirname(os.path.abspath(__file__)) + '/../inputs/E90-real-circles/per6_xy_poly.txt').read().replace('^', '**'))
Pp = sp.Poly(P, x)
d = Pp.degree()
cf = [sp.lambdify(y, c, 'numpy') for c in Pp.all_coeffs()]
disc = sp.discriminant(P, x)
lead = Pp.all_coeffs()[0]
facs = sp.factor_list(disc)[1] + sp.factor_list(lead)[1]
bps = []
for f, e in facs:
    for r in sp.Poly(f, y).nroots(n=40, maxsteps=300):
        r = complex(r)
        if all(abs(r-o) > 1e-10 for o in bps):
            bps.append(r)
y0 = complex(0.31, -0.47)          # generic base point
def roots(t):
    return np.roots([complex(c(t)) if callable(c) else complex(c) for c in cf])
def track(path, cur):
    cur = np.array(cur)
    for k in range(1, len(path)):
        a, b = path[k-1], path[k]
        sub = 1
        while True:
            tmp, ok = cur.copy(), True
            for j in range(1, sub+1):
                r = roots(a+(b-a)*j/sub)
                dist = np.abs(tmp[:, None]-r[None, :])
                idx = dist.argmin(1)
                s = np.sort(dist, 1)
                if len(set(idx)) < d or np.any(4*s[:, 0] > s[:, 1]):
                    ok = False
                    break
                tmp = r[idx]
            if ok:
                cur = tmp
                break
            sub *= 2
            if sub > 2**15:
                raise RuntimeError('tracking failed')
    return cur
x0 = roots(y0)
x0 = x0[np.lexsort((x0.imag, x0.real))]
def perm_of(x1):
    return [int(np.argmin(np.abs(x0-v))) for v in x1]
perms = []
for b in bps:
    others = [abs(b-o) for o in bps if o != b]
    rad = min(min(others)/3, 0.2) if others else 0.2
    st = b + rad*(y0-b)/abs(y0-b)
    m = 600
    line = [y0+(st-y0)*t/m for t in range(m+1)]
    th0 = cmath.phase(st-b)
    circ = [b+rad*cmath.exp(1j*(th0+2*math.pi*t/m)) for t in range(m+1)]
    perms.append(perm_of(track(line+circ[1:]+line[::-1][1:], x0)))
R = 4*max(abs(b) for b in bps) + 4
st = y0 + (R)*1
line = [y0+(st-y0)*t/800 for t in range(801)]
circ = [y0+R*cmath.exp(1j*2*math.pi*t/4000) for t in range(4001)]
pinf = perm_of(track(line+circ[1:]+line[::-1][1:], x0))
def ncyc(p):
    seen, c = set(), 0
    for i in range(d):
        if i not in seen:
            c += 1
            j = i
            while j not in seen:
                seen.add(j); j = p[j]
    return c
Rfin = sum(d-ncyc(p) for p in perms)
Rinf = d-ncyc(pinf)
g = (Rfin+Rinf-2*d+2)/2
# consistency: product of lassos in order of increasing argument about y0 equals the big loop
order = sorted(range(len(bps)), key=lambda k: cmath.phase(bps[k]-y0) % (2*math.pi))
comp = list(range(d))
for k in order:
    comp = [perms[k][comp[i]] for i in range(d)]
print('degree in x', d, '| distinct branch points (incl. nodes, leading-coefficient roots)', len(bps))
print('non-trivial lassos', sum(1 for p in perms if p != list(range(d))), '| ramification finite', Rfin,
      '| at infinity', Rinf, '| loop at infinity perm', pinf)
print('product of lassos (counterclockwise order about y0)', comp, '| equals loop at infinity:', comp == pinf)
print('Riemann-Hurwitz genus', g)
json.dump(dict(y0=[y0.real, y0.imag], x0=[[v.real, v.imag] for v in x0], R=R,
               bps=[[b.real, b.imag] for b in bps], perms=perms, pinf=pinf, genus=g),
          open('branch.json', 'w'), indent=1)
