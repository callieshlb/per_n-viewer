"""E95: certified vertex roots (python-flint acb_poly.roots, 512 bits) for trackmp curves; same as options_probe.py's froots."""
import flint
from gmpy2 import mpc, mpfr
USED = {}
def froots(C, y, prec=512):
    """escalates the working precision (coefficient evaluation and isolation, maxprec = 4 x) 512 -> 1024 -> 2048 -> 4096 bits when acb_poly.roots cannot isolate (needed near y = 0: 2048 bits at |y| = 0.0048); USED counts the precision used."""
    for p in (prec, 2*prec, 4*prec, 8*prec):
        try:
            r = _froots(C, y, p); USED[p] = USED.get(p, 0) + 1; flint.ctx.prec = prec; return r
        except ValueError:
            continue
    raise ValueError('flint roots failed up to %d bits' % (8*prec))
def _froots(C, y, prec):
    flint.ctx.prec = prec
    yy = flint.acb(complex(y).real, complex(y).imag); cs = []
    for row in C.coef_int:
        v = flint.acb(0)
        for a in reversed(row): v = v*yy + a
        cs.append(v)
    return flint.acb_poly(cs).roots(maxprec=4*prec)
def to_mpc(r):
    m = r.mid(); return mpc(mpfr(m.real.str(160, radius=False)), mpfr(m.imag.str(160, radius=False)))
