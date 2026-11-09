// E96 Per_n atlas: double-precision maps of every model, marked cycles, the free critical point, its fate, and the
// Weierstrass wp of the Per_5 torus. Shared by the shell, the plane/torus pane and node tests.
// Every map is f(z) = (A z^2 + B z + C)/(D z^2 + E z + F) with complex coefficients; z = null is infinity.
// Families (normal forms; see REGISTRATION.md):
//   quad  (n = 1): z^2 + c, parameter c, marked cycle {inf}, free critical point 0.
//   per2  (n = 2): f_t(z) = t(2z + 1)/(8 z^2), marked 0 -> inf -> 0, free critical point -1 (Milnor a/(z^2+2z): a = 8/t).
//   cor   (n >= 3): Caroline's chart f(z) = (z - x)(z - w)/z^2, w = (x + y - 1)/(x - 1), marked 0 -> inf -> 1 -> y -> ... -> x.
//         n = 3: parameter x, y = 0, w = 1.  n = 4: parameter y, x = y.  n = 5: torus coordinate u.  n >= 6: (x, y) on the curve.
'use strict';
var PM = (function () {
  function cmul(a, b) { return [a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]]; }
  function cdiv(a, b) { var d = b[0]*b[0]+b[1]*b[1]; return [(a[0]*b[0]+a[1]*b[1])/d, (a[1]*b[0]-a[0]*b[1])/d]; }
  function cadd(a, b) { return [a[0]+b[0], a[1]+b[1]]; }
  function csub(a, b) { return [a[0]-b[0], a[1]-b[1]]; }
  function cabs(a) { return Math.hypot(a[0], a[1]); }
  function cscale(a, s) { return [a[0]*s, a[1]*s]; }
  function chord(a, b) {
    if (a === null && b === null) return 0;
    if (a === null) { a = b; b = null; }
    if (b === null) return 1/Math.sqrt(1+a[0]*a[0]+a[1]*a[1]);
    var dr = a[0]-b[0], di = a[1]-b[1];
    return Math.sqrt(dr*dr+di*di)/Math.sqrt((1+a[0]*a[0]+a[1]*a[1])*(1+b[0]*b[0]+b[1]*b[1]));
  }
  // ---- Per_5 torus: X = wp(u), r = wp'(u) for g2 = 11/4, g3 = -3/8 (rectangular lattice W1, W2 i)
  var T5 = { g2: 11/4, g3: -3/8, W1: 3.09415950710224035, W2: 2.74573911808975367 };
  // wp(u) - e1 = (pi/W1)^2 (theta3(0) theta4(0) theta2(v)/theta1(v))^2, v = pi u/W1, nome q = exp(-pi W2/W1) ~ 0.0616
  // (the duplication formula was tried first and dropped: ~1 digit lost per doubling near the pole; log.md).
  // wp' = +-sqrt(4 wp^3 - g2 wp - g3), sign from a central difference of wp.
  T5.e1 = 0.75; T5.q = Math.exp(-Math.PI*T5.W2/T5.W1);
  function csin(z) { return [Math.sin(z[0])*Math.cosh(z[1]), Math.cos(z[0])*Math.sinh(z[1])]; }
  function ccos(z) { return [Math.cos(z[0])*Math.cosh(z[1]), -Math.sin(z[0])*Math.sinh(z[1])]; }
  function wpX(u) {
    var W1 = T5.W1, W2 = T5.W2, q = T5.q;
    var a = u[0] - W1*Math.round(u[0]/W1), b = u[1] - W2*Math.round(u[1]/W2);
    if (a === 0 && b === 0) return null;
    var v = [Math.PI*a/W1, Math.PI*b/W1], t1 = [0, 0], t2 = [0, 0], t3 = 1, t4 = 1;
    for (var n = 0; n < 6; n++) {
      var qq = Math.pow(q, (n + .5)*(n + .5)), arg = [(2*n+1)*v[0], (2*n+1)*v[1]];
      t1 = cadd(t1, cscale(csin(arg), 2*qq*(n % 2 ? -1 : 1))); t2 = cadd(t2, cscale(ccos(arg), 2*qq));
      if (n > 0) { var qn = Math.pow(q, n*n); t3 += 2*qn; t4 += 2*qn*(n % 2 ? -1 : 1); }
    }
    var R = cscale(cdiv(t2, t1), t3*t4), k = Math.PI/W1;
    return cadd([T5.e1, 0], cscale(cmul(R, R), k*k));
  }
  function csqrt(z) { var r = Math.hypot(z[0], z[1]), t = Math.atan2(z[1], z[0])/2, s = Math.sqrt(r); return [s*Math.cos(t), s*Math.sin(t)]; }
  function wp(u) {
    var X = wpX(u); if (X === null) return null;
    var r = csqrt(csub(csub(cscale(cmul(X, cmul(X, X)), 4), cscale(X, T5.g2)), [T5.g3, 0]));
    var h = 1e-5*(1 + cabs(u)), Xp = wpX([u[0]+h, u[1]]), Xm = wpX([u[0]-h, u[1]]);
    if (Xp && Xm) { var fd = cscale(csub(Xp, Xm), 1/(2*h)); if (cabs(csub(fd, r)) > cabs(cadd(fd, r))) r = [-r[0], -r[1]]; }
    return [X, r];
  }
  function torusXY(u) {                                     // u -> (x, y) on P5 = x^2 A - x B + C
    var w = wp(u);
    if (w === null) return { x: [1, 0], y: [1, 0], X: null, r: null };
    var X = w[0], r = w[1], s = cadd(X, [1.25, 0]);         // y = 1 - 1/(X + 5/4)
    if (cabs(s) < 1e-300) return { x: null, y: null, X: X, r: r };
    var y = csub([1, 0], cdiv([1, 0], s)), y1 = csub(y, [1, 0]), v = cmul(r, cmul(y1, y1));
    var y2 = cmul(y, y), A = cadd(csub(y2, [1, 0]), y), B = csub(cadd(cscale(y2, 2), y), [1, 0]), C = csub(cscale(y2, 2), y);
    var Bp = cadd(B, v), Bm = csub(B, v), x;
    if (cabs(Bp) >= cabs(Bm)) x = cabs(A) === 0 ? null : cdiv(Bp, cscale(A, 2)); else x = cdiv(cscale(C, 2), Bm);
    return { x: x, y: y, X: X, r: r };
  }
  // ---- maps
  function corMap(x, y, n) {
    var w = cdiv([x[0]+y[0]-1, x[1]+y[1]], [x[0]-1, x[1]]);
    var s = cadd(x, w), p = cmul(x, w);
    var M = { n: n, A: [1, 0], B: [-s[0], -s[1]], C: p, D: [1, 0], E: [0, 0], F: [0, 0], x: x, y: y, w: w };
    M.omega = cdiv(cscale(p, 2), s);
    M.marked = [[0, 0], null, [1, 0]];
    if (n >= 4) M.marked.push(y);
    for (var k = 4; k < n-1; k++) M.marked.push(fz(M, M.marked[k-1]));
    if (n >= 5) M.marked.push(x);
    return M;
  }
  // from a parameter message: family, n, and p (plane / torus coordinate) or x, y (curve point)
  function mapFrom(family, n, q) {
    var M;
    if (family === 'quad') { var c = q.p; M = { n: 1, A: [1, 0], B: [0, 0], C: c, D: [0, 0], E: [0, 0], F: [1, 0], omega: [0, 0], marked: [null] }; }
    else if (family === 'per2') { var t = q.p; M = { n: 2, A: [0, 0], B: cscale(t, 0.25), C: cscale(t, 0.125), D: [1, 0], E: [0, 0], F: [0, 0], omega: [-1, 0], marked: [[0, 0], null] }; }
    else if (family === 'cor3') M = corMap(q.p, [0, 0], 3);
    else if (family === 'cor4') M = corMap(q.p, q.p, 4);
    else if (family === 'torus5') { var P = torusXY(q.p); if (!P.x || !P.y) return null; M = corMap(P.x, P.y, 5); }
    else if (family === 'cor') M = corMap(q.x, q.y, n);
    else throw new Error('unknown family ' + family);
    M.scale = M.marked.map(function (a, j) {
      var m = 1; for (var i = 0; i < M.marked.length; i++) if (i !== j) m = Math.min(m, chord(a, M.marked[i])); return m; });
    return M;
  }
  function fz(M, z) {
    var num, den;
    if (z === null || cabs(z) > 1) {
      var w = z === null ? [0, 0] : cdiv([1, 0], z), w2 = cmul(w, w);
      num = cadd(cadd(M.A, cmul(M.B, w)), cmul(M.C, w2)); den = cadd(cadd(M.D, cmul(M.E, w)), cmul(M.F, w2));
    } else {
      var z2 = cmul(z, z);
      num = cadd(cadd(cmul(M.A, z2), cmul(M.B, z)), M.C); den = cadd(cadd(cmul(M.D, z2), cmul(M.E, z)), M.F);
    }
    if (cabs(den) <= 1e-300*cabs(num) || (cabs(den) === 0)) return null;
    var r = cdiv(num, den);
    return cabs(r) > 1e150 ? null : r;
  }
  // fate of z0: {kind: 'capture', j, t, phase} | {kind: 'cycle', period, t} | {kind: 'undecided'}  (as CORE.fate, rel = true)
  function fate(M, z, maxIt, delta, tolP) {
    var n = M.marked.length, ring = new Array(32), th = M.scale.map(function (s) { return delta*s; });
    for (var it = 0; it < maxIt; it++) {
      z = fz(M, z);
      for (var j = 0; j < n; j++) if (chord(z, M.marked[j]) < th[j]) return { kind: 'capture', j: j, t: it+1, phase: ((j-(it+1)) % n + n) % n };
      ring[it % 32] = z;
    }
    for (var per = 1; per <= 24; per++) {
      var zp = ring[(maxIt-1-per+3200) % 32], zl = ring[(maxIt-1) % 32];
      if (zp !== undefined && chord(zl, zp) < tolP) return { kind: 'cycle', period: per, t: maxIt };
    }
    return { kind: 'undecided' };
  }
  function orbit(M, z, N) { var out = [z]; for (var i = 0; i < N; i++) { z = fz(M, z); out.push(z); } return out; }
  function closure(M) { var a = M.marked[0], z = a; for (var i = 0; i < M.n; i++) z = fz(M, z); return chord(z, a); }
  function fmt(z, d) { d = d || 6; if (z === null) return '∞'; var f = function (v) { return (Math.abs(v) < 1e-12 ? 0 : v).toPrecision(d); };
    return f(z[0]) + (z[1] < 0 ? ' − ' : ' + ') + f(Math.abs(z[1])) + 'i'; }
  function fateText(f) { return f.kind === 'capture' ? 'captured by the marked cycle, phase ' + f.phase + ', time ' + f.t
    : f.kind === 'cycle' ? 'attracting cycle, period ' + f.period : 'undecided (bifurcation locus)'; }
  return { cmul: cmul, cdiv: cdiv, cadd: cadd, csub: csub, cabs: cabs, chord: chord, T5: T5, wp: wp, torusXY: torusXY,
           corMap: corMap, mapFrom: mapFrom, fz: fz, fate: fate, orbit: orbit, closure: closure, fmt: fmt, fateText: fateText };
})();
if (typeof module !== 'undefined') module.exports = PM;
