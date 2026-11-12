// E96 Per_n atlas: double-precision maps of every model, marked cycles, the free critical point, its fate, and the
// Weierstrass wp of the Per_5 torus. Shared by the shell, the plane/torus pane and node tests.
// Every map is f(z) = (A z^2 + B z + C)/(D z^2 + E z + F) with complex coefficients; z = null is infinity.
// Families (normal forms; see REGISTRATION.md):
//   quad  (n = 1): z^2 + c, parameter c, marked cycle {inf}, free critical point 0.
//   per2c (n = 2, Caroline 2026-11-09): f(z) = (z^2 - c)/(z^2 - 1), marked inf -> 1 -> inf (inf critical), free critical point 0.
//   cor3c (n = 3, Caroline): f(z) = (z - 1)(z - c/(2 - c))/z^2, parameter c = the free critical point (x = c/(2 - c)).
//   cor4c (n = 4, Caroline): y = 0.1/c + (sqrt 5 - 1)/2, f(z) = (z - y)(z - (2y - 1)/(y - 1))/z^2.
//   per2  (n = 2, superseded on the site): f_t(z) = t(2z + 1)/(8 z^2), marked 0 -> inf -> 0, free critical point -1.
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
  // (x, y) on P5 -> u (inverse of torusXY): X = 1/(1 - y) - 5/4, r = v/(y - 1)^2, v = 2(y^2 + y - 1)x - (2y^2 + y - 1);
  // Newton on wpX(u) = X from a grid of starts, then u or -u by the sign of wp' (as torus5.py's locate)
  function torusLocate(x, y) {
    var X = cadd(cdiv([1, 0], csub([1, 0], y)), [-1.25, 0]), y2 = cmul(y, y);
    var v = csub(cmul(cscale(cadd(csub(y2, [1, 0]), y), 2), x), csub(cadd(cscale(y2, 2), y), [1, 0])), y1 = csub(y, [1, 0]), r = cdiv(v, cmul(y1, y1));
    var best = null, bd = Infinity;
    for (var a = 0; a < 10; a++) for (var b = 0; b < 10; b++) {
      var u = [(a + .5)/10*T5.W1, (b + .5)/10*T5.W2];
      for (var it = 0; it < 60; it++) { var w = wp(u); if (!w) break; var du = cdiv(csub(w[0], X), w[1]); u = csub(u, du); if (cabs(du) < 1e-15) break; }
      var w2 = wp(u); if (!w2) continue; var e = cabs(csub(w2[0], X))/(1 + cabs(X));
      if (e < 1e-10) { var su = cabs(csub(w2[1], r)) <= cabs(cadd(w2[1], r)) ? u : [-u[0], -u[1]];
        su = [su[0] - T5.W1*Math.floor(su[0]/T5.W1), su[1] - T5.W2*Math.floor(su[1]/T5.W2)];
        var P = torusXY(su), d = (P.x ? cabs(csub(P.x, x)) : 1e9) + cabs(csub(P.y, y)); if (d < bd) { bd = d; best = su; } }
    }
    return best ? { u: best, err: bd } : null;
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

  // ---- external rays for z^2 + c (Per_1). Angle theta = p/q (BigInt, exact doubling).
  // Parameter ray: points c with Phi_M(c) = rho e^(2 pi i theta), Phi_M(c) ~ (f_c^n(c))^(1/2^n): Newton on f_c^n(c) = w.
  // Dynamic ray of f_c: points z with phi_c(z) = rho e^(2 pi i theta): Newton on f_c^n(z) = w.
  // Potential levels log rho_k = log(rho0) 2^(-k/S); n = ceil(k/S) keeps |w| in [rho0, rho0^2]. Stops when consecutive points
  // are closer than 1e-10 (landed, to the precision of the ray's approach) or at depth 50.
  function parseAngle(s) {
    s = String(s).trim(); var m;
    if ((m = s.match(/^(\d+)\s*\/\s*(\d+)$/))) { var p = BigInt(m[1]), q = BigInt(m[2]); if (q === 0n) return null; return { p: p % q, q: q, text: s }; }
    if ((m = s.match(/^0?\.(\d+)$/))) { var q2 = 10n ** BigInt(m[1].length); return { p: BigInt(m[1]) % q2, q: q2, text: s }; }
    return null;
  }
  function angAt(a, n) { var t = (a.p * (2n ** BigInt(n))) % a.q; return Number(t)/Number(a.q); }
  function traceRay(kind, a, c, S, kmax) {
    S = S || 8; kmax = kmax || 50*S;
    var rho0 = 100, L0 = Math.log(rho0), t0 = 2*Math.PI*angAt(a, 0);
    var v = [rho0*Math.cos(t0), rho0*Math.sin(t0)], out = [v], landed = false;
    for (var k = 1; k <= kmax; k++) {
      var logr = L0*Math.pow(2, -k/S), n = Math.ceil(k/S), th = 2*Math.PI*angAt(a, n), mod = Math.exp(Math.pow(2, n)*logr);
      var w = [mod*Math.cos(th), mod*Math.sin(th)], ok = false;
      for (var it = 0; it < 30; it++) {
        var z = v, dz = [1, 0];
        for (var j = 0; j < n; j++) {
          dz = kind === 'param' ? cadd(cscale(cmul(z, dz), 2), [1, 0]) : cscale(cmul(z, dz), 2);
          z = cadd(cmul(z, z), kind === 'param' ? v : c);
        }
        var step = cdiv(csub(z, w), dz); v = csub(v, step);
        if (!isFinite(v[0]) || !isFinite(v[1])) return { pts: out, landed: false, failed: true };
        if (cabs(step) < 1e-14*(1 + cabs(v))) { ok = true; break; }
      }
      out.push(v);
      if (!ok) return { pts: out, landed: false, failed: true };
      if (k > 4*S && cabs(csub(v, out[out.length - 2])) < 1e-10) { landed = true; break; }
    }
    return { pts: out, landed: landed, end: out[out.length - 1], depth: Math.ceil(out.length/S) };
  }

  // exact landing points for rational angles (the traced ray reaches a parabolic point only very slowly).
  // period of the angle under doubling: p >= 1 with 2^p p == p (mod q); preperiodic part -> null period.
  function anglePeriod(a) { var q = a.q; if (q % 2n === 0n) return null; var t = (a.p*2n) % q, p = 1; while (t !== a.p % q && p < 4096) { t = (t*2n) % q; p++; } return p; }
  function fpDerivs(z, c, p) {     // f_c^p(z) and d/dz, d/dc, d2/dz2, d2/dzdc
    var a = [1, 0], b = [0, 0], s2 = [0, 0], m = [0, 0];
    for (var j = 0; j < p; j++) {
      var z2 = cscale(z, 2);
      s2 = cadd(cscale(cmul(a, a), 2), cmul(z2, s2)); m = cadd(cscale(cmul(b, a), 2), cmul(z2, m));
      a = cmul(z2, a); b = cadd(cmul(z2, b), [1, 0]); z = cadd(cmul(z, z), c);
    }
    return { z: z, a: a, b: b, s: s2, m: m };
  }
  function dynLanding(c, z0, p) {                 // periodic point of period dividing p near z0 (Newton on f^p(z) = z)
    var z = z0;
    for (var it = 0; it < 100; it++) { var d = fpDerivs(z, c, p), st = cdiv(csub(d.z, z), csub(d.a, [1, 0])); z = csub(z, st); if (cabs(st) < 1e-15*(1 + cabs(z))) break; }
    var d2 = fpDerivs(z, c, p); return { z: z, res: cabs(csub(d2.z, z)), mult: d2.a };
  }
  function paramLanding(cE, a, p) {               // root: the p-cycle met by the dynamic ray, multiplier continued to 1
    var zr = traceRay('dyn', a, cE).pts, z = dynLanding(cE, zr[zr.length - 1], p).z, c = cE;
    var lam0 = fpDerivs(z, c, p).a, N = 200;
    for (var i = 1; i <= N; i++) {
      var t = i/N, tt = 1 - (1 - t)*(1 - t)*(1 - t), lam = cadd(cscale(lam0, 1 - tt), [tt*(1 - 1e-12), 0]);
      for (var it = 0; it < 60; it++) {
        var d = fpDerivs(z, c, p), F1 = csub(d.z, z), F2 = csub(d.a, lam);
        var J11 = csub(d.a, [1, 0]), J12 = d.b, J21 = d.s, J22 = d.m, det = csub(cmul(J11, J22), cmul(J12, J21));
        var dz = cdiv(csub(cmul(J22, F1), cmul(J12, F2)), det), dc = cdiv(csub(cmul(J11, F2), cmul(J21, F1)), det);
        z = csub(z, dz); c = csub(c, dc);
        if (cabs(dz) + cabs(dc) < 1e-15) break;
      }
    }
    var e = fpDerivs(z, c, p); return { c: c, z: z, res: cabs(csub(e.z, z)), mult: e.a };
  }
  function rayWithLanding(kind, a, c) {
    var r = traceRay(kind, a, c), p = anglePeriod(a), end = r.pts[r.pts.length - 1];
    if (p !== null && !r.landed) {
      if (kind === 'dyn') { var L = dynLanding(c, end, p); r.landing = L.z; r.how = 'Newton on f^' + p + '(z) = z, multiplier ' + fmt(L.mult, 4); }
      else { var P = paramLanding(end, a, p); r.landing = P.c; r.how = 'root: period-' + p + ' cycle with multiplier continued to 1 (' + fmt(P.mult, 6) + ')'; }
      r.jump = cabs(csub(r.landing, end));
    } else { r.landing = end; r.how = r.landed ? 'ray converged (|step| < 1e-10)' : 'end of traced ray (depth 50)'; }
    r.period = p; return r;
  }

  // ---- Hubbard tree of a center z^2 + c (period n): union of the regulated arcs [0, P^j(0)], built recursively from the
  // kneading sequence nu (nu[i-1] = side of P^i(0), i = 1..n-1, w.r.t. the two preimage rays of the characteristic angle):
  //   [0, P^j]   = the lift through P^j of [c, P^(j+1)] (f is injective on it: 0 is an end point);
  //   [P^i, P^j] = the lift through P^i of [P^(i+1), P^(j+1)] if nu_i == nu_j, else [P^i, 0] + [0, P^j];
  // depth 0 = straight segment. Indices mod n (P^n = P^0 = 0). Memoised; ~depth 14.
  function hubbardTree(c, n, nu, depth) {
    depth = depth === undefined ? Math.min(160, 14*n + 10) : depth;      // direction error halves once per period: depth ~ 14 n
    var P = [[0, 0]]; for (var i = 1; i < n; i++) { var q = P[i-1]; P.push(cadd(cmul(q, q), c)); }
    function side(i) { return nu[(i % n) - 1]; }
    function csq(z) { var r = Math.hypot(z[0], z[1]), t = Math.atan2(z[1], z[0])/2, s = Math.sqrt(r); return [s*Math.cos(t), s*Math.sin(t)]; }
    function resample(pts, N) {
      var L = [0]; for (var i = 1; i < pts.length; i++) L.push(L[i-1] + cabs(csub(pts[i], pts[i-1])));
      var T = L[L.length-1], out = [pts[0]], j = 1; if (T === 0) return [pts[0], pts[pts.length-1]];
      for (var k = 1; k < N; k++) { var t = T*k/N; while (j < pts.length - 1 && L[j] < t) j++;
        var u = (t - L[j-1])/((L[j] - L[j-1]) || 1); out.push(cadd(pts[j-1], cscale(csub(pts[j], pts[j-1]), u))); }
      out.push(pts[pts.length-1]); return out;
    }
    function liftFrom(z0, w) {                 // path w with w[0] = f(z0): continuous preimage starting at z0
      var out = [z0], prev = z0;
      for (var i = 1; i < w.length; i++) { var r = csq(csub(w[i], c)); if (cabs(csub(r, prev)) > cabs(cadd(r, prev))) r = [-r[0], -r[1]]; out.push(r); prev = r; }
      return out;
    }
    var memo = {};
    function arc(i, j, d) {                    // regulated arc from P^i to P^j (points), indices mod n
      i %= n; j %= n; var key = i + ',' + j + ',' + d; if (memo[key]) return memo[key];
      var res;
      if (i === j) res = [P[i]];
      else if (d === 0) { res = []; for (var t = 0; t <= 16; t++) res.push(cadd(P[i], cscale(csub(P[j], P[i]), t/16))); }
      else if (i === 0) res = liftFrom(P[j], arc(j + 1, 1, d - 1)).reverse();           // [0, P^j]
      else if (j === 0) res = liftFrom(P[i], arc(i + 1, 1, d - 1));                      // [P^i, 0]
      else if (side(i) === side(j)) res = liftFrom(P[i], arc(i + 1, j + 1, d - 1));
      else res = arc(i, 0, d).concat(arc(0, j, d).slice(1));
      res = resample(res, Math.min(800, Math.max(300, res.length)));
      res[0] = P[i]; res[res.length - 1] = P[j];                                           // exact end points
      memo[key] = res; return res;
    }
    var edges = []; for (var j = 1; j < n; j++) edges.push(arc(0, j, depth));
    return { P: P, edges: edges };
  }
  // Moebius map sending (z1, z2, z3) to (w1, w2, w3); points may be null (infinity). Returns z -> w (null = infinity).
  function moebius3(z, w) {
    function hv(p) { return p === null ? [[1, 0], [0, 0]] : [p, [1, 0]]; }
    function toStd(p) {                     // matrix sending p1 -> 0, p2 -> 1, p3 -> inf
      var a = hv(p[0]), b = hv(p[1]), cc = hv(p[2]);
      // row vectors: M = [[ (b x c) ... ]] built from the classical formula with homogeneous points
      var det = function (u, v) { return csub(cmul(u[0], v[1]), cmul(u[1], v[0])); };
      var k1 = det(b, cc), k2 = det(b, a);
      // z -> (z - p1)(p2 - p3)/((z - p3)(p2 - p1)) in homogeneous form: rows [a1', -a0'] scaled
      return [[cmul([a[1][0], a[1][1]], k1), cscale(cmul(a[0], k1), -1)], [cmul([cc[1][0], cc[1][1]], k2), cscale(cmul(cc[0], k2), -1)]];
    }
    var A = toStd(z), B = toStd(w);
    var dB = csub(cmul(B[0][0], B[1][1]), cmul(B[0][1], B[1][0]));
    var Bi = [[cdiv(B[1][1], dB), cdiv(cscale(B[0][1], -1), dB)], [cdiv(cscale(B[1][0], -1), dB), cdiv(B[0][0], dB)]];
    var Mx = [[cadd(cmul(Bi[0][0], A[0][0]), cmul(Bi[0][1], A[1][0])), cadd(cmul(Bi[0][0], A[0][1]), cmul(Bi[0][1], A[1][1]))],
              [cadd(cmul(Bi[1][0], A[0][0]), cmul(Bi[1][1], A[1][0])), cadd(cmul(Bi[1][0], A[0][1]), cmul(Bi[1][1], A[1][1]))]];
    return function (q) {
      var num = q === null ? Mx[0][0] : cadd(cmul(Mx[0][0], q), Mx[0][1]), den = q === null ? Mx[1][0] : cadd(cmul(Mx[1][0], q), Mx[1][1]);
      if (cabs(den) <= 1e-300*cabs(num) || cabs(den) === 0) return null; var r = cdiv(num, den); return cabs(r) > 1e150 ? null : r;
    };
  }
  // ---- maps
  function corMap(x, y, n) {
    var w = (n === 3 && y[0] === 0 && y[1] === 0) ? [1, 0] : cdiv([x[0]+y[0]-1, x[1]+y[1]], [x[0]-1, x[1]]);   // n = 3: w = 1 identically (0/0 at x = 1, the bitransitive center c = 1)
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
    else if (family === 'per2c') { var cc = q.p; M = { n: 2, A: [1, 0], B: [0, 0], C: [-cc[0], -cc[1]], D: [1, 0], E: [0, 0], F: [-1, 0], omega: [0, 0], marked: [null, [1, 0]] }; }
    else if (family === 'cor3c') { M = corMap(cdiv(q.p, csub([2, 0], q.p)), [0, 0], 3); }
    else if (family === 'cor4c') { var y4 = cadd(cdiv([0.1, 0], q.p), [(Math.sqrt(5) - 1)/2, 0]); M = corMap(y4, y4, 4); }
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
      for (var j = 0; j < n; j++) if (chord(z, M.marked[j]) < th[j]) return { kind: 'capture', j: j, t: it+1, phase: ((j-1-(it+1)) % n + n) % n };   // phase convention (Caroline 2026-11-10): the orbit reaches the critical value a_1 at times = -phase mod n
      ring[it % 32] = z;
    }
    for (var per = 1; per <= 24; per++) {
      var zp = ring[(maxIt-1-per+3200) % 32], zl = ring[(maxIt-1) % 32];
      if (zp !== undefined && chord(zl, zp) < tolP) return { kind: 'cycle', period: per, t: maxIt };
    }
    return { kind: 'undecided' };
  }
  // the center of the capture component containing parameter p (plane / torus families): smallest t' <= t such that Newton in
  // the parameter on G(q) = f_q^t'(omega_q) - a_j'(q) (j' = j - (t - t'), phase-consistent) converges near p to a center whose
  // fate is a capture of the same phase. Returns {preperiod t', landing a_j', center} or null.
  function captureCentre(family, n, p, f) {
    if (!f || f.kind !== 'capture') return null;
    for (var tp = 0; tp <= f.t; tp++) {
      var jp = ((f.j - (f.t - tp)) % n + n) % n, q = p.slice(), ok = false;
      var G = function (q) { var M = mapFrom(family, n, { p: q }); if (!M) return null; var z = M.omega;
        for (var i = 0; i < tp; i++) { z = fz(M, z); if (z === null) break; }
        var a = M.marked[jp]; if (z === null || a === null) return (z === null && a === null) ? [0, 0] : (a === null ? cdiv([1, 0], z || [1e300, 0]) : cdiv([1, 0], a)); return csub(z, a); };
      for (var it = 0; it < 60; it++) {
        var g0 = G(q); if (!g0) break; var h = 1e-7*(1 + cabs(q)), g1 = G(cadd(q, [h, 0])); if (!g1) break;
        var d = cdiv(cscale(g0, h), csub(g1, g0)); q = csub(q, d); if (!isFinite(q[0]) || cabs(csub(q, p)) > 2*(1 + cabs(p))) break;
        if (cabs(d) < 1e-13*(1 + cabs(q))) { ok = true; break; }
      }
      if (!ok) continue;
      var Mc = mapFrom(family, n, { p: q }); if (!Mc) continue;
      var fc = fate(Mc, Mc.omega, 400, 1e-6, 1e-9);
      if (fc.kind === 'capture' && fc.phase === f.phase && fc.t <= tp + 1) return { preperiod: tp, j: jp, center: q, phase: f.phase };
    }
    return null;
  }
  // ---- "View in Per_m" (REGISTRATION addendum 2026-11-10 (c)) ----
  function csolve(A, b) {                    // complex linear solve, Gaussian elimination with pivoting; A: k x k of [re, im]
    var k = b.length, i, j, r; A = A.map(function (row) { return row.slice(); }); b = b.slice();
    for (i = 0; i < k; i++) { var pv = i; for (r = i + 1; r < k; r++) if (cabs(A[r][i]) > cabs(A[pv][i])) pv = r;
      if (cabs(A[pv][i]) === 0) return null; var tr = A[i]; A[i] = A[pv]; A[pv] = tr; var tb = b[i]; b[i] = b[pv]; b[pv] = tb;
      for (r = i + 1; r < k; r++) { var f = cdiv(A[r][i], A[i][i]); for (j = i; j < k; j++) A[r][j] = csub(A[r][j], cmul(f, A[i][j])); b[r] = csub(b[r], cmul(f, b[i])); } }
    var x = new Array(k);
    for (i = k - 1; i >= 0; i--) { var s = b[i]; for (j = i + 1; j < k; j++) s = csub(s, cmul(A[i][j], x[j])); x[i] = cdiv(s, A[i][i]); }
    return x;
  }
  function iterN(M, z, m) { for (var i = 0; i < m; i++) { z = fz(M, z); if (z === null) return null; } return z; }
  // center of the component where omega is attracted to an m-cycle. sel = {p} (plane / torus families) or {x, y} (family 'cor', n >= 5).
  function cycleCentre(family, n, sel, m) {
    var two = family === 'cor', v0 = two ? [sel.x, sel.y] : [sel.p];
    var mk = function (v) { try { return two ? mapFrom('cor', n, { x: v[0], y: v[1] }) : mapFrom(family, n, { p: v[0] }); } catch (e) { return null; } };
    var M0 = mk(v0); if (!M0) return null;
    var z = M0.omega; for (var i = 0; i < 4000; i++) { z = fz(M0, z); if (z === null) z = fz(M0, null); }
    var cyc = [z]; for (i = 1; i < m; i++) cyc.push(fz(M0, cyc[i - 1]));
    var z0 = null; cyc.forEach(function (w) { if (w !== null && (z0 === null || cabs(w) < cabs(z0))) z0 = w; }); if (z0 === null) return null;
    var back = iterN(M0, z0, m); if (back === null || chord(back, z0) > 1e-6) return null;                 // not (yet) on an m-cycle
    var dF = function (M, w) { var h = 1e-6*(1 + cabs(w)), a = iterN(M, w, m), b = iterN(M, [w[0] + h, w[1]], m); if (!a || !b) return null; return cscale(csub(b, a), 1/h); };
    var rho0 = dF(M0, z0); if (!rho0 || !(cabs(rho0) < 1)) return null;
    var G, run = function (ts) { var u = v0.concat([z0]);
    var K = v0.length + 1;
    G = function (u, t) { var M = mk(u); if (!M) return null; var w = u[K - 1], a = iterN(M, w, m), d = dF(M, w); if (!a || !d) return null;
      var out = [csub(a, w), csub(d, cscale(rho0, t))];
      if (two) { var e = iterN(M, u[1], n - 4); if (!e) return null; out.unshift(csub(e, u[0])); }
      return out.every(function (q) { return isFinite(q[0]) && isFinite(q[1]); }) ? out : null; };
    for (var ti = 0; ti < ts.length; ti++) { var okk = false;
      for (var it = 0; it < 40; it++) { var g0 = G(u, ts[ti]); if (!g0) return null; var J = g0.map(function () { return []; });
        for (var c = 0; c < K; c++) { var h = 1e-7*(1 + cabs(u[c])), u2 = u.slice(); u2[c] = [u[c][0] + h, u[c][1]]; var g1 = G(u2, ts[ti]); if (!g1) return null;
          for (var r = 0; r < K; r++) J[r][c] = cscale(csub(g1[r], g0[r]), 1/h); }
        var st = csolve(J, g0); if (!st) return null; var big = 0;
        for (c = 0; c < K; c++) { u[c] = csub(u[c], st[c]); big = Math.max(big, cabs(st[c])/(1 + cabs(u[c]))); }
        if (!u.every(function (q) { return isFinite(q[0]) && isFinite(q[1]); })) return null;
        if (big < 1e-8) { okk = true; break; } }                                  // finite-difference (f^m)' limits this stage; exact polish below
      if (!okk) { if (typeof process !== 'undefined' && process.env.CCDBG) console.log('fail at t', ts[ti], 'u', u, 'G', G(u, ts[ti])); return null; } }
    var v = u.slice(0, K - 1), KV = K - 1, pol = false;                      // polish: f^m(omega) = omega exactly, in the parameter only
    var H = function (v) { var M = mk(v); if (!M) return null; var a = iterN(M, M.omega, m); if (!a) return null; var out = [csub(a, M.omega)];
      if (two) { var e = iterN(M, v[1], n - 4); if (!e) return null; out.unshift(csub(e, v[0])); } return out; };
    for (var it2 = 0; it2 < 40; it2++) { var h0 = H(v); if (!h0) return null; var J2 = h0.map(function () { return []; });
      for (var c2 = 0; c2 < KV; c2++) { var hh = 1e-7*(1 + cabs(v[c2])), v2 = v.slice(); v2[c2] = [v[c2][0] + hh, v[c2][1]]; var h1 = H(v2); if (!h1) return null;
        for (var r2 = 0; r2 < KV; r2++) J2[r2][c2] = cscale(csub(h1[r2], h0[r2]), 1/hh); }
      var s2 = csolve(J2, h0); if (!s2) return null; var b2 = 0; for (c2 = 0; c2 < KV; c2++) { v[c2] = csub(v[c2], s2[c2]); b2 = Math.max(b2, cabs(s2[c2])/(1 + cabs(v[c2]))); }
      if (b2 < 1e-13) { pol = true; break; } }
    if (!pol) return null; u = v;
    var Mc = mk(u); if (!Mc) return null; var bk = iterN(Mc, Mc.omega, m); if (bk === null || chord(bk, Mc.omega) > 1e-8) return null;
    for (i = 1; i < m; i++) { var wj = iterN(Mc, Mc.omega, i); if (wj !== null && chord(wj, Mc.omega) < 1e-8) return null; }   // exact period m
    for (i = 0; i < Mc.marked.length; i++) for (var j2 = i + 1; j2 < Mc.marked.length; j2++) if (chord(Mc.marked[i], Mc.marked[j2]) < 1e-7) return null;   // marked cycle of exact period n
    return two ? { x: u[0], y: u[1], M: Mc, rho0: rho0 } : { p: u[0], M: Mc, rho0: rho0 };
    };
    var coarse = [0.85, 0.7, 0.55, 0.4, 0.25, 0.12, 0.05, 0], fine = []; for (i = 59; i >= 0; i--) fine.push(i/60);
    return run(coarse) || run(fine);
  }
  // the map M (omega periodic of exact period m) written in Per_m(0): marked cycle = orbit of omega, free critical point = M.marked[0].
  // Returns {n: m, family, p} or {n: m, family: 'cor' | 'torus5', x, y} (as MATE.toChart), or {n: 1, family: 'quad', p: c}.
  function transfer(M, m) {
    var b = [M.omega], a0 = M.marked[0], R = function (z) { return fz(M, z); };
    for (var i = 1; i < m; i++) b.push(fz(M, b[i - 1]));
    if (m === 1) { var T = moebius3([b[0], a0, fz(M, a0)], [null, [0, 0], [1, 0]]);          // omega -> inf, a0 -> 0, f(a0) -> 1: g(z) = (z^2 + k)/(1 + k)... solved below
      // g = T f T^-1 is a quadratic polynomial A z^2 + B with g(0) = 1: B = 1, A = g(1) - 1; conjugate to z^2 + A B
      var Ti = moebius3([null, [0, 0], [1, 0]], [b[0], a0, fz(M, a0)]), g1 = T(fz(M, Ti([1, 0]))); if (g1 === null) return null;
      return { n: 1, family: 'quad', p: csub(g1, [1, 0]) }; }
    if (m === 2) { var T2 = moebius3([b[0], b[1], a0], [null, [1, 0], [0, 0]]); return { n: 2, family: 'per2c', p: T2(R(a0)) }; }
    var T3 = moebius3([b[0], b[1], b[2]], [[0, 0], null, [1, 0]]);
    var y = m >= 4 ? T3(b[3]) : [0, 0], x = T3(b[m - 1]), w = T3(a0), out = { n: m, x: x, y: y, omega: w };
    if (m === 3) { out.family = 'cor3c'; out.p = w; }
    else if (m === 4) { out.family = 'cor4c'; out.p = cdiv([0.1, 0], csub(y, [(Math.sqrt(5) - 1)/2, 0])); }
    else out.family = m === 5 ? 'torus5' : 'cor';
    return out;
  }
  // ---- bitransitive coordinate Phi (REGISTRATION addendum 2026-11-11 (j)). M: map; j: omega is in the immediate basin component of a_j (1 <= j <= n-1).
  // Returns Phi = Boettcher coordinate (of f^n at a_0) of v = f^(n-j)(omega), or null.
  function bitransFactors(M, j, ku) {
    var n = M.marked.length, T = M.marked[0] === null || M.marked[1] !== null ? moebius3([M.marked[0], M.marked[1], n >= 3 ? M.marked[2] : M.omega], [[0, 0], null, [1, 0]]) : function (z) { return z; };
    var Ti = M.marked[0] === null || M.marked[1] !== null ? moebius3([[0, 0], null, [1, 0]], [M.marked[0], M.marked[1], n >= 3 ? M.marked[2] : M.omega]) : function (z) { return z; };
    var F = function (w) { var z = Ti(w); for (var i = 0; i < n; i++) z = fz(M, z); return T(z); };
    // c_2 (F(w) = c_2 w^2 + ...) by a 16-point Cauchy sum on a circle that F maps well inside itself (addendum (k): the earlier central difference with
    // h = 1e-4 carried ~1e-8 relative rounding noise, which stopped the Newton continuation along internal rays; FAILURE kept in log.md)
    var c2 = null, h = 0.1, NC = 16, CM = null, m;
    for (var tr = 0; tr < 12 && !c2; tr++, h /= 2) { var S = [], ok = true; for (m = 0; m < 10; m++) S.push([0, 0]);
      for (var q = 0; q < NC; q++) { var v = F([h*Math.cos(2*Math.PI*q/NC), h*Math.sin(2*Math.PI*q/NC)]); if (v === null || !(cabs(v) < 0.5*h)) { ok = false; break; }
        for (m = 0; m < 10; m++) { var cm = Math.cos(2*Math.PI*q*m/NC), sm = Math.sin(2*Math.PI*q*m/NC); S[m][0] += v[0]*cm + v[1]*sm; S[m][1] += v[1]*cm - v[0]*sm; } }
      if (ok) { CM = S.map(function (x, mm) { return cscale(x, 1/(NC*Math.pow(h, mm))); }); c2 = CM[2]; } }
    if (!c2) return null;
    var hT = h*2/40, RT = CM.slice(2).map(function (x) { return cdiv(x, c2); });      // h was halved once more by the loop header; Taylor evaluation below |w| < hT
    var rTay = function (w) { var acc = RT[7]; for (var mm = 6; mm >= 0; mm--) acc = cadd(cmul(acc, w), RT[mm]); return acc; };     // F(w)/(c_2 w^2) = sum_(m>=2) (c_m/c_2) w^(m-2)
    if (!(cabs(c2) > 0) || !isFinite(c2[0] + c2[1])) return null;
    var z = M.omega; for (var i = 0; i < n - j; i++) z = fz(M, z); var w = T(z); if (w === null) return null;
    // phi(w) = u_0 prod_k r_k^(1/2^(k+1)), r_k = u_(k+1)/u_k^2. Returned raw: u0, the first KU factors r_k (their arguments must be chosen by
    // continuity in the PARAMETER, starting from 0 at the center: BTRI does that), and the principal product of the remaining factors.
    // FAILURES kept in log.md: principal arguments throughout (wrong when r_k has turned past the negative axis); unwrapping along the segment s w in the
    // dynamical plane (wrong for some parameters); iterating below |u| ~ 1e-9 (rounding).
    var KU = ku || 4, u = cmul(c2, w), u0 = u, rs = [], tail = [1, 0], e = 0.5, k;
    // addendum (k): near 0, F is evaluated by its Taylor polynomial (relative precision kept), so the product can run to |u| < 1e-13 without rounding noise
    for (k = 0; k < 70; k++) { if (!isFinite(u[0] + u[1]) || cabs(u) > 1e6) return null; if (cabs(u) < 1e-13) break;
      var r, u1;
      if (cabs(w) < hT) { r = rTay(w); u1 = cmul(r, cmul(u, u)); w = cdiv(u1, c2); }
      else { w = F(w); if (w === null) return null; u1 = cmul(c2, w); r = cdiv(u1, cmul(u, u)); }
      if (!(cabs(r) > 0) || !isFinite(r[0] + r[1])) return null;
      if (k < KU) rs.push(r); else { var mg = Math.pow(cabs(r), e), tg = Math.atan2(r[1], r[0])*e; tail = cmul(tail, [mg*Math.cos(tg), mg*Math.sin(tg)]); }
      u = u1; e /= 2; }
    if (cabs(u) >= 1e-3) return null; while (rs.length < KU) rs.push([1, 0]);
    return { u0: u0, r: rs, tail: tail };
  }
  // Phi from the raw factors and chosen arguments A[k] of r_k
  function bitransPhiFrom(Fc, A) { var phi = cmul(Fc.u0, Fc.tail), e = 0.5; for (var k = 0; k < Fc.r.length; k++) { var m = Math.pow(cabs(Fc.r[k]), e), t = A[k]*e; phi = cmul(phi, [m*Math.cos(t), m*Math.sin(t)]); e /= 2; } return phi; }
  // (addendum (k)) (T(F^(a+q)(v)) - T(F^a(v))) / T(F^a(v)), v = f^(n-j)(omega), F = f^n, in the coordinate a_0 = 0, a_1 = inf; null if an orbit point is at infinity there
  function bitransDiff(M, j, a, q) {
    var n = M.marked.length, id = !(M.marked[0] === null || M.marked[1] !== null), T = id ? function (z) { return z; } : moebius3([M.marked[0], M.marked[1], n >= 3 ? M.marked[2] : M.omega], [[0, 0], null, [1, 0]]);
    var z = M.omega, i; for (i = 0; i < n - j; i++) z = fz(M, z); for (i = 0; i < n*a; i++) z = fz(M, z); var z1 = z; for (i = 0; i < n*q; i++) z1 = fz(M, z1);
    var w0 = T(z), w1 = T(z1); if (w0 === null || w1 === null || !isFinite(w0[0] + w0[1] + w1[0] + w1[1]) || !(cabs(w0) > 0)) return null; return cdiv(csub(w1, w0), w0);   // divided by w0: removes the solution F^a(v) = a_0 (the center itself for a = 0)
  }
  // ---- puncture cards (REGISTRATION addendum 2026-11-12 (t)): at a map M NEAR a puncture: which marked points have (nearly) collided, where the free
  // critical point and its value sit, and the fixed points with their multipliers. Clusters: single linkage on chordal distances below thr, where thr is
  // placed in the largest gap (ratio) of the sorted pairwise distances, provided the gap ratio is >= 20 and the lower side is < 0.1; otherwise no cluster.
  function punctureInfo(M) {
    var n = M.marked.length, names = M.marked.map(function (_, i) { return 'a' + i; }), ds = [], i, j;
    for (i = 0; i < n; i++) for (j = i + 1; j < n; j++) ds.push([chord(M.marked[i], M.marked[j]), i, j]);
    ds.sort(function (u, v) { return u[0] - v[0]; }); var cut = -1, best = 20;
    for (i = 0; i + 1 < ds.length; i++) { if (!(ds[i][0] < 0.1)) break; var r = ds[i + 1][0]/Math.max(ds[i][0], 1e-300); if (r >= best) { best = r; cut = i; } }
    if (ds.length && cut < 0 && ds[ds.length - 1][0] < 0.1) cut = ds.length - 1;
    var par = names.map(function (_, k) { return k; }), find = function (k) { while (par[k] !== k) k = par[k]; return k; };
    for (i = 0; i <= cut; i++) par[find(ds[i][1])] = find(ds[i][2]);
    var G = {}; for (i = 0; i < n; i++) { var g = find(i); (G[g] = G[g] || []).push(i); }
    var clusters = Object.keys(G).map(function (g) { return G[g]; }).filter(function (c) { return c.length > 1; });
    var diam = clusters.map(function (c) { var d = 0; for (var u = 0; u < c.length; u++) for (var v = u + 1; v < c.length; v++) d = Math.max(d, chord(M.marked[c[u]], M.marked[c[v]])); return d; });
    var near = function (z) { var bd = Infinity, bi = -1; for (var k = 0; k < n; k++) { var d = chord(z, M.marked[k]); if (d < bd) { bd = d; bi = k; } } return { i: bi, d: bd }; };
    var out = { clusters: clusters, diam: diam, gap: cut >= 0 ? best : null, minOther: cut + 1 < ds.length ? ds[cut + 1][0] : null, omega: near(M.omega), value: near(fz(M, M.omega)), fixed: null };
    // fixed points: D z^3 + (E - A) z^2 + (F - B) z - C = 0 (plus infinity if D = 0); multipliers by a central difference
    if (cabs(M.D) > 0) { var c = [cscale(M.C, -1), csub(M.F, M.B), csub(M.E, M.A), M.D], z = [[0.4, 0.9], [-0.7, 0.3], [0.5, -0.8]], it, k;
      for (it = 0; it < 200; it++) { var mv = 0; for (k = 0; k < 3; k++) { var x = z[k], pv = cadd(cmul(cadd(cmul(cadd(cmul(c[3], x), c[2]), x), c[1]), x), c[0]), den = c[3];
          for (j = 0; j < 3; j++) if (j !== k) den = cmul(den, csub(x, z[j])); if (!(cabs(den) > 0)) continue; var st = cdiv(pv, den); z[k] = csub(x, st); mv = Math.max(mv, cabs(st)/(1 + cabs(z[k]))); } if (mv < 1e-15) break; }
      out.fixed = z.map(function (w) { var h = 1e-6*(1 + cabs(w)), a = fz(M, [w[0] + h, w[1]]), b = fz(M, [w[0] - h, w[1]]), fw = fz(M, w); return { z: w, mult: a && b ? cscale(csub(a, b), 1/(2*h)) : null, res: fw ? chord(fw, w) : null }; }); }
    return out;
  }
  function orbit(M, z, N) { var out = [z]; for (var i = 0; i < N; i++) { z = fz(M, z); out.push(z); } return out; }
  function closure(M) { var a = M.marked[0], z = a; for (var i = 0; i < M.n; i++) z = fz(M, z); return chord(z, a); }
  function fmt(z, d) { d = d || 6; if (z === null) return '∞'; var f = function (v) { return (Math.abs(v) < 1e-12 ? 0 : v).toPrecision(d); };
    return f(z[0]) + (z[1] < 0 ? ' − ' : ' + ') + f(Math.abs(z[1])) + 'i'; }
  function fateText(f) { return f.kind === 'capture' ? 'captured by the marked cycle, phase ' + f.phase + ', time ' + f.t
    : f.kind === 'cycle' ? 'attracting cycle, period ' + f.period : 'undecided (bifurcation locus)'; }
  return { cmul: cmul, cdiv: cdiv, cadd: cadd, csub: csub, cabs: cabs, chord: chord, T5: T5, wp: wp, torusXY: torusXY, torusLocate: torusLocate,
           corMap: corMap, mapFrom: mapFrom, captureCentre: captureCentre, bitransFactors: bitransFactors, bitransPhiFrom: bitransPhiFrom, bitransDiff: bitransDiff, punctureInfo: punctureInfo, cycleCentre: cycleCentre, transfer: transfer, iterN: iterN, parseAngle: parseAngle, traceRay: traceRay, rayWithLanding: rayWithLanding, hubbardTree: hubbardTree, moebius3: moebius3, anglePeriod: anglePeriod, fz: fz, fate: fate, orbit: orbit, closure: closure, fmt: fmt, fateText: fateText };
})();
if (typeof module !== 'undefined') module.exports = PM;
