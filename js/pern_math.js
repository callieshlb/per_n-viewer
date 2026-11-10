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

  // ---- Hubbard tree of a centre z^2 + c (period n): union of the regulated arcs [0, P^j(0)], built recursively from the
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
  // the centre of the capture component containing parameter p (plane / torus families): smallest t' <= t such that Newton in
  // the parameter on G(q) = f_q^t'(omega_q) - a_j'(q) (j' = j - (t - t'), phase-consistent) converges near p to a centre whose
  // fate is a capture of the same phase. Returns {preperiod t', landing a_j', centre} or null.
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
      if (fc.kind === 'capture' && fc.phase === f.phase && fc.t <= tp + 1) return { preperiod: tp, j: jp, centre: q, phase: f.phase };
    }
    return null;
  }
  function orbit(M, z, N) { var out = [z]; for (var i = 0; i < N; i++) { z = fz(M, z); out.push(z); } return out; }
  function closure(M) { var a = M.marked[0], z = a; for (var i = 0; i < M.n; i++) z = fz(M, z); return chord(z, a); }
  function fmt(z, d) { d = d || 6; if (z === null) return '∞'; var f = function (v) { return (Math.abs(v) < 1e-12 ? 0 : v).toPrecision(d); };
    return f(z[0]) + (z[1] < 0 ? ' − ' : ' + ') + f(Math.abs(z[1])) + 'i'; }
  function fateText(f) { return f.kind === 'capture' ? 'captured by the marked cycle, phase ' + f.phase + ', time ' + f.t
    : f.kind === 'cycle' ? 'attracting cycle, period ' + f.period : 'undecided (bifurcation locus)'; }
  return { cmul: cmul, cdiv: cdiv, cadd: cadd, csub: csub, cabs: cabs, chord: chord, T5: T5, wp: wp, torusXY: torusXY, torusLocate: torusLocate,
           corMap: corMap, mapFrom: mapFrom, captureCentre: captureCentre, parseAngle: parseAngle, traceRay: traceRay, rayWithLanding: rayWithLanding, hubbardTree: hubbardTree, moebius3: moebius3, anglePeriod: anglePeriod, fz: fz, fate: fate, orbit: orbit, closure: closure, fmt: fmt, fateText: fateText };
})();
if (typeof module !== 'undefined') module.exports = PM;
