// E96: the three triangles of a bitransitive component (REGISTRATION addendum 2026-11-11 (j)). Depends on pern_math.js (PM).
'use strict';
var BTRI = (function () {
  // mapOf(p) -> M or null; b: center (parameter point); j: slot. opts: N (grid), half (box half-size; found by doubling if absent), offs: [dx, dy] unused
  function fillBox(mapOf, b, j, half, N) {
    var W = new Float64Array(2*N*N).fill(NaN), AR = new Float64Array(4*N*N), st = new Uint8Array(N*N), dx = 2*half/(N - 1), x0 = b[0] - half, y0 = b[1] - half, ci = (N - 1) >> 1;
    var facAt = function (ix, iy) { var M = mapOf([x0 + ix*dx, y0 + iy*dx]); if (!M) return null; return PM.bitransFactors(M, j); };
    var roots = function (P) { var m = Math.cbrt(Math.hypot(P[0], P[1])), t0 = Math.atan2(P[1], P[0]), t = (t0 < -Math.PI + 1e-6 ? t0 + 2*Math.PI : t0)/3, o = [];   /* arg Phi = -pi (real-symmetric centers) is read as +pi: deterministic labels */ for (var k = 0; k < 3; k++) o.push([m*Math.cos(t + 2*Math.PI*k/3), m*Math.sin(t + 2*Math.PI*k/3)]); return o; };
    var idx = function (ix, iy) { return iy*N + ix; }, q = [], count = 0, rejDisc = 0, rejOut = 0, sup = 0, edge = false, DX = [1, -1, 0, 0, 1, 1, -1, -1], DY = [0, 0, 1, -1, 1, -1, 1, -1];
    var c0 = idx(ci, ci); st[c0] = 1; W[2*c0] = 0; W[2*c0 + 1] = 0; count++;            // the center: all arguments 0
    var assign = function (ix, iy, par) {                               // par: index of an assigned neighbor (arguments continued from it); returns true if accepted
      var id = idx(ix, iy), Fc = facAt(ix, iy); if (!Fc) { st[id] = 2; rejOut++; return false; }
      var A = [0, 0, 0, 0]; for (var k = 0; k < 4; k++) { var a = Math.atan2(Fc.r[k][1], Fc.r[k][0]), ref = AR[4*par + k]; A[k] = a + 2*Math.PI*Math.round((ref - a)/(2*Math.PI)); }
      var P = PM.bitransPhiFrom(Fc, A); if (!isFinite(P[0] + P[1]) || Math.hypot(P[0], P[1]) >= 1) { st[id] = 2; rejOut++; return false; }
      var R = roots(P), best = R[0], px = (ix - ci)*dx, py = (iy - ci)*dx, pred = kappa ? [kappa[0]*px - kappa[1]*py, kappa[0]*py + kappa[1]*px] : null;
      if (pred && Math.hypot(pred[0], pred[1]) < 0.3) { var bp = Infinity; for (k = 0; k < 3; k++) { var dq = Math.hypot(R[k][0] - pred[0], R[k][1] - pred[1]); if (dq < bp) { bp = dq; best = R[k]; } } }   // near the center W ~ kappa (p - b): the nearest-to-parent rule is ill-conditioned there (FAILURE kept in log.md)
      else if (par !== c0) { var rx = W[2*par], ry = W[2*par + 1], bd = Infinity; for (k = 0; k < 3; k++) { var dd = Math.hypot(R[k][0] - rx, R[k][1] - ry); if (dd < bd) { bd = dd; best = R[k]; } }
        var mr = Math.hypot(rx, ry), mb = Math.hypot(best[0], best[1]); if (bd > 0.3 || (mr > 0.8 && mb < mr - 0.07)) { st[id] = 3; rejDisc++; return false; } }
      st[id] = 1; W[2*id] = best[0]; W[2*id + 1] = best[1]; for (k = 0; k < 4; k++) AR[4*id + k] = A[k]; count++; sup = Math.max(sup, Math.hypot(best[0], best[1])); q.push(ix, iy); return true; };
    var kappa = null, F0 = facAt(ci + 1, ci); if (F0) { var P0 = PM.bitransPhiFrom(F0, F0.r.map(function (r) { return Math.atan2(r[1], r[0]); })), k0 = roots(P0)[0]; kappa = [k0[0]/dx, k0[1]/dx]; }
    if (!assign(ci + 1, ci, c0)) return { W: W, N: N, x0: x0, y0: y0, dx: dx, count: count, fail: 'no value next to the center' };
    var head = 0;
    while (head < q.length) { var ix = q[head++], iy = q[head++], me = idx(ix, iy);
      for (var d = 0; d < 8; d++) { var nx = ix + DX[d], ny = iy + DY[d]; if (nx < 0 || ny < 0 || nx >= N || ny >= N) { edge = true; continue; } if (st[idx(nx, ny)]) continue; assign(nx, ny, me); } }
    return { W: W, st: st, N: N, x0: x0, y0: y0, dx: dx, kappa: kappa, count: count, rejDisc: rejDisc, rejOut: rejOut, sup: sup, edge: edge, half: half };
  }
  function fill(mapOf, b, j, opts) {
    opts = opts || {}; var half = opts.half || 1e-3*(1 + Math.hypot(b[0], b[1])), R = null;
    for (var it = 0; it < 16; it++) { R = fillBox(mapOf, b, j, half, 41); if (R.fail) return R; if (!R.edge) break; half *= 2; }
    return fillBox(mapOf, b, j, half*1.15, opts.N || 241);
  }
  // triangle index (0, 1, 2) and internal angle in [0, 1) of a W value
  function tri(wx, wy) { var a = Math.atan2(wy, wx); if (a < 0) a += 2*Math.PI; var s = a/(2*Math.PI/3), k = Math.min(2, Math.floor(s)); return [k, s - k]; }
  // ---- internal ray (REGISTRATION addendum 2026-11-11 (k)). k: triangle (0, 1, 2); t: internal angle in [0, 1). opts.kappa, opts.dx: from fill (same labels as the overlay).
  function ray(mapOf, b, j, k, t, opts) {
    opts = opts || {}; var KU = 48, TWO = 2*Math.PI, th = TWO*(k + t)/3, dir = [Math.cos(th), Math.sin(th)];
    var fac = function (p) { var M = mapOf(p); return M ? PM.bitransFactors(M, j, KU) : null; };
    var roots = function (P) { var m = Math.cbrt(Math.hypot(P[0], P[1])), a0 = Math.atan2(P[1], P[0]), a = (a0 < -Math.PI + 1e-6 ? a0 + TWO : a0)/3, o = []; for (var i = 0; i < 3; i++) o.push([m*Math.cos(a + TWO*i/3), m*Math.sin(a + TWO*i/3)]); return o; };
    // W at p: arguments continued from Aref (null: principal), cube root nearest to Wpred
    var ev = function (p, Aref, Wpred) { var Fc = fac(p); if (!Fc) return null; var A = [], jump = 0;
      for (var i = 0; i < KU; i++) { var a = Math.atan2(Fc.r[i][1], Fc.r[i][0]); if (Aref) { var a2 = a + TWO*Math.round((Aref[i] - a)/TWO); jump = Math.max(jump, Math.abs(a2 - Aref[i])); a = a2; } A.push(a); }
      var P = PM.bitransPhiFrom(Fc, A); if (!isFinite(P[0] + P[1])) return null; var R = roots(P), best = R[0], bd = Infinity;
      for (i = 0; i < 3; i++) { var d = Math.hypot(R[i][0] - Wpred[0], R[i][1] - Wpred[1]); if (d < bd) { bd = d; best = R[i]; } }
      return { W: best, A: A, jump: jump, off: bd }; };
    // reference point b + dx (principal arguments, principal cube root: the fill's convention) and kappa
    var dx = opts.dx, kappa = opts.kappa, e0;
    if (!dx) { dx = 1e-3*(1 + Math.hypot(b[0], b[1])); for (var it = 0; it < 30; it++) { var F0 = fac([b[0] + dx, b[1]]); if (F0) { var P0 = PM.bitransPhiFrom(F0, F0.r.map(function (r) { return Math.atan2(r[1], r[0]); })); if (Math.hypot(P0[0], P0[1]) < 1e-6) break; } dx /= 4; } kappa = null; }
    var Fr = fac([b[0] + dx, b[1]]); if (!Fr) return { fail: 'no value next to the center' };
    var Ar = Fr.r.map(function (r) { return Math.atan2(r[1], r[0]); }), Pr = PM.bitransPhiFrom(Fr, Ar), Wr = roots(Pr)[0]; if (!kappa) kappa = [Wr[0]/dx, Wr[1]/dx];
    var kap2 = kappa[0]*kappa[0] + kappa[1]*kappa[1], lin = function (W) { return [(W[0]*kappa[0] + W[1]*kappa[1])/kap2, (W[1]*kappa[0] - W[0]*kappa[1])/kap2]; };   // W / kappa
    // walk the circle |p - b| = dx from direction 0 to the ray's direction
    var d0 = lin(dir), al = Math.atan2(d0[1], d0[0]), ns = Math.max(2, Math.ceil(Math.abs(al)/(TWO/180))), A = Ar, W = Wr, p = [b[0] + dx, b[1]];
    for (var s = 1; s <= ns; s++) { var be = al*s/ns; p = [b[0] + dx*Math.cos(be), b[1] + dx*Math.sin(be)]; var pred = [kappa[0]*(p[0] - b[0]) - kappa[1]*(p[1] - b[1]), kappa[0]*(p[1] - b[1]) + kappa[1]*(p[0] - b[0])];
      e0 = ev(p, A, pred); if (!e0) return { fail: 'no value on the small circle' }; A = e0.A; W = e0.W; }
    // march in rho
    var rho = Math.hypot(W[0], W[1]), pts = [], hs = 1e-6/Math.sqrt(kap2), step = 0.02, nev = 0, why = 'reached 1 - rho < 1e-8';
    var newton = function (p0, A0, Wt) { var q = p0.slice(), cur = ev(q, A0, Wt); nev++; if (!cur) return null;
      for (var i = 0; i < 12; i++) { var r = [cur.W[0] - Wt[0], cur.W[1] - Wt[1]], rm = Math.hypot(r[0], r[1]); if (rm < 1e-11) return { p: q, A: cur.A, W: cur.W, res: rm };
        var h = hs*Math.max(1e-3, Math.min(1, 1 - Math.hypot(Wt[0], Wt[1]))), e1 = ev([q[0] + h, q[1]], cur.A, cur.W); nev++; if (!e1) return null;
        var D = [(e1.W[0] - cur.W[0])/h, (e1.W[1] - cur.W[1])/h], dd = D[0]*D[0] + D[1]*D[1]; if (!(dd > 0)) return null;
        var dp = [(r[0]*D[0] + r[1]*D[1])/dd, (r[1]*D[0] - r[0]*D[1])/dd]; q = [q[0] - dp[0], q[1] - dp[1]];
        var nx = ev(q, cur.A, Wt); nev++; if (!nx || nx.jump > 1.5) return null; if (i >= 3 && Math.hypot(nx.W[0] - Wt[0], nx.W[1] - Wt[1]) > 0.5*rm && rm > 1e-9) return null; cur = nx; }
      r = Math.hypot(cur.W[0] - Wt[0], cur.W[1] - Wt[1]); return r < 1e-8 ? { p: q, A: cur.A, W: cur.W, res: r } : null; };
    // first: correct the starting point onto the ray at its own rho
    var st = newton(p, A, [rho*dir[0], rho*dir[1]]); if (!st) return { fail: 'could not get onto the ray near the center' }; p = st.p; A = st.A; pts.push([rho, p[0], p[1]]);
    for (var guard = 0; guard < 4000; guard++) { if (1 - rho < 1e-8) break; var dr = Math.min(step, (1 - rho)/4), r1 = rho + dr, sN = newton(p, A, [r1*dir[0], r1*dir[1]]);
      if (!sN || Math.hypot(sN.p[0] - p[0], sN.p[1] - p[1]) > 0.5/Math.sqrt(kap2)) { step = dr/2; if (step < 1e-13) { why = 'continuation stopped (step < 1e-13)'; break; } continue; }
      p = sN.p; A = sN.A; rho = r1; pts.push([rho, p[0], p[1]]); step = Math.min(0.02, step*1.5); }
    if (guard >= 4000) why = 'step limit';
    var drift = null; for (var i = pts.length - 1; i >= 0; i--) if (1 - pts[i][0] >= 10*(1 - rho)) { drift = Math.hypot(pts[i][1] - p[0], pts[i][2] - p[1]); break; }
    return { pts: pts, end: p, rho: rho, drift: drift, why: why, nev: nev, kappa: kappa, dx: dx };
  }
  // preperiod a and period q of the angle {p, q} (BigInt) under doubling
  function prePer(ang) { var seen = new Map(), x = ang.p % ang.q, i = 0; while (!seen.has(x)) { seen.set(x, i++); x = (2n*x) % ang.q; if (i > 4000) return null; } var a = seen.get(x); return [a, i - a]; }
  // Newton on F^(a+q)(v) = F^a(v) from p0
  function polish(mapOf, j, a, q, p0) {
    var g = function (p) { var M = mapOf(p); return M ? PM.bitransDiff(M, j, a, q) : null; }, p = p0.slice(), sc = 1e-7*(1 + Math.hypot(p0[0], p0[1])), r = null;
    for (var i = 0; i < 40; i++) { var g0 = g(p); if (!g0) return null; r = Math.hypot(g0[0], g0[1]); if (r < 1e-13) break; var g1 = g([p[0] + sc, p[1]]); if (!g1) return null;
      var D = [(g1[0] - g0[0])/sc, (g1[1] - g0[1])/sc], dd = D[0]*D[0] + D[1]*D[1]; if (!(dd > 0)) return null; var dp = [(g0[0]*D[0] + g0[1]*D[1])/dd, (g0[1]*D[0] - g0[0]*D[1])/dd];
      p = [p[0] - dp[0], p[1] - dp[1]]; if (!isFinite(p[0] + p[1])) return null; if (Math.hypot(dp[0], dp[1]) < 1e-15*(1 + Math.hypot(p[0], p[1]))) { var gl = g(p); r = gl ? Math.hypot(gl[0], gl[1]) : NaN; break; } }
    return { p: p, res: r, it: i };
  }
  // periodic internal angle (period q): the landing point is NOT a critical-orbit relation (a critical point cannot be periodic and repelling); the equation
  // registered instead (addendum (k2)) is: a point z of period n q under f with multiplier 1. Unknowns (z, p); start z = f^(n-j)(omega) at the ray end.
  function polishPar(mapOf, j, q, p0) {
    var M0 = mapOf(p0); if (!M0) return null; var n = M0.marked.length, z = M0.omega, i; for (i = 0; i < n - j; i++) z = PM.fz(M0, z); if (!z) return null;
    var G = function (M, w) { for (var i = 0; i < n*q; i++) { w = PM.fz(M, w); if (!w || !isFinite(w[0] + w[1])) return null; } return w; };
    var H = function (z, p) { var M = mapOf(p); if (!M) return null; var h = 1e-3*Math.hypot(z[0], z[1]), g0 = G(M, z), a = G(M, [z[0] + h, z[1]]), b = G(M, [z[0] - h, z[1]]), c = G(M, [z[0] + 2*h, z[1]]), d = G(M, [z[0] - 2*h, z[1]]);
      if (!g0 || !a || !b || !c || !d) return null; var D = [(8*(a[0] - b[0]) - (c[0] - d[0]))/(12*h), (8*(a[1] - b[1]) - (c[1] - d[1]))/(12*h)]; return [g0[0] - z[0], g0[1] - z[1], D[0] - 1, D[1] - 1]; };
    var cm = function (a, b) { return [a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0]]; }, cdv = function (a, b) { var d = b[0]*b[0] + b[1]*b[1]; return [(a[0]*b[0] + a[1]*b[1])/d, (a[1]*b[0] - a[0]*b[1])/d]; };
    var p = p0.slice(), res = null, sz = 1e-6*Math.hypot(z[0], z[1]), sp = 1e-7*(1 + Math.hypot(p[0], p[1]));
    for (var it = 0; it < 60; it++) { var h0 = H(z, p); if (!h0) return null; res = Math.hypot(h0[0], h0[1])/Math.hypot(z[0], z[1]) + Math.hypot(h0[2], h0[3]); if (res < 1e-11) break;
      var hz = H([z[0] + sz, z[1]], p), hp = H(z, [p[0] + sp, p[1]]); if (!hz || !hp) return null;
      var a11 = [(hz[0] - h0[0])/sz, (hz[1] - h0[1])/sz], a21 = [(hz[2] - h0[2])/sz, (hz[3] - h0[3])/sz], a12 = [(hp[0] - h0[0])/sp, (hp[1] - h0[1])/sp], a22 = [(hp[2] - h0[2])/sp, (hp[3] - h0[3])/sp];
      var det = [cm(a11, a22)[0] - cm(a12, a21)[0], cm(a11, a22)[1] - cm(a12, a21)[1]]; if (!(Math.hypot(det[0], det[1]) > 0)) return null;
      var r1 = [h0[0], h0[1]], r2 = [h0[2], h0[3]], dz = cdv([cm(r1, a22)[0] - cm(a12, r2)[0], cm(r1, a22)[1] - cm(a12, r2)[1]], det), dp = cdv([cm(a11, r2)[0] - cm(a21, r1)[0], cm(a11, r2)[1] - cm(a21, r1)[1]], det);
      z = [z[0] - dz[0], z[1] - dz[1]]; p = [p[0] - dp[0], p[1] - dp[1]]; if (!isFinite(z[0] + z[1] + p[0] + p[1])) return null; }
    var hf = H(z, p); return hf ? { p: p, z: z, res: res, it: it, mult: [hf[2] + 1, hf[3]] } : null;
  }
  // ray + landing point. ang: PM.parseAngle result. Exact landing point only for strictly preperiodic angles (Misiurewicz relation, checked: residual, repelling,
  // close to the ray end). Periodic angles: the ray end only (polishPar is NOT validated and not used: log.md 2026-11-11 (d)).
  function boundaryPoint(mapOf, b, j, k, ang, opts) {
    var t = Number(ang.p)/Number(ang.q), r = ray(mapOf, b, j, k, t, opts); if (r.fail) return r; var pp = prePer(ang); r.t = t; r.k = k; r.pp = pp; r.exact = false; r.point = r.end;
    if (pp && pp[0] > 0) { var P = polish(mapOf, j, pp[0], pp[1], r.end);
      if (P && P.res < 1e-10 && Math.hypot(P.p[0] - r.end[0], P.p[1] - r.end[1]) <= 100*(r.drift || 0) + 1e-9) { var M = mapOf(P.p), n = M.marked.length, z = M.omega, i; for (i = 0; i < n - j + n*pp[0]; i++) z = PM.fz(M, z);
        var Fq = function (w) { for (var i = 0; i < n*pp[1] && w; i++) w = PM.fz(M, w); return w; }, h = 1e-6*(Math.hypot(z[0], z[1]) || 1), z2 = Fq([z[0] + h, z[1]]), z3 = Fq([z[0] - h, z[1]]);
        var mult = z2 && z3 ? Math.hypot(z2[0] - z3[0], z2[1] - z3[1])/(2*h) : NaN; if (mult > 1) { r.exact = true; r.point = P.p; r.mult = mult; r.res = P.res; r.gap = Math.hypot(P.p[0] - r.end[0], P.p[1] - r.end[1]); } } }
    return r;
  }
  return { fill: fill, fillBox: fillBox, tri: tri, ray: ray, prePer: prePer, polish: polish, polishPar: polishPar, boundaryPoint: boundaryPoint };
})();
if (typeof module !== 'undefined') module.exports = BTRI;
