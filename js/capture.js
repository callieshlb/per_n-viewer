// E96 captures (REGISTRATION addendum 2026-11-10 (g)). Depends on pern_math.js (PM) and mating.js (MATE.toChart).
'use strict';
var CAP = (function () {
  var cmul = PM.cmul, cdiv = PM.cdiv, cadd = PM.cadd, csub = PM.csub, cabs = PM.cabs;
  function csqrt(z) { var r = Math.hypot(z[0], z[1]), t = Math.atan2(z[1], z[0])/2, s = Math.sqrt(r); return [s*Math.cos(t), s*Math.sin(t)]; }
  var fc = function (z, c) { return cadd(cmul(z, z), c); };
  function cycle(c, n) { var o = [[0, 0]]; for (var i = 1; i < n; i++) o.push(fc(o[i - 1], c)); return o; }
  // exact landing point of the dynamic ray of angle a (periodic or preperiodic): ray end + Newton on f^(l+m)(z) = f^l(z)
  function landing(c, a, deep) {
    var r = PM.traceRay('dyn', a, c, 8, deep ? 2400 : undefined), z = r.pts[r.pts.length - 1], q = a.q, l = 0; while (q % 2n === 0n) { q /= 2n; l++; }
    var m = 1; while ((2n**BigInt(m) - 1n) % q !== 0n && m < 60) m++;
    for (var it = 0; it < 60; it++) { var w = z, dw = [1, 0], wl = null, dl = null;
      for (var i = 0; i < l + m; i++) { if (i === l) { wl = w; dl = dw; } dw = cmul([2*w[0], 2*w[1]], dw); w = fc(w, c); }
      if (l === 0) { wl = z; dl = [1, 0]; } var F = csub(w, wl), dF = csub(dw, dl), d = cdiv(F, dF); z = csub(z, d); if (cabs(d) < 1e-14) break; }
    return { z: z, ray: r.pts, moved: cabs(csub(z, r.pts[r.pts.length - 1])), pre: l, per: m };
  }
  // center of the Fatou component of w (basin of the n-cycle of z^2 + c), with the path from w to it; null if w escapes / undecided
  function centerOf(c, n, cyc, w, rho) {
    var z = w, t, k = -1;
    for (t = 0; t <= 400; t++) { if (cabs(z) > 50) return null; for (var j = 0; j < n; j++) if (cabs(csub(z, cyc[j])) < rho) { k = j; break; } if (k >= 0) break; z = fc(z, c); }
    if (k < 0) return null;
    // orbit w_0 .. w_t; lift the segment w_t -> cyc[k] back one iterate at a time (z -> +-sqrt(z - c), by continuity from w_(i-1))
    var orb = [w], i, it, u, du, d; for (i = 0; i < t; i++) orb.push(fc(orb[i], c));
    var path = [], NP = 48; for (i = 0; i <= NP; i++) path.push(cadd(cyc[k], [(1 - i/NP)*(orb[t][0] - cyc[k][0]), (1 - i/NP)*(orb[t][1] - cyc[k][1])]));
    for (var lev = t - 1; lev >= 0; lev--) { var prev = orb[lev], np = [];
      var prevSrc = null, lift = function (pt, dep) {                    // adaptive: insert midpoints (in the path being lifted) while the nearest root is not clear
        var r = csqrt(csub(pt, c)), dm = cabs(csub(r, prev)), dp = cabs(cadd(r, prev)); if (dp < dm) { r = [-r[0], -r[1]]; var tt = dm; dm = dp; dp = tt; }
        if (prevSrc && dp > 0 && dm > 0.3*dp && dep < 40) { var mid = [(prevSrc[0] + pt[0])/2, (prevSrc[1] + pt[1])/2]; lift(mid, dep + 1); lift(pt, dep + 1); return; }
        np.push(r); prev = r; prevSrc = pt; };
      for (i = 0; i < path.length; i++) lift(path[i], 0);
      if (np.length > 20000) return null;
      var eL = np[np.length - 1];                                       // snap the end point to an exact cycle point (or its negative) when it is one: rounding is amplified by sqrt at c otherwise
      for (var js = 0; js < n; js++) { if (cabs(csub(eL, cyc[js])) < 1e-6) np[np.length - 1] = cyc[js].slice(); else if (cabs(cadd(eL, cyc[js])) < 1e-6 && js !== 0) np[np.length - 1] = [-cyc[js][0], -cyc[js][1]]; }
      path = np; }
    var x = path[path.length - 1];
    var chk = x; for (i = 0; i < t; i++) chk = fc(chk, c); if (cabs(csub(chk, cyc[k])) > 1e-6) return null;
    // polish: the first time the center's orbit meets the cycle (f^t' is univalent near x; the lift above is degenerate when the orbit crosses U_0)
    var y = x, t1 = -1, k1 = -1; for (i = 0; i <= t && t1 < 0; i++) { for (var j2 = 0; j2 < n; j2++) if (cabs(csub(y, cyc[j2])) < 1e-3) { t1 = i; k1 = j2; } y = fc(y, c); }
    if (t1 < 0) return null;
    if (t1 === 0) x = cyc[k1].slice(); else for (it = 0; it < 60; it++) { u = x; du = [1, 0]; for (i = 0; i < t1; i++) { du = cmul([2*u[0], 2*u[1]], du); u = fc(u, c); } d = cdiv(csub(u, cyc[k1]), du); x = csub(x, d); if (cabs(d) < 1e-15) break; }
    path[path.length - 1] = x;
    return { x: x, t: t1, k: k1, path: path };
  }
  // Fatou components at the landing point z: [{x, t, k, frac, path (from z into the component to x)}], best first
  function componentsAt(c, n, z) {
    var cyc = cycle(c, n), dmin = Infinity, i, j; for (i = 0; i < n; i++) for (j = i + 1; j < n; j++) dmin = Math.min(dmin, cabs(csub(cyc[i], cyc[j])));
    var rho = 0.02*dmin, found = [], NS = 96, radii = [3e-3, 1e-3];
    radii.forEach(function (r, ri) { for (var s = 0; s < NS; s++) { var w = [z[0] + r*Math.cos(2*Math.PI*s/NS), z[1] + r*Math.sin(2*Math.PI*s/NS)], R = centerOf(c, n, cyc, w, rho); if (!R) continue;
        var e = null; for (var q = 0; q < found.length; q++) if (cabs(csub(found[q].x, R.x)) < 1e-7) e = found[q];
        if (!e) { e = { x: R.x, t: R.t, k: R.k, cnt: [0, 0], best: [null, null] }; found.push(e); } e.cnt[ri]++; e.best[ri] = e.best[ri] || []; e.best[ri].push(R); } });
    var out = found.filter(function (e) { return e.cnt[0] >= 0.08*NS && e.cnt[1] >= 0.08*NS; }).map(function (e) { var L = e.best[1], mid = L[Math.floor(L.length/2)];
      var per = -1; for (var q = 0; q < n; q++) if (cabs(csub(e.x, cyc[q])) < 1e-8) per = q;
      return { x: e.x, t: e.t, k: e.k, frac: Math.min(e.cnt[0], e.cnt[1])/NS, path: [z].concat(mid.path), cycIndex: per }; });
    out.sort(function (u, v) { return v.frac - u.frac; }); return out;
  }
  // ---- the pullback machinery. Critical points 0, inf. R(z) = (v2 z^2 + kappa)/(z^2 + lambda), kappa/lambda = v1 (v1 = inf: lambda = 0).
  // state: v1, v2 (points or null), tg[p] (targets, point or null), q (normalising index), Zq. roots(): new positions nearest to Z (ambiguity reported)
  function solve(Z, st) {
    var wq = st.tg[st.q], zq2 = cmul(Z[st.q], Z[st.q]), mu;
    if (st.v1 === null) mu = wq === null ? null : cmul(zq2, csub(wq, st.v2));
    else mu = wq === null ? [-zq2[0], -zq2[1]] : cdiv(cmul(zq2, csub(wq, st.v2)), csub(st.v1, wq));
    if (mu === null || !isFinite(mu[0] + mu[1])) return null;
    var kappa = st.v1 === null ? mu : cmul(mu, st.v1), lambda = st.v1 === null ? [0, 0] : mu, out = [], amb = 0;
    for (var p = 0; p < Z.length; p++) { var w = st.tg[p], z2;
      if (st.fixed[p] !== undefined) { out.push(st.fixed[p]); continue; }
      if (w === null) z2 = [-lambda[0], -lambda[1]]; else { var den = csub(w, st.v2); if (cabs(den) < 1e-300) return null; z2 = cdiv(csub(kappa, cmul(lambda, w)), den); }
      var r = csqrt(z2), old = Z[p], dm = cabs(csub(r, old)), dp = cabs(cadd(r, old)); if (dp < dm) { r = [-r[0], -r[1]]; var tt = dm; dm = dp; dp = tt; }
      amb = Math.max(amb, dp > 0 ? dm/dp : 1); out.push(r); }
    return { Z: out, amb: amb, mu: mu, kappa: kappa, lambda: lambda };
  }
  var lerp = function (a, b, t) { if (a === null || b === null) return t < 1 ? a : b; return [a[0] + t*(b[0] - a[0]), a[1] + t*(b[1] - a[1])]; };
  function mix(A, B, t) { return { v1: lerp(A.v1, B.v1, t), v2: lerp(A.v2, B.v2, t), tg: A.tg.map(function (w, i) { return lerp(w, B.tg[i], t); }), q: B.q, fixed: B.fixed }; }
  // continue the roots from state A (roots Z) to state B, subdividing while the nearest root is not clear
  function go(Z, A, B, depth, stat) {
    var S = solve(Z, B); stat.solves++;
    if (S && S.amb < 0.25) return S;
    if (depth > 26) { stat.forced++; return S; }
    var Mid = mix(A, B, 0.5), H = go(Z, A, Mid, depth + 1, stat); if (!H) return null; return go(H.Z, Mid, B, depth + 1, stat);
  }
  // capture of z^2 + c (period n) at the angle a (BigInt p/q), optional choice among several components. opts: R0 (start radius), sub (path refinement)
  function capture(c, n, a, opts) {
    opts = opts || {}; var cyc = cycle(c, n), L = landing(c, a); if (L.moved > 1e-3) L = landing(c, a, true); if (L.moved > 0.05) return { ok: false, why: 'the dynamic ray did not land well (Newton moved ' + L.moved.toExponential(1) + ')' };
    if (opts.stopAtLanding) {                                           // addendum (j), T3: the path stops at the landing point: critical value on the Julia set
      var th0 = 2*Math.PI*Number(a.p)/Number(a.q), bt = []; for (var e0 = Math.log10(opts.R0 || 1e8); e0 > 2.01; e0 -= 0.25) bt.push([Math.pow(10, e0)*Math.cos(th0), Math.pow(10, e0)*Math.sin(th0)]);
      return pull(c, n, cyc, { x: L.z, t: -1, k: -1, cycIndex: -1 }, false, bt.concat(L.ray, [L.z]), opts, [], L); }
    var comps = componentsAt(c, n, L.z); if (!comps.length) return { ok: false, why: 'the ray lands at a point that is not on the boundary of a Fatou component (no component fills an arc around the landing point)', landing: L.z };
    var V = comps[Math.min(opts.choice || 0, comps.length - 1)];
    if (opts.wantCenter) { V = null; comps.forEach(function (v) { if (cabs(csub(v.x, opts.wantCenter)) < 1e-6) V = v; }); if (!V) return { ok: false, why: 'the ray does not land on the boundary of the chosen component', comps: comps, landing: L.z }; }
    if (opts.wantCyc !== undefined) { V = null; comps.forEach(function (v) { if (v.cycIndex === opts.wantCyc) V = v; }); if (!V) return { ok: false, why: 'the ray does not land on the boundary of U_' + opts.wantCyc, comps: comps, landing: L.z }; }
    var typeII = V.cycIndex >= 0, i, j;
    if (typeII && V.cycIndex === 1) return { ok: false, why: 'the ray lands on the critical value component U_1 of f itself: no capture', landing: L.z, comps: comps };
    if (typeII && V.cycIndex === 0 && n === 1) return { ok: false, why: 'period 1' };
    var th = 2*Math.PI*Number(a.p)/Number(a.q), R0 = opts.R0 || 1e8, beta = [];
    for (var e = Math.log10(R0); e > 2.01; e -= 0.25) beta.push([Math.pow(10, e)*Math.cos(th), Math.pow(10, e)*Math.sin(th)]);
    beta = beta.concat(L.ray, V.path); var sub = opts.sub || 1;
    if (sub > 1) { var fine = [beta[0]]; for (i = 1; i < beta.length; i++) for (j = 1; j <= sub; j++) fine.push(lerp(beta[i - 1], beta[i], j/sub)); beta = fine; }
    if (sub < 1) { var coarse = []; for (i = 0; i < beta.length; i += Math.round(1/sub)) coarse.push(beta[i]); coarse.push(beta[beta.length - 1]); beta = coarse; }
    return pull(c, n, cyc, V, typeII, beta, opts, comps, L);
  }
  function pull(c, n, cyc, V, typeII, beta, opts, comps, L) {
    var i, j;
    // marked points: names + positions in the plane of f, images
    var P = [], img = [], crit2 = -1;                                    // P[i] = position (null = inf); img[i] = index of G(P[i]); crit2 = index of the second critical point (inf)
    if (!typeII) { for (i = 0; i < n; i++) { P.push(cyc[i]); img.push((i + 1) % n); }
      var x = V.x, first = P.length, guard = 0; while (guard++ < 500) { var hit = -1; for (j = 0; j < P.length; j++) if (P[j] && cabs(csub(P[j], x)) < 1e-8) hit = j; if (hit >= 0) { if (P.length > first) img[P.length - 1] = hit; else first = hit; break; } P.push(x); img.push(P.length); x = fc(x, c); }
      P.push(null); img.push(first); crit2 = P.length - 1;
    } else { var slot = (V.cycIndex - 1 + n) % n;                          // the second critical point takes the place of a_slot; its value is a_(slot+1) = x
      if (slot === 0) return { ok: false, why: 'this would put both critical points at a_0' };
      for (i = 0; i < n; i++) { P.push(i === slot ? null : cyc[i]); img.push((i + 1) % n); } crit2 = slot; }
    var N = P.length, stat = { solves: 0, forced: 0 };
    var q = -1; for (i = 1; i < N && q < 0; i++) if (P[i] !== null && i !== 0) q = i; if (q < 0) return { ok: false, why: 'no finite marked point to normalize with (period 2, capture at U_0: the orbifold point)' };
    var fixed = {}; fixed[0] = [0, 0]; fixed[crit2] = null; fixed[q] = P[q];
    var Z = P.map(function (p) { return p; });
    // ---- first pullback: targets = positions in the plane of f; v2 runs along beta; (type II) the target of the point before the slot runs back along zeta
    var before = typeII ? (crit2 - 1 + n) % n : -1, zeta = null;
    if (typeII) { zeta = []; var zprev = null;                             // zeta(s): the preimage of beta(s) continued from inf (far out: sqrt branch nearest the ray direction th/2 or th/2 + pi chosen at the END)
      // build from the end: zeta ends at cyc[crit2]; lift beta backwards from its end point by continuity
      var zz = cyc[crit2]; for (i = beta.length - 1; i >= 0; i--) { var rt = csqrt(csub(beta[i], c)); if (cabs(csub(rt, zz)) > cabs(cadd(rt, zz))) rt = [-rt[0], -rt[1]]; zeta[i] = rt; zz = rt; } }
    var mkState = function (s) { var tg = P.map(function (_, p) { var w = P[img[p]]; return w; });
      if (typeII) { tg[before] = s >= beta.length - 1 ? null : zeta[beta.length - 1 - s]; tg[crit2] = beta[s]; } else tg[crit2] = beta[s];
      return { v1: tg[0], v2: beta[s], tg: tg, q: q, fixed: fixed }; };
    var A = mkState(0), S0 = solve(Z, A); if (!S0) return { ok: false, why: 'start failed' }; Z = S0.Z;
    for (var s = 1; s < beta.length; s++) { var B = mkState(s), S = go(Z, A, B, 0, stat); if (!S) return { ok: false, why: 'the continuation along the path broke down at step ' + s + ' of ' + beta.length, comps: comps, landing: L.z }; Z = S.Z; A = B; }
    // ---- further pullbacks
    var stateOf = function (Zc) { var tg = Zc.map(function (_, p) { return Zc[img[p]]; }); return { v1: tg[0], v2: tg[crit2], tg: tg, q: q, fixed: fixed }; };
    var chordZ = function (U, W) { var m = 0; for (var p = 0; p < N; p++) m = Math.max(m, PM.chord(U[p], W[p])); return m; };
    var it, change = Infinity, last = null, dminHist = [];
    for (it = 0; it < (opts.maxIt || 3000); it++) { var B2 = stateOf(Z), S2 = go(Z, A, B2, 0, stat); if (!S2) return { ok: false, why: 'pullback broke down at iteration ' + it, comps: comps };
      change = chordZ(S2.Z, Z); Z = S2.Z; A = B2; last = S2;
      var dmin = Infinity; for (i = 0; i < N; i++) for (j = i + 1; j < N; j++) dmin = Math.min(dmin, PM.chord(Z[i], Z[j]));
      if (dmin < 1e-9) return { ok: false, why: 'marked points collide (degenerating: obstruction?)', iters: it, dmin: dmin, comps: comps, landing: L.z, typeII: typeII, component: V };
      if (change < 1e-15) break; }
    var fin = stateOf(Z), Sf = solve(Z, fin); if (!Sf) return { ok: false, why: 'final map degenerate' };
    var Rm = function (z) { if (z === null) return fin.v2; var z2 = cmul(z, z), den = cadd(z2, Sf.lambda); if (cabs(den) === 0) return null; return cdiv(cadd(cmul(fin.v2 === null ? [1, 0] : fin.v2, z2), Sf.kappa), den); };
    var res = 0; for (i = 0; i < N; i++) res = Math.max(res, PM.chord(Rm(Z[i]), Z[img[i]]));
    var Zc = []; for (i = 0; i < n; i++) Zc.push(Z[i]);
    var m = { ok: change < 1e-12 && res < 1e-10, iters: it, change: change, residual: res, nf: n, Zf: Zc, R: Rm, Z: Z, img: img, typeII: typeII, component: V, comps: comps, landing: L.z, stat: stat,
      predicted: typeII ? { slot: crit2 } : { t: V.t + 1, k: V.k } };
    if (!m.ok) m.why = 'did not converge (change ' + change.toExponential(1) + ', residual ' + res.toExponential(1) + ')';
    return m;
  }
  // the bounded Fatou component containing the point w: {x, t, k, cycIndex} or null
  function componentOf(c, n, w) { var cyc = cycle(c, n), dmin = Infinity, i, j; for (i = 0; i < n; i++) for (j = i + 1; j < n; j++) dmin = Math.min(dmin, cabs(csub(cyc[i], cyc[j])));
    var R = centerOf(c, n, cyc, w, 0.02*dmin); if (!R) return null; var per = -1; for (i = 0; i < n; i++) if (cabs(csub(R.x, cyc[i])) < 1e-8) per = i; return { x: R.x, t: R.t, k: R.k, cycIndex: per }; }
  // the two rays landing at the root of the component V (center V.x, reaching the cycle at a_k after t steps). ks = [k-, k+] numerators over 2^n - 1.
  // Returns {minus: {p, q, text}, plus: {...}, root: point} (BigInt angles), or {why}
  function rootRays(c, n, ks, V) {
    var Q = 2n**BigInt(n) - 1n, k = V.k === 0 ? n : V.k, out = {}, orb = [V.x], i; for (i = 0; i < V.t; i++) orb.push(fc(orb[i], c));
    var sides = [['minus', BigInt(Math.min(ks[0], ks[1]))], ['plus', BigInt(Math.max(ks[0], ks[1]))]];
    for (var si = 0; si < 2; si++) { var a = { p: (sides[si][1]*(2n**BigInt(k - 1))) % Q, q: Q }, L = landing(c, a), r = L.z; if (L.moved > 1e-3) { L = landing(c, a, true); r = L.z; } if (L.moved > 0.05) return { why: 'a root ray did not land well' };
      for (var lev = V.t - 1; lev >= 0; lev--) {                          // lift the root and the angle one step back along the orbit of the center
        var rt = csqrt(csub(r, c)); if (cabs(csub(rt, orb[lev])) > cabs(cadd(rt, orb[lev]))) rt = [-rt[0], -rt[1]]; r = rt;
        var h1 = { p: a.p, q: a.q*2n }, h2 = { p: a.p + a.q, q: a.q*2n }, L1 = landing(c, h1); if (L1.moved > 1e-3) L1 = landing(c, h1, true);
        var d1 = cabs(csub(L1.z, r)), d2 = cabs(cadd(L1.z, r)); if (Math.min(d1, d2) > 1e-4) return { why: 'could not follow the root ray back (landing point ' + Math.min(d1, d2).toExponential(1) + ' from the lifted root)' };
        a = d1 <= d2 ? h1 : h2; }
      var g = function (x, y) { return y ? g(y, x % y) : x; }, gg = g(a.p, a.q); a = { p: a.p/gg, q: a.q/gg }; a.text = a.p + '/' + a.q; out[sides[si][0]] = a; out.root = r; }
    return out;
  }
  return { capture: capture, componentsAt: componentsAt, componentOf: componentOf, rootRays: rootRays, landing: landing, cycle: cycle };
})();
if (typeof module !== 'undefined') module.exports = CAP;
