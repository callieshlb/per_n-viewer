// E96 mating calculator (REGISTRATION addendum 2026-11-10). Depends on pern_math.js (PM).
// spiderCentre: Hubbard-Schleicher spider algorithm (port of E90 pern.spider) + Newton polish: the center of period n with angle k/(2^n - 1).
// misiurewicz: the landing point of a preperiodic parameter ray (PM.traceRay) + Newton on f^(l+p)(0) = f^l(0).
// mate: Thurston pullback on the formal mating f ⊔ g (rays theta of f glued to -theta of g), critical points of f, g at 0, inf.
'use strict';
var MATE = (function () {
  var cmul = PM.cmul, cdiv = PM.cdiv, cadd = PM.cadd, csub = PM.csub, cabs = PM.cabs;
  function cexp2pi(t) { return [Math.cos(2*Math.PI*t), Math.sin(2*Math.PI*t)]; }
  function csqrt(z) { var r = Math.hypot(z[0], z[1]), t = Math.atan2(z[1], z[0])/2, s = Math.sqrt(r); return [s*Math.cos(t), s*Math.sin(t)]; }
  function newtonCentre(c, n) {
    for (var it = 0; it < 80; it++) { var z = [0, 0], dz = [0, 0];
      for (var i = 0; i < n; i++) { dz = cadd(cmul([2*z[0], 2*z[1]], dz), [1, 0]); z = cadd(cmul(z, z), c); }
      var d = cdiv(z, dz); c = csub(c, d); if (cabs(d) < 1e-16*(1 + cabs(c))) break; }
    return c;
  }
  function spiderCentre(k, n) {
    if (n === 1) return [0, 0];
    var D = Math.pow(2, n) - 1, th = [], R = 1e4, J = 12, rad = [];
    for (var j = 0; j < n; j++) th.push(((k*Math.pow(2, j)) % D)/D);
    for (var i = 0; i <= J; i++) rad.push(Math.pow(R, Math.pow(0.5, i)));
    var legs = [];
    for (j = 0; j < n; j++) { var u = cexp2pi(th[j]); legs.push({ pts: rad.map(function (r) { return [r*u[0], r*u[1]]; }), foot: j !== n - 1 ? u : [0, 0] }); }
    var cOld = legs[0].foot;
    for (var it = 0; it < 600; it++) {
      var c = legs[0].foot, nw = [];
      for (j = 0; j < n; j++) {
        var src = legs[(j + 1) % n], u2 = cexp2pi(th[j]), seq = src.pts.map(function (p) { return csub(p, c); }), last = src.pts[src.pts.length - 1];
        for (var t = 1; t <= 8; t++) seq.push(csub(cadd(last, [(src.foot[0] - last[0])*t/8, (src.foot[1] - last[1])*t/8]), c));
        var z = csqrt(seq[0]), target = [Math.sqrt(R)*u2[0], Math.sqrt(R)*u2[1]];
        if (cabs(csub(z, target)) > cabs(cadd(z, target))) z = [-z[0], -z[1]];
        var lifted = [z];
        for (var s = 1; s < seq.length; s++) { var y = csqrt(seq[s]); if (cabs(csub(y, z)) > cabs(cadd(y, z))) y = [-y[0], -y[1]]; lifted.push(y); z = y; }
        nw.push({ pts: [[R*u2[0], R*u2[1]]].concat(lifted.slice(0, J)), foot: j === n - 1 ? [0, 0] : lifted[lifted.length - 1] });
      }
      legs = nw; c = legs[0].foot;
      if (it > 8 && cabs(csub(c, cOld)) < 1e-11) break;
      cOld = c;
    }
    return newtonCentre(legs[0].foot, n);
  }
  // angle p/q (BigInt) -> PCF parameter: {c, kind: 'centre' | 'misiurewicz', period, preperiod}
  function pcfFromAngle(a) {
    var q = a.q, p = a.p % q, l = 0;
    while (q % 2n === 0n && q > 1n) { q /= 2n; l++; }                // q = 2^l * odd
    if (q === 1n) { if (p === 0n) return { c: [0, 0], kind: 'centre', period: 1, preperiod: 0 }; }
    var per = PM.anglePeriod({ p: p % a.q, q: a.q }) || null;
    if (l === 0) {
      // periodic angle k/q: write k/(2^n - 1)
      var n = 1, t = 1n; while ((2n**BigInt(n) - 1n) % q !== 0n && n < 60) n++;
      var k = Number(p*((2n**BigInt(n) - 1n)/q));
      return { c: spiderCentre(k, n), kind: 'centre', period: n, preperiod: 0, k: k };
    }
    var r = PM.traceRay('param', a, null), c = r.pts[r.pts.length - 1]; var c0 = c;
    // period m of the eventual cycle of the angle, preperiod l (exact in the angle; the critical value's preperiod is l)
    var odd = q, m = 1; while ((2n**BigInt(m) - 1n) % odd !== 0n && m < 60) m++;
    for (var it2 = 0; it2 < 60; it2++) {                            // Newton on f^(l+1+m)(0) - f^(l+1)(0)
      var z = [0, 0], dz = [0, 0], zl = null, dzl = null;
      for (var i = 0; i < l + 1 + m; i++) { dz = cadd(cmul([2*z[0], 2*z[1]], dz), [1, 0]); z = cadd(cmul(z, z), c); if (i === l) { zl = z; dzl = dz; } }   // FIX 2026-11-10: the critical VALUE has preperiod l, so f^(l+1+m)(0) = f^(l+1)(0)
      var F = csub(z, zl), dF = csub(dz, dzl), d = cdiv(F, dF); c = csub(c, d); if (cabs(d) < 1e-15) break;
    }
    return { c: c, kind: 'misiurewicz', period: m, preperiod: l, rayLanded: r.landed, rayEnd: c0, moved: cabs(csub(c, c0)) };
  }
  // critical orbit data: list of points (index 0 = critical point 0, then c, ...), F[i] = index of the image
  function orbitOf(c, maxLen) {
    var pts = [[0, 0]], z = [0, 0];
    for (var i = 0; i < (maxLen || 80); i++) {
      z = cadd(cmul(z, z), c);
      for (var j = 0; j < pts.length; j++) if (cabs(csub(z, pts[j])) < 1e-8*(1 + cabs(z))) {
        var F = pts.map(function (_, i2) { return i2 + 1; }); F[pts.length - 1] = j; return { pts: pts, F: F };
      }
      pts.push(z);
    }
    return null;                                                    // not PCF to this length
  }
  // Thurston pullback. f: center c_f; g: PCF c_g. Returns {ok, Z: {f: [...], g: [...]}, R: {v1, v2, mu}, iters, residual, ...}
  function mate(cf, cg, opts) {
    opts = opts || {}; var S = opts.S || 200, maxIt = opts.maxIt || 4000;
    var Of = orbitOf(cf), Og = orbitOf(cg); if (!Of || !Og) return { ok: false, why: 'not post-critically finite' };
    if (Of.F[Of.pts.length - 1] !== 0) return { ok: false, why: 'f: critical point not periodic' };
    var nf = Of.pts.length, ng = Og.pts.length, s = cdiv([1, 0], cf);
    // positions; f-side index i, g-side index i (g[0] = critical point at infinity = null)
    var Zf = Of.pts.map(function (x) { return cmul(x, s); }), Zg = Og.pts.map(function (y, i) { return i === 0 ? null : cmul(cdiv([S, 0], y), s); });
    var hist = [], it, change = Infinity, mu = null, v1, v2;
    function pre(w, mu, v1, v2) {             // z^2 for R(z) = w; w, v2 may be null (infinity)
      if (w === null) return v2 === null ? null : [-mu[0], -mu[1]];
      if (v2 === null) return cmul(mu, csub(w, v1));                 // R(z) = (z^2 + mu v1)/mu, i.e. z^2 = mu (w - v1)
      var den = csub(v2, w); if (cabs(den) === 0) return null; return cmul(mu, cdiv(csub(w, v1), den));
    }
    function root(z2, old) { if (z2 === null) return null; var r = csqrt(z2); if (old !== null && cabs(csub(r, old)) > cabs(cadd(r, old))) r = [-r[0], -r[1]]; return r; }
    for (it = 0; it < maxIt; it++) {
      v1 = Zf[Of.F[0]]; v2 = Zg[Og.F[0]];                            // positions of the critical values
      // normalisation: new position of f's critical value (index 1) is 1
      var w1 = Of.F[1] === 0 ? Zf[0] : Zf[Of.F[1]];
      var base = pre(w1, [1, 0], v1, v2);                            // z^2 at mu = 1
      if (base === null || cabs(base) === 0) return { ok: false, why: 'normalisation degenerate', iters: it };
      mu = cdiv([1, 0], base);
      var nZf = Zf.map(function (_, i) { if (i === 0) return [0, 0]; if (i === 1) return [1, 0]; return root(pre(Zf[Of.F[i]], mu, v1, v2), Zf[i]); });
      var nZg = Zg.map(function (_, i) { if (i === 0) return null; var t = Og.F[i] === 0 ? null : Zg[Og.F[i]]; return root(pre(t, mu, v1, v2), Zg[i]); });
      change = 0;
      for (var i3 = 0; i3 < nf; i3++) change = Math.max(change, PM.chord(nZf[i3], Zf[i3]));
      for (i3 = 0; i3 < ng; i3++) change = Math.max(change, PM.chord(nZg[i3], Zg[i3]));
      Zf = nZf; Zg = nZg;
      // degeneration: two marked points collide
      var all = Zf.concat(Zg), dmin = Infinity;
      for (var a = 0; a < all.length; a++) for (var b = a + 1; b < all.length; b++) dmin = Math.min(dmin, PM.chord(all[a], all[b]));
      if (dmin < 1e-9) return { ok: false, why: 'marked points collide (degenerating: obstruction?)', iters: it, dmin: dmin, Zf: Zf, Zg: Zg };
      if (change < 1e-15) break;
    }
    v1 = Zf[Of.F[0]]; v2 = Zg[Og.F[0]];
    var Rm = function (z) {                                          // the mating as a map of the sphere
      if (z === null) return v2;
      var z2 = cmul(z, z); if (v2 === null) return cadd(cdiv(z2, mu), v1);
      var den = cadd(z2, mu); if (cabs(den) === 0) return null; return cdiv(cadd(cmul(v2, z2), cmul(mu, v1)), den); };
    var res = 0;
    Zf.forEach(function (z, i) { res = Math.max(res, PM.chord(Rm(z), Zf[Of.F[i]])); });
    Zg.forEach(function (z, i) { res = Math.max(res, PM.chord(Rm(z), Og.F[i] === 0 ? null : Zg[Og.F[i]])); });
    return { ok: change < 1e-12 && res < 1e-10, iters: it, change: change, residual: res, Zf: Zf, Zg: Zg, nf: nf, ng: ng, Of: Of, Og: Og, mu: mu, v1: v1, v2: v2, R: Rm };
  }
  // the mating in Per_n coordinates (n = period of f): chart points a_k = Zf[k]; Moebius T: a0 -> 0, a1 -> inf, a2 -> 1.
  function toChart(m) {
    var n = m.nf, a = m.Zf, om = null;                               // g's critical point (the free one) is at infinity
    if (n === 1) return null;
    if (n === 2) {                                                   // Caroline's Per_2 form: a0 -> inf, a1 -> 1, omega -> 0; c = f(0)
      var T2 = PM.moebius3([a[0], a[1], om], [null, [1, 0], [0, 0]]);
      return { n: 2, family: 'per2c', p: T2(m.R(om)) };
    }
    var T = PM.moebius3([a[0], a[1], a[2]], [[0, 0], null, [1, 0]]);
    var y = n >= 4 ? T(a[3]) : [0, 0], x = T(a[n - 1]), w = T(om);
    var out = { n: n, x: x, y: y, omega: w };
    if (n === 3) { out.family = 'cor3c'; out.p = w; }                // c = the free critical point
    else if (n === 4) { out.family = 'cor4c'; out.p = cdiv([0.1, 0], csub(y, [(Math.sqrt(5) - 1)/2, 0])); }
    else out.family = n === 5 ? 'torus5' : 'cor';
    return out;
  }
  // Lavaurs' algorithm: the pairs of characteristic angles [k0, k1] / (2^n - 1) of all period-n components, n = 2..P (cached, extended on demand).
  // Smallest unused angle of period n is joined to the next unused angle of period n such that the chord crosses no earlier chord.
  var LAV = { P: 1, pairs: {}, ends: [], part: {} };                    // ends: sorted angles (float) of all chords so far; part[angle] = partner
  function lavaursOld(P) {
    for (var n = LAV.P + 1; n <= P; n++) { var D = Math.pow(2, n) - 1, ang = [], k, d;
      for (k = 1; k < D; k++) { var per = n; for (d = 1; d < n; d++) if (n % d === 0 && (k*(Math.pow(2, d) - 1)) % D === 0) { per = d; break; } if (per === n) ang.push(k); }
      var used = new Uint8Array(ang.length), out = [], lo = function (x) { var a = 0, b = LAV.ends.length; while (a < b) { var m = (a + b) >> 1; if (LAV.ends[m] < x) a = m + 1; else b = m; } return a; };
      for (var i = 0; i < ang.length; i++) { if (used[i]) continue; var A = ang[i]/D;
        for (var j = i + 1; j < ang.length; j++) { if (used[j]) continue; var B = ang[j]/D, ok = true;
          for (var e = lo(A); e < LAV.ends.length && LAV.ends[e] < B; e++) { var q = LAV.part[LAV.ends[e]]; if (q < A || q > B) { ok = false; break; } }
          if (!ok) continue;
          used[i] = used[j] = 1; out.push([ang[i], ang[j]]); LAV.ends.splice(lo(A), 0, A); LAV.ends.splice(lo(B), 0, B); LAV.part[A] = B; LAV.part[B] = A; break; } }
      LAV.pairs[n] = out; LAV.P = n; }
    return LAV.pairs;
  }
  // Lavaurs' algorithm, linear time per period (2026-11-12; the first version was quadratic and stopped being usable at period 16): the chords of lower
  // period cut the disk into gaps; the period-n angles of one gap, in increasing order, are joined 1st-2nd, 3rd-4th, ... A sweep over all endpoints with a
  // stack of open chords gives each new angle its gap (the innermost open chord). Checked against the old version for n <= 14 (tests/m10_lavaurs_fast.js).
  var LAV2 = { P: 1, pairs: { 1: [] }, ends: [] };          // ends: [angle, +id (opening) or -id (closing)], sorted by angle
  function lavaurs(P) {
    for (var n = LAV2.P + 1; n <= P; n++) { var D = Math.pow(2, n) - 1, ang = [], k, d, divs = [];
      for (d = 1; d < n; d++) if (n % d === 0) divs.push(Math.pow(2, d) - 1);
      for (k = 1; k < D; k++) { var prim = true; for (d = 0; d < divs.length; d++) if ((k*divs[d]) % D === 0) { prim = false; break; } if (prim) ang.push(k); }
      var E = LAV2.ends, ei = 0, stack = [], last = {}, out = [], nid = E.length/2 + 1;
      for (var i = 0; i < ang.length; i++) { var A = ang[i]/D;
        while (ei < E.length && E[ei][0] < A) { if (E[ei][1] > 0) stack.push(E[ei][1]); else stack.pop(); ei++; }
        var g = stack.length ? stack[stack.length - 1] : 0;
        if (last[g] === undefined) last[g] = ang[i]; else { out.push([last[g], ang[i]]); last[g] = undefined; } }
      out.sort(function (u, v) { return u[0] - v[0]; });
      var add = []; out.forEach(function (pr) { var id = nid++; add.push([pr[0]/D, id], [pr[1]/D, -id]); }); add.sort(function (u, v) { return u[0] - v[0]; });
      var M = [], x = 0, y = 0; while (x < E.length || y < add.length) { if (y >= add.length || (x < E.length && E[x][0] < add[y][0])) M.push(E[x++]); else M.push(add[y++]); }
      LAV2.ends = M; LAV2.pairs[n] = out; LAV2.P = n; }
    return LAV2.pairs;
  }
  // Vein of theta = a.p/a.q (REGISTRATION addendum 2026-11-10 (b), M7; extended (d), M9): the centers of period 2..P whose wake (k0/Q, k1/Q) strictly
  // contains theta, i.e. the hyperbolic components the combinatorial arc from 0 to the landing point of theta passes through. The K of lowest period
  // (ties: wider wake first) are kept and returned in order along the arc (widest wake first). polys[n] = poly_info_n<n>.json (n <= 8, exact c);
  // for n > 8 the center comes from the spider algorithm (c: null until centerOf is called).
  function veinList(a, polys, K, P) {
    var out = [], th = Number(a.p)/Number(a.q), pairs = lavaurs(P || 8);
    for (var n = 2; n <= (P || 8); n++) { var Q = Math.pow(2, n) - 1, Qb = BigInt(Q);
      pairs[n].forEach(function (pr) { var k0 = BigInt(pr[0]), k1 = BigInt(pr[1]);
        if (k0*a.q < a.p*Qb && a.p*Qb < k1*a.q) out.push({ n: n, k0: pr[0], k1: pr[1], Q: Q, width: (pr[1] - pr[0])/Q, c: null }); }); }
    out.sort(function (u, v) { return u.n - v.n || v.width - u.width; });
    if (K !== undefined && K !== null) out = out.slice(0, K);
    out.sort(function (u, v) { return v.width - u.width; });
    out.forEach(function (v) { var D = polys && polys[v.n], Pn = D && D.find(function (q) { return q.ks.indexOf(v.k0) >= 0; });
      if (Pn) { v.c = [Number(Pn.c_str[0]), Number(Pn.c_str[1])]; v.label = Pn.label; } });
    return out;
  }
  function centerOf(v) {                     // center of a vein entry (spider for n > 8); null if the result does not have exact period n
    if (v.c) return v.c; var c = spiderCentre(v.k0, v.n), o = orbitOf(c, v.n + 2);
    if (!o || o.pts.length !== v.n || o.F[v.n - 1] !== 0) return null; v.c = c; return c;
  }
  // period of the free critical point (infinity) under the mated map m.R (chordal tolerance), or 0
  function omegaPeriod(m, maxP) { var z = null;
    for (var i = 1; i <= (maxP || 40); i++) { z = m.R(z); if (z === null || PM.chord(z, [1e9, 0]) < 1e-7) return i; } return 0; }
  return { lavaurs: lavaurs, lavaursOld: lavaursOld, centerOf: centerOf, veinList: veinList, omegaPeriod: omegaPeriod, spiderCentre: spiderCentre, pcfFromAngle: pcfFromAngle, orbitOf: orbitOf, mate: mate, toChart: toChart, newtonCentre: newtonCentre };
})();
if (typeof module !== 'undefined') module.exports = MATE;
