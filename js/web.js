// E96 features 6 / 7 (REGISTRATION addendum 2026-11-12 (u)), beta: V_f (veins from a polynomial f to its periodic captures) and T'_f (the image in
// Per_n, under mating, of the angles landing on the ringed, pruned Hubbard tree of f). Depends on PM, MATE, CAP. Pure functions (node-testable).
'use strict';
var WEBS = (function () {
  var gcd = function (a, b) { while (b) { var t = a % b; a = b; b = t; } return a; };
  var frac = function (p, q) { var g = gcd(p < 0n ? -p : p, q); return { p: p/g, q: q/g, text: (p/g) + '/' + (q/g) }; };
  var num = function (a) { return Number(a.p)/Number(a.q); };
  // binary digits of t = p/q in [0, 1]: { pre: '...', per: '...' } (t = 1 is given as pre '', per '1')
  function digits(a) { if (a.p === a.q) return { pre: '', per: '1' }; var x = a.p % a.q, seen = new Map(), bits = '', i = 0;
    while (!seen.has(x)) { seen.set(x, i++); x = 2n*x; if (x >= a.q) { bits += '1'; x -= a.q; } else bits += '0'; } var s = seen.get(x); return { pre: bits.slice(0, s), per: bits.slice(s) }; }
  function fromDigits(pre, per) { var a = BigInt(pre.length), b = BigInt(per.length), A = pre.length ? BigInt('0b' + pre) : 0n, B = BigInt('0b' + per); return frac(A*((1n << b) - 1n) + B, (1n << a)*((1n << b) - 1n)); }
  // angles landing on the boundary of the critical VALUE component U_1 of the center P (period n): tau(t) = t with digit 0 -> B-, 1 -> B+ (Douady tuning)
  function tau(P, n, t) { var ks = P.ks.slice().sort(function (u, v) { return u - v; }), Bm = ks[0].toString(2).padStart(n, '0'), Bp = ks[1].toString(2).padStart(n, '0'), d = digits(t), sub = function (s) { return s.split('').map(function (c) { return c === '0' ? Bm : Bp; }).join(''); };
    return fromDigits(sub(d.pre), sub(d.per)); }
  function dbl(a, k) { var p = (a.p*(1n << BigInt(k))) % a.q; return frac(p, a.q); }
  // inner angles t of preperiod l', period p', l' + p' <= depth (plus 0 and 1)
  function innerAngles(depth, perDepth) { var out = [{ p: 0n, q: 1n, text: '0' }, { p: 1n, q: 1n, text: '1' }], seen = new Set(['0/1', '1/1']);
    for (var l = 0; l <= depth - 1; l++) for (var p = 1; p <= (l === 0 ? Math.max(depth, perDepth || 0) : depth - l); p++) { var q = (1n << BigInt(l))*((1n << BigInt(p)) - 1n); for (var k = 1n; k < q; k++) { var a = frac(k, q), key = a.p + '/' + a.q; if (!seen.has(key)) { seen.add(key); out.push(a); } } }
    out.sort(function (u, v) { return num(u) - num(v); }); return out; }
  // ring: for k = 1..n (k = n is U_0), the angles 2^(k-1) tau(t)
  function ringAngles(P, n, depth, perDepth) { var T = innerAngles(depth, perDepth), out = [];
    for (var k = 1; k <= n; k++) T.forEach(function (t) { out.push({ comp: k % n, t: t, th: dbl(tau(P, n, t), k - 1) }); }); return out; }
  // critical value alpha-sector: the wake (k0/Q', k0+1/Q') of the p/q-bulb of the main cardioid containing f's angle (the limb of f); returned as numbers
  function sector(P, n) { var th = Math.min.apply(null, P.ks)/(Math.pow(2, n) - 1), th2 = Math.max.apply(null, P.ks)/(Math.pow(2, n) - 1);
    for (var q = 2; q <= n; q++) { var Q = Math.pow(2, q) - 1, k0 = Math.floor(th*Q + 1e-12); if (gcd(k0, Q) >= 1 && k0/Q <= th + 1e-12 && (k0 + 1)/Q >= th2 - 1e-12) { // the pair (k0, k0+1)/Q must be a period-q ray pair: rotation check
        var orb = [], x = k0; for (var i = 0; i < q; i++) { orb.push(x); x = (2*x) % Q; } if (x === k0 && orb.indexOf(k0 + 1) >= 0) return [k0/Q, (k0 + 1)/Q, q]; } }
    return null; }
  // tree-interior candidates: (pre)periodic angles with l + p <= L whose landing point has >= 2 of them and separates two periodic components
  function treeAngles(P, n, L) { var c = [Number(P.c_str[0]), Number(P.c_str[1])], Q = Math.pow(2, n) - 1, groups = [], out = [], roots = [];
    for (var k = 0; k < n; k++) roots.push([P.orbit_minus[k]/Q, P.orbit_plus[k]/Q]);
    for (var l = 0; l <= L - 1; l++) for (var p = 1; p <= L - l; p++) { var q = (1n << BigInt(l))*((1n << BigInt(p)) - 1n);
        for (var kk = 1n; kk < q; kk++) { var a = frac(kk, q), d = digits(a); if (d.pre.length !== l || d.per.length !== p) continue; var Ld = CAP.landing(c, a); if (!(Ld.moved < 1e-2) || !isFinite(Ld.z[0] + Ld.z[1])) continue;
          var g = null; for (var i = 0; i < groups.length; i++) if (Math.hypot(groups[i].z[0] - Ld.z[0], groups[i].z[1] - Ld.z[1]) < 1e-7) { g = groups[i]; break; } if (!g) { g = { z: Ld.z, a: [] }; groups.push(g); } g.a.push(a); } }
    groups.forEach(function (g) { if (g.a.length < 2) return; var A = g.a.map(num).sort(function (u, v) { return u - v; }), arc = function (t) { for (var i = 0; i < A.length; i++) if (t < A[i] - 1e-12) return i; return 0; }, used = {}, onRoot = false;
      roots.forEach(function (r) { if (A.some(function (t) { return Math.abs(t - r[0]) < 1e-12 || Math.abs(t - r[1]) < 1e-12; })) { onRoot = true; return; } var i0 = arc(r[0]), i1 = arc(r[1]); if (i0 === i1) used[i0] = true; });
      if (onRoot || Object.keys(used).length < 2) return; g.a.forEach(function (a) { out.push({ th: a, z: g.z, rays: g.a.length }); }); });
    return out; }
  // the point of Per_n for the angle th of f: the mating f ⊔ g, g of angle 1 - th (convention of the site; checked in tests/vt_proto.js)
  function place(cf, th) { var ao = frac((th.q - th.p) % th.q, th.q), g = MATE.pcfFromAngle(ao); if (!g || !g.c || !isFinite(g.c[0] + g.c[1]) || (g.kind === 'misiurewicz' && !(g.moved < 0.05))) return { ok: false, why: 'g not found' };
    // which points of g's critical orbit does the mating identify? orbit point 1 + i of g has g-angle 2^i (1 - th), i.e. f-angle 2^i th: group by landing point on J_f
    var d = digits(th), L = d.pre.length + d.per.length, lands = [], rep = {}, merged = 0;
    for (var i = 0; i < L; i++) { var ai = dbl(th, i), Ld = CAP.landing(cf, ai); lands.push(Ld.moved < 0.3 && isFinite(Ld.z[0] + Ld.z[1]) ? Ld.z : null); }   // the ray end can be far from a slowly approached point; Newton polishes it
    for (i = 0; i < L; i++) { rep[1 + i] = 1 + i; if (!lands[i]) continue; for (var j = 0; j < i; j++) if (lands[j] && Math.hypot(lands[i][0] - lands[j][0], lands[i][1] - lands[j][1]) < 1e-7) { rep[1 + i] = rep[1 + j]; merged++; break; } }
    // NOT YET HANDLED (log.md 2026-11-12): when several points of g's critical orbit are identified by the mating (merged > 0), the pullback on the formal
    // mating makes them collide; tracking them as one point (MATE.mate opts.merge) was tried and still collides, so it is not used. Such angles are reported.
    if (g.kind !== 'misiurewicz') merged = 0;        // a center g has its critical orbit in the Fatou set: nothing is identified
    if (merged && !(typeof WEBS_TRYMERGE !== 'undefined' && WEBS_TRYMERGE)) return { ok: false, why: 'the mating identifies ' + (merged + 1) + ' points of the critical orbit of g (not handled yet)', merged: merged, ident: true };
    var m = MATE.mate(cf, g.c, merged ? { merge: rep } : {});
    if (!m.ok) { // glued kind by its outcome (U4): only points of g's critical orbit have collided with each other, f's marked points stay apart
      var glued = false; if (m.Zf && m.Zg && g.kind === 'misiurewicz') { var gg = 0, other = 0, A = m.Zf.map(function (z) { return ['f', z]; }).concat(m.Zg.map(function (z) { return ['g', z]; }));
        for (var u = 0; u < A.length; u++) for (var v = u + 1; v < A.length; v++) if (PM.chord(A[u][1], A[v][1]) < 1e-6) { if (A[u][0] === 'g' && A[v][0] === 'g') gg++; else other++; } glued = gg > 0 && other === 0; }
      return { ok: false, why: glued ? 'the mating identifies points of the critical orbit of g (not handled yet)' : m.why, merged: merged, ident: glued }; }
    var ch = MATE.toChart(m); if (!ch) return { ok: false, why: 'no chart point' }; return { ok: true, x: ch.x, y: ch.y, p: ch.p, kind: g.kind, residual: m.residual, g: ao, merged: merged, m: m }; }
  function inSector(sec, th) { var t = num(th); return sec && t > sec[0] + 1e-12 && t < sec[1] - 1e-12; }
  function buildT(P, n, opts) { opts = opts || {}; var cf = [Number(P.c_str[0]), Number(P.c_str[1])], sec = sector(P, n), R = ringAngles(P, n, opts.ring === undefined ? 3 : opts.ring, opts.ringPer || Math.max(3, Math.floor(20/n))), TR   /* periodic inner angles up to period 20/n: g is then a center of period <= 20 (spider) */ = treeAngles(P, n, opts.tree || 7), dots = [], st = { ring: R.length, tree: TR.length, pruned: 0, obstructed: 0, identified: 0, placed: 0 };
    R.forEach(function (r) { if (inSector(sec, r.th)) { st.pruned++; return; } var q = place(cf, r.th); if (!q.ok) { if (q.ident) st.identified++; else st.obstructed++; return; } st.placed++; dots.push({ x: q.x, y: q.y, p: q.p, part: 'ring', comp: r.comp, t: num(r.t), th: r.th.text, periodic: q.kind === 'centre' }); });
    TR.forEach(function (r) { if (inSector(sec, r.th)) { st.pruned++; return; } var q = place(cf, r.th); if (!q.ok) { if (q.ident) st.identified++; else st.obstructed++; return; } st.placed++; dots.push({ x: q.x, y: q.y, p: q.p, part: 'tree', th: r.th.text, periodic: q.kind === 'centre', rays: r.rays }); });
    return { dots: dots, sector: sec, stats: st }; }
  // V_f: for every periodic capture of f (table entries), the chain f, matings along the vein of 1 - theta, f ⊔ g_(1-theta), the capture's center
  function buildV(P, n, table, polys, mis, K) { var cf = [Number(P.c_str[0]), Number(P.c_str[1])], paths = [], st = { captures: 0, edges: 0, obstructedEnd: 0, veinPts: 0 }, f0 = MATE.toChart(MATE.mate(cf, [0, 0]));
    table.forEach(function (b) { if (!b.x || !b.y) return; b.src.forEach(function (s) { if (s.f !== P.label) return; st.captures++; var a = PM.parseAngle(s.angle), ao = frac(a.q - a.p, a.q), V = [], pts = [];
        for (var Pt = 8; Pt <= 14; Pt += 2) { V = MATE.veinList(ao, polys, K || 10, Pt, mis); if (V.length >= (K || 10)) break; }
        if (f0) pts.push({ x: f0.x, y: f0.y, p: f0.p });
        V.forEach(function (v) { var cg = MATE.centerOf(v); if (!cg) return; var m = MATE.mate(cf, cg); if (!m.ok) return; var ch = MATE.toChart(m); if (ch) { pts.push({ x: ch.x, y: ch.y, p: ch.p }); st.veinPts++; } });
        var e = place(cf, a); if (e.ok) pts.push({ x: e.x, y: e.y, p: e.p, end: true }); else st.obstructedEnd++;
        pts.push({ x: b.x, y: b.y, p: b.p, capture: true }); paths.push({ pts: pts, comp: s.i - 1, side: s.side, angle: s.angle, slot: b.slot }); st.edges++; }); });
    return { paths: paths, stats: st }; }
  return { digits: digits, tau: tau, ringAngles: ringAngles, treeAngles: treeAngles, sector: sector, place: place, buildT: buildT, buildV: buildV, innerAngles: innerAngles };
})();
if (typeof module !== 'undefined') module.exports = WEBS;
