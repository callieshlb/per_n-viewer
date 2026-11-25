// PP62 viewer core (adapted from PP24 pp24_viewer_core.js): the CORE interface that E93's live page expects (init, mapData, fate, newton,
// coeffs, roots, Grel, cabs, csub), backed by E96 pern_math.js (PM.fate) and ../js/pp62_core.js (PP62.mapYX).
// Chart (ruling R2): y = 1/(s - 3/2), s = a + b, x = 64 y^6 M(s) a, M(s) = s (s+1)^2 (s^2+4s+1); P_m(x, y) is the monic quintic in x.
// Ruling R1 (x exact): the root of P_m(., y) nearest the chart x is computed through the SAME root in a: the quintic G(a, s - a) = sum c_k(s) a^k
// (P_m(64 y^6 M(s) a, y) = G(a, s - a) times a nonzero function of y at non-special y; c_k have coefficients <= 2618, P_m up to 2.8e25),
// all 5 roots by Aberth iteration in double, nearest to the chart a, mapped back to x.
var CORE = (function () {
  var cm = PM.cmul, ca = PM.cadd, cs = PM.csub, cdv = PM.cdiv, cab = PM.cabs;
  var Q = [[1, 16, 112, 452, 1166, 2013, 2364, 1878, 981, 317, 56, 4], [0, -8, -96, -490, -1391, -2408, -2618, -1770, -707, -148, -12], [0, 4, 58, 319, 899, 1433, 1319, 676, 172, 16], [0, -2, -29, -150, -377, -500, -344, -110, -12], [0, 0, 7, 41, 90, 90, 38, 5], [0, 0, -1, -6, -10, -6, -1]];          // c_k(s) ascending in s, k = 0..5 (viewer_quintic.json)
  var GT = [[1, 10, 1], [2, 10, 0], [10, 9, 2], [32, 9, 1], [22, 9, 0], [47, 8, 3], [208, 8, 2], [275, 8, 1], [114, 8, 0], [136, 7, 4], [770, 7, 3], [1492, 7, 2], [1214, 7, 1], [356, 7, 0], [267, 6, 5], [1842, 6, 4], [4676, 6, 3], [5605, 6, 2], [3216, 6, 1], [712, 6, 0], [370, 5, 6], [3008, 5, 5], [9392, 5, 4], [14720, 5, 3], [12358, 5, 2], [5304, 5, 1], [912, 5, 0], [365, 4, 7], [3412, 4, 6], [12572, 4, 5], [24095, 4, 4], [26233, 4, 3], [16331, 4, 2], [5400, 4, 1], [732, 4, 0], [252, 3, 8], [2658, 3, 7], [11220, 3, 6], [25180, 3, 5], [33248, 3, 4], [26646, 3, 3], [12712, 3, 2], [3308, 3, 1], [360, 3, 0], [116, 2, 9], [1360, 2, 8], [6432, 2, 7], [16397, 2, 6], [25163, 2, 5], [24319, 2, 4], [14885, 2, 3], [5584, 2, 2], [1168, 2, 1], [104, 2, 0], [32, 1, 10], [412, 1, 9], [2146, 1, 8], [6078, 1, 7], [10528, 1, 6], [11776, 1, 5], [8674, 1, 4], [4174, 1, 3], [1260, 1, 2], [216, 1, 1], [16, 1, 0], [4, 0, 11], [56, 0, 10], [317, 0, 9], [981, 0, 8], [1878, 0, 7], [2364, 0, 6], [2013, 0, 5], [1166, 0, 4], [452, 0, 3], [112, 0, 2], [16, 0, 1], [1, 0, 0]];        // terms [coef, i, j] of G(a, b) (G_terms.json)
  function polyS(co, s) { var v = [0, 0]; for (var i = co.length - 1; i >= 0; i--) v = ca(cm(v, s), [co[i], 0]); return v; }
  function sOf(y) { return ca([1.5, 0], cdv([1, 0], y)); }
  function Mof(s) { var s1 = ca(s, [1, 0]), q = ca(ca(cm(s, s), cm([4, 0], s)), [1, 0]); return cm(cm(s, cm(s1, s1)), q); }
  function scale(y) { var y2 = cm(y, y), y6 = cm(cm(y2, y2), y2); return cm([64, 0], cm(y6, Mof(sOf(y)))); }   // x = scale(y) a
  function coeffs(y) { var s = sOf(y); return { y: y, s: s, c: Q.map(function (co) { return polyS(co, s); }), k: scale(y) }; }
  function pev(c, z) { var p = [0, 0], d = [0, 0]; for (var i = c.length - 1; i >= 0; i--) { d = ca(cm(d, z), p); p = ca(cm(p, z), c[i]); } return [p, d]; }
  function rootsA(q) {          // Aberth on the quintic in a (double)
    var n = 5, c = q.c, lead = c[5], r = [];
    var R = 0; for (var k = 0; k < 5; k++) R = Math.max(R, Math.pow(cab(cdv(c[k], lead)), 1/(5 - k)));
    for (var k = 0; k < n; k++) r.push([2*R*Math.cos(2*Math.PI*k/n + 0.4), 2*R*Math.sin(2*Math.PI*k/n + 0.4)]);
    for (var it = 0; it < 500; it++) { var mx = 0;
      for (var k = 0; k < n; k++) { var pd = pev(c, r[k]); if (cab(pd[1]) === 0) continue; var w = cdv(pd[0], pd[1]), sum = [0, 0];
        for (var j = 0; j < n; j++) if (j !== k) sum = ca(sum, cdv([1, 0], cs(r[k], r[j])));
        var den = cs([1, 0], cm(w, sum)), dz = cdv(w, den); r[k] = cs(r[k], dz); mx = Math.max(mx, cab(dz)/(1 + cab(r[k]))); }
      if (mx < 1e-15) break; }
    return r; }
  function roots(y) { var q = coeffs(y); return rootsA(q).map(function (a) { return cm(q.k, a); }); }
  // 'newton' kept as the name E93's page calls: returns [exact root nearest x, 0, distance to the next root]
  function newton(x, q) { var r = rootsA(q).map(function (a) { return cm(q.k, a); });
    var d = r.map(function (z) { return cab(cs(z, x)); }), o = d.map(function (v, i) { return i; }).sort(function (i, j) { return d[i] - d[j]; });
    return [r[o[0]], 0, cab(cs(r[o[1]], r[o[0]]))]; }
  function Grel(a, b) { var s = [0, 0], tot = 0; GT.forEach(function (t) { var v = [t[0], 0]; for (var i = 0; i < t[1]; i++) v = cm(v, a); for (var j = 0; j < t[2]; j++) v = cm(v, b); s = ca(s, v); tot += cab(v); }); return cab(s)/tot; }
  return { init: function () {}, mapData: function (x, y) { return PP62.mapYX(y, x); }, fate: function (M, z, it, d, tol) { return PM.fate(M, z, it, d, tol); },
           newton: newton, coeffs: coeffs, roots: roots, Grel: Grel, cabs: PM.cabs, csub: PM.csub, Q: Q };
})();
