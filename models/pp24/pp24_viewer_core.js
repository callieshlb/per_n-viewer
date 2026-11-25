// PP24 viewer core: the CORE interface that E93's live page expects (init, mapData, fate, newton, coeffs, cabs, csub), backed by
// E96 pern_math.js (PM.fate) and ../js/pp24_core.js (PP24.mapYX). Curve in the chart: P(x, y) = x^2 - B(y) x + C(y),
// B = y(y + 2)(2y + 1), C = (y + 1)^3 (2y + 1)(y^2 + 3y + 1); y = s = a + b, x = (s + 1)^2 a.
// Ruling 2026-11-25: x is computed EXACTLY as the root of the monic P(., y) nearest the chart's x (the chart x only selects the sheet).
var CORE = (function () {
  var cm = PM.cmul, ca = PM.cadd, cs = PM.csub, cdv = PM.cdiv;
  function csqrt(z) { var r = Math.hypot(z[0], z[1]), t = Math.atan2(z[1], z[0])/2, s = Math.sqrt(r); return [s*Math.cos(t), s*Math.sin(t)]; }
  function coeffs(y) { var y1 = ca(y, [1, 0]), y2 = ca(cm([2, 0], y), [1, 0]), B = cm(cm(y, ca(y, [2, 0])), y2),
      C = cm(cm(cm(cm(y1, y1), y1), y2), ca(ca(cm(y, y), cm([3, 0], y)), [1, 0])); return { B: B, C: C }; }
  function roots(y) { var q = coeffs(y), d = csqrt(cs(cm(q.B, q.B), cm([4, 0], q.C))); return [cm([0.5, 0], ca(q.B, d)), cm([0.5, 0], cs(q.B, d))]; }
  // 'newton' kept as the name E93's page calls: returns [exact root nearest x, 0 iterations, distance to the other root]
  function newton(x, q) { var r = roots(q.y); var d0 = PM.cabs(cs(r[0], x)), d1 = PM.cabs(cs(r[1], x)); return d0 <= d1 ? [r[0], 0, d1] : [r[1], 0, d0]; }
  function coeffsY(y) { return { y: y }; }
  function Grel(a, b) {          // |G(a, b)| / sum |terms| of the PrePer_{2,4} quartic (double)
    var T = [[1, 4, 0], [4, 3, 1], [6, 3, 0], [7, 2, 2], [19, 2, 1], [11, 2, 0], [6, 1, 3], [22, 1, 2], [22, 1, 1], [6, 1, 0], [2, 0, 4], [9, 0, 3], [12, 0, 2], [6, 0, 1], [1, 0, 0]];
    var s = [0, 0], tot = 0; T.forEach(function (t) { var v = [t[0], 0]; for (var i = 0; i < t[1]; i++) v = cm(v, a); for (var j = 0; j < t[2]; j++) v = cm(v, b); s = ca(s, v); tot += PM.cabs(v); });
    return PM.cabs(s)/tot; }
  return { init: function () {}, mapData: function (x, y) { return PP24.mapYX(y, x); }, fate: function (M, z, it, d, tol) { return PM.fate(M, z, it, d, tol); },
           newton: newton, coeffs: coeffsY, roots: roots, Grel: Grel, cabs: PM.cabs, csub: PM.csub };
})();
