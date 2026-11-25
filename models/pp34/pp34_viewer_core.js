// PP34 viewer core (adapted from PP24 pp24_viewer_core.js): the CORE interface E93's live page expects, backed by E96 pern_math.js
// (PM.fate) and ../js/pp34_core.js (PP34.mapYX). Curve in the chart: monic quartic P(x, y) = x^4 + c3 x^3 + c2 x^2 + c1 x + c0 with
// c3 = -2y(y+1)(4y^2+9y+3), c2 = y(y+1)^2(2y^2+4y+1)(15y^3+40y^2+20y+2), c1 = -y^2(y+1)^3(2y+1)(2y^2+4y+1)^2(7y^2+19y+6),
// c0 = y^2(y+1)^4(2y+1)(3y+1)(y^2+3y+1)(2y^2+4y+1)^3; y = s = a + b, x = s(2s^2+4s+1) a.
// Ruling 2026-11-25: x is computed EXACTLY (double) as the root of the monic P(., y) nearest the chart's x (chart x selects the sheet only).
var CORE = (function () {
  var cm = PM.cmul, ca = PM.cadd, cs = PM.csub, cdv = PM.cdiv;
  function pw(z, k) { var r = [1, 0]; for (var i = 0; i < k; i++) r = cm(r, z); return r; }
  function poly(y, cf) { var r = [0, 0]; for (var i = 0; i < cf.length; i++) r = ca(cm(r, y), [cf[i], 0]); return r; }   // cf highest first
  function coeffs(y) {
    var y1 = ca(y, [1, 0]), q = poly(y, [2, 4, 1]);
    var c3 = cm([-2, 0], cm(cm(y, y1), poly(y, [4, 9, 3])));
    var c2 = cm(cm(cm(y, pw(y1, 2)), q), poly(y, [15, 40, 20, 2]));
    var c1 = cm([-1, 0], cm(cm(cm(cm(pw(y, 2), pw(y1, 3)), poly(y, [2, 1])), pw(q, 2)), poly(y, [7, 19, 6])));
    var c0 = cm(cm(cm(cm(cm(cm(pw(y, 2), pw(y1, 4)), poly(y, [2, 1])), poly(y, [3, 1])), poly(y, [1, 3, 1])), pw(q, 3)), [1, 0]);
    return [[1, 0], c3, c2, c1, c0];
  }
  function peval(c, x) { var r = [0, 0]; for (var i = 0; i < c.length; i++) r = ca(cm(r, x), c[i]); return r; }
  function roots(y) {          // Durand-Kerner on the monic quartic, then Newton polish (double)
    var c = coeffs(y), sc = 0; for (var i = 1; i < 5; i++) sc = Math.max(sc, Math.pow(PM.cabs(c[i]), 1/i));
    sc = sc || 1; var z = []; for (var k = 0; k < 4; k++) z.push([sc*Math.cos(0.4 + 2*Math.PI*k/4), sc*Math.sin(0.4 + 2*Math.PI*k/4)]);
    for (var it = 0; it < 500; it++) { var mv = 0;
      for (var k = 0; k < 4; k++) { var den = [1, 0]; for (var j = 0; j < 4; j++) if (j !== k) den = cm(den, cs(z[k], z[j]));
        var d = cdv(peval(c, z[k]), den); z[k] = cs(z[k], d); mv = Math.max(mv, PM.cabs(d)/(1 + PM.cabs(z[k]))); }
      if (mv < 1e-15) break; }
    var dc = [[4, 0], cm([3, 0], c[1]), cm([2, 0], c[2]), c[3]];
    for (var k = 0; k < 4; k++) for (var it = 0; it < 3; it++) { var dv = peval(dc, z[k]); if (PM.cabs(dv) > 0) z[k] = cs(z[k], cdv(peval(c, z[k]), dv)); }
    return z;
  }
  // 'newton' kept as the name E93's page calls: returns [exact root nearest x, 0, distance from x to the second-nearest root]
  function newton(x, q) { var r = roots(q.y).map(function (z) { return [PM.cabs(cs(z, x)), z]; }).sort(function (u, v) { return u[0] - v[0]; }); return [r[0][1], 0, r[1][0]]; }
  function coeffsY(y) { return { y: y }; }
  var T = [[1, 8, 0], [8, 7, 1], [8, 7, 0], [31, 6, 2], [65, 6, 1], [34, 6, 0], [74, 5, 3], [238, 5, 2], [250, 5, 1], [86, 5, 0], [117, 4, 4], [505, 4, 3], [789, 4, 2], [532, 4, 1], [131, 4, 0], [124, 3, 5], [664, 3, 4], [1356, 3, 3], [1332, 3, 2], [632, 3, 1], [116, 3, 0], [85, 2, 6], [535, 2, 5], [1325, 2, 4], [1670, 2, 3], [1134, 2, 2], [393, 2, 1], [54, 2, 0], [34, 1, 7], [242, 1, 6], [692, 1, 5], [1042, 1, 4], [896, 1, 3], [440, 1, 2], [114, 1, 1], [12, 1, 0], [6, 0, 8], [47, 0, 7], [150, 0, 6], [258, 0, 5], [263, 0, 4], [163, 0, 3], [60, 0, 2], [12, 0, 1], [1, 0, 0]];
  function Grel(a, b) {          // |G(a, b)| / sum |terms| of the PrePer_{3,4} octic (double)
    var s = [0, 0], tot = 0; T.forEach(function (t) { var v = cm(cm([t[0], 0], pw(a, t[1])), pw(b, t[2])); s = ca(s, v); tot += PM.cabs(v); });
    return PM.cabs(s)/tot; }
  return { init: function () {}, mapData: function (x, y) { return PP34.mapYX(y, x); }, fate: function (M, z, it, d, tol) { return PM.fate(M, z, it, d, tol); },
           newton: newton, coeffs: coeffsY, roots: roots, Grel: Grel, cabs: PM.cabs, csub: PM.csub };
})();
if (typeof module !== 'undefined') module.exports = CORE;
