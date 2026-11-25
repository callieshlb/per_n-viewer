// PP43 viewer core (adapted from PP24 pp24_viewer_core.js): the CORE interface that E93's live page expects (init, mapData, fate, newton,
// coeffs, cabs, csub), backed by E96 pern_math.js (PM.fate) and ../js/pp43_core.js (PP43.mapYX). Curve in the chart: the monic quartic
// P(x, y) = x^4 + c3 x^3 + c2 x^2 + c1 x + c0 with c3 = -y(2y+1)(3y+5), c2 = y(y+1)^2(y+2)(2y+1)(7y+1), c1 = -y^2(y+1)^4(2y+1)(7y^2+19y+6),
// c0 = y^2(y+1)^6(2y+1)(3y+1)(y^2+3y+1); y = s = a + b, x = s(s+1)a.
// Ruling 2026-11-25: x is computed EXACTLY as the root of the monic P(., y) nearest the chart's x (the chart x only selects the sheet):
// all four roots by Durand-Kerner (double), each polished by Newton.
var CORE = (function () {
  var cm = PM.cmul, ca = PM.cadd, cs = PM.csub, cdv = PM.cdiv;
  function sc(k, z) { return [k*z[0], k*z[1]]; }
  function coeffs4(y) {
    var y1 = ca(y, [1, 0]), y2 = ca(sc(2, y), [1, 0]), yy = cm(y, y), y11 = cm(y1, y1), y14 = cm(y11, y11);
    var c3 = sc(-1, cm(cm(y, y2), ca(sc(3, y), [5, 0])));
    var c2 = cm(cm(cm(cm(y, y11), ca(y, [2, 0])), y2), ca(sc(7, y), [1, 0]));
    var c1 = sc(-1, cm(cm(cm(yy, y14), y2), ca(ca(sc(7, yy), sc(19, y)), [6, 0])));
    var c0 = cm(cm(cm(cm(yy, cm(y14, y11)), y2), ca(sc(3, y), [1, 0])), ca(ca(yy, sc(3, y)), [1, 0]));
    return [c0, c1, c2, c3, [1, 0]];
  }
  function pv(c, x) { var p = [0, 0], d = [0, 0]; for (var k = c.length - 1; k >= 0; k--) { d = ca(cm(d, x), p); p = ca(cm(p, x), c[k]); } return [p, d]; }
  function roots(y) {
    var c = coeffs4(y), R = 1; for (var k = 0; k < 4; k++) R = Math.max(R, Math.pow(PM.cabs(c[k]), 1/(4 - k)));
    var z = [0, 1, 2, 3].map(function (k) { var t = 0.4 + 2*Math.PI*k/4; return [1.3*R*Math.cos(t), 1.3*R*Math.sin(t)]; });
    for (var it = 0; it < 500; it++) {
      var mx = 0;
      for (var i = 0; i < 4; i++) { var den = [1, 0]; for (var j = 0; j < 4; j++) if (j !== i) den = cm(den, cs(z[i], z[j]));
        var st = cdv(pv(c, z[i])[0], den); z[i] = cs(z[i], st); mx = Math.max(mx, PM.cabs(st)/(1 + PM.cabs(z[i]))); }
      if (mx < 1e-16) break;
    }
    return z.map(function (x) { for (var k = 0; k < 3; k++) { var q = pv(c, x); if (PM.cabs(q[1]) === 0) break; x = cs(x, cdv(q[0], q[1])); } return x; });
  }
  // 'newton' kept as the name E93's page calls: returns [exact root nearest x, 0 iterations, distance from x to the next-nearest root]
  function newton(x, q) { var r = roots(q.y).map(function (z) { return [PM.cabs(cs(z, x)), z]; }).sort(function (u, v) { return u[0] - v[0]; });
    return [r[0][1], 0, r[1][0]]; }
  function coeffsY(y) { return { y: y }; }
  var T = [[1,7,0],[7,6,1],[7,6,0],[23,5,2],[48,5,1],[25,5,0],[47,4,3],[153,4,2],[163,4,1],[57,4,0],[64,3,4],[281,3,3],[443,3,2],[299,3,1],[73,3,0],[56,2,5],[301,2,4],[606,2,3],[576,2,2],[258,2,1],[43,2,0],[28,1,6],[173,1,5],[410,1,4],[483,1,3],[299,1,2],[92,1,1],[11,1,0],[6,0,7],[41,0,6],[109,0,5],[149,0,4],[114,0,3],[49,0,2],[11,0,1],[1,0,0]];
  function Grel(a, b) {          // |G(a, b)| / sum |terms| of the PrePer_{4,3} septic (double)
    var s = [0, 0], tot = 0; T.forEach(function (t) { var v = [t[0], 0]; for (var i = 0; i < t[1]; i++) v = cm(v, a); for (var j = 0; j < t[2]; j++) v = cm(v, b); s = ca(s, v); tot += PM.cabs(v); });
    return PM.cabs(s)/tot; }
  return { init: function () {}, mapData: function (x, y) { return PP43.mapYX(y, x); }, fate: function (M, z, it, d, tol) { return PM.fate(M, z, it, d, tol); },
           newton: newton, coeffs: coeffsY, roots: roots, Grel: Grel, cabs: PM.cabs, csub: PM.csub };
})();
if (typeof module !== 'undefined') module.exports = CORE;
