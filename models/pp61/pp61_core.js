// PP61 (adapted from PP24 pp24_core.js): the PrePer_{6,1} map from a site-chart point (a, b), built exactly as E96 PM.ppMap does for
// (k = 6, n = 1): f(z) = 1 + b/z + a/z^2 (A = 1, B = b, C = a, D = 1, E = F = 0), omega' = -2a/b, w_0 = 0, w_1 = inf, w_2 = 1,
// w_j = f(w_(j-1)); tail = w_0 .. w_5; marked = the landing fixed point w_6; scale as PM.mapFrom (one marked point: 1).
// Chart of the conformal model (projection m, PP61 REGISTRATION): y = 1/(s - 3/2), s = a + b, x' = 64 s (s + 1)^4 a y^6, so
// s = 3/2 + 1/y, a = x'/(64 y^6 s (s + 1)^4), b = s - a.
var PM = (typeof module !== 'undefined') ? require('./pern_math.js') : window.PM;
var PP61 = (function () {
  function mapAB(a, b) {
    var M = { n: 1, k: 6, A: [1, 0], B: b, C: a, D: [1, 0], E: [0, 0], F: [0, 0], a: a, b: b };
    M.omega = PM.cabs(b) === 0 ? null : PM.cdiv(PM.cmul([-2, 0], a), b);
    var w = [[0, 0], null, [1, 0]]; for (var j = 3; j < 7; j++) w.push(PM.fz(M, w[j - 1]));
    M.tail = w.slice(0, 6); M.marked = w.slice(6, 7);
    M.scale = M.marked.map(function (p, j) { var m = 1; for (var i = 0; i < M.marked.length; i++) if (i !== j) m = Math.min(m, PM.chord(p, M.marked[i])); return m; });
    return M;
  }
  function abFromYX(y, x) {
    var s = PM.cadd([1.5, 0], PM.cdiv([1, 0], y)), y2 = PM.cmul(y, y), y6 = PM.cmul(PM.cmul(y2, y2), y2);
    var s1 = PM.cadd(s, [1, 0]), s12 = PM.cmul(s1, s1), M4 = PM.cmul(s, PM.cmul(s12, s12));
    var a = PM.cdiv(x, PM.cmul([64, 0], PM.cmul(y6, M4))); return [a, PM.csub(s, a)];
  }
  function mapYX(y, x) { var ab = abFromYX(y, x); return mapAB(ab[0], ab[1]); }
  return { mapAB: mapAB, mapYX: mapYX, abFromYX: abFromYX };
})();
if (typeof module !== 'undefined') module.exports = PP61;
