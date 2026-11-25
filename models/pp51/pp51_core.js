// PP51: the PrePer_{5,1} map from a site-chart point (a, b), built exactly as E96 PM.ppMap does for a family (k = 5, n = 1):
// f(z) = 1 + b/z + a/z^2 (PM coefficients A = 1, B = b, C = a, D = 1, E = F = 0), omega' = -2a/b, w_0 = 0, w_1 = inf, w_2 = 1,
// w_j = f(w_(j-1)); tail = w_0 .. w_4; marked = the landing fixed point w_5; scale as PM.mapFrom (1 for a single marked point).
// Chart of the conformal model (PP51 ruling): y = 1/(s - 1), x = s (s + 1)^2 a y^4, s = a + b; so s = 1 + 1/y, a = x/(y^4 s (s + 1)^2), b = s - a.
var PM = (typeof module !== 'undefined') ? require('./pern_math.js') : window.PM;
var PP51 = (function () {
  var K = 5, N = 1;
  function mapAB(a, b) {
    var M = { n: N, k: K, A: [1, 0], B: b, C: a, D: [1, 0], E: [0, 0], F: [0, 0], a: a, b: b };
    M.omega = PM.cabs(b) === 0 ? null : PM.cdiv(PM.cmul([-2, 0], a), b);
    var w = [[0, 0], null, [1, 0]]; for (var j = 3; j < K + N; j++) w.push(PM.fz(M, w[j - 1]));
    M.tail = w.slice(0, K); M.marked = w.slice(K, K + N);
    M.scale = M.marked.map(function (p, j) { var m = 1; for (var i = 0; i < M.marked.length; i++) if (i !== j) m = Math.min(m, PM.chord(p, M.marked[i])); return m; });
    return M;
  }
  function abFromYX(y, x) {
    var s = PM.cadd([1, 0], PM.cdiv([1, 0], y)), s1 = PM.cadd(s, [1, 0]), y2 = PM.cmul(y, y);
    var den = PM.cmul(PM.cmul(y2, y2), PM.cmul(s, PM.cmul(s1, s1))), a = PM.cdiv(x, den);
    return [a, PM.csub(s, a)];
  }
  function mapYX(y, x) { var ab = abFromYX(y, x); return mapAB(ab[0], ab[1]); }
  return { mapAB: mapAB, mapYX: mapYX, abFromYX: abFromYX };
})();
if (typeof module !== 'undefined') module.exports = PP51;
