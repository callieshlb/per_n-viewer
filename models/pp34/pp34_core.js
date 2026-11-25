// PP34 (adapted from PP24 pp24_core.js): the PrePer_{3,4} map from a site-chart point (a, b), built exactly as E96 PM.ppMap (k = 3, n = 4):
// f(z) = 1 + b/z + a/z^2 (A = 1, B = b, C = a, D = 1, E = F = 0), omega' = -2a/b, w_0 = 0, w_1 = inf, w_2 = 1, w_j = f(w_(j-1));
// tail = w_0 .. w_2; marked = the landing cycle w_3 .. w_6; scale as PM.mapFrom. Dynamics: PM.fz, PM.fate (E96).
// Chart of the conformal model: y = s = a + b, x = L a with L = s(2s^2 + 4s + 1) (PP34 REGISTRATION), so a = x/L, b = y - a.
var PM = (typeof module !== 'undefined') ? require('./pern_math.js') : window.PM;
var PP34 = (function () {
  var K = 3, N = 4;
  function mapAB(a, b) {
    var M = { n: N, k: K, A: [1, 0], B: b, C: a, D: [1, 0], E: [0, 0], F: [0, 0], a: a, b: b };
    M.omega = PM.cabs(b) === 0 ? null : PM.cdiv(PM.cmul([-2, 0], a), b);
    var w = [[0, 0], null, [1, 0]]; for (var j = 3; j < K + N; j++) w.push(PM.fz(M, w[j - 1]));
    M.tail = w.slice(0, K); M.marked = w.slice(K, K + N);
    M.scale = M.marked.map(function (p, j) { var m = 1; for (var i = 0; i < M.marked.length; i++) if (i !== j) m = Math.min(m, PM.chord(p, M.marked[i])); return m; });
    return M;
  }
  function Lof(y) { var y2 = PM.cmul(y, y); return PM.cmul(y, [2*y2[0] + 4*y[0] + 1, 2*y2[1] + 4*y[1]]); }
  function abFromYX(y, x) { var a = PM.cdiv(x, Lof(y)); return [a, PM.csub(y, a)]; }
  function mapYX(y, x) { var ab = abFromYX(y, x); return mapAB(ab[0], ab[1]); }
  return { mapAB: mapAB, mapYX: mapYX, abFromYX: abFromYX, K: K, N: N };
})();
if (typeof module !== 'undefined') module.exports = PP34;
