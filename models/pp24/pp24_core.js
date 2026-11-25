// PP24: the PrePer_{2,4} map from a site-chart point (a, b), built exactly as E96 PM.ppMap does for a family (k = 2, n = 4):
// f(z) = 1 + b/z + a/z^2 (PM coefficients A = 1, B = b, C = a, D = 1, E = F = 0), omega' = -2a/b, w_0 = 0, w_1 = inf, w_2 = 1,
// w_j = f(w_(j-1)); tail = w_0, w_1; marked = the landing cycle w_2 .. w_5; scale as PM.mapFrom. Dynamics: PM.fz, PM.fate (E96).
// Chart of the conformal model: y = s = a + b, x = (s + 1)^2 a (PP24 REGISTRATION), so a = x/(y + 1)^2, b = y - a.
var PM = (typeof module !== 'undefined') ? require('./pern_math.js') : window.PM;
var PP24 = (function () {
  function mapAB(a, b) {
    var M = { n: 4, k: 2, A: [1, 0], B: b, C: a, D: [1, 0], E: [0, 0], F: [0, 0], a: a, b: b };
    M.omega = PM.cabs(b) === 0 ? null : PM.cdiv(PM.cmul([-2, 0], a), b);
    var w = [[0, 0], null, [1, 0]]; for (var j = 3; j < 6; j++) w.push(PM.fz(M, w[j - 1]));
    M.tail = w.slice(0, 2); M.marked = w.slice(2, 6);
    M.scale = M.marked.map(function (p, j) { var m = 1; for (var i = 0; i < M.marked.length; i++) if (i !== j) m = Math.min(m, PM.chord(p, M.marked[i])); return m; });
    return M;
  }
  function abFromYX(y, x) { var e = [y[0] + 1, y[1]], a = PM.cdiv(x, PM.cmul(e, e)); return [a, PM.csub(y, a)]; }
  function mapYX(y, x) { var ab = abFromYX(y, x); return mapAB(ab[0], ab[1]); }
  return { mapAB: mapAB, mapYX: mapYX, abFromYX: abFromYX };
})();
if (typeof module !== 'undefined') module.exports = PP24;
