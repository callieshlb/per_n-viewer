// PP62: the PrePer_{6,2} map from a site-chart point (a, b), built exactly as E96 PM.ppMap does for a family (k = 6, n = 2) (copied
// from PP24 pp24_core.js): f(z) = 1 + b/z + a/z^2 (PM coefficients A = 1, B = b, C = a, D = 1, E = F = 0), omega' = -2a/b, w_0 = 0,
// w_1 = inf, w_2 = 1, w_j = f(w_(j-1)); tail = w_0 .. w_5; marked = the landing cycle w_6, w_7; scale as PM.mapFrom.
// Chart of the conformal model (PP62 REGISTRATION, ruling R2): y = 1/(s - 3/2), s = a + b, x = 64 X y^6 with X = s (s+1)^2 (s^2+4s+1) a,
// so s = 3/2 + 1/y, X = x / (64 y^6), a = X / (s (s+1)^2 (s^2+4s+1)), b = s - a.
var PM = (typeof module !== 'undefined') ? require('./pern_math.js') : window.PM;
var PP62 = (function () {
  function mapAB(a, b) {
    var M = { n: 2, k: 6, A: [1, 0], B: b, C: a, D: [1, 0], E: [0, 0], F: [0, 0], a: a, b: b };
    M.omega = PM.cabs(b) === 0 ? null : PM.cdiv(PM.cmul([-2, 0], a), b);
    var w = [[0, 0], null, [1, 0]]; for (var j = 3; j < 8; j++) w.push(PM.fz(M, w[j - 1]));
    M.tail = w.slice(0, 6); M.marked = w.slice(6, 8);
    M.scale = M.marked.map(function (p, j) { var m = 1; for (var i = 0; i < M.marked.length; i++) if (i !== j) m = Math.min(m, PM.chord(p, M.marked[i])); return m; });
    return M;
  }
  function abFromYX(y, x) {
    var s = PM.cadd([1.5, 0], PM.cdiv([1, 0], y));
    var y2 = PM.cmul(y, y), y6 = PM.cmul(PM.cmul(y2, y2), y2);
    var X = PM.cdiv(x, PM.cmul([64, 0], y6));
    var s1 = PM.cadd(s, [1, 0]), q = PM.cadd(PM.cadd(PM.cmul(s, s), PM.cmul([4, 0], s)), [1, 0]);
    var Ms = PM.cmul(PM.cmul(s, PM.cmul(s1, s1)), q);
    var a = PM.cdiv(X, Ms); return [a, PM.csub(s, a)];
  }
  function mapYX(y, x) { var ab = abFromYX(y, x); return mapAB(ab[0], ab[1]); }
  return { mapAB: mapAB, mapYX: mapYX, abFromYX: abFromYX };
})();
if (typeof module !== 'undefined') module.exports = PP62;
