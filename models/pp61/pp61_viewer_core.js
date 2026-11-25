// PP61 viewer core (adapted from PP24 pp24_viewer_core.js): the CORE interface that E93's live page expects (init, mapData, fate, newton,
// coeffs, cabs, csub), backed by E96 pern_math.js (PM.fate) and ../js/pp61_core.js (PP61.mapYX).
// Chart (projection m): y = 1/(s - 3/2), s = a + b, x' = 64 s (s + 1)^4 a y^6 (P_m monic of degree 8 in x', coefficients up to 1.5e40).
// Ruling 2026-11-25 (x exact from y, chart x selects the sheet): for degree 8 the exact root is computed by Newton's method on
// g(a) = G(a, s - a), evaluated term by term in double-double arithmetic (double alone cannot resolve the root near the places
// at [-1:1:0], where the terms cancel to 1e-16 of their size), started from the chart's a = x'/(64 s (s + 1)^4 y^6)
// (registered: "Newton iterations started at the chart x"). Global root finding in double (Aberth on the expanded g or on P_s) was
// tried on the 31 polynomial points and rejected: cancellation cost up to 30% accuracy. aberth/cofs are kept for diagnostics only.
var CORE = (function () {
  var cm = PM.cmul, ca = PM.cadd, cs = PM.csub, cdv = PM.cdiv, cab = PM.cabs;
  var GT = [[1, 16, 0], [16, 15, 1], [16, 15, 0], [124, 14, 2], [252, 14, 1], [128, 14, 0], [618, 13, 3], [1910, 13, 2], [1960, 13, 1], [668, 13, 0], [2216, 12, 4], [9235, 12, 3], [14325, 12, 2], [9810, 12, 1], [2504, 12, 0], [6056, 11, 5], [31812, 11, 4], [66116, 11, 3], [68036, 11, 2], [34692, 11, 1], [7016, 11, 0], [13028, 10, 6], [82540, 10, 5], [214788, 10, 4], [294249, 10, 3], [224020, 10, 2], [89913, 10, 1], [14866, 10, 0], [22460, 9, 7], [166256, 9, 6], [518206, 9, 5], [882936, 9, 4], [888982, 9, 3], [529202, 9, 2], [172488, 9, 1], [23742, 9, 0], [31285, 8, 8], [264050, 8, 7], [954950, 8, 6], [1935812, 8, 5], [2408175, 8, 4], [1883649, 8, 3], [904838, 8, 2], [243999, 8, 1], [28262, 8, 0], [35232, 7, 9], [332492, 7, 8], [1361844, 7, 7], [3182272, 7, 6], [4680242, 7, 5], [4495536, 7, 4], [2820720, 7, 3], [1114584, 7, 2], [251514, 7, 1], [24668, 7, 0], [31880, 6, 10], [331036, 6, 9], [1506406, 6, 8], [3962037, 6, 7], [6677274, 6, 6], [7539583, 6, 5], [5777685, 6, 4], [2966404, 6, 3], [975927, 6, 2], [185576, 6, 1], [15464, 6, 0], [22866, 5, 11], [257770, 5, 10], [1283076, 5, 9], [3727970, 5, 8], [7033174, 5, 7], [9052924, 5, 6], [8114632, 5, 5], [5064104, 5, 4], [2154932, 5, 3], [594842, 5, 2], [95720, 5, 1], [6790, 5, 0], [12706, 4, 12], [153759, 4, 11], [826483, 4, 10], [2613171, 4, 9], [5419393, 4, 8], [7772222, 4, 7], [7906102, 4, 6], [5746428, 4, 5], [2959989, 4, 4], [1052635, 4, 3], [244949, 4, 2], [33425, 4, 1], [2018, 4, 0], [5272, 3, 13], [67840, 3, 12], [389600, 3, 11], [1324144, 3, 10], [2975254, 3, 9], [4670844, 3, 8], [5272446, 3, 7], [4331940, 3, 6], [2588858, 3, 5], [1110164, 3, 4], [331482, 3, 3], [65116, 3, 2], [7520, 3, 1], [384, 3, 0], [1536, 2, 14], [20852, 2, 13], [126814, 2, 12], [458642, 2, 11], [1103464, 2, 10], [1869846, 2, 9], [2302025, 2, 8], [2091236, 2, 7], [1407514, 2, 6], [697491, 2, 5], [249985, 2, 4], [62648, 2, 3], [10332, 2, 2], [997, 2, 1], [42, 2, 0], [280, 1, 15], [3984, 1, 14], [25470, 1, 13], [97206, 1, 12], [248016, 1, 11], [448502, 1, 10], [594002, 1, 9], [586494, 1, 8], [434810, 1, 7], [241600, 1, 6], [99504, 1, 5], [29698, 1, 4], [6174, 1, 3], [834, 1, 2], [64, 1, 1], [2, 1, 0], [24, 0, 16], [356, 0, 15], [2378, 0, 14], [9511, 0, 13], [25530, 0, 12], [48811, 0, 11], [68772, 0, 10], [72797, 0, 9], [58422, 0, 8], [35572, 0, 7], [16308, 0, 6], [5531, 0, 5], [1344, 0, 4], [221, 0, 3], [22, 0, 2], [1, 0, 1]];    // [coefficient, i, j]: c a^i b^j (pp61_poly.txt)
  var PS = [[0, 0, 0, 0, 0, 1, 50, 1215, 19124, 219248, 1951934, 14049733, 84049018, 426331907, 1861514066, 7079033920, 23664779416, 70064287696, 184840009512, 436676649988, 927596549912, 1777626342378, 3081599762388, 4842836770194, 6910868423880, 8966093000880, 10584273852180, 11373255445830, 11124653941980, 9901744197030, 8013736958916, 5890498119780, 3926133009520, 2367904738544, 1288845378408, 631031071124, 276800006168, 108245329253, 37511291962, 11434051591, 3037465700, 695069216, 134962846, 21796957, 2848378, 289287, 21418, 1028, 24], [0, 0, 0, 1, 44, 927, 12446, 119497, 872064, 5013082, 23167582, 86758554, 261324718, 609552578, 954392958, 181072810, -5118819522, -22330782594, -63477325410, -142714338900, -269359832610, -438356188560, -623994781470, -783530169300, -872308596330, -863616532470, -761457203550, -598079459850, -418139869974, -259754640066, -142979646038, -69465028954, -29634799878, -11027642326, -3548764598, -976702653, -226577186, -43446287, -6701352, -798917, -69070, -3852, -104, 0, 0, 0, 0, 0, 0], [0, 0, 0, -8, -286, -4874, -52655, -404257, -2341392, -10582668, -37982344, -108498034, -241660710, -387435290, -293329925, 638573265, 3331033880, 8651409980, 16699321350, 26191586260, 34566136270, 39040151010, 38068937255, 32189218865, 23633863440, 15054470140, 8295674860, 3935157390, 1595433182, 547366274, 156766381, 36784575, 6886568, 988948, 102262, 6776, 216, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 2, 82, 1462, 15561, 112456, 589340, 2316966, 6900169, 15209962, 22133664, 8039280, -65078686, -239397264, -525657952, -869008244, -1154728230, -1267115434, -1162682086, -896659818, -581373531, -315804200, -142681356, -53002514, -15913987, -3765954, -676028, -86492, -7024, -272, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -7, -175, -1955, -12856, -54295, -144777, -178394, 352616, 2566482, 7761330, 16089086, 25322648, 31451498, 31355382, 25249920, 16411000, 8548853, 3520653, 1120669, 265920, 44269, 4611, 226, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, -14, -114, -834, -5242, -24750, -84826, -213664, -403128, -577984, -634502, -533490, -340926, -162694, -56142, -13228, -1903, -126, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 16, 117, 539, 1763, 4235, 7438, 9398, 8364, 5083, 1999, 457, 46, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-2, -15, -48, -87, -92, -49, -10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]];     // P_s(x, s) = sum_k x^k sum_j PS[k][j] s^j: monic of degree 8 in x, x = s (s + 1)^4 a (s-projection, pp61_P_s.txt)
  function cofs(s) { return PS.map(function (row) { var acc = [0, 0]; for (var j = row.length - 1; j >= 0; j--) acc = ca(cm(acc, s), [row[j], 0]); return acc; }); }
  function hornerD(c, z) { var p = [0, 0], dp = [0, 0]; for (var k = c.length - 1; k >= 0; k--) { dp = ca(cm(dp, z), p); p = ca(cm(p, z), c[k]); } return [p, dp]; }
  function aberth(c) {     // all roots of the monic polynomial with coefficients c[0..n] (c[n] = 1), Aberth-Ehrlich iteration
    var n = c.length - 1, R = 0; for (var k = 0; k < n; k++) R = Math.max(R, Math.pow(cab(c[k]), 1/(n - k)));
    var z = []; for (var k = 0; k < n; k++) z.push([R*Math.cos(2*Math.PI*k/n + 0.4), R*Math.sin(2*Math.PI*k/n + 0.4)]);
    for (var it = 0; it < 800; it++) { var mx = 0;
      for (var k = 0; k < n; k++) { var h = hornerD(c, z[k]); if (cab(h[1]) === 0) continue; var r = cdv(h[0], h[1]), sm = [0, 0];
        for (var j = 0; j < n; j++) if (j !== k) sm = ca(sm, cdv([1, 0], cs(z[k], z[j])));
        var w = cdv(r, cs([1, 0], cm(r, sm))); z[k] = cs(z[k], w); mx = Math.max(mx, cab(w)/(1 + cab(z[k]))); }
      if (mx < 1e-16) break; }
    return z;
  }
  function gdirect(a, s) {   // g(a) = G(a, s - a) and dg/da = G_a - G_b, evaluated term by term (well conditioned)
    var b = cs(s, a), g = [0, 0], dg = [0, 0], ap = [[1, 0]], bp = [[1, 0]];
    for (var k = 1; k <= 16; k++) { ap.push(cm(ap[k - 1], a)); bp.push(cm(bp[k - 1], b)); }
    GT.forEach(function (t) { var c = t[0], i = t[1], j = t[2]; g = ca(g, cm([c, 0], cm(ap[i], bp[j])));
      if (i > 0) dg = ca(dg, cm([c*i, 0], cm(ap[i - 1], bp[j]))); if (j > 0) dg = cs(dg, cm([c*j, 0], cm(ap[i], bp[j - 1]))); });
    return [g, dg];
  }
  function aberthG(z, s) {  // Aberth-Ehrlich on g(a) = G(a, s - a) (degree 8) from the starting values z, term-by-term evaluation
    z = z.slice(); var n = z.length;
    for (var it = 0; it < 400; it++) { var mx = 0;
      for (var k = 0; k < n; k++) { var h = gdirect(z[k], s); if (cab(h[1]) === 0) continue; var r = cdv(h[0], h[1]), sm = [0, 0];
        for (var j = 0; j < n; j++) if (j !== k) sm = ca(sm, cdv([1, 0], cs(z[k], z[j])));
        var w = cdv(r, cs([1, 0], cm(r, sm))); z[k] = cs(z[k], w); mx = Math.max(mx, cab(w)/(1e-300 + cab(z[k]))); }
      if (mx < 1e-15) break; }
    return z;
  }
  // ---- double-double complex arithmetic (Dekker/Knuth; ~32 significant digits) for the Newton step on g(a) = G(a, s - a) ----
  var SPL = 134217729;
  function twoSum(a, b) { var s = a + b, bb = s - a; return [s, (a - (s - bb)) + (b - bb)]; }
  function split(a) { var t = SPL*a, hi = t - (t - a); return [hi, a - hi]; }
  function twoProd(a, b) { var p = a*b, A = split(a), B = split(b); return [p, ((A[0]*B[0] - p) + A[0]*B[1] + A[1]*B[0]) + A[1]*B[1]]; }
  function ddAdd(x, y) { var s = twoSum(x[0], y[0]), e = s[1] + x[1] + y[1], r = twoSum(s[0], e); return r; }
  function ddMul(x, y) { var p = twoProd(x[0], y[0]); p[1] += x[0]*y[1] + x[1]*y[0]; return twoSum(p[0], p[1]); }
  function ddNeg(x) { return [-x[0], -x[1]]; }
  function ddDiv(x, y) { var q1 = x[0]/y[0], r = ddAdd(x, ddNeg(ddMul(y, [q1, 0]))), q2 = r[0]/y[0]; r = ddAdd(r, ddNeg(ddMul(y, [q2, 0]))); var q3 = r[0]/y[0]; return ddAdd(ddAdd([q1, 0], [q2, 0]), [q3, 0]); }
  function D(v) { return [v, 0]; }
  function zAdd(z, w) { return [ddAdd(z[0], w[0]), ddAdd(z[1], w[1])]; }
  function zSub(z, w) { return [ddAdd(z[0], ddNeg(w[0])), ddAdd(z[1], ddNeg(w[1]))]; }
  function zMul(z, w) { return [ddAdd(ddMul(z[0], w[0]), ddNeg(ddMul(z[1], w[1]))), ddAdd(ddMul(z[0], w[1]), ddMul(z[1], w[0]))]; }
  function zScale(z, c) { return [ddMul(z[0], D(c)), ddMul(z[1], D(c))]; }
  function zDiv(z, w) { var den = ddAdd(ddMul(w[0], w[0]), ddMul(w[1], w[1])), n = zMul(z, [w[0], ddNeg(w[1])]); return [ddDiv(n[0], den), ddDiv(n[1], den)]; }
  function zOf(c) { return [D(c[0]), D(c[1])]; }
  function zTo(z) { return [z[0][0] + z[0][1], z[1][0] + z[1][1]]; }
  function gdirectDD(a, s) {
    var b = zSub(s, a), g = zOf([0, 0]), dg = zOf([0, 0]), ap = [zOf([1, 0])], bp = [zOf([1, 0])];
    for (var k = 1; k <= 16; k++) { ap.push(zMul(ap[k - 1], a)); bp.push(zMul(bp[k - 1], b)); }
    GT.forEach(function (t) { var c = t[0], i = t[1], j = t[2]; g = zAdd(g, zScale(zMul(ap[i], bp[j]), c));
      if (i > 0) dg = zAdd(dg, zScale(zMul(ap[i - 1], bp[j]), c*i)); if (j > 0) dg = zSub(dg, zScale(zMul(ap[i], bp[j - 1]), c*j)); });
    return [g, dg];
  }
  function sOf(y) { return ca([1.5, 0], cdv([1, 0], y)); }
  function y6(y) { var y2 = cm(y, y); return cm(cm(y2, y2), y2); }
  function scaleOf(y) { var s = sOf(y), s1 = ca(s, [1, 0]), s12 = cm(s1, s1); return cm([64, 0], cm(y6(y), cm(s, cm(s12, s12)))); }
  // 'newton' kept as the name E93's page calls: returns [exact root x' nearest the chart x', 0, distance (in x') to the next root]
  function newton(x, q) {   // [x' of the Newton root from the chart value (double-double), iterations (-1: no convergence), |a moved|/|a|]
    var y = q.y, M = scaleOf(y), s = sOf(y), a0 = cdv(x, M), a = zOf(a0), sz = zOf(s), it;
    for (it = 0; it < 60; it++) { var h = gdirectDD(a, sz), d = zDiv(h[0], h[1]); a = zSub(a, d);
      if (cab(zTo(d)) <= 1e-26*(1 + cab(zTo(a)))) break; }
    if (it >= 60) it = -1;
    var af = zTo(a); return [cm(M, af), it, cab(cs(af, a0))/(1e-300 + cab(a0))];
  }
  function coeffsY(y) { return { y: y }; }
  function Grel(a, b) {
    var s = [0, 0], tot = 0; GT.forEach(function (t) { var v = [t[0], 0]; for (var i = 0; i < t[1]; i++) v = cm(v, a); for (var j = 0; j < t[2]; j++) v = cm(v, b); s = ca(s, v); tot += cab(v); });
    return cab(s)/tot; }
  return { init: function () {}, mapData: function (x, y) { return PP61.mapYX(y, x); }, fate: function (M, z, it, d, tol) { return PM.fate(M, z, it, d, tol); },
           newton: newton, coeffs: coeffsY, Grel: Grel, aberth: aberth, cofs: cofs, cabs: PM.cabs, csub: PM.csub };
})();
if (typeof module !== 'undefined') module.exports = CORE;
