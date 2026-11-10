// E96 puncture explorer, part 1 (REGISTRATION addendum 2026-11-10 (f)). Depends on pern_math.js (PM).
'use strict';
var SPH = (function () {
  var ch = PM.chord;
  // M: map; eps: neck threshold. Returns {labels, pts, necks, compatible, spheres: [{dirs: [[idx...]], image, consistent, period, T, Ti}], ...}
  function analyze(M, eps) {
    eps = eps || 1e-2; var n = M.marked.length, pts = M.marked.concat([M.omega, PM.fz(M, M.omega)]), labels = M.marked.map(function (_, j) { return 'a' + j; }).concat(['ω', 'f(ω)']), N = pts.length, i, j, k, l;
    var img = function (q) { return q < n ? (q + 1) % n : q === n ? n + 1 : -1; };                 // a_j -> a_(j+1), omega -> f(omega); f(omega) has no image in X
    for (i = N - 1; i >= n; i--) for (j = 0; j < i; j++) if (ch(pts[i], pts[j]) < 1e-9) { labels[j] += ' = ' + labels[i]; pts.splice(i, 1); labels.splice(i, 1); break; }   // coincident extra points (omega = a_j, f(omega) = omega at centers)
    N = pts.length; var nOm = labels.findIndex(function (t) { return t.split(' = ')[0] === 'ω' || t.indexOf('= ω') >= 0; }), nFo = labels.findIndex(function (t) { return t.split(' = ')[0] === 'f(ω)'; });
    img = function (q) { return q < n ? (q + 1) % n : (q === nOm && nFo >= 0) ? nFo : -1; };
    var d = pts.map(function (p) { return pts.map(function (q) { return ch(p, q); }); });
    var sep = function (i, j, k, l) { var num = d[i][j]*d[k][l]; return Math.max(num/Math.max(d[i][k]*d[j][l], 1e-300), num/Math.max(d[i][l]*d[j][k], 1e-300)); };
    var necks = [];                                                   // side not containing point 1 (a_1), as a bitmask
    for (var mask = 1; mask < (1 << N); mask++) { if (mask & 2) continue; var A = [], B = []; for (i = 0; i < N; i++) (mask >> i & 1 ? A : B).push(i);
      if (A.length < 2 || B.length < 2) continue; var worst = 0;
      for (i = 0; i < A.length && worst < eps; i++) for (j = i + 1; j < A.length && worst < eps; j++) for (k = 0; k < B.length && worst < eps; k++) for (l = k + 1; l < B.length; l++) { worst = Math.max(worst, sep(A[i], A[j], B[k], B[l])); if (worst >= eps) break; }
      if (worst < eps) necks.push({ mask: mask, A: A, worst: worst }); }
    var compatible = true;                                            // laminar: any two are nested or disjoint
    for (i = 0; i < necks.length; i++) for (j = i + 1; j < necks.length; j++) { var x = necks[i].mask, y = necks[j].mask, c = x & y; if (c !== 0 && c !== x && c !== y) compatible = false; }
    if (!compatible) return { labels: labels, pts: pts, necks: necks, compatible: false, spheres: [], eps: eps };
    var full = (1 << N) - 1, nodes = [{ mask: full, parent: -1 }].concat(necks.slice().sort(function (u, v) { return v.A.length - u.A.length; }).map(function (q) { return { mask: q.mask, parent: -1, worst: q.worst }; }));
    for (i = 1; i < nodes.length; i++) { var best = 0; for (j = 1; j < nodes.length; j++) if (j !== i && (nodes[i].mask & nodes[j].mask) === nodes[i].mask && nodes[j].mask !== nodes[i].mask) { if (best === 0 || (nodes[j].mask & nodes[best].mask) === nodes[j].mask) best = j; } nodes[i].parent = best; }
    var bits = function (m) { var o = []; for (var q = 0; q < N; q++) if (m >> q & 1) o.push(q); return o; };
    var spheres = nodes.map(function (nd, idx) { var kids = [], cover = 0; nodes.forEach(function (c, ci) { if (ci !== idx && c.parent === idx && ci !== 0) { kids.push(bits(c.mask)); cover |= c.mask; } });
      var dirs = kids.concat(bits(nd.mask & ~cover).map(function (q) { return [q]; })); if (idx !== 0) dirs.unshift(bits(full & ~nd.mask));      // first direction: toward the root
      return { dirs: dirs, neck: nd.worst }; });
    // median vertex of three leaves: the sphere in which they lie in three different directions
    var dirOf = function (S, q) { for (var t = 0; t < S.dirs.length; t++) if (S.dirs[t].indexOf(q) >= 0) return t; return -1; };
    var median = function (p, q, r) { for (var s = 0; s < spheres.length; s++) { var a = dirOf(spheres[s], p), b = dirOf(spheres[s], q), c2 = dirOf(spheres[s], r); if (a !== b && a !== c2 && b !== c2) return s; } return -1; };
    spheres.forEach(function (S) { var md = S.dirs.map(function (D) { return D.filter(function (q) { return img(q) >= 0; }); }).filter(function (D) { return D.length; }), imgs = {};
      S.markedDirs = md.length;
      if (md.length < 3) { S.image = -1; S.consistent = null; return; }
      for (var a = 0; a < md.length; a++) for (var b = a + 1; b < md.length; b++) for (var c2 = b + 1; c2 < md.length; c2++)
        md[a].forEach(function (p) { md[b].forEach(function (q) { md[c2].forEach(function (r) { var p1 = img(p), q1 = img(q), r1 = img(r); imgs[(p1 === q1 || p1 === r1 || q1 === r1) ? -1 : median(p1, q1, r1)] = 1; }); }); });
      var ks = Object.keys(imgs).map(Number); S.consistent = ks.length === 1; S.image = ks.length === 1 ? ks[0] : -2; S.images = ks; });
    spheres.forEach(function (S, s) { var t = s, per = 0; for (var it = 1; it <= spheres.length + 1; it++) { t = spheres[t] ? spheres[t].image : -1; if (t < 0) break; if (t === s) { per = it; break; } } S.period = per; });
    // coordinate of each sphere: representatives of three directions -> -1, 1, inf (the direction with the most points goes to inf)
    spheres.forEach(function (S) { var order = S.dirs.map(function (D, t) { return t; }).sort(function (u, v) { return S.dirs[v].length - S.dirs[u].length; });
      var rep = function (t) { return pts[S.dirs[t][0]]; }, src = [rep(order[1]), rep(order[2]), rep(order[0])];
      S.T = PM.moebius3(src, [[-1, 0], [1, 0], null]); S.Ti = PM.moebius3([[-1, 0], [1, 0], null], src); S.atInf = order[0]; });
    var name = function (D) { return D.map(function (q) { return labels[q]; }).join(' '); };
    spheres.forEach(function (S, s) { S.text = 'S' + s + ': ' + S.dirs.map(function (D) { return '{' + name(D) + '}'; }).join(' | '); });
    return { labels: labels, pts: pts, necks: necks, compatible: true, spheres: spheres, eps: eps };
  }
  function summary(R) {
    if (!R.compatible) return 'necks are not compatible at this threshold (no tree): ' + R.necks.length + ' necks';
    return R.spheres.map(function (S, s) { return S.text + '  →  ' + (S.image >= 0 ? 'S' + S.image : S.image === -1 ? (S.markedDirs < 3 ? '? (fewer than 3 directions with points whose image is tracked)' : 'collapses') : 'inconsistent ' + JSON.stringify(S.images)) +
      (S.period ? '  [periodic, period ' + S.period + ']' : ''); }).join('\n');
  }
  return { analyze: analyze, summary: summary };
})();
if (typeof module !== 'undefined') module.exports = SPH;
