// E96 (REGISTRATION addendum (x), 2026-11-12): core entropy of z^2 + c_theta from the angle theta = num/den, by Thurston's algorithm on pairs of
// postcritical angles. Pure functions, node-testable.
'use strict';
var ENT = (function () {
  // spectral radius of the nonnegative matrix given by successor lists, as the max over strongly connected components
  function rho(succ) {
    var N = succ.length, idx = new Array(N).fill(-1), low = new Array(N).fill(0), on = new Array(N).fill(false), st = [], comp = new Array(N).fill(-1), nc = 0, c = 0;
    for (var s = 0; s < N; s++) if (idx[s] < 0) {           // iterative Tarjan
      var stack = [[s, 0]]; idx[s] = low[s] = c++; st.push(s); on[s] = true;
      while (stack.length) { var top = stack[stack.length - 1], v = top[0];
        if (top[1] < succ[v].length) { var w = succ[v][top[1]++]; if (idx[w] < 0) { idx[w] = low[w] = c++; st.push(w); on[w] = true; stack.push([w, 0]); } else if (on[w]) low[v] = Math.min(low[v], idx[w]); }
        else { stack.pop(); if (stack.length) { var u = stack[stack.length - 1][0]; low[u] = Math.min(low[u], low[v]); }
          if (low[v] === idx[v]) { var x; do { x = st.pop(); on[x] = false; comp[x] = nc; } while (x !== v); nc++; } } } }
    var best = 0;
    for (var k = 0; k < nc; k++) { var nodes = []; for (var i = 0; i < N; i++) if (comp[i] === k) nodes.push(i);
      var edges = 0; nodes.forEach(function (i) { succ[i].forEach(function (j) { if (comp[j] === k) edges++; }); }); if (!edges) continue;
      if (edges === nodes.length) { best = Math.max(best, 1); continue; }             // a single cycle
      var pos = {}; nodes.forEach(function (i, a) { pos[i] = a; }); var vec = nodes.map(function () { return 1; }), lam = 1;
      for (var it = 0; it < 20000; it++) { var nw = vec.slice();                       // power iteration on I + A (primitive on a component)
        nodes.forEach(function (i, a) { succ[i].forEach(function (j) { if (comp[j] === k) nw[a] += vec[pos[j]]; }); });
        var m = 0; for (var a = 0; a < nw.length; a++) m = Math.max(m, nw[a]); for (a = 0; a < nw.length; a++) nw[a] /= m; var d = Math.abs(m - lam); lam = m; vec = nw; if (d < 1e-15*m && it > 50) break; }
      best = Math.max(best, lam - 1); }
    return best;
  }
  // core entropy of the angle num/den (integers, 0 < num < den). Returns { h, lambda, orbit: number of postcritical angles, pairs }
  function core(num, den) {
    var D = 2*den, orb = [], seen = {}, a = ((num % den) + den) % den;
    while (!seen[a]) { seen[a] = true; orb.push(a); a = (2*a) % den; }
    var e1 = num % D, e2 = (num + den) % D, lo = Math.min(e1, e2), hi = Math.max(e1, e2);
    var side = function (t) { var p = (2*t) % D; return p === lo || p === hi ? 0 : (p > lo && p < hi ? 1 : -1); };   // positions doubled: the diameter ends are lo, hi
    var id = {}, pairs = [], key = function (p, q) { return p < q ? p + ',' + q : q + ',' + p; };
    for (var i = 0; i < orb.length; i++) for (var j = i + 1; j < orb.length; j++) { id[key(orb[i], orb[j])] = pairs.length; pairs.push([orb[i], orb[j]]); }
    var th = ((num % den) + den) % den, succ = pairs.map(function (pq) { var p = pq[0], q = pq[1], p2 = (2*p) % den, q2 = (2*q) % den, out = [], add = function (u, v) { if (u !== v) out.push(id[key(u, v)]); };
      if (side(p)*side(q) < 0) { add(p2, th); add(th, q2); } else add(p2, q2); return out; });
    var r = rho(succ); return { h: r > 1 ? Math.log(r) : 0, lambda: Math.max(1, r), orbit: orb.length, pairs: pairs.length };
  }
  return { core: core, rho: rho };
})();
if (typeof module !== 'undefined') module.exports = ENT;
