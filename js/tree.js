// E96 Hubbard tree off-center (REGISTRATION addendum 2026-11-10 (e)). Depends on pern_math.js (PM).
// Coordinates: w with a_0 = 0, omega = inf, a_1 = c0. psi: w -> chart point (Moebius). A tree is a list of arcs (polylines of [re, im]).
'use strict';
var HT = (function () {
  var cmul = PM.cmul, cdiv = PM.cdiv, cadd = PM.cadd, csub = PM.csub, cabs = PM.cabs;
  function csqrt(z) { var r = Math.hypot(z[0], z[1]), t = Math.atan2(z[1], z[0])/2, s = Math.sqrt(r); return [s*Math.cos(t), s*Math.sin(t)]; }
  var dist = function (a, b) { return Math.hypot(a[0] - b[0], a[1] - b[1]); };
  // both chart preimages of the chart point u (null = infinity) under M
  function preimages(M, u) {
    var a, b, c;
    if (u === null) { a = M.D; b = M.E; c = M.F; } else { a = csub(M.A, cmul(u, M.D)); b = csub(M.B, cmul(u, M.E)); c = csub(M.C, cmul(u, M.F)); }
    var sc = Math.max(cabs(a), cabs(b), cabs(c));
    if (cabs(a) < 1e-13*sc) return [null, cabs(b) < 1e-13*sc ? null : cdiv([-c[0], -c[1]], b)];
    var disc = csqrt(csub(cmul(b, b), cmul([4*a[0], 4*a[1]], c))), q = (b[0]*disc[0] + b[1]*disc[1] >= 0) ? cadd(b, disc) : csub(b, disc);
    q = [-q[0]/2, -q[1]/2];
    if (cabs(q) < 1e-300) return [[0, 0], [0, 0]];
    return [cdiv(q, a), cdiv(c, q)];
  }
  // setup for the map M with center parameter c0: psi, its inverse, marked points in w, the w-preimage function
  function frame(M, c0) {
    var psi = PM.moebius3([[0, 0], c0, null], [M.marked[0], M.marked[1], M.omega]), inv = PM.moebius3([M.marked[0], M.marked[1], M.omega], [[0, 0], c0, null]);
    return { psi: psi, inv: inv, marked: M.marked.map(function (a) { return inv(a); }),
      pre: function (w) { var r = preimages(M, psi(w)); return [inv(r[0]), inv(r[1])]; }, f: function (w) { return inv(PM.fz(M, psi(w))); } };
  }
  function toGraph(arcs, tol) {
    var nodes = [], edges = [], ends = [];
    var endId = function (p) { for (var i = 0; i < ends.length; i++) if (dist(nodes[ends[i]], p) < tol) return ends[i]; nodes.push(p.slice()); ends.push(nodes.length - 1); return nodes.length - 1; };
    arcs.forEach(function (arc) { if (arc.length < 2) return; var prev = endId(arc[0]);
      for (var i = 1; i < arc.length; i++) { var id; if (i === arc.length - 1) id = endId(arc[i]); else { nodes.push(arc[i].slice()); id = nodes.length - 1; } if (id !== prev) edges.push([prev, id]); prev = id; } });
    return { nodes: nodes, edges: edges };
  }
  function nearest(nodes, p) { var b = -1, bd = Infinity; for (var i = 0; i < nodes.length; i++) { if (!nodes[i]) continue; var d = dist(nodes[i], p); if (d < bd) { bd = d; b = i; } } return [b, bd]; }
  function resample(path, h) {
    var L = 0, i; for (i = 1; i < path.length; i++) L += dist(path[i], path[i - 1]);
    var k = Math.max(1, Math.ceil(L/h)), out = [path[0]], acc = 0, j = 1, step = L/k, target = step;
    for (i = 1; i < path.length && j < k; i++) { var d = dist(path[i], path[i - 1]);
      while (j < k && acc + d >= target) { var t = d > 0 ? (target - acc)/d : 0; out.push([path[i - 1][0] + t*(path[i][0] - path[i - 1][0]), path[i - 1][1] + t*(path[i][1] - path[i - 1][1])]); j++; target = j*step; }
      acc += d; }
    out.push(path[path.length - 1]); return out;
  }
  // one pullback step. arcs: current tree (w coords). Returns {arcs, why?}
  function step(arcs, F, h, mpos) {                                  // mpos[j]: where a_j sits on the CURRENT tree (default: its exact position)
    mpos = mpos || F.marked; var n = F.marked.length, G = toGraph(arcs, h*1e-3), N = G.nodes.length, i, j;
    for (i = 0; i < N; i++) if (!isFinite(G.nodes[i][0] + G.nodes[i][1])) return { why: 'tree reaches infinity' };
    var mk = mpos.map(function (a) { return a === null ? [-1, Infinity] : nearest(G.nodes, a); });
    for (j = 0; j < n; j++) if (mk[j][1] > 5*h) return { why: 'marked point a' + j + ' is off the tree' };
    for (j = 0; j < n; j++) if (F.marked[j] === null) return { why: 'a marked point is at infinity in the tree coordinate' };
    var P = new Array(2*N), pre = G.nodes.map(function (z) { return F.pre(z); });
    for (i = 0; i < N; i++) { if (pre[i][0] === null || pre[i][1] === null) return { why: 'a preimage of the tree is at infinity (the tree meets f(omega))' }; P[2*i] = pre[i][0]; P[2*i + 1] = pre[i][1]; }
    var par = new Int32Array(2*N); for (i = 0; i < 2*N; i++) par[i] = i; var find = function (x) { while (par[x] !== x) { par[x] = par[par[x]]; x = par[x]; } return x; };
    var a1 = mk[1][0]; par[2*a1 + 1] = 2*a1; P[2*a1] = [0, 0];                               // the two preimages of a_1 are a_0
    var adj = {}, add = function (u, v) { u = find(u); v = find(v); if (u === v) return; (adj[u] = adj[u] || new Set()).add(v); (adj[v] = adj[v] || new Set()).add(u); };
    G.edges.forEach(function (e) { var u = e[0], v = e[1];
      if (dist(P[2*u], P[2*v]) + dist(P[2*u + 1], P[2*v + 1]) <= dist(P[2*u], P[2*v + 1]) + dist(P[2*u + 1], P[2*v])) { add(2*u, 2*v); add(2*u + 1, 2*v + 1); } else { add(2*u, 2*v + 1); add(2*u + 1, 2*v); } });
    var keep = new Set();                                                                    // new marked nodes: the preimage of a_(j+1) nearest to a_j
    for (j = 0; j < n; j++) { var src = mk[(j + 1) % n][0], c0 = find(2*src), c1 = find(2*src + 1), t = F.marked[j];
      var pick = dist(P[c0], t) <= dist(P[c1], t) ? c0 : c1; if (dist(P[pick], t) > 1e-6*(1 + cabs(t)) && !(mpos !== F.marked)) return { why: 'the preimage tree misses a' + j }; P[pick] = t.slice(); keep.add(pick); }
    var leaves = Object.keys(adj).map(Number).filter(function (u) { return adj[u].size === 1 && !keep.has(u); });
    while (leaves.length) { var u = leaves.pop(); if (!adj[u] || adj[u].size !== 1 || keep.has(u)) continue; var v = adj[u].values().next().value; delete adj[u]; adj[v].delete(u); if (adj[v].size === 1 && !keep.has(v)) leaves.push(v); if (adj[v].size === 0) delete adj[v]; }
    var root = keep.values().next().value, seen = new Set([root]), stack = [root], ne = 0;      // connected component of the marked points
    while (stack.length) { var x = stack.pop(); if (!adj[x]) continue; adj[x].forEach(function (y) { ne++; if (!seen.has(y)) { seen.add(y); stack.push(y); } }); }
    var missing = []; keep.forEach(function (k) { if (!seen.has(k)) missing.push(k); }); if (missing.length) return { why: 'the preimage tree is disconnected (' + missing.length + ' marked points cut off)' };
    if (ne/2 !== seen.size - 1) return { why: 'the preimage is not a tree (a loop appeared)' };
    var special = function (u) { return keep.has(u) || adj[u].size !== 2; }, out = [], done = new Set();
    seen.forEach(function (u) { if (!special(u)) return; adj[u].forEach(function (v) { if (done.has(u + ',' + v)) return; var path = [P[u]], prev = u, curN = v;
        while (!special(curN)) { path.push(P[curN]); var nx = null; adj[curN].forEach(function (y) { if (y !== prev) nx = y; }); prev = curN; curN = nx; }
        path.push(P[curN]); done.add(u + ',' + v); done.add(curN + ',' + prev); out.push(resample(path, h)); }); });
    var tot = 0; out.forEach(function (q) { tot += q.length; }); if (tot > 12000) return { why: 'the tree length blows up (more than 12000 nodes)' };
    return { arcs: out };
  }
  function hausdorff(A, B) {                 // max over nodes of A of the distance to the nodes of B, and conversely
    var fa = [].concat.apply([], A), fb = [].concat.apply([], B), one = function (X, Y) { var m = 0; for (var i = 0; i < X.length; i++) { var b = Infinity; for (var j = 0; j < Y.length; j++) { var d = dist(X[i], Y[j]); if (d < b) b = d; } if (b > m) m = b; } return m; };
    return Math.max(one(fa, fb), one(fb, fa));
  }
  function shape(arcs, h) { var G = toGraph(arcs, h*1e-3), deg = {}; G.edges.forEach(function (e) { deg[e[0]] = (deg[e[0]] || 0) + 1; deg[e[1]] = (deg[e[1]] || 0) + 1; });
    return Object.keys(deg).map(function (k) { return deg[k]; }).filter(function (d) { return d >= 3; }).sort().join(','); }
  // deform: iterate from arcs0. Returns {arcs (w coords), steps, defect, h, diam, why?, shape}
  function deform(M, c0, arcs0, opts) {
    opts = opts || {}; var F = frame(M, c0), fl = [].concat.apply([], arcs0), x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    fl.forEach(function (p) { x0 = Math.min(x0, p[0]); x1 = Math.max(x1, p[0]); y0 = Math.min(y0, p[1]); y1 = Math.max(y1, p[1]); });
    var diam = Math.hypot(x1 - x0, y1 - y0), h = diam/(opts.res || 500), arcs = arcs0.map(function (a) { return resample(a, h); }), defect = Infinity, k, why = null;
    var mpos = opts.mpos || null;
    for (k = 0; k < (opts.maxSteps || 80); k++) { var r = step(arcs, F, h, mpos); mpos = null; if (r.why) { why = r.why; break; } defect = hausdorff(r.arcs, arcs); arcs = r.arcs; if (defect < (opts.tol || 0.45)*h) { k++; break; } }
    return { arcs: arcs, steps: k, defect: defect, h: h, diam: diam, why: why, shape: shape(arcs, h), F: F, mpos: F.marked };
  }
  return { deform: deform, frame: frame, step: step, hausdorff: hausdorff, preimages: preimages, resample: resample, shape: shape };
})();
if (typeof module !== 'undefined') module.exports = HT;
