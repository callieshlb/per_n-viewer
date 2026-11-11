// E96: the three triangles of a bitransitive component (REGISTRATION addendum 2026-11-11 (j)). Depends on pern_math.js (PM).
'use strict';
var BTRI = (function () {
  // mapOf(p) -> M or null; b: center (parameter point); j: slot. opts: N (grid), half (box half-size; found by doubling if absent), offs: [dx, dy] unused
  function fillBox(mapOf, b, j, half, N) {
    var W = new Float64Array(2*N*N).fill(NaN), AR = new Float64Array(4*N*N), st = new Uint8Array(N*N), dx = 2*half/(N - 1), x0 = b[0] - half, y0 = b[1] - half, ci = (N - 1) >> 1;
    var facAt = function (ix, iy) { var M = mapOf([x0 + ix*dx, y0 + iy*dx]); if (!M) return null; return PM.bitransFactors(M, j); };
    var roots = function (P) { var m = Math.cbrt(Math.hypot(P[0], P[1])), t = Math.atan2(P[1], P[0])/3, o = []; for (var k = 0; k < 3; k++) o.push([m*Math.cos(t + 2*Math.PI*k/3), m*Math.sin(t + 2*Math.PI*k/3)]); return o; };
    var idx = function (ix, iy) { return iy*N + ix; }, q = [], count = 0, rejDisc = 0, rejOut = 0, sup = 0, edge = false, DX = [1, -1, 0, 0, 1, 1, -1, -1], DY = [0, 0, 1, -1, 1, -1, 1, -1];
    var c0 = idx(ci, ci); st[c0] = 1; W[2*c0] = 0; W[2*c0 + 1] = 0; count++;            // the center: all arguments 0
    var assign = function (ix, iy, par) {                               // par: index of an assigned neighbor (arguments continued from it); returns true if accepted
      var id = idx(ix, iy), Fc = facAt(ix, iy); if (!Fc) { st[id] = 2; rejOut++; return false; }
      var A = [0, 0, 0, 0]; for (var k = 0; k < 4; k++) { var a = Math.atan2(Fc.r[k][1], Fc.r[k][0]), ref = AR[4*par + k]; A[k] = a + 2*Math.PI*Math.round((ref - a)/(2*Math.PI)); }
      var P = PM.bitransPhiFrom(Fc, A); if (!isFinite(P[0] + P[1]) || Math.hypot(P[0], P[1]) >= 1) { st[id] = 2; rejOut++; return false; }
      var R = roots(P), best = R[0], px = (ix - ci)*dx, py = (iy - ci)*dx, pred = kappa ? [kappa[0]*px - kappa[1]*py, kappa[0]*py + kappa[1]*px] : null;
      if (pred && Math.hypot(pred[0], pred[1]) < 0.3) { var bp = Infinity; for (k = 0; k < 3; k++) { var dq = Math.hypot(R[k][0] - pred[0], R[k][1] - pred[1]); if (dq < bp) { bp = dq; best = R[k]; } } }   // near the center W ~ kappa (p - b): the nearest-to-parent rule is ill-conditioned there (FAILURE kept in log.md)
      else if (par !== c0) { var rx = W[2*par], ry = W[2*par + 1], bd = Infinity; for (k = 0; k < 3; k++) { var dd = Math.hypot(R[k][0] - rx, R[k][1] - ry); if (dd < bd) { bd = dd; best = R[k]; } }
        var mr = Math.hypot(rx, ry), mb = Math.hypot(best[0], best[1]); if (bd > 0.3 || (mr > 0.8 && mb < mr - 0.07)) { st[id] = 3; rejDisc++; return false; } }
      st[id] = 1; W[2*id] = best[0]; W[2*id + 1] = best[1]; for (k = 0; k < 4; k++) AR[4*id + k] = A[k]; count++; sup = Math.max(sup, Math.hypot(best[0], best[1])); q.push(ix, iy); return true; };
    var kappa = null, F0 = facAt(ci + 1, ci); if (F0) { var P0 = PM.bitransPhiFrom(F0, F0.r.map(function (r) { return Math.atan2(r[1], r[0]); })), k0 = roots(P0)[0]; kappa = [k0[0]/dx, k0[1]/dx]; }
    if (!assign(ci + 1, ci, c0)) return { W: W, N: N, x0: x0, y0: y0, dx: dx, count: count, fail: 'no value next to the center' };
    var head = 0;
    while (head < q.length) { var ix = q[head++], iy = q[head++], me = idx(ix, iy);
      for (var d = 0; d < 8; d++) { var nx = ix + DX[d], ny = iy + DY[d]; if (nx < 0 || ny < 0 || nx >= N || ny >= N) { edge = true; continue; } if (st[idx(nx, ny)]) continue; assign(nx, ny, me); } }
    return { W: W, st: st, N: N, x0: x0, y0: y0, dx: dx, count: count, rejDisc: rejDisc, rejOut: rejOut, sup: sup, edge: edge, half: half };
  }
  function fill(mapOf, b, j, opts) {
    opts = opts || {}; var half = opts.half || 1e-3*(1 + Math.hypot(b[0], b[1])), R = null;
    for (var it = 0; it < 16; it++) { R = fillBox(mapOf, b, j, half, 41); if (R.fail) return R; if (!R.edge) break; half *= 2; }
    return fillBox(mapOf, b, j, half*1.15, opts.N || 241);
  }
  // triangle index (0, 1, 2) and internal angle in [0, 1) of a W value
  function tri(wx, wy) { var a = Math.atan2(wy, wx); if (a < 0) a += 2*Math.PI; var s = a/(2*Math.PI/3), k = Math.min(2, Math.floor(s)); return [k, s - k]; }
  return { fill: fill, fillBox: fillBox, tri: tri };
})();
if (typeof module !== 'undefined') module.exports = BTRI;
