// Per_6(0) explorer core (E92). Double precision. Shared by the page, its web workers, and node tests.
// Curve: P7(x, y) = sum_k x^k sum_j coef[k][j] y^j = 0 in Caroline's chart
//   f(z) = (z - x)(z - w)/z^2,  w = (x + y - 1)/(x - 1);  y = a_3, x = a_5; marked orbit 0 -> inf -> 1 -> y -> f(y) -> x -> 0.
// Sheets: the 5 roots x(y), labeled by continuation from the base point y0 along the straight segment,
// with a detour around the node at y = 0 (trivial monodromy) on the segment's own side.
'use strict';
var CORE = (function () {
  var D = null;            // data: coef, y0, x0, bps, perms, centers, punctures
  var NODE_R = 0.12;
  var specials = [];       // branch points (incl. the node) for step-size control

  function init(data) {
    D = data;
    specials = data.bps.map(function (b) { return [b[0], b[1]]; });
  }

  // ---- complex helpers on [re, im]
  function cmul(a, b) { return [a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]]; }
  function cdiv(a, b) { var d = b[0]*b[0]+b[1]*b[1]; return [(a[0]*b[0]+a[1]*b[1])/d, (a[1]*b[0]-a[0]*b[1])/d]; }
  function cabs(a) { return Math.hypot(a[0], a[1]); }
  function csub(a, b) { return [a[0]-b[0], a[1]-b[1]]; }
  function cadd(a, b) { return [a[0]+b[0], a[1]+b[1]]; }
  function chord(a, b) {     // chordal distance; b may be null = infinity
    if (b === null) return 1/Math.sqrt(1+a[0]*a[0]+a[1]*a[1]);
    var dr = a[0]-b[0], di = a[1]-b[1];
    return Math.sqrt(dr*dr+di*di)/Math.sqrt((1+a[0]*a[0]+a[1]*a[1])*(1+b[0]*b[0]+b[1]*b[1]));
  }

  // ---- coefficients c_k(y) and dc_k/dy
  function coeffs(y) {
    var n = D.coef.length, c = [], dc = [];
    for (var k = 0; k < n; k++) {
      var row = D.coef[k], m = row.length;
      var pr = 0, pi = 0, qr = 0, qi = 0;     // Horner for value and derivative
      for (var j = m-1; j >= 0; j--) {
        var tr = qr*y[0]-qi*y[1]+pr, ti = qr*y[1]+qi*y[0]+pi;   // q = q*y + p
        qr = tr; qi = ti;
        var ur = pr*y[0]-pi*y[1]+row[j], ui = pr*y[1]+pi*y[0];   // p = p*y + a_j
        pr = ur; pi = ui;
      }
      c.push([pr, pi]); dc.push([qr, qi]);
    }
    // leading coefficient is exactly (y-1)^2; avoid cancellation near y = 1
    var e = [y[0]-1, y[1]];
    c[n-1] = cmul(e, e); dc[n-1] = [2*e[0], 2*e[1]];
    return [c, dc];
  }
  // P, dP/dx, dP/dy at (x, y) given coefficients
  function pEval(x, cc) {
    var c = cc[0], dc = cc[1], n = c.length;
    var p = [0, 0], px = [0, 0], py = [0, 0];
    for (var k = n-1; k >= 0; k--) {
      px = cadd(cmul(px, x), p);
      p = cadd(cmul(p, x), c[k]);
      py = cadd(cmul(py, x), dc[k]);
    }
    return [p, px, py];
  }
  function newton(x, cc, maxit) {
    for (var it = 0; it < maxit; it++) {
      var e = pEval(x, cc);
      var dx = cdiv(e[0], e[1]);
      x = csub(x, dx);
      if (cabs(dx) <= 1e-11*(1+cabs(x))) return [x, it+1];   // 1e-13 is unreachable next to the near-double roots at the node y = 0
    }
    return [x, -1];
  }
  function distSpecial(y) {
    var m = Infinity;
    for (var i = 0; i < specials.length; i++) { var d = Math.hypot(y[0]-specials[i][0], y[1]-specials[i][1]); if (d < m) m = d; }
    return m;
  }
  // waypoints from a to b, detouring around the node at 0 on the segment's own side
  function arcTo(r, f1, f2, n) {        // points on |y| = r from angle f1 to f2 (shorter way); excludes f1
    var df = ((f2-f1+Math.PI) % (2*Math.PI) + 2*Math.PI) % (2*Math.PI) - Math.PI, out = [];
    for (var i = 1; i <= n; i++) { var f = f1+df*i/n; out.push([r*Math.cos(f), r*Math.sin(f)]); }
    return out;
  }
  function waypoints(a, b) {
    var v = csub(b, a), L2 = v[0]*v[0]+v[1]*v[1];
    if (L2 === 0) return [a, b];
    var ra = cabs(a), rb = cabs(b);
    // endpoints inside the node disk: the only special point there is the node (trivial monodromy),
    // so any path in the disk is equivalent; never approach 0 closer than the endpoints do.
    if (rb < NODE_R || ra < NODE_R) {
      var tt = Math.min(1, Math.max(0, -(a[0]*v[0]+a[1]*v[1])/L2));
      var cl = cabs([a[0]+tt*v[0], a[1]+tt*v[1]]);
      if (cl >= 0.5*Math.min(ra, rb)) return [a, b];
      var r = Math.max(Math.min(Math.max(ra, rb), NODE_R), 1e-12);
      var A = [r*a[0]/ra, r*a[1]/ra], Bp = [r*b[0]/rb, r*b[1]/rb];
      var path = [a];
      if (ra > r) path = waypoints(a, A); else if (ra < r) path.push(A);
      path = path.concat(arcTo(r, Math.atan2(a[1], a[0]), Math.atan2(b[1], b[0]), 24));
      if (rb !== r) path.push(b);
      return path;
    }
    var t = -(a[0]*v[0]+a[1]*v[1])/L2;
    var cl = [a[0]+t*v[0], a[1]+t*v[1]];
    if (!(t > 0 && t < 1 && cabs(cl) < NODE_R && cabs(a) > NODE_R && cabs(b) > NODE_R)) return [a, b];
    var A = L2, Bq = 2*(a[0]*v[0]+a[1]*v[1]), C = a[0]*a[0]+a[1]*a[1]-NODE_R*NODE_R;
    var disc = Math.sqrt(Math.max(Bq*Bq-4*A*C, 0));
    var t1 = (-Bq-disc)/(2*A), t2 = (-Bq+disc)/(2*A);
    var p1 = [a[0]+t1*v[0], a[1]+t1*v[1]], p2 = [a[0]+t2*v[0], a[1]+t2*v[1]];
    var f1 = Math.atan2(p1[1], p1[0]), f2 = Math.atan2(p2[1], p2[0]);
    var df = ((f2-f1+Math.PI) % (2*Math.PI) + 2*Math.PI) % (2*Math.PI) - Math.PI;
    var out = [a, p1];
    for (var i = 1; i <= 24; i++) { var f = f1+df*i/24; out.push([NODE_R*Math.cos(f), NODE_R*Math.sin(f)]); }
    out.push(b);
    return out;
  }
  // track the roots xs (array) along the straight segment a -> b with adaptive steps;
  // a step is accepted when every root's Newton converges from its tangent predictor and roots stay distinct.
  function trackSegment(xs, a, b) {
    var v = csub(b, a), L = cabs(v);
    if (L === 0) return xs;
    var t = 0, cur = xs.slice(), cc = coeffs(a);
    var dxdy = cur.map(function (x) { var e = pEval(x, cc); return cdiv([-e[2][0], -e[2][1]], e[1]); });
    var guard = 0;
    while (t < 1 && guard < 200000) {
      guard++;
      var ycur = [a[0]+t*v[0], a[1]+t*v[1]];
      var h = Math.min(0.5, Math.max(1e-7, 0.15*distSpecial(ycur)));
      var dt = Math.min(1-t, h/L);
      var ok = false;
      for (var tries = 0; tries < 40 && !ok; tries++) {
        var t2 = t+dt, y2 = [a[0]+t2*v[0], a[1]+t2*v[1]], dy = [v[0]*dt, v[1]*dt];
        var cc2 = coeffs(y2), nx = [], good = true;
        for (var r = 0; r < cur.length && good; r++) {
          var pred = cadd(cur[r], cmul(dxdy[r], dy));
          var res = newton(pred, cc2, 12);
          if (res[1] < 0 || cabs(csub(res[0], pred)) > 0.3*(cabs(cmul(dxdy[r], dy))+1e-9*(1+cabs(cur[r])))+1e-12) good = false;
          nx.push(res[0]);
        }
        if (good && nx.length > 1) {
          for (var i = 0; i < nx.length && good; i++) for (var j = i+1; j < nx.length; j++)
            if (chord(nx[i], nx[j]) < 1e-10) { good = false; break; }
        }
        if (good) {
          ok = true; t = t2; cur = nx; cc = cc2;
          dxdy = cur.map(function (x) { var e = pEval(x, cc); return cdiv([-e[2][0], -e[2][1]], e[1]); });
        } else dt /= 2;
      }
      if (!ok) return null;          // could not continue (very close to a branch point)
    }
    return cur;
  }
  function track(xs, a, b) {
    var w = waypoints(a, b), cur = xs;
    for (var i = 1; i < w.length && cur; i++) cur = trackSegment(cur, w[i-1], w[i]);
    return cur;
  }
  // the 5 globally labeled roots at y (atlas sheets)
  function globalRoots(y) { return track(D.x0.slice(), D.y0, y); }
  function sheetOf(y, x) {
    var g = globalRoots(y);
    if (!g) return -1;
    var best = 0, bd = Infinity;
    for (var k = 0; k < g.length; k++) { var d = chord(g[k], x); if (d < bd) { bd = d; best = k; } }
    return best;
  }

  // ---- dynamics: f(z) = 1 - s/z + p/z^2 with s = x + w, p = x w.  z = null means infinity.
  function mapData(x, y) {
    var w = cdiv([x[0]+y[0]-1, x[1]+y[1]], [x[0]-1, x[1]]);
    var s = cadd(x, w), p = cmul(x, w);
    var om = cdiv([2*p[0], 2*p[1]], s);
    var M = { x: x, y: y, w: w, s: s, p: p, omega: om };
    M.a4 = fz(M, y);
    M.marked = [[0, 0], null, [1, 0], y, M.a4, x];
    // local scale of each marked point: chordal distance to the nearest other marked point
    M.scale = M.marked.map(function (a, j) {
      var m = 1;
      for (var i = 0; i < 6; i++) if (i !== j) m = Math.min(m, chordN(a, M.marked[i]));
      return m;
    });
    // chart of the census: F = 1 + b/z + c/z^2 with b = -s, c = p
    M.b = [-s[0], -s[1]]; M.c = p;
    return M;
  }
  function fz(M, z) {
    if (z === null) return [1, 0];
    var r2 = z[0]*z[0]+z[1]*z[1];
    if (r2 < 1e-300) return null;
    if (r2 > 1e300) return [1, 0];
    var iz = [z[0]/r2, -z[1]/r2];
    var t1 = cmul(M.s, iz), iz2 = cmul(iz, iz), t2 = cmul(M.p, iz2);
    return [1-t1[0]+t2[0], -t1[1]+t2[1]];
  }
  function chordN(a, b) {   // chordal distance allowing nulls (infinity)
    if (a === null && b === null) return 0;
    if (a === null) return chord(b, null);
    return chord(a, b);
  }
  // fate of a starting point z under f.
  // returns {kind: 'capture', j, t, phase} | {kind: 'cycle', period, t, index} | {kind: 'undecided'}
  // rel = true: capture radius at a_j is delta * (distance from a_j to the nearest other marked point)
  function fate(M, z, maxIt, delta, tolP, cycle, rel) {
    var ring = new Array(32);
    var th = M.scale.map(function (s) { return rel ? delta*s : delta; });
    for (var it = 0; it < maxIt; it++) {
      z = fz(M, z);
      for (var j = 0; j < 6; j++) {
        if (chordN(z, M.marked[j]) < th[j]) return { kind: 'capture', j: j, t: it+1, phase: ((j-1-(it+1)) % 6 + 6) % 6 };
      }
      if (cycle) {
        for (var q = 0; q < cycle.pts.length; q++)
          if (chordN(z, cycle.pts[q]) < delta) return { kind: 'cycle', period: cycle.pts.length, t: it+1, index: ((q-(it+1)) % cycle.pts.length + cycle.pts.length) % cycle.pts.length };
      }
      ring[it % 32] = z;
    }
    if (!cycle) {
      for (var per = 1; per <= 24; per++) {
        var zp = ring[(maxIt-1-per+3200) % 32], zl = ring[(maxIt-1) % 32];
        if (zp !== undefined && chordN(zl, zp) < tolP) return { kind: 'cycle', period: per, t: maxIt };
      }
    }
    return { kind: 'undecided' };
  }
  // attracting cycle of the free critical point, if any (for the dynamical plane)
  function freeCycle(M, maxIt, tolP) {
    var z = M.omega, seq = [];
    var n = Math.max(maxIt, 200);
    for (var it = 0; it < n; it++) {
      z = fz(M, z);
      for (var j = 0; j < 6; j++) if (chordN(z, M.marked[j]) < 1e-9) return null;     // captured: no free cycle
    }
    var z0 = z;
    for (var per = 1; per <= 48; per++) {
      z = fz(M, z);
      seq.push(z);
      if (chordN(z, z0) < tolP) return { pts: seq };
    }
    return null;
  }

  return { init: init, coeffs: coeffs, pEval: pEval, newton: newton, track: track, waypoints: waypoints,
           globalRoots: globalRoots, sheetOf: sheetOf, mapData: mapData, fz: fz, fate: fate, freeCycle: freeCycle,
           chord: chord, chordN: chordN, cabs: cabs, cmul: cmul, cdiv: cdiv, csub: csub, cadd: cadd,
           get D() { return D; } };
})();
if (typeof module !== 'undefined') module.exports = CORE;
