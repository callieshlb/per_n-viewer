// E96 (2026-11-18): angled internal address of a hyperbolic component of the Mandelbrot set from one of its characteristic angles, via wakes.
// The components whose wakes contain theta lie on the vein from 0 to c(theta) and are nested; W_0 = main cardioid, W_(k+1) = the component of lowest
// period among those inside the wake of W_k (that is the definition of the internal address). The arrow W_k -> W_(k+1) gets the angle p/q of the limb
// of W_k containing theta: the p/q-limb of W (angles w-, w+ of period s) is the wake of the satellite W ⊥ R_p/q, whose angles are obtained from the
// rabbit's angles by Douady's substitution 0 -> w-, 1 -> w+. Needs MATE.lavaurs (pairs of every period) and RESC.rabbit. Exact (BigInt) comparisons.
'use strict';
var ADDR = (function () {
  var cmpF = function (a, b) { var l = a[0]*b[1], r = b[0]*a[1]; return l < r ? -1 : l > r ? 1 : 0; };
  var bits = function (k, len) { var t = k.toString(2); while (t.length < len) t = '0' + t; return t; };
  var val = function (per) { return [BigInt('0b' + per), (1n << BigInt(per.length)) - 1n]; };
  // all components of period <= P whose closed wake contains theta = [num, den] (BigInt), outer to inner
  function vein(th, P) {
    var pairs = MATE.lavaurs(P), out = [{ n: 1, k0: 0, k1: 1, Q: 1, w: 1 }];
    for (var n = 2; n <= P; n++) { var Q = BigInt(Math.pow(2, n) - 1);
      pairs[n].forEach(function (pr) { if (cmpF([BigInt(pr[0]), Q], th) <= 0 && cmpF(th, [BigInt(pr[1]), Q]) <= 0) out.push({ n: n, k0: pr[0], k1: pr[1], Q: Number(Q), w: (pr[1] - pr[0])/Number(Q) }); }); }
    out.sort(function (a, b) { return b.w - a.w; }); return out;
  }
  // limb p/q of W (period s, angles k0, k1 over 2^s - 1; s = 1: the main cardioid) containing theta
  function limb(W, th, qmax) {
    var wm = W.n === 1 ? '0' : bits(W.k0, W.n), wp = W.n === 1 ? '1' : bits(W.k1, W.n);
    for (var q = 2; q <= qmax; q++) for (var p = 1; p < q; p++) { var g = function (a, b) { return b ? g(b, a % b) : a; }; if (g(p, q) !== 1) continue;
      var R = RESC.rabbit(p, q); if (!R) continue; var sub = function (k) { var b = bits(k, q), s = ''; for (var i = 0; i < q; i++) s += b[i] === '0' ? wm : wp; return s; };
      var lo = val(sub(R[0])), hi = val(sub(R[1])); if (cmpF(lo, th) <= 0 && cmpF(th, hi) <= 0) return [p, q]; }
    return null;
  }
  // angled internal address of the center with characteristic angle k/(2^n - 1)
  function of(n, k) {
    var th = [BigInt(k), (1n << BigInt(n)) - 1n], V = vein(th, n), cur = V[0], seq = [cur], i = 0;
    while (cur.n !== n) { var best = null, bi = -1; for (var j = i + 1; j < V.length; j++) if (!best || V[j].n < best.n) { best = V[j]; bi = j; } if (!best) break; seq.push(best); cur = best; i = bi; }
    var addr = seq.map(function (w) { return w.n; }), ang = [];
    for (var t = 0; t + 1 < seq.length; t++) { var q2 = Math.floor(seq[t + 1].n/seq[t].n) + 2, L = limb(seq[t], th, Math.max(2, q2)); ang.push(L); }
    return { address: addr, angled: ang };
  }
  // internal address only, from the kneading sequence (Lau-Schleicher): nu_j = 1 iff 2^(j-1) theta lies in (theta/2, (theta+1)/2), nu_n = * ;
  // rho(r) = min { j > r : nu_j != nu_(j-r) }, address = 1, rho(1), rho(rho(1)), ... up to n. Fast; used to filter before computing angles.
  function kneadAddr(n, k) {
    var Q = (1n << BigInt(n)) - 1n, D = 2n*Q, t = BigInt(k), lo = t, hi = t + Q, nu = [null];       // positions doubled: theta/2 <-> t over 2Q
    for (var j = 1; j <= n; j++) { var x = (2n*((t*(1n << BigInt(j - 1))) % Q)) % D; nu.push(j === n ? '*' : (x > lo && x < hi ? 1 : 0)); }
    var A = [1], r = 1; while (r < n) { var j2 = r + 1; while (j2 <= n && nu[j2] === nu[j2 - r]) j2++; if (j2 > n) break; A.push(j2); r = j2; }
    return A;
  }
  function fmt(A) { var s = String(A.address[0]); for (var i = 1; i < A.address.length; i++) { var a = A.angled[i - 1]; s += (a ? ' –' + a[0] + '/' + a[1] + '→ ' : ' → ') + A.address[i]; } return s; }
  // query: "1 3 6", "1-3-6", "1 → 3 → 6", or with angles after a period: "1 1/3 3 1/2 6", "1 -1/3-> 3 -> 6". Returns { address, angles: [[p, q] | null per arrow] } or null
  function parse(text) {
    var T = String(text).replace(/[→>]/g, ' ').replace(/(\d)\s*-\s*(?=\d+\s*\/)/g, '$1 ').replace(/-/g, ' ').replace(/[(),]/g, ' ').trim().split(/\s+/).filter(Boolean), addr = [], ang = [];
    for (var i = 0; i < T.length; i++) { var m = T[i].match(/^(\d+)\/(\d+)$/);
      if (m) { if (!addr.length || ang.length !== addr.length - 1) return null; ang.push([+m[1], +m[2]]); }
      else if (/^\d+$/.test(T[i])) { if (ang.length < addr.length) ang.push(null); addr.push(+T[i]); } else return null; }
    if (!addr.length) return null; if (addr[0] !== 1) { addr.unshift(1); ang.unshift(null); }
    while (ang.length > addr.length - 1) ang.pop();
    for (var k = 1; k < addr.length; k++) if (!(addr[k] > addr[k - 1])) return null;
    return { address: addr, angles: ang };
  }
  return { vein: vein, limb: limb, of: of, fmt: fmt, parse: parse, kneadAddr: kneadAddr };
})();
if (typeof module !== 'undefined') module.exports = ADDR;
