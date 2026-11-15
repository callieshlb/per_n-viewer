// E96 (REGISTRATION addendum (y), 2026-11-15): active rescaling limit of a period-n center f of the p/q-limb, from its angles.
// f = R_{p/q} ⊥ g  ->  g ⊔ C;   f in the decoration of R_{p/q} at the dyadic l/2^k  ->  t ⊔ C, t = Misiurewicz polynomial of angle l/2^k. Pure, BigInt.
'use strict';
var RESC = (function () {
  // characteristic angles of the rabbit R_{p/q} as numerators over 2^q - 1 (q >= 2): the period-q orbit of rotation number p/q, shortest sector
  function rabbit(p, q) {
    var Q = (1 << q) - 1;
    for (var k = 1; k < Q; k++) { var orb = [k], x = (2*k) % Q; while (x !== k && orb.length <= q) { orb.push(x); x = (2*x) % Q; } if (orb.length !== q) continue;
      var s = orb.slice().sort(function (a, b) { return a - b; }), ok = true;
      for (var i = 0; i < q && ok; i++) if (s[(s.indexOf(orb[i]) + p) % q] !== orb[(i + 1) % q]) ok = false;      // doubling moves each ray p places counterclockwise
      if (!ok) continue; var best = -1, bw = Infinity; for (i = 0; i < q; i++) { var w = (s[(i + 1) % q] - s[i] + Q) % Q; if (w < bw) { bw = w; best = i; } }
      return [s[best], s[(best + 1) % q]]; }
    return null;
  }
  var bits = function (k, len) { var t = k.toString(2); while (t.length < len) t = '0' + t; return t; };
  var cmp = function (a, b) { var l = a[0]*b[1], r = b[0]*a[1]; return l < r ? -1 : l > r ? 1 : 0; };   // fractions [num, den] of BigInt
  // value of the binary expansion .pre (per)^inf as a fraction
  var val = function (pre, per) { var L = BigInt(pre.length), P = (1n << BigInt(per.length)) - 1n; return [(pre ? BigInt('0b' + pre) : 0n)*P + BigInt('0b' + per), (1n << L)*P]; };
  // k = numerator of an angle of f over 2^n - 1, limb p/q. Returns { kind: 'tuned', m, inner: numerator over 2^m - 1 } | { kind: 'dyadic', l, k, e } | { kind: 'outside' } | { kind: 'undecided' }
  function active(n, k, p, q) {
    var R = rabbit(p, q); if (!R) return { kind: 'undecided', why: 'no rabbit' };
    var wm = bits(R[0], q), wp = bits(R[1], q), th = [BigInt(k), (1n << BigInt(n)) - 1n], W = function (e) { var s = ''; for (var i = 0; i < e.length; i++) s += e[i] === '0' ? wm : wp; return s; };
    if (cmp(th, val('', wm)) < 0 || cmp(th, val('', wp)) > 0) return { kind: 'outside', rabbit: R };
    var e = '';
    for (var it = 0; it < 4*n + 8; it++) {
      if (n % q === 0 && e.length === n/q) { var m = n/q; if (cmp(th, val('', W(e))) === 0) return { kind: 'tuned', m: m, inner: parseInt(e, 2), rabbit: R }; }
      var lo = val(W(e) + wm, wp), hi = val(W(e) + wp, wm);
      if (cmp(th, lo) <= 0) e += '0'; else if (cmp(th, hi) >= 0) e += '1'; else return { kind: 'dyadic', l: 2*parseInt(e || '0', 2) + 1, k: e.length + 1, e: e, rabbit: R };
    }
    return { kind: 'undecided', why: 'no gap found to depth ' + e.length, rabbit: R };
  }
  return { rabbit: rabbit, active: active };
})();
if (typeof module !== 'undefined') module.exports = RESC;
