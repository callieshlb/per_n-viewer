// E99 copy of E92 3-D step 2 (uses per6_core_mesh.js; logs min root separation): roots per vertex (atlas labels), sheet maps along every edge, groups at branch vertices.
const C = require('../explorer/per6_core_mesh.js'); C.init(require('../explorer/data.json'));
const M = require('./sphere_mesh.json');
const fs = require('fs');
const nV = M.y.length, EPS = 1e-4;
const nearest = (x, xs) => { let b = -1, bd = Infinity; xs.forEach((y, i) => { const d = C.chord(x, y); if (d < bd) { bd = d; b = i; } }); return [b, bd]; };
const isB = v => M.bp[v] >= 0;
let retried = 0, retryDisagree = 0;
const toU = y => { const r2 = y[0]*y[0]+y[1]*y[1]; return [2*y[0]/(r2+1), 2*y[1]/(r2+1), (r2-1)/(r2+1)]; };
const toY = u => [u[0]/(1-u[2]), u[1]/(1-u[2])];
function sphereTrack(xs, a, c) {          // follow the great-circle arc a -> c on the y-sphere (needed when |y| is large)
  const ua = toU(a), uc = toU(c);
  const dot = Math.max(-1, Math.min(1, ua[0]*uc[0]+ua[1]*uc[1]+ua[2]*uc[2])), om = Math.acos(dot), K = 16;
  let cur = xs, prev = a;
  for (let k = 1; k <= K && cur; k++) {
    const t = k/K, s1 = Math.sin((1-t)*om)/Math.sin(om), s2 = Math.sin(t*om)/Math.sin(om);
    const u = [s1*ua[0]+s2*uc[0], s1*ua[1]+s2*uc[1], s1*ua[2]+s2*uc[2]];
    const q = k === K ? c : toY(u);
    cur = trackStraight(cur, prev, q); prev = q;
  }
  return cur;
}
function trackEdge(xs, a, c) {
  if (Math.max(Math.hypot(a[0], a[1]), Math.hypot(c[0], c[1])) > 1) return sphereTrack(xs, a, c);
  return trackStraight(xs, a, c);
}
function trackStraight(xs, a, c) {          // straight; on failure via a midpoint pushed off by 20% either side
  let r = C.track(xs, a, c);
  if (r) return r;
  retried++;
  const v = [c[0]-a[0], c[1]-a[1]], mid = [(a[0]+c[0])/2, (a[1]+c[1])/2], res = [];
  for (const s of [0.2, -0.2]) {
    const m = [mid[0] - s*v[1], mid[1] + s*v[0]];
    const r1 = C.track(xs, a, m); const r2 = r1 && C.track(r1, m, c);
    if (r2) res.push(r2);
  }
  if (res.length === 2 && res[0].some((x, i) => C.chord(x, res[1][i]) > 1e-8)) retryDisagree++;
  return res.length ? res[0] : null;
}
const R = new Array(nV), group = {}, off = {};
let fails = 0, worst = 0;
const t0 = Date.now(); let minSep = Infinity, minSepV = -1;
for (let v = 0; v < nV; v++) {
  if (isB(v)) continue;
  R[v] = C.globalRoots(M.y[v]); if (!R[v]) { fails++; continue; }
  for (let i = 0; i < 5; i++) for (let j = i+1; j < 5; j++) { const d = C.chord(R[v][i], R[v][j]); if (d < minSep) { minSep = d; minSepV = v; } }
}
console.log(`normal vertices: min chordal root separation ${minSep.toExponential(3)} at vertex ${minSepV} y = ${M.y[minSepV]}; global-root failures ${fails}; ${(Date.now()-t0)/1000} s`);
// branch vertices: roots at an offset point, groups = cycles of the small loop around b
for (let v = 0; v < nV; v++) {
  if (!isB(v)) continue;
  const b = M.y[v];
  let found = null;
  for (const eps of [1e-4, 3e-4, 1e-3, 3e-3]) {
    for (const th of [0.3, 1.9, 3.4, 5.0]) {
      const o = [b[0] + eps*Math.cos(th), b[1] + eps*Math.sin(th)];
      const xs = C.globalRoots(o); if (!xs) continue;
      let cur = xs, prev = o;
      for (let i = 1; i <= 64 && cur; i++) { const q = [b[0] + eps*Math.cos(th + 2*Math.PI*i/64), b[1] + eps*Math.sin(th + 2*Math.PI*i/64)]; cur = C.track(cur, prev, q); prev = q; }
      if (!cur) continue;
      found = { o, xs, eps, th, cur }; break;
    }
    if (found) break;
  }
  if (!found) { fails++; console.log('branch vertex failed', v, b); continue; }
  off[v] = found;
  const xs = found.xs, cur = found.cur;
  const p = cur.map(x => nearest(x, xs)[0]);
  const g = new Array(5).fill(-1); let ng = 0;
  for (let i = 0; i < 5; i++) { if (g[i] >= 0) continue; let j = i; while (g[j] < 0) { g[j] = ng; j = p[j]; } ng++; }
  group[v] = { g, ng, perm: p };
}
// edges
const out = [];
let bb = 0, badPerm = 0;
for (const [p, q] of M.E) {
  if (isB(p) && isB(q)) { bb++; out.push(null); continue; }
  const [a, c] = isB(p) ? [q, p] : [p, q];           // a normal
  if (!isB(c)) {
    if (!R[a] || !R[c]) { fails++; out.push(null); continue; }
    const tr = trackEdge(R[a], M.y[a], M.y[c]);
    if (!tr) { fails++; out.push(null); continue; }
    const m = tr.map(x => { const [i, d] = nearest(x, R[c]); worst = Math.max(worst, d); return i; });
    if (new Set(m).size !== 5) badPerm++;
    out.push({ a, c, m });                            // root i at a -> root m[i] at c
  } else {
    const b = M.y[c], dir = [M.y[a][0] - b[0], M.y[a][1] - b[1]], L = Math.hypot(dir[0], dir[1]);
    if (!off[c] || !R[a]) { fails++; out.push(null); continue; }
    const E2 = off[c].eps, s = [b[0] + E2*dir[0]/L, b[1] + E2*dir[1]/L];
    const tr = trackEdge(R[a], M.y[a], s);            // to the offset point on this edge
    if (!tr) { fails++; out.push(null); continue; }
    // carry the offset-point roots of b round the small circle to s, to read off groups
    const o = off[c].o, a0 = Math.atan2(o[1]-b[1], o[0]-b[0]), a1 = Math.atan2(s[1]-b[1], s[0]-b[0]);
    let da = ((a1 - a0 + Math.PI) % (2*Math.PI) + 2*Math.PI) % (2*Math.PI) - Math.PI;
    let cur = off[c].xs, prev = o;
    for (let i = 1; i <= 32 && cur; i++) { const qq = [b[0] + E2*Math.cos(a0 + da*i/32), b[1] + E2*Math.sin(a0 + da*i/32)]; cur = C.track(cur, prev, qq); prev = qq; }
    if (!cur) { fails++; out.push(null); continue; }
    const m = tr.map(x => { const [i, d] = nearest(x, cur); worst = Math.max(worst, d); return group[c].g[i]; });
    out.push({ a, c, m, toGroup: true });             // root i at a -> group m[i] at branch vertex c
  }
}
console.log(`retried edges ${retried} (both sides disagree: ${retryDisagree}); vertices ${nV}, edges ${M.E.length}; tracking failures ${fails}; non-permutation edge maps ${badPerm}; ` +
            `edges between two branch vertices ${bb}; worst chordal match ${worst.toExponential(2)}; ${(Date.now()-t0)/1000} s`);
console.log('branch-vertex groups (number of surface points over each):', Object.values(group).map(g => g.ng).join(' '), '| offsets used', Object.values(off).map(o => o.eps).join(' '));
fs.writeFileSync('edge_perms.json', JSON.stringify({ R: R.map(r => r || null), group, edges: out }));
