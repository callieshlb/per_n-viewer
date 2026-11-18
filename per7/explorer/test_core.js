// E93 A2 (dynamics part) and core self-check: per7_core.js on explorer/data.json (atlas_points.py).
const CORE = require('./per7_core.js');
const D = require('./data.json');
CORE.init(D);
let bad = 0, maxP = 0, sheetMis = 0, sheetFail = 0;
const fates = {};
for (const c of D.centers) {
  const e = CORE.pEval(c.x, CORE.coeffs(c.y)); maxP = Math.max(maxP, CORE.cabs(e[0]));
  const M = CORE.mapData(c.x, c.y), f = CORE.fate(M, M.omega, 600, 1e-4, 1e-7, null);
  const key = f.kind + (f.period ? ':' + f.period : '') + (f.j !== undefined ? ':' + f.j : '');
  fates[key] = (fates[key] || 0) + 1;
  if (!(f.kind === 'cycle' && f.period === 1)) { bad++; console.log('centre fate', c.label, JSON.stringify(f)); }
  const s = CORE.sheetOf(c.y, c.x);
  if (s < 0) sheetFail++; else if (s !== c.sheet) { sheetMis++; console.log('centre sheet (double core) differs', c.label, s, c.sheet, c.y); }
}
console.log('A2 dynamics: centres', D.centers.length, '| fates', JSON.stringify(fates), '| not attracting period 1:', bad,
            '| max |P| at centres (double)', maxP.toExponential(2));
console.log('double-precision core vs 256-bit atlas sheets at centres: mismatches', sheetMis, '| continuation failures', sheetFail);
let pm = 0, pf = 0;
for (const p of D.punctures) { const s = CORE.sheetOf(p.y, p.x); if (s < 0) pf++; else if (s !== p.sheet) { pm++; console.log('puncture sheet differs', p.name, 'js', s, 'mp', p.sheet, 'y', p.y); } }
console.log('punctures: double core vs 256-bit sheets: mismatches', pm, '| continuation failures (expected near y = 0, 1)', pf, 'of', D.punctures.length);
// marked cycle consistency at a centre: f(a_6) = 0 and a_6 = f(a_5)
const c0 = D.centers[0], M0 = CORE.mapData(c0.x, c0.y);
console.log('marked points at centre', c0.label, M0.marked.length, '| |f(a_5) - a_6|', CORE.cabs(CORE.csub(CORE.fz(M0, M0.marked[5]), M0.marked[6])).toExponential(2),
            '| |f(a_6)|', CORE.cabs(CORE.fz(M0, M0.marked[6])).toExponential(2));
