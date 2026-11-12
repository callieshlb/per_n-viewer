// E96 shell: model menu, the parameter pane (iframe, one page per manifest entry), the dynamical plane, the status line.
// Parameter panes post {e96: 'param', model, p | x, y}; the map is built here by PM.mapFrom(family, n, message).
'use strict';
(async function () {
  const $ = s => document.getElementById(s);
  const man = await fetch('models/manifest.json?v=20261112204321').then(r => r.json());
  let cur = null, M = null, sel = null, frac = 0.55, FACTS = '';
  // ---- layout
  function layout() { const W = $('main').clientWidth, x = Math.round(W*frac);
    $('left').style.width = x + 'px'; $('split').style.left = x + 'px'; $('right').style.left = (x + 4) + 'px'; }
  window.addEventListener('resize', layout); layout();
  { let drag = false; $('split').addEventListener('pointerdown', e => { drag = true; $('split').setPointerCapture(e.pointerId); $('frame').style.pointerEvents = 'none'; });
    $('split').addEventListener('pointermove', e => { if (!drag) return; const r = $('main').getBoundingClientRect();
      frac = Math.min(0.85, Math.max(0.15, (e.clientX - r.left)/r.width)); layout(); });
    $('split').addEventListener('pointerup', () => { drag = false; $('frame').style.pointerEvents = ''; }); }
  // ---- model menu
  for (const m of man.models) {
    const b = document.createElement('button'); b.innerHTML = m.short; b.id = 'b_' + m.id; b.title = m.titleText || m.id;
    if (m.placeholder) { b.disabled = true; b.title = m.placeholder; } else b.onclick = () => choose(m.id);
    $('models').appendChild(b);
  }
  function choose(id) {
    cur = man.models.find(m => m.id === id); if (!cur) return;
    for (const m of man.models) $('b_' + m.id).classList.toggle('on', m.id === id);
    const I = cur.info || {};
    $('kind').textContent = ({ plane: 'parameter plane', torus: 'flat torus', hyperbolic: 'hyperbolic model' }[cur.kind] || '') +
      ` · genus ${I.genus} · ${I.centers} polynomial centers · ${I.punctures} punctures · ${I.circles} real circle${I.circles === 1 ? '' : 's'}` + (I.sheets ? ` · ${I.sheets} sheets over the y-line` : '');
    hideCard(); $('dover').style.display = 'none'; $('vin').style.display = 'none'; TREEP = null;
    $('frame').src = cur.page + (cur.page.includes('?') ? '&' : '?') + 'v=20261112204321'; history.replaceState(null, '', '#' + id);
    M = null; sel = null; DRAYS = []; RAYA = []; $('dhint').style.display = ''; dv.redraw();
    FACTS = `Polynomials z² + c are labeled by the even numerator k of their characteristic angles k/(2^n − 1).\n${cur.titleText}: genus ${I.genus}, ${I.punctures} puncture${I.punctures === 1 ? '' : 's'}, ${I.centers} polynomial center${I.centers === 1 ? '' : 's'}, ${I.circles} real circle${I.circles === 1 ? '' : 's'}${I.sheets ? ', ' + I.sheets + ' sheets over the y-line' : ''}.\n${I.tiling || ''}`;
    $('status').textContent = FACTS + `\n${cur.kind === 'hyperbolic' ? 'The model loads in the left pane (tens of MB; the first load can take a while), then click a point.' : 'Click a parameter.'}`;
  }
  // ---- dynamical plane
  const opt = { mk: true, ml: true, orb: false, orbN: 60, omode: false, ok: 20, inv: false, tree: true }; let ORB0 = null, TREE = null;
  // polynomial view (Caroline 2026-11-10): w = T(z), T Moebius with T(a0) = 0, T(omega) = inf, scaled so that T f T^-1 =
  // (w^2 + B)/(C w^2 + 1) (= w^2 + c at a center, where omega is fixed). T = s T0, T0(a0) = 0, T0(a1) = 1, T0(omega) = inf.
  let TV = null;     // { T: z -> w, Ti: w -> z, mat: inverse as 2x2 }
  function polyView() {
    TV = null; if (!M) return;
    if (SVIEW) { TV = makeTV(SVIEW.T, SVIEW.Ti, [1, 0]); return; }
    if (!opt.inv) return;
    const T0 = PM.moebius3([M.marked[0], M.marked[1], M.omega], [[0, 0], [1, 0], null]), T0i = PM.moebius3([[0, 0], [1, 0], null], [M.marked[0], M.marked[1], M.omega]);
    const g0 = w => T0(PM.fz(M, T0i(w))), b0 = g0([0, 0]), aI = g0(null), g1 = g0([1, 0]);
    let A; if (aI === null) A = PM.csub(g1, b0); else { const C = PM.cdiv(PM.csub(g1, b0), PM.csub(aI, g1)); A = PM.cmul(aI, C); }
    const s = A, si = PM.cdiv([1, 0], s);
    const T = z => { const w = T0(z); return w === null ? null : PM.cmul(s, w); }, Ti = w => T0i(w === null ? null : PM.cmul(si, w));
    TV = makeTV(T, Ti, s);
  }
  let SVIEW = null;                                    // sphere coordinate (puncture explorer), overrides the polynomial view
  function makeTV(T, Ti, s) {
    // the inverse as a matrix for the shader: images of 0, 1, inf under Ti give (mA w + mB)/(mC w + mD)
    const z0 = Ti([0, 0]), z1 = Ti([1, 0]), zi = Ti(null), hv = z => z === null ? [[1, 0], [0, 0]] : [z, [1, 0]];
    // matrix with columns (zi-vector * a, z0-vector * b), a + b normalized so that column sum gives z1
    const P = hv(zi), Q = hv(z0), R = hv(z1), det = PM.csub(PM.cmul(P[0], Q[1]), PM.cmul(P[1], Q[0]));
    const a = PM.cdiv(PM.csub(PM.cmul(R[0], Q[1]), PM.cmul(R[1], Q[0])), det), b = PM.cdiv(PM.csub(PM.cmul(P[0], R[1]), PM.cmul(P[1], R[0])), det);
    return { T, Ti, s, m: [PM.cmul(P[0], a), PM.cmul(Q[0], b), PM.cmul(P[1], a), PM.cmul(Q[1], b)] };
  }
  const J = z => !TV ? z : TV.T(z);                                                       // plane point of z
  const SC = z => { const w = J(z); return w === null ? null : dv.toScreen(w); };
  const dv = new PlaneView($('dyn'), {
    ctr: [0, 0], upp: 0.01, maxIt: 500, delta: 1e-3,
    setUniforms: (gl, U) => {
      if (!M) { gl.uniform1i(U('fam'), 1); gl.uniform1i(U('maxIt'), 0); return; }
      gl.uniform1i(U('fam'), 0); gl.uniform1i(U('invz'), TV ? 1 : 0); if (TV) { gl.uniform2f(U('mA'), ...TV.m[0]); gl.uniform2f(U('mB'), ...TV.m[1]); gl.uniform2f(U('mC'), ...TV.m[2]); gl.uniform2f(U('mD'), ...TV.m[3]); } GLSL.setPalette(gl, U, cur && cur.palette, cur && cur.flat);
      ['A', 'B', 'C', 'D', 'E', 'F'].forEach(k => gl.uniform2f(U('u' + k), M[k][0], M[k][1]));
      const mk = new Float32Array(24), mi = new Int32Array(12);
      M.marked.forEach((a, j) => { if (a === null) mi[j] = 1; else { mk[2*j] = a[0]; mk[2*j + 1] = a[1]; } });
      gl.uniform2fv(U('uMk'), mk); gl.uniform1iv(U('uMkI'), mi); gl.uniform1i(U('nM'), M.marked.length);
    },
    drawOverlay: (g, W, H) => {
      if (!M) { g.fillStyle = '#111'; g.fillRect(0, 0, W, H); return; }
      if (opt.orb) {
        const o = PM.orbit(M, M.omega, opt.orbN); g.strokeStyle = 'rgba(255,230,80,.75)'; g.lineWidth = 1; g.beginPath(); let on = false;
        for (const z of o) { if (z === null) { on = false; continue; } const s = SC(z); if (!s) { on = false; continue; } if (on) g.lineTo(...s); else g.moveTo(...s); on = true; }
        g.stroke(); g.fillStyle = '#ffe650';
        for (const z of o) if (z) { const s = SC(z); if (s) g.fillRect(s[0] - 1.5, s[1] - 1.5, 3, 3); }
      }
      const RC = ['#ffe650', '#7cf', '#f9a', '#9f9', '#fc6', '#c9f', '#6ff', '#f66'];
      DRAYS.forEach((r, i) => { const col = RC[i % RC.length];
        g.beginPath(); r.pts.forEach((q, j) => { const s = SC(q); if (s) (j ? g.lineTo(...s) : g.moveTo(...s)); });
        g.strokeStyle = 'rgba(0,0,0,.7)'; g.lineWidth = 3.5; g.stroke(); g.strokeStyle = col; g.lineWidth = 1.8; g.stroke();
        const s = SC(r.landing) || [-99, -99]; g.beginPath(); g.arc(s[0], s[1], 4, 0, 2*Math.PI); g.fillStyle = col; g.fill();
        g.font = '600 12px system-ui'; g.lineWidth = 3; g.strokeStyle = 'rgba(0,0,0,.85)'; g.strokeText(r.a.text, s[0] + 6, s[1] - 8); g.fillStyle = col; g.fillText(r.a.text, s[0] + 6, s[1] - 8); });
      if (opt.omode && ORB0) {                     // orbit mode: z0 -> z1 -> ... -> zk, straight segments
        const o = PM.orbit(M, ORB0, opt.ok); g.lineWidth = 1.6; let prev = null;
        o.forEach((z, i) => { const s = SC(z); if (!s) { prev = null; return; }
          if (prev) { g.beginPath(); g.moveTo(...prev); g.lineTo(...s); g.strokeStyle = 'rgba(0,0,0,.6)'; g.lineWidth = 3.2; g.stroke(); g.strokeStyle = '#ffffff'; g.lineWidth = 1.4; g.stroke(); }
          g.beginPath(); g.arc(s[0], s[1], i === 0 ? 4.5 : 2.6, 0, 2*Math.PI); g.fillStyle = i === 0 ? '#ff3b6b' : '#ffffff'; g.fill();
          if (opt.ml && (i === 0 || i === o.length - 1 || o.length <= 12)) { g.font = '600 11px system-ui'; g.lineWidth = 3; g.strokeStyle = 'rgba(0,0,0,.85)'; g.strokeText('z' + i, s[0] + 5, s[1] - 7); g.fillStyle = '#fff'; g.fillText('z' + i, s[0] + 5, s[1] - 7); }
          prev = s; });
      }
      if (opt.tree && TREE) {                    // Hubbard tree of the clicked center, moved into this chart
        for (const arc of TREE) { g.beginPath(); let on = false, last = null;
          for (const z of arc) { const s = SC(z); if (!s || (last && Math.hypot(s[0] - last[0], s[1] - last[1]) > W)) { on = false; last = s; continue; }
            if (on) g.lineTo(...s); else g.moveTo(...s); on = true; last = s; }
          g.strokeStyle = 'rgba(0,0,0,.8)'; g.lineWidth = 4; g.stroke(); g.strokeStyle = '#ffffff'; g.lineWidth = 2; g.stroke(); }
      }
      g.font = '600 12px system-ui'; g.textBaseline = 'middle';
      const dot = (z, fill, lab) => { const s = SC(z); if (!s) return; if (s[0] < -20 || s[1] < -20 || s[0] > W + 20 || s[1] > H + 20) return;
        g.beginPath(); g.arc(s[0], s[1], 5, 0, 2*Math.PI); g.fillStyle = fill; g.fill(); g.strokeStyle = '#141414'; g.lineWidth = 1.2; g.stroke();
        if (opt.ml && lab) { g.lineWidth = 3; g.strokeStyle = 'rgba(0,0,0,.85)'; g.strokeText(lab, s[0] + 8, s[1]); g.fillStyle = fill; g.fillText(lab, s[0] + 8, s[1]); } };
      if (opt.mk) M.marked.forEach((a, j) => { if (a || opt.inv) dot(a, '#ffffff', 'a' + j); });
      dot(M.omega, '#ffe650', 'ω');
    },
    onPick: p => {
      if (!M) return; p = TV ? (TV.Ti(p) || [1e300, 0]) : p; if (opt.omode) { ORB0 = p; dv.redraw(); }
      const f = PM.fate(M, p, dv.maxIt, dv.delta, 1e-7);
      $('status').textContent = statusText() + `\nclicked z = ${PM.fmt(p)}: ${PM.fateText(f)}`;
    } });
  function fit() {
    if (!M) return; if (SVIEW) { const [W0, H0] = dv.size(); dv.ctr = [0, 0]; dv.upp = 7/Math.max(50, Math.min(W0, H0)); dv.redraw(); return; }
    if (!M) return; const pts = M.marked.concat([M.omega, [0, 0]]).map(J).filter(a => a && PM.cabs(a) < 1e3);
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (const z of pts) { x0 = Math.min(x0, z[0]); x1 = Math.max(x1, z[0]); y0 = Math.min(y0, z[1]); y1 = Math.max(y1, z[1]); }
    const [W, H] = dv.size(), span = Math.max(x1 - x0, y1 - y0, 1.5)*1.8;
    dv.ctr = [(x0 + x1)/2, (y0 + y1)/2]; dv.upp = span/Math.max(50, Math.min(W, H)); dv.redraw();
  }
  function statusText() {
    if (!M) return '';
    const par = (cur.kind !== 'plane') && sel.x && sel.y ? `x = a${M.n - 1} = ${PM.fmt(sel.x)}   y = a3 = ${PM.fmt(sel.y)}${sel.p ? '   (' + cur.param + ' = ' + PM.fmt(sel.p) + ')' : ''}` : `${cur.param} = ${PM.fmt(sel.p)}`;
    const f = PM.fate(M, M.omega, dv.maxIt, dv.delta, 1e-7);
    let comp = '';
    if (f.kind === 'capture' && cur.kind !== 'hyperbolic' && sel.p) { const C = PM.captureCentre(cur.family, cur.n, sel.p, f);
      comp = C ? `capture component: at its center (${cur.param} = ${PM.fmt(C.center, 6)}) ${C.preperiod === 0 ? 'the free critical point IS a' + C.j : 'the free critical point lands on a' + C.j + ' after ' + C.preperiod + ' step' + (C.preperiod === 1 ? '' : 's')}: preperiod ${C.preperiod}, period ${cur.n} (phase ${C.phase})`
               : `captured (phase ${f.phase}); center not found by Newton from here`; }
    else if (f.kind === 'capture') comp = `captured by the marked cycle (phase ${f.phase})`;
    if (f.kind === 'capture' && cur.family === 'quad') comp = 'escapes to ∞ (outside the Mandelbrot set)';
    else if (f.kind === 'cycle') comp = f.period === 1 ? 'hyperbolic component with an attracting fixed point (polynomial-like: the free critical point is attracted to a fixed point)' : `hyperbolic component of period ${f.period} (the free critical point is attracted to a ${f.period}-cycle)`;
    else comp = 'undecided: bifurcation locus (or beyond the iteration count)';
    VIEWM = f.kind === 'cycle' && f.period <= 9 ? f.period : 0;   // 9: Per_9 wired in 2026-11-12
    $('vin').style.display = VIEWM ? '' : 'none';
    if (VIEWM) $('vin').textContent = VIEWM === cur.n ? `View the involution image in Per_${VIEWM}  (v)` : `View this component's center in Per_${VIEWM}  (v)`;
    return `${cur.titleText}${sel.label ? '  ·  ' + sel.label : ''}   ${par}\n` +
      `f(z) = ((${PM.fmt(M.A, 5)}) z² + (${PM.fmt(M.B, 5)}) z + (${PM.fmt(M.C, 5)})) / ((${PM.fmt(M.D, 5)}) z² + (${PM.fmt(M.E, 5)}) z + (${PM.fmt(M.F, 5)}))\n` +
      `marked cycle: ${M.marked.map((a, j) => 'a' + j + ' = ' + PM.fmt(a, 5)).join(',  ')}   (closes to ${PM.closure(M).toExponential(1)})\n` +
      (sel.amb ? `NOTE: two sheets of the curve are close here (nearest / second-nearest root distance from the chart value = ${sel.amb.toFixed(2)}); the sheet shown is the nearest one and may not be the one of the picture.\n` : '') +
      (sel.raw ? `NOTE: x here is the disk viewer's own chart value, not re-solved on the curve (the re-solve did not converge, or moved x by more than the guard allows); the map may not be on the curve at all: check "closes to" above.\n` : '') +
      `free critical point ω = ${PM.fmt(M.omega, 6)}: ${comp}   |   iterates ${dv.maxIt}, δ ${dv.delta.toExponential(0)}`;
  }
  // ---- Per_1: the clicked hyperbolic component (REGISTRATION addendum 2026-11-12 (s)): its center, the two parameter rays landing at its root, and a
  // button to its center with Hubbard tree and polynomial card. Angles: poly_info (period <= 9) or the tables poly/centres_p{10..14}.json.
  let COMP = null; const CENT = {};
  function loadCen(m) { if (m < 10 || m > 14) return Promise.resolve(null); if (!CENT[m]) CENT[m] = fetch(`poly/centres_p${m}.json`).then(r => r.ok ? r.json() : null).catch(() => null); return CENT[m]; }
  async function compInfo() {
    $('cen').style.display = 'none'; COMP = null; if (!M || !cur || cur.family !== 'quad' || !sel.p) return;
    const Mc = M, f = PM.fate(M, M.omega, dv.maxIt, dv.delta, 1e-7); if (f.kind !== 'cycle') return; const m = f.period;
    const C = PM.cycleCentre('quad', 1, { p: sel.p }, m);
    if (!C || !C.p) { $('status').textContent += `\ncomponent of period ${m}: center not found from here (click nearer its middle)`; return; }
    const c0 = C.p, Q = 2**m - 1; let ent = null, why = '';
    if (m === 1) ent = { ks: [0, 0], label: 0, c: [0, 0] };
    else if (m <= 9) { const D = await loadPoly(m); if (D) { let bd = Infinity; for (const P of D) { const d = Math.hypot(P.c[0] - c0[0], P.c[1] - c0[1]); if (d < bd) { bd = d; ent = P; } } if (!(bd < 1e-8)) { ent = null; why = `no tabulated center within 1e-8 (nearest ${bd.toExponential(1)})`; } } else why = 'table not available'; }
    else if (m <= 14) { const T = await loadCen(m); if (T) { let bd = Infinity, bi = -1; for (let i = 0; i < T.length; i += 4) { const d = Math.hypot(T[i + 2] - c0[0], T[i + 3] - c0[1]); if (d < bd) { bd = d; bi = i; } }
        if (bd < 1e-8) ent = { ks: [T[bi], T[bi + 1]], label: T[bi] % 2 === 0 ? T[bi] : T[bi + 1], c: [T[bi + 2], T[bi + 3]] }; else why = `no tabulated center within 1e-8 (nearest ${bd.toExponential(1)})`; } else why = 'table not available'; }
    else why = 'angles are tabulated up to period 14 only';
    if (M !== Mc) return;
    let t = `\ncomponent of period ${m}: center c = ${PM.fmt(c0, 12)}`;
    if (ent && m >= 2) { t += `; parameter rays landing at its root: ${ent.ks[0]}/${Q} and ${ent.ks[1]}/${Q}`;
      for (const u of ent.tuning || []) { const Qd = 2**u.d - 1, Qm = 2**u.inner_period - 1; t += `\n  tuning: (${u.base_ks[0]}/${Qd}, ${u.base_ks[1]}/${Qd})${u.base_name ? ' ' + u.base_name : ''} ⊥ (${u.inner_ks[0]}/${Qm}, ${u.inner_ks[1]}/${Qm})${u.inner_name ? ' ' + u.inner_name : ''}`; } }
    else if (m >= 2) t += `; root rays: ${why}`;
    $('status').textContent += t; COMP = { m, c: c0, ent };
    try { $('frame').contentWindow.postMessage({ e96: 'compRays', n: m, ks: ent && m >= 2 ? ent.ks : null, c: c0 }, '*'); } catch (e) {}
    $('cen').textContent = `Go to this component's center${m >= 2 && m <= 9 && ent ? ': Hubbard tree and polynomial card' : ''}  (h)`; $('cen').style.display = '';
  }
  function goCenter() {
    if (!COMP) return; const { m, c, ent } = COMP, Q = 2**m - 1; let marker = null;
    if (ent && m >= 2) marker = { kind: 'centre', label: String(ent.label), n: m, ks: ent.ks, title: `center ${ent.label} of period ${m} (angles ${ent.ks[0]}/${Q}, ${ent.ks[1]}/${Q})`,
      lines: [`z^2 + c,  c = ${ent.c_str ? ent.c_str[0] + (ent.c_str[1].startsWith('-') ? ' - ' : ' + ') + ent.c_str[1].replace('-', '') + 'i' : PM.fmt(c, 14)}`, `angles ${ent.ks[0]}/${Q}, ${ent.ks[1]}/${Q}`]
        .concat(ent.kneading ? [`kneading ${ent.kneading}`, `internal address ${ent.address.join(' -> ')}`] : [])
        .concat((ent.tuning || []).map(u => `TUNING: (${u.base_ks[0]}/${2**u.d - 1}, ${u.base_ks[1]}/${2**u.d - 1})${u.base_name ? ' ' + u.base_name : ''} ⊥ (${u.inner_ks[0]}/${2**u.inner_period - 1}, ${u.inner_ks[1]}/${2**u.inner_period - 1})${u.inner_name ? ' ' + u.inner_name : ''}`)) };
    handleParam({ e96: 'param', p: c, marker, label: `center of the period-${m} component` });
    try { $('frame').contentWindow.postMessage({ e96: 'compRays', n: m, ks: ent && m >= 2 ? ent.ks : null, c, go: true }, '*'); } catch (e) {}
  }
  $('cen').onclick = goCenter;
  window.E96TEST_tree = () => TREE ? (TREE.arcs ? TREE.arcs.length : (TREE.length !== undefined ? TREE.length : 'set')) : 0;
  window.addEventListener('keydown', e => { if (e.key === 'h' && COMP && e.target.tagName !== 'INPUT' && e.target.tagName !== 'SELECT' && !e.metaKey && !e.ctrlKey) { e.preventDefault(); goCenter(); } });
  window.addEventListener('message', ev => { if (ev.data && ev.data.e96 === 'key' && ev.data.key === 'h') goCenter(); });
  // ---- "View in Per_m" (REGISTRATION addendum 2026-11-10 (c)): the center of the clicked component, as a point of Per_m (m = period of omega's cycle)
  let VIEWM = 0, VSEND = null;
  function viewM() {
    if (!M || !VIEWM || !cur) return; const m = VIEWM, n0 = cur.n, from = cur.titleText, hyp = cur.kind === 'hyperbolic';
    const C = PM.cycleCentre(hyp ? 'cor' : cur.family, n0, hyp ? { x: sel.x, y: sel.y } : { p: sel.p }, m);
    if (!C) { $('status').textContent += `\nView in Per_${m}: the center of this component was not found from here (click nearer its middle).`; return; }
    const T = PM.transfer(C.M, m); if (!T) { $('status').textContent += `\nView in Per_${m}: transfer failed.`; return; }
    const label = `≅ ${from.replace('(0)', '')}${sel.label ? ' ' + sel.label : ''}`;
    if (T.n !== n0) choose('per' + T.n);
    let msg;
    if (cur.kind === 'plane') msg = { p: T.p };
    else if (cur.kind === 'torus') { const L = PM.torusLocate(T.x, T.y); if (!L) { $('status').textContent += '\ntorus point not found'; return; } msg = { p: L.u }; }
    else msg = { x: T.x, y: T.y };
    const note = `${m === n0 ? 'Involution image' : 'The same rational map'}: center of the component clicked in ${from}` +
      (C.p ? ` (there: ${PM.fmt(C.p, 8)})` : ` (there: x = ${PM.fmt(C.x, 8)}, y = ${PM.fmt(C.y, 8)})`) + `, now with the other critical cycle marked.`;
    handleParam(Object.assign({ e96: 'param', label, keepCard: true }, msg)); $('status').textContent += '\n' + note;
    clearInterval(VSEND && VSEND.t); let tries = 0;
    const post = () => { tries++; if (tries > 150) { clearInterval(VSEND.t); VSEND = null; return; }
      try { $('frame').contentWindow.postMessage(Object.assign({ e96: 'mark', label, go: true }, msg), '*'); } catch (e) {} };
    VSEND = { label, t: setInterval(post, 1200) }; post();
  }
  window.addEventListener('message', ev => { const d = ev.data; if (d && d.e96 === 'markDone' && VSEND && d.label === VSEND.label) { clearInterval(VSEND.t); VSEND = null;
    if (!d.ok) $('status').textContent += `\nnot located in the disk (best chart residual ${d.res != null ? d.res.toExponential(1) : '?'})`; } });
  // hotkeys for the toggles (Caroline 2026-11-10). Dynamical pane (focus here): t tree, i polynomial view, b bands, l labels, k marked cycle, +/- iterates.
  // Parameter-pane toggles typed here are forwarded to the pane: r real circles, g G_n edges, e violet arrows, p period colors.
  const DKEY = { t: 'dtree', i: 'dinv', b: 'dbands', l: 'dml', k: 'dmk' };
  function dynKey(key) {
    if (DKEY[key]) { const el = $(DKEY[key]); el.checked = !el.checked; el.dispatchEvent(new Event('change')); return true; }
    if (key === '+' || key === '=' || key === '-') { const el = $('dit'); el.value = Math.max(+el.min, Math.min(+el.max, +el.value + (key === '-' ? -1 : 1)*100)); el.dispatchEvent(new Event('input')); return true; }
    return false;
  }
  window.addEventListener('keydown', e => { if (e.target.tagName === 'INPUT' && e.target.type !== 'checkbox' && e.target.type !== 'range' || e.target.tagName === 'SELECT' || e.metaKey || e.ctrlKey || e.altKey) return;
    if (dynKey(e.key)) { e.preventDefault(); return; }
    if ('rgepja'.includes(e.key) && e.key.length === 1) { try { $('frame').contentWindow.postMessage({ e96: 'key', key: e.key }, '*'); } catch (err) {} } });
  window.addEventListener('message', ev => { const d = ev.data; if (d && d.e96 === 'key' && (d.key === 't' || d.key === 'i')) dynKey(d.key); });
  $('vin').onclick = viewM;
  window.addEventListener('keydown', e => { if (e.key === 'v' && e.target.tagName !== 'INPUT' && e.target.tagName !== 'SELECT') { e.preventDefault(); viewM(); } });
  window.addEventListener('message', ev => { if (ev.data && ev.data.e96 === 'key' && ev.data.key === 'v') viewM(); });
  let RAYA = [], DRAYS = [];
  function dynRays() { DRAYS = (M && cur && cur.family === 'quad') ? RAYA.map(t => { const a = PM.parseAngle(t); return { a, ...PM.rayWithLanding('dyn', a, M.C) }; }) : []; }
  window.addEventListener('message', ev => {
    const d0 = ev.data; if (d0 && d0.e96 === 'rays') { RAYA = d0.angles || []; dynRays(); dv.redraw(); return; }
    if (d0 && d0.e96 === 'param' && cur) handleParam(d0);
  });
  function handleParam(d) {
    const q = { p: d.p, x: d.x, y: d.y };
    const fam = cur.kind === 'hyperbolic' ? 'cor' : cur.family;
    const m = PM.mapFrom(fam, cur.n, q);
    const finite = m && [m.A, m.B, m.C, m.D, m.E, m.F, m.omega].every(z => isFinite(z[0]) && isFinite(z[1]));
    if (!finite) { $('status').textContent = `${cur.titleText}: the map degenerates at this point (puncture); the dynamical plane keeps the previous map`;
      const mk0 = parseMarker(d); if (mk0) showMarker(mk0); TREE = null; dv.redraw(); return; }
    SVIEW = null; if (SPHR) { SPHR = null; $('sphl').textContent = 'new map: press find (u) again'; }
    M = m; sel = { p: d.p || null, x: d.x || (m.x || null), y: d.y || (m.y || null), label: d.label, raw: !!d.raw, amb: d.amb || 0 };
    if (cur.family === 'quad' || cur.family === 'per2') { sel.x = null; sel.y = null; }
    $('dhint').style.display = 'none'; dynRays(); polyView(); fit(); $('status').textContent = statusText(); ORB0 = null;
    TREE = null; const mk = parseMarker(d); if (mk) { showMarker(mk); if (mk.kind === 'centre') showTree(mk); } else if (!d.keepCard) hideCard();
    if (!mk || mk.kind !== 'centre') setTimeout(followTree, 30);
    if (cur.family === 'quad' && !mk) compInfo(); else { $('cen').style.display = 'none'; COMP = null; }
  }
  window.SHELL = { get M() { return M; }, get sel() { return sel; }, dv, choose, viewM, handleParam };
  // ---- controls
  $('dit').oninput = e => { dv.maxIt = +e.target.value; $('ditv').textContent = dv.maxIt; dv.redraw(); if (M) $('status').textContent = statusText(); };
  $('ddel').oninput = e => { dv.delta = Math.pow(10, +e.target.value); $('ddelv').textContent = dv.delta.toExponential(0); dv.redraw(); };
  $('dmk').onchange = e => { opt.mk = e.target.checked; dv.redraw(); };
  $('dml').onchange = e => { opt.ml = e.target.checked; dv.redraw(); };
  $('dorb').onchange = e => { opt.orb = e.target.checked; dv.redraw(); };
  $('dorbn').oninput = e => { opt.orbN = Math.max(1, Math.min(2000, +e.target.value || 1)); dv.redraw(); };
  $('dhome').onclick = fit;
  // f of a mating: a label (even numerator k), or either characteristic angle p/q of a period-n center (Caroline 2026-11-10)
  async function resolveF(text, n) {
    const D = n >= 2 ? await loadPoly(n) : null; if (!D) return null; const t = String(text).trim(), Q = 2**n - 1; let P = null, m;
    if ((m = t.match(/^(\d+)\s*\/\s*(\d+)$/))) { const k = +m[1]*Q/+m[2]; if (Number.isInteger(k)) P = D.find(q => q.ks.includes(k)); }
    else if (/^\d+$/.test(t)) P = D.find(q => String(q.label) === t) || D.find(q => q.ks.includes(+t));
    return P || null;
  }
  async function mateFromText(ft, gt, K, useMis) {
    if (!cur) return; const n = cur.n;
    if (n < 2) { $('status').textContent = 'Mating: choose Per_n with n ≥ 2 (in Per_1, z² ⊔ g = g).'; return; }
    const P = await resolveF(ft, n);
    if (!P) { $('status').textContent = `Mating: "${ft}" is not a period-${n} center (give its label k, or either angle k/${2**n - 1})`; return; }
    const cf = [Number(P.c_str[0]), Number(P.c_str[1])]; let vtxt = '';
    if (K > 0) vtxt = await doVein(cf, n, P.label, gt, K, useMis !== false);
    doMate(cf, n, P.label, gt, null);
    if (vtxt) { const o = MATEOUT && MATEOUT.out; if (MATEOUT) MATEOUT.base += '\n' + vtxt; $('status').textContent += '\n' + vtxt; }
  }
  // mating vein (REGISTRATION addendum 2026-11-10 (b), M7): f ⊔ g_i for the centers g_i (period ≤ 8) on the arc from 0 to g
  let MISTAB = null;
  async function doVein(cf, n, lab, text, K, useMis) {
    const a = PM.parseAngle(text); if (!a) return '';
    const polys = {}; for (let m = 2; m <= 8; m++) polys[m] = await loadPoly(m);
    if (useMis && !MISTAB) MISTAB = await fetch('models/misiurewicz.json').then(r => r.ok ? r.json() : []).catch(() => []);
    const mis = useMis ? MISTAB : null;
    let V = [], P = 8; for (const Pt of [8, 10, 12, 14, 16, 18, 20]) { P = Pt; V = MATE.veinList(a, polys, K, Pt, mis); if (V.length >= K) break; }       // periods up to 20 (linear-time Lavaurs + spider), M9, M10
    const W = $('frame').contentWindow; let ok = 0; const obs = [];
    for (const v of V) { const cg = MATE.centerOf(v), nm = v.mis ? v.name : `${v.k0}/${v.Q}`; if (!cg) { obs.push(nm + ' (center not found)'); continue; }
      const m = MATE.mate(cf, cg); if (!m.ok) { obs.push(nm); continue; }
      const ch = MATE.toChart(m); if (!ch) continue; let msg;
      if (cur.kind === 'plane') msg = { p: ch.p };
      else if (cur.kind === 'torus') { const L = PM.torusLocate(ch.x, ch.y); if (!L) { obs.push(nm + ' (not located)'); continue; } msg = { p: L.u }; }
      else msg = { x: ch.x, y: ch.y };
      W.postMessage(Object.assign({ e96: 'mark', label: nm, vein: true, quiet: true }, msg), '*'); ok++; }
    const nm_ = V.filter(v => v.mis).length;
    return `vein of ${a.text}: ${V.length} points on the arc from z² to g (${V.length - nm_} centers of period ≤ ${P}${useMis ? `, ${nm_} Misiurewicz points with preperiod + period ≤ ${Math.min(P, 12)}` : ''}; asked K = ${K}${V.length < K ? `; there are no more up to period 20${useMis ? ' / preperiod + period 12' : ''}` : ''}): ${ok} mated (small yellow dots, labeled by the lower angle of g_i; an even denominator means preperiodic)` +
      (obs.length ? `; ${obs.length} obstructed or not placed: ${obs.slice(0, 12).join(', ')}${obs.length > 12 ? ', …' : ''}` : '');
  }
  function keyM() {                                   // 'm': mating prompt: f angle, g angle, Enter
    const dl = $('mdlg'); dl.style.display = 'block'; $('mdn').textContent = cur ? cur.n : 'n';
    const cardLab = (document.querySelector('#dct b') || {}).textContent || '', m = cardLab.match(/^cent(?:er|re) (\d+)/);
    $('mdf').value = ''; $('mdg').value = ''; $('mdf').placeholder = m ? m[1] + ' (this card)' : 'angle or label'; $('mdf').dataset.def = m ? m[1] : '';
    $('mdf').focus();          // always: type f, Tab (or Enter), type g, Enter (Caroline 2026-11-10)
  }
  $('mdf').addEventListener('keydown', e => { if (e.key === 'Enter') { $('mdg').focus(); $('mdg').select(); } });
  $('mdg').addEventListener('keydown', e => { if (e.key === 'Enter') { $('mdlg').style.display = 'none'; mateFromText($('mdf').value.trim() || $('mdf').dataset.def || '', $('mdg').value, $('mdv').checked ? Math.max(0, Math.min(400, +$('mdk').value || 0)) : 0, $('mdm').checked); } });
  $('mdk').addEventListener('keydown', e => { if (e.key === 'Enter') { $('mdg').focus(); } });
  $('mdc').onclick = () => { $('frame').contentWindow.postMessage({ e96: 'clearMarks' }, '*'); };
  $('mdlg').addEventListener('keydown', e => { if (e.key === 'Escape') $('mdlg').style.display = 'none'; });
  $('mdx').onclick = () => { $('mdlg').style.display = 'none'; };
  window.addEventListener('keydown', e => { if (e.key === 'm' && e.target.tagName !== 'INPUT' && e.target.tagName !== 'SELECT') { e.preventDefault(); keyM(); } });
  window.addEventListener('message', ev => { if (ev.data && ev.data.e96 === 'key' && ev.data.key === 'm') { window.focus(); keyM(); } });
  window.addEventListener('message', ev => {                       // mating box in a model's side panel
    const d = ev.data; if (!d || d.e96 !== 'mateReq' || !cur) return; mateFromText(d.label, d.angle, Math.max(0, Math.min(400, +d.K || 0)), d.mis !== false);
  });
  window.addEventListener('keydown', e => {          // 'o': orbit of the point under the cursor (dynamical plane)
    if (e.target.tagName === 'INPUT' || !M || e.key !== 'o' || !dv.mouse) return;
    ORB0 = TV ? (TV.Ti(dv.mouse) || [1e300, 0]) : dv.mouse; opt.omode = true; $('omode').checked = true; dv.redraw(); });
  $('dbands').onchange = e => { dv.bands = e.target.checked; dv.redraw(); };
  $('dinv').onchange = e => { opt.inv = e.target.checked; polyView(); fit(); };
  $('dtree').onchange = e => { opt.tree = e.target.checked; dv.redraw(); };
  $('omode').onchange = e => { opt.omode = e.target.checked; if (!opt.omode) ORB0 = null; dv.redraw(); };
  $('ok').oninput = e => { opt.ok = +e.target.value; $('okv').textContent = opt.ok; dv.redraw(); };
  // ---- marker card in the dynamical window: census data, polynomial card (Hubbard tree, rays), Thurston obstruction
  function hideCard() { $('dcard').style.display = 'none'; }
  $('dcx').onclick = hideCard;
  $('dovx').onclick = () => { $('dover').style.display = 'none'; $('dovf').src = 'about:blank'; };
  $('dovm').onclick = () => { const m = $('dover').classList.toggle('min'); $('dovm').textContent = m ? 'Restore' : 'Minimise'; };
  $('dcm').onclick = () => { const m = $('dcard').classList.toggle('min'); $('dcm').textContent = m ? '+' : '–'; };
  function overlay(title, tabs) {
    $('dovt').textContent = title; $('dovtabs').innerHTML = '';
    tabs.forEach((t, i) => { const b = document.createElement('button'); b.className = 'sm'; b.textContent = t[0];
      b.onclick = () => { $('dovf').src = t[1]; [...$('dovtabs').children].forEach(x => x.style.fontWeight = ''); b.style.fontWeight = '700'; };
      $('dovtabs').appendChild(b); if (i === 0) b.onclick(); });
    $('dover').classList.remove('min'); $('dovm').textContent = 'Minimise'; $('dover').style.display = 'block';
  }
  function parseMarker(d) {        // structured (plane / torus pages) or the hyperbolic viewers' info text
    if (d.marker) return d.p && !d.marker.p ? Object.assign({ p: d.p }, d.marker) : d.marker;
    if (!d.info) return null;
    const L = d.info.split('\n'), t = L[0].trim(), m = t.match(/^(center|center|puncture) (\S+)/);
    if (!m) return null;
    const mk = { kind: m[1] === 'puncture' ? 'puncture' : 'centre', label: m[2], title: t, lines: L.slice(1) };   // titles may say center or center
    if (m[1] === 'puncture') { const sh = d.info.match(/shared by centers: ([^\n]+)/); mk.labels = sh ? [...sh[1].matchAll(/(?:^|, )(\d+)/g)].map(x => +x[1]) : []; }
    return mk;
  }
  const POLY = {};
  function loadPoly(n) { if (n > 9) return Promise.resolve(null);   /* polynomial cards exist for periods 2..9 */
    if (!POLY[n]) POLY[n] = fetch(`poly/poly_info_n${n}.json`).then(r => r.ok ? r.json() : null).catch(() => null); return POLY[n]; }
  async function showTree(mk) {                 // the center's Hubbard tree (polynomial coordinates, E92/E96 poly_info), via psi:
    TREE = null; const n = mk.n || cur.n, Mc = M; if (n < 2 || !Mc) return;     // 0 -> a0, c -> a1, inf -> omega (Moebius)
    const lab = mk.ks ? (mk.ks.find(k => k % 2 === 0) ?? mk.ks[0]) : +mk.label, D = await loadPoly(n); if (!D || M !== Mc) return;
    const P = D.find(q => q.label === lab); if (!P) return;
    const pc = [Number(P.c_str[0]), Number(P.c_str[1])], psi = cur.family === 'quad' ? (q => q) : PM.moebius3([[0, 0], pc, null], [Mc.marked[0], Mc.marked[1], Mc.omega]);
    let err = 0, z = [0, 0]; if (cur.family !== 'quad') for (let k = 0; k < n; k++) { err = Math.max(err, PM.chord(psi(z), Mc.marked[k])); z = PM.cadd(PM.cmul(z, z), pc); }
    TREE = P.tree.map(arc => arc.map(q => psi(q))); TREEERR = err; dv.redraw();
    TREEP = { model: cur.id, n, lab, P, c0: pc, arcsW: P.tree, mpos: null };
    $('status').textContent = statusText() + `\nHubbard tree of polynomial ${lab} (white), moved into this chart: max chordal |psi(P^k(0)) - a_k| = ${err.toExponential(1)}`;
  }
  // ---- captures (REGISTRATION addendum 2026-11-10 (g)); 'c': f ^ theta
  function keyC() { if (!cur) return; $('cdlg').style.display = 'block'; $('cdn').textContent = cur.n; $('cdf').value = ''; $('cdg').value = '';
    const cardLab = (document.querySelector('#dct b') || {}).textContent || '', m = cardLab.match(/^cent(?:er|re) (\d+)/); $('cdf').placeholder = m ? m[1] + ' (this card)' : 'angle or label'; $('cdf').dataset.def = m ? m[1] : ''; $('cdf').focus(); }
  async function doCapture(ft, gt, choice) {
    const n = cur.n, st = t => { $('status').textContent = t; };
    if (n < 2) { st('Capture: choose Per_n with n ≥ 2.'); return; }
    const P = await resolveF(ft, n); if (!P) { st(`Capture: "${ft}" is not a period-${n} center (give its label k, or either angle k/${2**n - 1})`); return; }
    const a = PM.parseAngle(gt); if (!a) { st('Capture: enter the angle θ as p/q'); return; }
    const cf = [Number(P.c_str[0]), Number(P.c_str[1])], m = CAP.capture(cf, n, a, { choice: Math.max(0, (choice || 1) - 1) }), label = `${P.label} ^ ${a.text}`;
    const cl = m.comps ? `\nFatou components at the landing point of the ray ${a.text} (z = ${PM.fmt(m.landing, 6)}): ` + m.comps.map((v, i) => `#${i + 1}: center ${PM.fmt(v.x, 5)}${v.cycIndex >= 0 ? ' = x' + (v.cycIndex === 0 ? n : v.cycIndex) + ' (periodic)' : ', reaches the cycle after ' + v.t + ' step' + (v.t === 1 ? '' : 's')}`).join(';  ') : '';
    if (!m.ok) { st(`${label}: no capture: ${m.why}${m.iters !== undefined ? ' at pullback ' + m.iters : ''}` + cl); return; }
    const ch = MATE.toChart(m), kind = m.typeII ? `bitransitive (the free critical point sits on the marked cycle in place of x${m.predicted.slot})` : `capture (type C): the free critical point reaches the marked cycle after ${m.predicted.t} steps`;
    const base = `${label}: ${kind}; ${m.iters} pullbacks, residual ${m.residual.toExponential(1)}` + cl;
    let msg;
    if (cur.kind === 'plane') { if (!ch.p || !isFinite(ch.p[0] + ch.p[1])) { st(base + `\nThis map is at ∞ of the Per_${n} chart shown here, so it cannot be marked in the plane.`); return; } msg = { p: ch.p }; }
    else if (cur.kind === 'torus') { const L = PM.torusLocate(ch.x, ch.y); if (!L) { st(base + '\ntorus point not found'); return; } msg = { p: L.u }; }
    else msg = { x: ch.x, y: ch.y };
    $('frame').contentWindow.postMessage(Object.assign({ e96: 'mark', label, go: cur.kind !== 'hyperbolic' }, msg), '*');
    handleParam(Object.assign({ e96: 'param', label, keepCard: true }, msg));
    st(base + (cur.kind === 'hyperbolic' ? ' (locating it in the disk…)' : ' (yellow star)')); if (cur.kind === 'hyperbolic') MATEOUT = { out: { set textContent(t) { $('status').textContent = t; } }, base, label };
  }
  // (a) captures requested by a click in the polynomial card (REGISTRATION addendum 2026-11-11 (i))
  const PEND = new Map();                               // label -> interval: marks waiting for the pane (after a model switch)
  function placeMark(n, ch, label) {                    // returns a short text; switches to Per_n if needed
    if (!cur || cur.n !== n) choose('per' + n);
    let msg; if (cur.kind === 'plane') { if (!ch.p || !isFinite(ch.p[0] + ch.p[1])) return 'at ∞ of the Per_' + n + ' chart (not drawn)'; msg = { p: ch.p }; }
    else if (cur.kind === 'torus') { const L = PM.torusLocate(ch.x, ch.y); if (!L) return 'torus point not found'; msg = { p: L.u }; } else msg = { x: ch.x, y: ch.y };
    handleParam(Object.assign({ e96: 'param', label, keepCard: true }, msg)); clearInterval(PEND.get(label)); let tries = 0;
    const post = () => { if (++tries > 150) { clearInterval(PEND.get(label)); PEND.delete(label); return; } try { $('frame').contentWindow.postMessage(Object.assign({ e96: 'mark', label, go: true }, msg), '*'); } catch (e) {} };
    PEND.set(label, setInterval(post, 1200)); post(); return 'yellow star';
  }
  window.addEventListener('message', ev => { const d = ev.data; if (d && d.e96 === 'markDone' && PEND.has(d.label)) { clearInterval(PEND.get(d.label)); PEND.delete(d.label);
    if (!d.ok) $('status').textContent += `\n${d.label}: not located in the disk (best chart residual ${d.res != null ? d.res.toExponential(1) : '?'})`; } });
  window.addEventListener('message', async ev => { const d = ev.data; if (!d || d.e96 !== 'capReq') return;
    const n = d.n, D = await loadPoly(n), P = D && D.find(q => String(q.label) === String(d.label)); if (!P) return; const c = [Number(P.c_str[0]), Number(P.c_str[1])], lines = [];
    const kindOf = m => m.typeII ? `bitransitive (ω in place of x${m.predicted.slot})` : `capture, ω reaches the cycle after ${m.predicted.t} steps`;
    if (d.angle) { const a = PM.parseAngle(d.angle); if (d.comp === 1) { $('status').textContent = `${P.label} ^ ${d.angle}: this ray lands at the root of U1, the critical value component of f itself: no capture`; return; }
      const m = CAP.capture(c, n, a, { wantCyc: d.comp % n }), label = `${P.label} ^ ${d.angle}`;
      if (!m.ok) { $('status').textContent = `${label} (ray at the root of U${d.comp}): no capture: ${m.why}`; return; }
      const where = placeMark(n, MATE.toChart(m), label); $('status').textContent = `${label} (ray at the root of U${d.comp}): ${kindOf(m)}; ${m.iters} pullbacks, residual ${m.residual.toExponential(1)}; ${where}`; return; }
    const V = CAP.componentOf(c, n, d.point), zt = PM.fmt(d.point, 5); if (!V) { $('status').textContent = `Polynomial ${P.label}: the point ${zt} is not in a bounded Fatou component (click inside one, or on a ray)`; return; }
    const vname = V.cycIndex >= 0 ? `U${V.cycIndex === 0 ? n : V.cycIndex} (periodic)` : `the component with center ${PM.fmt(V.x, 5)} (reaches the cycle after ${V.t} step${V.t === 1 ? '' : 's'})`;
    if (V.cycIndex === 1) { $('status').textContent = `Polynomial ${P.label}: ${zt} is in U1, the critical value component of f itself: no capture`; return; }
    const R = CAP.rootRays(c, n, P.ks, V); if (R.why) { $('status').textContent = `Polynomial ${P.label}, ${vname}: ${R.why}`; return; }
    const res = {}; for (const side of ['minus', 'plus']) { const m = CAP.capture(c, n, R[side], { wantCenter: V.x }); res[side] = m.ok ? { m, ch: MATE.toChart(m) } : { why: m.why }; }
    const keyOf = r => r.ch.x && n >= 5 ? [r.ch.x, r.ch.y] : [r.ch.p], same = res.minus.ch && res.plus.ch && keyOf(res.minus).every((v, i) => PM.chord(v, keyOf(res.plus)[i]) < 1e-7);
    const head = `Polynomial ${P.label}, click at ${zt}: ${vname}; its root rays are ${R.minus.text} (−) and ${R.plus.text} (+).`;
    if (same) { const label = `${P.label} ^ ${R.minus.text}`, where = placeMark(n, res.minus.ch, label); $('status').textContent = `${head}\nBoth sides give the same map: ${kindOf(res.minus.m)}; residual ${res.minus.m.residual.toExponential(1)}; ${where} "${label}"`; return; }
    for (const side of ['minus', 'plus']) { const r = res[side], sg = side === 'minus' ? '−' : '+'; if (!r.ch) { lines.push(`side ${sg} (${R[side].text}): no capture: ${r.why}`); continue; }
      const label = `${P.label} ^ ${R[side].text}`, where = placeMark(n, r.ch, label); lines.push(`side ${sg}: ${kindOf(r.m)}; residual ${r.m.residual.toExponential(1)}; ${where} "${label}"`); }
    $('status').textContent = head + '\n' + lines.join('\n');
  });
  $('cdf').addEventListener('keydown', e => { if (e.key === 'Enter') $('cdg').focus(); });
  $('cdg').addEventListener('keydown', e => { if (e.key === 'Enter') { $('cdlg').style.display = 'none'; doCapture($('cdf').value.trim() || $('cdf').dataset.def || '', $('cdg').value, +$('cdc').value); } });
  $('cdlg').addEventListener('keydown', e => { if (e.key === 'Escape') $('cdlg').style.display = 'none'; });
  $('cdx').onclick = () => { $('cdlg').style.display = 'none'; };
  window.addEventListener('keydown', e => { if (e.key === 'c' && !e.metaKey && !e.ctrlKey && e.target.tagName !== 'INPUT' && e.target.tagName !== 'SELECT') { e.preventDefault(); keyC(); } });
  window.addEventListener('message', ev => { if (ev.data && ev.data.e96 === 'key' && ev.data.key === 'c') { window.focus(); keyC(); } });
  // ---- puncture explorer, part 1 (REGISTRATION addendum 2026-11-10 (f)): spheres of the marked configuration at the selected map
  let SPHR = null;
  function spheres() {
    const box = $('sphl'); SVIEW = null; if (!M) { box.textContent = 'click a parameter first (near a puncture)'; return; }
    const eps = Math.pow(10, +$('sphe').value), R = SPH.analyze(M, eps); SPHR = R; box.innerHTML = '';
    const add = (txt, fn, title) => { const d = document.createElement('div'); if (fn) { const b = document.createElement('button'); b.className = 'sm'; b.textContent = 'view'; b.title = title || ''; b.onclick = fn; d.appendChild(b); d.appendChild(document.createTextNode(' ')); }
      d.appendChild(document.createTextNode(txt)); box.appendChild(d); };
    if (!R.compatible) { add(`no tree at ε = ${eps.toExponential(0)}: the ${R.necks.length} necks cross (try a smaller ε)`); polyView(); fit(); return; }
    add(`${R.spheres.length} sphere${R.spheres.length === 1 ? '' : 's'} at ε = ${eps.toExponential(0)} (points: ${R.labels.join(', ')})`);
    add('whole plane (chart coordinate)', () => { SVIEW = null; polyView(); fit(); });
    SPH.summary(R).split('\n').forEach((line, s) => add(line, () => { SVIEW = R.spheres[s]; polyView(); fit(); }, 'the dynamical plane in the coordinate of this sphere: its largest direction at ∞, two others at −1 and 1'));
    add('→ is the sphere containing the images of three points from different directions; "inconsistent" = different choices give different spheres.');
    polyView(); fit();
  }
  $('sphb').onclick = spheres;
  $('sphe').oninput = e => { $('sphev').textContent = Math.pow(10, +e.target.value).toExponential(0); if (SPHR) spheres(); };
  window.addEventListener('keydown', e => { if (e.key === 'u' && e.target.tagName !== 'INPUT' && e.target.tagName !== 'SELECT') { $('dside').classList.remove('min'); spheres(); } });
  window.addEventListener('message', ev => { if (ev.data && ev.data.e96 === 'key' && ev.data.key === 'u') { $('dside').classList.remove('min'); spheres(); } });
  // ---- Hubbard tree off-center (REGISTRATION addendum 2026-11-10 (e)): the tree of the last polynomial, continued to the selected map
  let TREEP = null, TREEJOB = 0;
  async function followTree() {
    if (!M || !cur || cur.n < 2 || cur.family === 'quad' || !opt.tree) return; const Mc = M, job = ++TREEJOB, n = cur.n, hyp = cur.kind === 'hyperbolic';
    const f = PM.fate(Mc, Mc.omega, 2000, 1e-6, 1e-9); let note = '';
    if (f.kind === 'cycle' && f.period === 1) {                     // inside a polynomial component: whose?
      const C = PM.cycleCentre(hyp ? 'cor' : cur.family, n, hyp ? { x: sel.x, y: sel.y } : { p: sel.p }, 1), T = C && PM.transfer(C.M, 1), D = await loadPoly(n);
      if (job !== TREEJOB || M !== Mc) return;
      const P = T && D && D.find(q => PM.cabs(PM.csub([Number(q.c_str[0]), Number(q.c_str[1])], T.p)) < 1e-6);
      if (P && (!TREEP || TREEP.P !== P || TREEP.model !== cur.id)) TREEP = { model: cur.id, n, lab: P.label, P, c0: [Number(P.c_str[0]), Number(P.c_str[1])], arcsW: P.tree, mpos: null };
    }
    if (!TREEP || TREEP.model !== cur.id) return;
    const orb0 = () => { const o = [[0, 0]]; for (let k = 1; k < n; k++) o.push(PM.cadd(PM.cmul(o[k - 1], o[k - 1]), TREEP.c0)); return o; };
    let R = HT.deform(Mc, TREEP.c0, TREEP.arcsW, { mpos: TREEP.mpos || orb0(), res: 400, maxSteps: 70 }), from = 'continued from the previous tree';
    if (R.why && TREEP.arcsW !== TREEP.P.tree) { const R2 = HT.deform(Mc, TREEP.c0, TREEP.P.tree, { mpos: orb0(), res: 400, maxSteps: 70 }); if (!R2.why) { R = R2; from = 'restarted from the center\'s tree'; } }
    if (M !== Mc || job !== TREEJOB) return;
    if (!R.why) { TREEP.arcsW = R.arcs; TREEP.mpos = R.mpos; }
    TREE = (R.why && R.steps === 0) ? null : R.arcs.map(arc => arc.map(q => R.F.psi(q))); dv.redraw();
    note = R.why ? `Hubbard tree of polynomial ${TREEP.lab}: the pullback stops here after ${R.steps} step${R.steps === 1 ? '' : 's'}: ${R.why}` + (TREE ? ' (last tree shown)' : '')
      : `Hubbard tree of polynomial ${TREEP.lab} (white) for THIS map: regulated hull of the marked cycle by pullback, ${from}; ${R.steps} steps, last movement ${(R.defect/R.h).toFixed(2)} node spacings${R.defect > 0.5*R.h ? ' (NOT yet settled: the drawn tree is approximate; click again nearby to continue it)' : ''}, branch degrees [${R.shape}]`;
    $('status').textContent = statusText() + '\n' + note;
  }
  let TREEERR = null;
  const GN = {};
  function loadG(n) { if (n > 8) return Promise.resolve(null);   /* no G_9 yet */
    if (!GN[n]) GN[n] = fetch(`models/Gn/G${n}.json`).then(r => r.ok ? r.json() : null).catch(() => null); return GN[n]; }
  async function evenPath(n, lab) {              // shortest path in G_n with strictly more factors of 2 at each step, to 2^(n-1)
    const G = await loadG(n); if (!G) return null; const L = G.labels, v2 = k => { let v = 0; while (k > 0 && k % 2 === 0) { k /= 2; v++; } return v; };
    const adj = L.map(() => new Set()); for (const e of G.edges) if (e.b !== null && e.b >= 0) { adj[e.a].add(e.b); adj[e.b].add(e.a); }
    const s = L.indexOf(lab), t = L.indexOf(2**(n - 1)); if (s < 0 || t < 0) return null;
    const prev = new Map([[s, -1]]), q = [s];
    while (q.length) { const i = q.shift(); if (i === t) break; for (const j of adj[i]) if (!prev.has(j) && v2(L[j]) > v2(L[i])) { prev.set(j, i); q.push(j); } }
    if (!prev.has(t)) return []; const path = []; for (let i = t; i !== -1; i = prev.get(i)) path.unshift(L[i]); return path;
  }
  function parseC(mk) {                         // the polynomial c of a center card (census lines 'c = re +/- im i')
    if (mk.pc) return mk.pc.map(Number);
    const t = (mk.lines || []).join('\n').match(/c = ([-+]?[\d.eE-]+) ([+-]) ([\d.eE-]+)i/); return t ? [Number(t[1]), (t[2] === '-' ? -1 : 1)*Number(t[3])] : null;
  }
  // ---- mating calculator (mating.js; REGISTRATION addendum 2026-11-10, checks M1-M5 in log.md)
  function doMate(cf, n, lab, text, out) {
    out = out || { set textContent(t) { $('status').textContent = t; } };
    const a = PM.parseAngle(text); if (!a || !cf) { out.textContent = 'enter an angle p/q'; return; }
    const g = MATE.pcfFromAngle(a), m = MATE.mate(cf, g.c), gname = `g = ${g.kind === 'centre' ? 'center' : g.kind} of angle ${a.text} (c = ${PM.fmt(g.c, 6)})`;
    if (g.kind === 'misiurewicz' && g.moved > 0.05) { out.textContent = `g not found: the parameter ray of ${a.text} did not land well enough (Newton moved ${g.moved.toExponential(1)} from the ray's end)`; return; }
    if (!m.ok) { out.textContent = `not converged: ${m.why} (${gname})`; return; }
    const ch = MATE.toChart(m), label = `${lab} ⊔ ${a.text}`;
    if (!ch || cur.n !== n) { out.textContent = `f ⊔ g lies in Per_${n}; open Per_${n} to place it (coming next)`; return; }
    let msg;
    if (cur.kind === 'plane') msg = { p: ch.p };
    else if (cur.kind === 'torus') { const L = PM.torusLocate(ch.x, ch.y); if (!L) { out.textContent = 'torus point not found'; return; } msg = { p: L.u }; }
    else msg = { x: ch.x, y: ch.y };
    $('frame').contentWindow.postMessage(Object.assign({ e96: 'mark', label }, msg), '*');
    handleParam(Object.assign({ e96: 'param', label, keepCard: true }, msg));
    const base = `${label}: converged in ${m.iters} steps, residual ${m.residual.toExponential(1)}; ${gname}`;
    out.textContent = base + (cur.kind === 'hyperbolic' ? ' (locating it in the disk…)' : ' (yellow marker)');
    if (cur.kind === 'hyperbolic') MATEOUT = { out, base, label };
  }
  let MATEOUT = null;
  window.addEventListener('message', ev => {         // the hyperbolic viewer's answer to 'mark' (chart inversion, M6)
    const d = ev.data; if (!d || d.e96 !== 'markDone') return;
    if (!d.ok && (!MATEOUT || d.label !== MATEOUT.label)) { $('status').textContent += `\nvein point ${d.label}: not located in the disk (chart residual ${d.res != null ? d.res.toExponential(1) : '?'})`; return; }
    if (!MATEOUT || d.label !== MATEOUT.label) return;
    MATEOUT.out.textContent = MATEOUT.base + (d.ok ? ` (yellow star in the disk; chart residual ${d.res.toExponential(1)})`
      : ` (not located in the disk${d.res != null ? ': best chart residual ' + d.res.toExponential(1) : ''})`);
    MATEOUT = null;
  });
  // jump to the polynomial with label lab in Per_n (retries until the model has loaded and answers with its card)
  let GOTO = null;
  function viewIn(n, lab) {
    choose('per' + n); clearInterval(GOTO && GOTO.t); let tries = 0;
    GOTO = { lab: String(lab), t: setInterval(() => { tries++; if (tries > 120) { clearInterval(GOTO.t); GOTO = null; return; }
      try { $('frame').contentWindow.postMessage({ e96: 'goto', label: lab }, '*'); } catch (e) {} }, 1500) };
  }
  // puncture card (REGISTRATION addendum 2026-11-12 (t)): what happens to the marked cycle near the puncture. Plane / torus models: the maps at the
  // puncture's chart point displaced by 1e-3 and 1e-5; disk models: the map at the census loop point (the marker's own x, y).
  function punctLines(mk) {
    const out = [], maps = [];
    try {
      if (cur.kind === 'hyperbolic') { if (M) maps.push(['at the census loop point (next to the puncture)', M]); }
      else if (mk.p) for (const e of [1e-3, 1e-5]) { let Mx = null; try { Mx = PM.mapFrom(cur.family, cur.n, { p: [mk.p[0] + e*Math.cos(0.7), mk.p[1] + e*Math.sin(0.7)] }); } catch (err) {} if (Mx && isFinite(Mx.omega[0] + Mx.omega[1])) maps.push([`at distance ${e} from the puncture (chart coordinate)`, Mx]); }
      for (const [nm, Mx] of maps) { const I = PM.punctureInfo(Mx);
        out.push(`— ${nm}:`);
        out.push(I.clusters.length ? '  marked points that have nearly collided: ' + I.clusters.map((c, i) => '{' + c.map(k => 'a' + k).join(', ') + '} (diameter ' + I.diam[i].toExponential(1) + ')').join(', ') + (I.minOther !== null ? `; all other marked distances ≥ ${I.minOther.toExponential(1)}` : '')
                                 : `  no marked points have collided (smallest marked distance ${I.minOther !== null ? I.minOther.toExponential(1) : '?'})`);
        out.push(`  free critical point ω: nearest marked point a${I.omega.i} at ${I.omega.d.toExponential(1)}; its value f(ω): nearest a${I.value.i} at ${I.value.d.toExponential(1)}   (a0 = marked critical point, a1 = marked critical value; distances are chordal)`);
        if (I.fixed) out.push('  fixed points: ' + I.fixed.map(f => `z = ${PM.fmt(f.z, 4)}, multiplier ${f.mult ? PM.fmt(f.mult, 4) + ' (|·| = ' + PM.cabs(f.mult).toExponential(2) + ')' : '?'}`).join(';  ')); }
    } catch (e) { out.push('(puncture data not computed: ' + e.message + ')'); }
    return out;
  }
  function showMarker(mk) {
    if (GOTO && mk.label === GOTO.lab) { clearInterval(GOTO.t); GOTO = null; }
    const n = mk.n || cur.n, D = 2**n - 1, esc = t => String(t).replace(/&/g, '&amp;').replace(/</g, '&lt;');
    $('dct').innerHTML = '<b>' + esc(mk.title) + '</b>\n' + esc(mk.lines.join('\n'));
    if (mk.kind === 'puncture') { const extra = punctLines(mk); if (extra.length) $('dct').innerHTML += '\n' + esc(extra.join('\n')); }
    const B = $('dcb'); B.innerHTML = '';
    const btn = (txt, fn) => { const b = document.createElement('button'); b.className = 'sm'; b.textContent = txt; b.onclick = fn; B.appendChild(b); };
    if (n >= 2) {
      if (mk.kind === 'centre') {
        const lab = mk.ks ? mk.ks.find(k => k % 2 === 0) ?? mk.ks[0] : +mk.label;
        if (n !== cur.n && man.models.some(q => q.id === 'per' + n && !q.placeholder))
          btn(`View in Per_${n}`, () => viewIn(n, lab));
        btn('Polynomial card: Hubbard tree, rays, addresses', () => overlay(`Per_${n}: polynomial ${lab}`, [[`${lab}/${D}`, `poly/index.html#n=${n}&label=${lab}`]]));
        btn('Obstruction Γ(θ) of its puncture', () => overlay(`Thurston obstruction Γ(θ), θ = ${lab}/${D}`, [[`θ = ${lab}/${D}`, `obstruction/index.html#p=${lab}&q=${D}`]]));
        const mrow = document.createElement('div'); mrow.style.cssText = 'width:100%;margin-top:4px';
        mrow.innerHTML = 'Mate with g of angle <input class="mang" size="9" placeholder="e.g. 1/7"> <button class="sm mbtn">f ⊔ g</button> <span class="mres"></span>';
        B.appendChild(mrow);
        const cf = parseC(mk);
        mrow.querySelector('.mbtn').onclick = () => doMate(cf, n, lab, mrow.querySelector('.mang').value, mrow.querySelector('.mres'));
        mrow.querySelector('.mang').addEventListener('keydown', e => { if (e.key === 'Enter') mrow.querySelector('.mbtn').click(); });
        if (n >= 5) evenPath(n, lab).then(path => { if (!path) return; const d = document.createElement('div'); d.style.cssText = 'width:100%;margin-top:4px';
          d.innerHTML = path.length ? `G<sub>${n}</sub> path to ${2**(n - 1)} (more factors of 2 at each step; click to go): ` : `no increasing G<sub>${n}</sub> path from ${lab}`;
          path.forEach((L, i) => { const a = document.createElement('a'); a.href = '#'; a.textContent = L; a.style.margin = '0 3px';
            a.onclick = ev => { ev.preventDefault(); $('frame').contentWindow.postMessage({ e96: 'goto', label: L }, '*'); }; if (i) d.append('→'); d.appendChild(a); });
          B.appendChild(d); });
      } else {
        let labs = mk.labels || (mk.thetas ? mk.thetas.filter(k => k % 2 === 0) : []);
        if (!labs.length && mk.thetas) labs = mk.thetas;
        const openObs = () => overlay(`Puncture ${mk.label}: Thurston obstruction Γ(θ) (θ = an angle of a center sharing it)`,
          labs.map(L => [`θ = ${L}/${D}`, `obstruction/index.html#p=${L}&q=${D}`]));
        btn("Thurston obstruction Γ(θ)", openObs);
        btn('Polynomial cards of the sharing centers', () => overlay(`Puncture ${mk.label}: polynomials sharing it`, labs.map(L => [`${L}`, `poly/index.html#n=${n}&label=${L}`])));
      }
    }
    $('dcard').style.display = 'block';
  }
  $('dmin').onclick = () => { $('dside').classList.toggle('min'); $('dmin').textContent = $('dside').classList.contains('min') ? '+' : '–'; };
  if ($('right').clientWidth < 900) $('dmin').onclick();
  window.addEventListener('hashchange', () => { const h2 = location.hash.slice(1); if (man.models.some(m => m.id === h2 && !m.placeholder) && (!cur || cur.id !== h2)) choose(h2); });
  const h = location.hash.slice(1); choose(man.models.some(m => m.id === h && !m.placeholder) ? h : 'per1');
})().catch(e => { document.getElementById('status').textContent = 'error: ' + e.message; console.error(e); });
