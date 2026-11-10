// E96 shell: model menu, the parameter pane (iframe, one page per manifest entry), the dynamical plane, the status line.
// Parameter panes post {e96: 'param', model, p | x, y}; the map is built here by PM.mapFrom(family, n, message).
'use strict';
(async function () {
  const $ = s => document.getElementById(s);
  const man = await fetch('models/manifest.json?v=20261110021841').then(r => r.json());
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
    hideCard(); $('dover').style.display = 'none';
    $('frame').src = cur.page + (cur.page.includes('?') ? '&' : '?') + 'v=20261110021841'; history.replaceState(null, '', '#' + id);
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
    TV = null; if (!opt.inv || !M) return;
    const T0 = PM.moebius3([M.marked[0], M.marked[1], M.omega], [[0, 0], [1, 0], null]), T0i = PM.moebius3([[0, 0], [1, 0], null], [M.marked[0], M.marked[1], M.omega]);
    const g0 = w => T0(PM.fz(M, T0i(w))), b0 = g0([0, 0]), aI = g0(null), g1 = g0([1, 0]);
    let A; if (aI === null) A = PM.csub(g1, b0); else { const C = PM.cdiv(PM.csub(g1, b0), PM.csub(aI, g1)); A = PM.cmul(aI, C); }
    const s = A, si = PM.cdiv([1, 0], s);
    const T = z => { const w = T0(z); return w === null ? null : PM.cmul(s, w); }, Ti = w => T0i(w === null ? null : PM.cmul(si, w));
    // the inverse as a matrix for the shader: images of 0, 1, inf under Ti give (mA w + mB)/(mC w + mD)
    const z0 = Ti([0, 0]), z1 = Ti([1, 0]), zi = Ti(null), hv = z => z === null ? [[1, 0], [0, 0]] : [z, [1, 0]];
    // matrix with columns (zi-vector * a, z0-vector * b), a + b normalized so that column sum gives z1
    const P = hv(zi), Q = hv(z0), R = hv(z1), det = PM.csub(PM.cmul(P[0], Q[1]), PM.cmul(P[1], Q[0]));
    const a = PM.cdiv(PM.csub(PM.cmul(R[0], Q[1]), PM.cmul(R[1], Q[0])), det), b = PM.cdiv(PM.csub(PM.cmul(P[0], R[1]), PM.cmul(P[1], R[0])), det);
    TV = { T, Ti, s, m: [PM.cmul(P[0], a), PM.cmul(Q[0], b), PM.cmul(P[1], a), PM.cmul(Q[1], b)] };
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
    return `${cur.titleText}${sel.label ? '  ·  ' + sel.label : ''}   ${par}\n` +
      `f(z) = ((${PM.fmt(M.A, 5)}) z² + (${PM.fmt(M.B, 5)}) z + (${PM.fmt(M.C, 5)})) / ((${PM.fmt(M.D, 5)}) z² + (${PM.fmt(M.E, 5)}) z + (${PM.fmt(M.F, 5)}))\n` +
      `marked cycle: ${M.marked.map((a, j) => 'a' + j + ' = ' + PM.fmt(a, 5)).join(',  ')}   (closes to ${PM.closure(M).toExponential(1)})\n` +
      `free critical point ω = ${PM.fmt(M.omega, 6)}: ${comp}   |   iterates ${dv.maxIt}, δ ${dv.delta.toExponential(0)}`;
  }
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
    M = m; sel = { p: d.p || null, x: d.x || (m.x || null), y: d.y || (m.y || null), label: d.label };
    if (cur.family === 'quad' || cur.family === 'per2') { sel.x = null; sel.y = null; }
    $('dhint').style.display = 'none'; dynRays(); polyView(); fit(); $('status').textContent = statusText(); ORB0 = null;
    TREE = null; const mk = parseMarker(d); if (mk) { showMarker(mk); if (mk.kind === 'centre') showTree(mk); } else if (!d.keepCard) hideCard();
  }
  window.SHELL = { get M() { return M; }, get sel() { return sel; }, dv, choose };
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
  async function mateFromText(ft, gt, K) {
    if (!cur) return; const n = cur.n;
    if (n < 2) { $('status').textContent = 'Mating: choose Per_n with n ≥ 2 (in Per_1, z² ⊔ g = g).'; return; }
    const P = await resolveF(ft, n);
    if (!P) { $('status').textContent = `Mating: "${ft}" is not a period-${n} center (give its label k, or either angle k/${2**n - 1})`; return; }
    const cf = [Number(P.c_str[0]), Number(P.c_str[1])]; let vtxt = '';
    if (K > 0) vtxt = await doVein(cf, n, P.label, gt, K);
    doMate(cf, n, P.label, gt, null);
    if (vtxt) { const o = MATEOUT && MATEOUT.out; if (MATEOUT) MATEOUT.base += '\n' + vtxt; $('status').textContent += '\n' + vtxt; }
  }
  // mating vein (REGISTRATION addendum 2026-11-10 (b), M7): f ⊔ g_i for the centers g_i (period ≤ 8) on the arc from 0 to g
  async function doVein(cf, n, lab, text, K) {
    const a = PM.parseAngle(text); if (!a) return '';
    const polys = {}; for (let m = 2; m <= 8; m++) polys[m] = await loadPoly(m);
    const V = MATE.veinList(a, polys, K), W = $('frame').contentWindow; let ok = 0; const obs = [];
    for (const v of V) { const m = MATE.mate(cf, v.c), nm = `${v.k0}/${v.Q}`; if (!m.ok) { obs.push(nm); continue; }
      const ch = MATE.toChart(m); if (!ch) continue; let msg;
      if (cur.kind === 'plane') msg = { p: ch.p };
      else if (cur.kind === 'torus') { const L = PM.torusLocate(ch.x, ch.y); if (!L) { obs.push(nm + ' (not located)'); continue; } msg = { p: L.u }; }
      else msg = { x: ch.x, y: ch.y };
      W.postMessage(Object.assign({ e96: 'mark', label: nm, vein: true, quiet: true }, msg), '*'); ok++; }
    return `vein of ${a.text}: ${V.length} centers of period ≤ 8 on the arc from z² to g (K = ${K}): ${ok} mated (small yellow dots, labeled by the lower angle of g_i)` +
      (obs.length ? `; ${obs.length} obstructed or not placed: ${obs.join(', ')}` : '');
  }
  function keyM() {                                   // 'm': mating prompt: f angle, g angle, Enter
    const dl = $('mdlg'); dl.style.display = 'block'; $('mdn').textContent = cur ? cur.n : 'n';
    const cardLab = (document.querySelector('#dct b') || {}).textContent || '', m = cardLab.match(/^cent(?:er|re) (\d+)/);
    if (m && !$('mdf').value) $('mdf').value = m[1];
    const first = $('mdf').value ? $('mdg') : $('mdf'); first.focus(); first.select();
  }
  $('mdf').addEventListener('keydown', e => { if (e.key === 'Enter') { $('mdg').focus(); $('mdg').select(); } });
  $('mdg').addEventListener('keydown', e => { if (e.key === 'Enter') { $('mdlg').style.display = 'none'; mateFromText($('mdf').value, $('mdg').value, $('mdv').checked ? Math.max(0, Math.min(200, +$('mdk').value || 0)) : 0); } });
  $('mdk').addEventListener('keydown', e => { if (e.key === 'Enter') { $('mdg').focus(); } });
  $('mdc').onclick = () => { $('frame').contentWindow.postMessage({ e96: 'clearMarks' }, '*'); };
  $('mdlg').addEventListener('keydown', e => { if (e.key === 'Escape') $('mdlg').style.display = 'none'; });
  $('mdx').onclick = () => { $('mdlg').style.display = 'none'; };
  window.addEventListener('keydown', e => { if (e.key === 'm' && e.target.tagName !== 'INPUT' && e.target.tagName !== 'SELECT') { e.preventDefault(); keyM(); } });
  window.addEventListener('message', ev => { if (ev.data && ev.data.e96 === 'key' && ev.data.key === 'm') { window.focus(); keyM(); } });
  window.addEventListener('message', ev => {                       // mating box in a model's side panel
    const d = ev.data; if (!d || d.e96 !== 'mateReq' || !cur) return; mateFromText(d.label, d.angle);
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
    if (d.marker) return d.marker;
    if (!d.info) return null;
    const L = d.info.split('\n'), t = L[0].trim(), m = t.match(/^(center|center|puncture) (\S+)/);
    if (!m) return null;
    const mk = { kind: m[1] === 'puncture' ? 'puncture' : 'centre', label: m[2], title: t, lines: L.slice(1) };   // titles may say center or center
    if (m[1] === 'puncture') { const sh = d.info.match(/shared by centers: ([^\n]+)/); mk.labels = sh ? [...sh[1].matchAll(/(?:^|, )(\d+)/g)].map(x => +x[1]) : []; }
    return mk;
  }
  const POLY = {};
  function loadPoly(n) { if (!POLY[n]) POLY[n] = fetch(`poly/poly_info_n${n}.json`).then(r => r.ok ? r.json() : null).catch(() => null); return POLY[n]; }
  async function showTree(mk) {                 // the center's Hubbard tree (polynomial coordinates, E92/E96 poly_info), via psi:
    TREE = null; const n = mk.n || cur.n, Mc = M; if (n < 2 || !Mc) return;     // 0 -> a0, c -> a1, inf -> omega (Moebius)
    const lab = mk.ks ? (mk.ks.find(k => k % 2 === 0) ?? mk.ks[0]) : +mk.label, D = await loadPoly(n); if (!D || M !== Mc) return;
    const P = D.find(q => q.label === lab); if (!P) return;
    const pc = [Number(P.c_str[0]), Number(P.c_str[1])], psi = cur.family === 'quad' ? (q => q) : PM.moebius3([[0, 0], pc, null], [Mc.marked[0], Mc.marked[1], Mc.omega]);
    let err = 0, z = [0, 0]; if (cur.family !== 'quad') for (let k = 0; k < n; k++) { err = Math.max(err, PM.chord(psi(z), Mc.marked[k])); z = PM.cadd(PM.cmul(z, z), pc); }
    TREE = P.tree.map(arc => arc.map(q => psi(q))); TREEERR = err; dv.redraw();
    $('status').textContent = statusText() + `\nHubbard tree of polynomial ${lab} (white), moved into this chart: max chordal |psi(P^k(0)) - a_k| = ${err.toExponential(1)}`;
  }
  let TREEERR = null;
  const GN = {};
  function loadG(n) { if (!GN[n]) GN[n] = fetch(`models/Gn/G${n}.json`).then(r => r.ok ? r.json() : null).catch(() => null); return GN[n]; }
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
  function showMarker(mk) {
    if (GOTO && mk.label === GOTO.lab) { clearInterval(GOTO.t); GOTO = null; }
    const n = mk.n || cur.n, D = 2**n - 1, esc = t => String(t).replace(/&/g, '&amp;').replace(/</g, '&lt;');
    $('dct').innerHTML = '<b>' + esc(mk.title) + '</b>\n' + esc(mk.lines.join('\n'));
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
