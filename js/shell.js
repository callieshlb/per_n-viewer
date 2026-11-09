// E96 shell: model menu, the parameter pane (iframe, one page per manifest entry), the dynamical plane, the status line.
// Parameter panes post {e96: 'param', model, p | x, y}; the map is built here by PM.mapFrom(family, n, message).
'use strict';
(async function () {
  const $ = s => document.getElementById(s);
  const man = await fetch('models/manifest.json').then(r => r.json());
  let cur = null, M = null, sel = null, frac = 0.55;
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
    $('kind').textContent = { plane: 'parameter plane', torus: 'flat torus (genus 1)', hyperbolic: 'hyperbolic model (uniformized compact curve)' }[cur.kind] || '';
    $('frame').src = cur.page; history.replaceState(null, '', '#' + id);
    M = null; sel = null; $('dhint').style.display = ''; dv.redraw();
    $('status').textContent = `${cur.titleText}: ${cur.kind === 'hyperbolic' ? 'loading the model (tens of MB; the first load can take a while)…' : 'click a parameter'}`;
  }
  // ---- dynamical plane
  const opt = { mk: true, ml: true, orb: false, orbN: 60 };
  const dv = new PlaneView($('dyn'), {
    ctr: [0, 0], upp: 0.01, maxIt: 500, delta: 1e-3,
    setUniforms: (gl, U) => {
      if (!M) { gl.uniform1i(U('fam'), 1); gl.uniform1i(U('maxIt'), 0); return; }
      gl.uniform1i(U('fam'), 0);
      ['A', 'B', 'C', 'D', 'E', 'F'].forEach(k => gl.uniform2f(U('u' + k), M[k][0], M[k][1]));
      const mk = new Float32Array(24), mi = new Int32Array(12);
      M.marked.forEach((a, j) => { if (a === null) mi[j] = 1; else { mk[2*j] = a[0]; mk[2*j + 1] = a[1]; } });
      gl.uniform2fv(U('uMk'), mk); gl.uniform1iv(U('uMkI'), mi); gl.uniform1i(U('nM'), M.marked.length);
    },
    drawOverlay: (g, W, H) => {
      if (!M) { g.fillStyle = '#111'; g.fillRect(0, 0, W, H); return; }
      if (opt.orb) {
        const o = PM.orbit(M, M.omega, opt.orbN); g.strokeStyle = 'rgba(255,230,80,.75)'; g.lineWidth = 1; g.beginPath(); let on = false;
        for (const z of o) { if (z === null) { on = false; continue; } const s = dv.toScreen(z); if (on) g.lineTo(...s); else g.moveTo(...s); on = true; }
        g.stroke(); g.fillStyle = '#ffe650';
        for (const z of o) if (z) { const s = dv.toScreen(z); g.fillRect(s[0] - 1.5, s[1] - 1.5, 3, 3); }
      }
      g.font = '600 12px system-ui'; g.textBaseline = 'middle';
      const dot = (z, fill, lab) => { const s = dv.toScreen(z); if (s[0] < -20 || s[1] < -20 || s[0] > W + 20 || s[1] > H + 20) return;
        g.beginPath(); g.arc(s[0], s[1], 5, 0, 2*Math.PI); g.fillStyle = fill; g.fill(); g.strokeStyle = '#141414'; g.lineWidth = 1.2; g.stroke();
        if (opt.ml && lab) { g.lineWidth = 3; g.strokeStyle = 'rgba(0,0,0,.85)'; g.strokeText(lab, s[0] + 8, s[1]); g.fillStyle = fill; g.fillText(lab, s[0] + 8, s[1]); } };
      if (opt.mk) M.marked.forEach((a, j) => { if (a) dot(a, '#ffffff', 'a' + j); });
      dot(M.omega, '#ffe650', 'ω');
    },
    onPick: p => {
      if (!M) return; const f = PM.fate(M, p, dv.maxIt, dv.delta, 1e-7);
      $('status').textContent = statusText() + `\nclicked z = ${PM.fmt(p)}: ${PM.fateText(f)}`;
    } });
  function fit() {
    if (!M) return; const pts = M.marked.filter(a => a && PM.cabs(a) < 1e3).concat([M.omega, [0, 0]]);
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (const z of pts) { x0 = Math.min(x0, z[0]); x1 = Math.max(x1, z[0]); y0 = Math.min(y0, z[1]); y1 = Math.max(y1, z[1]); }
    const [W, H] = dv.size(), span = Math.max(x1 - x0, y1 - y0, 1.5)*1.8;
    dv.ctr = [(x0 + x1)/2, (y0 + y1)/2]; dv.upp = span/Math.max(50, Math.min(W, H)); dv.redraw();
  }
  function statusText() {
    if (!M) return '';
    const par = (cur.kind !== 'plane') && sel.x && sel.y ? `x = a${M.n - 1} = ${PM.fmt(sel.x)}   y = a3 = ${PM.fmt(sel.y)}${sel.p ? '   (' + cur.param + ' = ' + PM.fmt(sel.p) + ')' : ''}` : `${cur.param} = ${PM.fmt(sel.p)}`;
    const f = PM.fate(M, M.omega, dv.maxIt, dv.delta, 1e-7);
    return `${cur.titleText}${sel.label ? '  ·  ' + sel.label : ''}   ${par}\n` +
      `f(z) = ((${PM.fmt(M.A, 5)}) z² + (${PM.fmt(M.B, 5)}) z + (${PM.fmt(M.C, 5)})) / ((${PM.fmt(M.D, 5)}) z² + (${PM.fmt(M.E, 5)}) z + (${PM.fmt(M.F, 5)}))\n` +
      `marked cycle: ${M.marked.map((a, j) => 'a' + j + ' = ' + PM.fmt(a, 5)).join(',  ')}   (closes to ${PM.closure(M).toExponential(1)})\n` +
      `free critical point ω = ${PM.fmt(M.omega, 6)}: ${PM.fateText(f)}   |   iterates ${dv.maxIt}, δ ${dv.delta.toExponential(0)}`;
  }
  window.addEventListener('message', ev => {
    const d = ev.data; if (!d || d.e96 !== 'param' || !cur) return;
    const q = { p: d.p, x: d.x, y: d.y };
    const fam = cur.kind === 'hyperbolic' ? 'cor' : cur.family;
    const m = PM.mapFrom(fam, cur.n, q);
    if (!m) { $('status').textContent = `${cur.titleText}: degenerate map at this point (puncture)`; return; }
    M = m; sel = { p: d.p || null, x: d.x || (m.x || null), y: d.y || (m.y || null), label: d.label };
    if (cur.family === 'quad' || cur.family === 'per2') { sel.x = null; sel.y = null; }
    $('dhint').style.display = 'none'; fit(); $('status').textContent = statusText();
  });
  window.SHELL = { get M() { return M; }, get sel() { return sel; }, dv, choose };
  // ---- controls
  $('dit').oninput = e => { dv.maxIt = +e.target.value; $('ditv').textContent = dv.maxIt; dv.redraw(); if (M) $('status').textContent = statusText(); };
  $('ddel').oninput = e => { dv.delta = Math.pow(10, +e.target.value); $('ddelv').textContent = dv.delta.toExponential(0); dv.redraw(); };
  $('dmk').onchange = e => { opt.mk = e.target.checked; dv.redraw(); };
  $('dml').onchange = e => { opt.ml = e.target.checked; dv.redraw(); };
  $('dorb').onchange = e => { opt.orb = e.target.checked; dv.redraw(); };
  $('dorbn').oninput = e => { opt.orbN = Math.max(1, Math.min(2000, +e.target.value || 1)); dv.redraw(); };
  $('dhome').onclick = fit;
  $('dmin').onclick = () => { $('dside').classList.toggle('min'); $('dmin').textContent = $('dside').classList.contains('min') ? '+' : '–'; };
  if ($('right').clientWidth < 900) $('dmin').onclick();
  const h = location.hash.slice(1); choose(man.models.some(m => m.id === h && !m.placeholder) ? h : 'per1');
})().catch(e => { document.getElementById('status').textContent = 'error: ' + e.message; console.error(e); });
