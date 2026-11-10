// E96: parameter pane for kind = plane (Per_1 .. Per_4) and kind = torus (Per_5). URL: plane.html?m=<model id>
'use strict';
(async function () {
  const id = new URLSearchParams(location.search).get('m') || 'per1';
  const man = await fetch('models/manifest.json').then(r => r.json());
  const E = man.models.find(m => m.id === id);
  const C = await fetch(E.census).then(r => r.json());
  const torus = E.kind === 'torus', T5 = PM.T5;
  const $ = s => document.getElementById(s);
  $('title').innerHTML = E.title; $('nf').innerHTML = E.normal; $('about').innerHTML = E.about || '';
  document.title = E.titleText || id;
  $('legend').innerHTML = GLSL.LEGEND(E.n, E.palette).map(([c, t]) => `<span class="sw" style="background:${c}"></span><span>${t}</span>`).join('');
  if (torus) ['gnl', 'gn', 'sml', 'sm'].forEach(k => { $(k).style.display = ''; });
  const opt = { mk: true, ml: true, gn: false, sm: true, rc: true };
  if (E.real) ['rcl', 'rc'].forEach(k => { $(k).style.display = ''; });
  ['mk', 'ml', 'gn', 'sm', 'rc'].forEach(k => { $(k).onchange = e => { opt[k] = e.target.checked; pv.redraw(); }; });
  $('minb').onclick = () => { $('side').classList.toggle('min'); $('minb').textContent = $('side').classList.contains('min') ? '+' : '–'; };
  let picked = null;
  function lattice(W, H) {           // visible lattice translates (torus)
    const a = pv.toPlane(0, H), b = pv.toPlane(W, 0), out = [];
    const m0 = Math.floor(a[0]/T5.W1) - 1, m1 = Math.ceil(b[0]/T5.W1) + 1, l0 = Math.floor(a[1]/T5.W2) - 1, l1 = Math.ceil(b[1]/T5.W2) + 1;
    if ((m1 - m0)*(l1 - l0) > 400) return [[0, 0]];
    for (let m = m0; m <= m1; m++) for (let l = l0; l <= l1; l++) out.push([m*T5.W1, l*T5.W2]);
    return out;
  }
  const markers = C.centres.map(c => ({ ...c, kind: 'centre' })).concat(C.punctures.filter(p => p.p).map(p => ({ ...p, kind: 'puncture' })));
  function drawOverlay(g, W, H) {
    const offs = torus ? lattice(W, H) : [[0, 0]];
    if (torus && opt.sm) {                    // seams: the fundamental rectangle's edges
      g.strokeStyle = 'rgba(255,255,255,.55)'; g.lineWidth = 1; g.setLineDash([5, 4]);
      for (const o of offs) { const p = pv.toScreen(o), q = pv.toScreen([o[0] + T5.W1, o[1] + T5.W2]);
        g.strokeRect(p[0], q[1], q[0] - p[0], p[1] - q[1]); }
      g.setLineDash([]);
    }
    if (E.real && opt.rc) {                   // real circles: where the normal form has real coefficients (log.md)
      const lines = [];
      if (torus) { const a = pv.toPlane(0, H), b = pv.toPlane(0, 0), h2 = T5.W2/2;
        for (let k = Math.floor(a[1]/h2) - 1; k <= Math.ceil(b[1]/h2) + 1; k++) lines.push(k*h2); }
      else lines.push(0);
      for (const im of lines) { const y = pv.toScreen([0, im])[1]; if (y < -5 || y > H + 5) continue;
        g.beginPath(); g.moveTo(0, y); g.lineTo(W, y); g.strokeStyle = 'rgba(0,0,0,.75)'; g.lineWidth = 4.5; g.stroke();
        g.strokeStyle = '#00e5ff'; g.lineWidth = 2.2; g.stroke(); }
    }
    if (torus && opt.gn && C.edges) {
      for (const o of offs) for (const e of C.edges) {
        const P = [e.P[0] + o[0], e.P[1] + o[1]], Q = [e.Q[0] + o[0], e.Q[1] + o[1]], a = pv.toScreen(P), b = pv.toScreen(Q);
        if (Math.max(a[0], b[0]) < 0 || Math.min(a[0], b[0]) > W || Math.max(a[1], b[1]) < 0 || Math.min(a[1], b[1]) > H) continue;
        g.beginPath(); g.moveTo(...a); g.lineTo(...b); g.strokeStyle = 'rgba(0,0,0,.7)'; g.lineWidth = 4; g.stroke();
        g.strokeStyle = '#ffffff'; g.lineWidth = 2; g.stroke();
      }
    }
    if (opt.mk) {
      g.font = '600 12px system-ui'; g.textBaseline = 'middle';
      for (const o of offs) for (const m of markers) {
        const s = pv.toScreen([m.p[0] + o[0], m.p[1] + o[1]]); if (s[0] < -10 || s[1] < -10 || s[0] > W + 10 || s[1] > H + 10) continue;
        g.beginPath(); g.arc(s[0], s[1], 5, 0, 2*Math.PI); g.fillStyle = m.kind === 'centre' ? '#fff' : '#e8173a'; g.fill();
        g.strokeStyle = '#141414'; g.lineWidth = 1.2; g.stroke();
        if (opt.ml) { g.lineWidth = 3; g.strokeStyle = 'rgba(0,0,0,.85)'; g.strokeText(m.label, s[0] + 8, s[1]);
          g.fillStyle = m.kind === 'centre' ? '#fff' : '#ffd0d8'; g.fillText(m.label, s[0] + 8, s[1]); }
      }
    }
    if (RAYS.length) drawRays(g, W, H);
    if (picked) { const s = pv.toScreen(picked); g.strokeStyle = '#fff'; g.lineWidth = 2;
      g.beginPath(); g.moveTo(s[0] - 9, s[1]); g.lineTo(s[0] + 9, s[1]); g.moveTo(s[0], s[1] - 9); g.lineTo(s[0], s[1] + 9); g.stroke(); }
  }
  // ---- external rays (manifest rays: true; z^2 + c only)
  let RAYS = [];
  const RAYC = ['#ffe650', '#7cf', '#f9a', '#9f9', '#fc6', '#c9f', '#6ff', '#f66'];
  function drawRays(g, W, H) {
    RAYS.forEach((r, i) => { const col = RAYC[i % RAYC.length];
      g.beginPath(); r.pts.forEach((q, j) => { const s = pv.toScreen(q); j ? g.lineTo(...s) : g.moveTo(...s); });
      g.strokeStyle = 'rgba(0,0,0,.7)'; g.lineWidth = 3.5; g.stroke(); g.strokeStyle = col; g.lineWidth = 1.8; g.stroke();
      if (r.jump !== undefined) { const a = pv.toScreen(r.pts[r.pts.length - 1]), b = pv.toScreen(r.landing);
        g.setLineDash([3, 3]); g.beginPath(); g.moveTo(...a); g.lineTo(...b); g.stroke(); g.setLineDash([]); }
      const s = pv.toScreen(r.landing); g.beginPath(); g.arc(s[0], s[1], 4, 0, 2*Math.PI); g.fillStyle = col; g.fill(); g.strokeStyle = '#141414'; g.lineWidth = 1; g.stroke();
      if (opt.ml) { g.font = '600 12px system-ui'; g.lineWidth = 3; g.strokeStyle = 'rgba(0,0,0,.85)'; g.strokeText(r.a.text, s[0] + 6, s[1] - 8); g.fillStyle = col; g.fillText(r.a.text, s[0] + 6, s[1] - 8); } });
  }
  function setRays(text) {
    const angs = text.split(',').map(t => PM.parseAngle(t)).filter(Boolean);
    RAYS = angs.map(a => ({ a, ...PM.rayWithLanding('param', a, null) }));
    $('rayout').textContent = RAYS.map(r => `${r.a.text}: lands at c = ${PM.fmt(r.landing, 10)}\n  (${r.how}${r.period ? ', angle period ' + r.period : ', preperiodic angle'})`).join('\n');
    parent.postMessage({ e96: 'rays', angles: angs.map(a => a.text) }, '*'); pv.redraw();
  }
  function showCard(title, lines) {
    const esc = t => String(t).replace(/&/g, '&amp;').replace(/</g, '&lt;');
    $('cardt').innerHTML = '<b>' + esc(title) + '</b>\n' + esc(lines.join('\n')); $('card').style.display = 'block';
  }
  $('cardx').onclick = () => { $('card').style.display = 'none'; };
  function pick(p, marker) {
    picked = p; pv.redraw();
    const M = PM.mapFrom(E.family, E.n, { p });
    const lines = [];
    let x = null, y = null;
    if (torus) { const P = PM.torusXY(p); x = P.x; y = P.y; lines.push(`u = ${PM.fmt(p)}  (flat coordinate)`, `y = a3 = ${PM.fmt(y)}`, `x = a4 = ${PM.fmt(x)}`); }
    else lines.push(`${E.param} = ${PM.fmt(p)}`);
    if (!M) lines.push('degenerate map here (puncture)');
    else { const f = PM.fate(M, M.omega, pv.maxIt, pv.delta, 1e-7); lines.push('free critical point: ' + PM.fateText(f)); }
    if (marker) lines.unshift(...marker.lines, '');
    showCard(marker ? marker.title : 'Parameter', lines);
    parent.postMessage({ e96: 'param', model: id, p, x, y, label: marker ? marker.title : null,
      marker: marker ? { kind: marker.kind, label: marker.label, title: marker.title, lines: marker.lines.concat(lines.slice(marker.lines.length + 1)),
                         ks: marker.ks || null, thetas: marker.thetas || null } : null }, '*');
  }
  const pv = new PlaneView($('view'), {
    ctr: E.view.ctr, upp: E.view.upp, maxIt: 600, delta: 1e-3,
    setUniforms: (gl, U) => { gl.uniform1i(U('fam'), E.famId); GLSL.setPalette(gl, U, E.palette, E.flat);
      gl.uniform1f(U('W1'), T5.W1); gl.uniform1f(U('W2'), T5.W2); gl.uniform1f(U('qn'), Math.exp(-Math.PI*T5.W2/T5.W1)); gl.uniform1f(U('e1'), 0.75); },
    drawOverlay,
    onPick: (p, s) => {
      let best = null, bd = 12;
      if (opt.mk) { const [W, H] = pv.size(); for (const o of torus ? lattice(W, H) : [[0, 0]]) for (const m of markers) {
        const q = pv.toScreen([m.p[0] + o[0], m.p[1] + o[1]]), d = Math.hypot(q[0] - s[0], q[1] - s[1]); if (d < bd) { bd = d; best = { m, o }; } } }
      if (best) pick([best.m.p[0] + best.o[0], best.m.p[1] + best.o[1]], best.m); else pick(p, null);
    } });
  window.PV = pv;
  $('plain').onchange = e => { pv.plain = e.target.value === '1'; pv.redraw(); };
  $('bands').onchange = e => { pv.bands = e.target.checked; pv.redraw(); };
  $('pit').oninput = e => { pv.maxIt = +e.target.value; $('pitv').textContent = pv.maxIt; pv.redraw(); };
  $('del').oninput = e => { pv.delta = Math.pow(10, +e.target.value); $('delv').textContent = pv.delta.toExponential(0); pv.redraw(); };
  $('home').onclick = () => { pv.ctr = E.view.ctr.slice(); pv.upp = E.view.upp; pv.redraw(); };
  function buildList() {
    const q = $('mfilter').value.trim().toLowerCase(), box = $('mlist'); box.innerHTML = '';
    const all = C.centres.map(c => ({ ...c, kind: 'centre' })).concat(C.punctures.map(p => ({ ...p, kind: 'puncture' })));
    for (const it of all) {
      const text = `${it.kind} ${it.label}: ${it.summary || ''}`; if (q && !text.toLowerCase().includes(q)) continue;
      const d = document.createElement('div'); d.textContent = text; d.title = it.lines.join('\n');
      d.onclick = () => { if (it.p) { pv.ctr = it.p.slice(); pv.redraw(); pick(it.p, it); } else showCard(it.title, it.lines); };
      box.appendChild(d);
    }
  }
  $('mfilter').oninput = buildList; buildList();
  window.addEventListener('message', ev => { const d = ev.data;
    if (d && d.e96 === 'goto' && d.label !== undefined) { const m = markers.find(q => q.kind === 'centre' && q.label === String(d.label)); if (m) { pv.ctr = m.p.slice(); pv.redraw(); pick(m.p, m); } return; }
    if (d && d.e96 === 'goto' && d.p) { pv.ctr = d.p.slice(); pv.redraw(); pick(d.p, null); } });
  if (E.rays) { $('raybox').style.display = ''; $('raygo').onclick = () => setRays($('rayin').value);
    $('rayin').addEventListener('keydown', e => { if (e.key === 'Enter') setRays(e.target.value); });
    $('rayclr').onclick = () => { RAYS = []; $('rayout').textContent = ''; parent.postMessage({ e96: 'rays', angles: [] }, '*'); pv.redraw(); }; }
  if (window.innerWidth < 650) $('minb').onclick();
  window.E96READY = true;
})().catch(e => { document.body.innerHTML = '<pre style="color:#fff;padding:20px">' + e.stack + '</pre>'; console.error(e); });
