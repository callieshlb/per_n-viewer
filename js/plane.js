// E96: parameter pane for kind = plane (Per_1 .. Per_4) and kind = torus (Per_5). URL: plane.html?m=<model id>
'use strict';
(async function () {
  const id = new URLSearchParams(location.search).get('m') || 'per1';
  const man = await fetch('models/manifest.json?v=20261110013643').then(r => r.json());
  const E = man.models.find(m => m.id === id);
  const C = await fetch(E.census + '?v=20261110013643').then(r => r.json());
  const torus = E.kind === 'torus', T5 = PM.T5;
  const $ = s => document.getElementById(s);
  $('title').innerHTML = E.title; $('nf').innerHTML = E.normal; $('about').innerHTML = E.about || '';
  document.title = E.titleText || id;
  $('legend').innerHTML = GLSL.LEGEND(E.n, E.palette).map(([c, t]) => `<span class="sw" style="background:${c}"></span><span>${t}</span>`).join('');
  if (torus) ['gnl', 'gn', 'vil', 'vi', 'sml', 'sm'].forEach(k => { $(k).style.display = ''; });
  const opt = { mk: true, ml: true, gn: false, vi: false, sm: true, rc: true };
  if (E.real) ['rcl', 'rc'].forEach(k => { $(k).style.display = ''; });
  ['mk', 'ml', 'gn', 'vi', 'sm', 'rc'].forEach(k => { $(k).onchange = e => { opt[k] = e.target.checked; pv.redraw(); }; });
  $('minb').onclick = () => { $('side').classList.toggle('min'); $('minb').textContent = $('side').classList.contains('min') ? '+' : '–'; };
  let picked = null; const MARKS = [];
  function lattice(W, H) {           // visible lattice translates (torus)
    const a = pv.toPlane(0, H), b = pv.toPlane(W, 0), out = [];
    const m0 = Math.floor(a[0]/T5.W1) - 1, m1 = Math.ceil(b[0]/T5.W1) + 1, l0 = Math.floor(a[1]/T5.W2) - 1, l1 = Math.ceil(b[1]/T5.W2) + 1;
    if ((m1 - m0)*(l1 - l0) > 400) return [[0, 0]];
    for (let m = m0; m <= m1; m++) for (let l = l0; l <= l1; l++) out.push([m*T5.W1, l*T5.W2]);
    return out;
  }
  const markers = C.centers.map(c => ({ ...c, kind: 'centre' })).concat(C.punctures.filter(p => p.p).map(p => ({ ...p, kind: 'puncture' })));
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
        const comp = torus ? ((Math.round(im/(T5.W2/2)) % 2) + 2) % 2 : 0;           // torus: Im u = 0 and Im u = W2/2 are the two real circles
        g.beginPath(); g.moveTo(0, y); g.lineTo(W, y); g.strokeStyle = 'rgba(0,0,0,.75)'; g.lineWidth = 4.5; g.stroke();
        g.strokeStyle = ['#00e5ff', '#ffdb1a'][comp]; g.lineWidth = 2.2; g.stroke(); }
    }
    if (torus && opt.gn && C.edges) {
      for (const o of offs) for (const e of C.edges) {
        const P = [e.P[0] + o[0], e.P[1] + o[1]], Q = [e.Q[0] + o[0], e.Q[1] + o[1]], a = pv.toScreen(P), b = pv.toScreen(Q);
        if (Math.max(a[0], b[0]) < 0 || Math.min(a[0], b[0]) > W || Math.max(a[1], b[1]) < 0 || Math.min(a[1], b[1]) > H) continue;
        g.beginPath(); g.moveTo(...a); g.lineTo(...b); g.strokeStyle = 'rgba(0,0,0,.7)'; g.lineWidth = 4; g.stroke();
        g.strokeStyle = '#ffffff'; g.lineWidth = 2; g.stroke();
      }
    }
    if (torus && opt.vi && C.edges) {          // violet: lo -> hi with v2(hi) > v2(lo) and v2(hi) maximal among lo's neighbors (ties kept)
      const v2 = k => { k = +k; let r = 0; while (k > 0 && k % 2 === 0) { k /= 2; r++; } return r; }, MV = {};
      for (const e of C.edges) { MV[e.a] = Math.max(MV[e.a] ?? -1, v2(e.b)); MV[e.b] = Math.max(MV[e.b] ?? -1, v2(e.a)); }
      for (const o of offs) for (const e of C.edges) {
        const d = v2(e.b) - v2(e.a); if (!d) continue;
        const lo = d > 0 ? e.a : e.b, hi = d > 0 ? e.b : e.a; if (v2(hi) < MV[lo]) continue;
        let P = [e.P[0] + o[0], e.P[1] + o[1]], Q = [e.Q[0] + o[0], e.Q[1] + o[1]]; if (d < 0) [P, Q] = [Q, P];
        const a = pv.toScreen(P), b = pv.toScreen(Q);
        if (Math.max(a[0], b[0]) < 0 || Math.min(a[0], b[0]) > W || Math.max(a[1], b[1]) < 0 || Math.min(a[1], b[1]) > H) continue;
        const L = Math.hypot(b[0] - a[0], b[1] - a[1]); if (L < 1) continue;
        const ux = (b[0] - a[0])/L, uy = (b[1] - a[1])/L, tip = [b[0] - 7*ux, b[1] - 7*uy];
        g.beginPath(); g.moveTo(...a); g.lineTo(...tip); g.strokeStyle = 'rgba(0,0,0,.7)'; g.lineWidth = 5; g.stroke();
        g.strokeStyle = '#b388ff'; g.lineWidth = 2.6; g.stroke();
        g.beginPath(); g.moveTo(...tip); g.lineTo(tip[0] - 11*ux + 5*uy, tip[1] - 11*uy - 5*ux); g.lineTo(tip[0] - 11*ux - 5*uy, tip[1] - 11*uy + 5*ux); g.closePath();
        g.fillStyle = '#b388ff'; g.fill(); g.strokeStyle = '#141414'; g.lineWidth = 1; g.stroke();
      }
    }
    if (opt.mk) {
      g.font = '600 12px system-ui'; g.textBaseline = 'middle';
      for (const o of offs) for (const m of markers) {
        const s = pv.toScreen([m.p[0] + o[0], m.p[1] + o[1]]); if (s[0] < -10 || s[1] < -10 || s[0] > W + 10 || s[1] > H + 10) continue;
        g.beginPath(); g.arc(s[0], s[1], 5, 0, 2*Math.PI); g.fillStyle = m.kind === 'centre' ? (E.n >= 5 && m.label === String(2**(E.n - 1)) ? '#b388ff' : '#fff') : '#e8173a'; g.fill();
        g.strokeStyle = '#141414'; g.lineWidth = 1.2; g.stroke();
        if (opt.ml) { g.lineWidth = 3; g.strokeStyle = 'rgba(0,0,0,.85)'; g.strokeText(m.label, s[0] + 8, s[1]);
          g.fillStyle = m.kind === 'centre' ? '#fff' : '#ffd0d8'; g.fillText(m.label, s[0] + 8, s[1]); }
      }
    }
    if (RAYS.length) drawRays(g, W, H);
    if (E.rays) drawPerN(g, W, H);
    for (const mk of MARKS) for (const o of torus ? lattice(W, H) : [[0, 0]]) {        // matings (yellow)
      const s = pv.toScreen([mk.p[0] + o[0], mk.p[1] + o[1]]); if (s[0] < -10 || s[1] < -10 || s[0] > W + 10 || s[1] > H + 10) continue;
      g.beginPath(); for (let k = 0; k < 10; k++) { const r = k % 2 ? 3 : 7, t = -Math.PI/2 + k*Math.PI/5; g.lineTo(s[0] + r*Math.cos(t), s[1] + r*Math.sin(t)); } g.closePath();
      g.fillStyle = '#ffe600'; g.fill(); g.strokeStyle = '#141414'; g.lineWidth = 1.2; g.stroke();
      if (opt.ml) { g.font = '600 12px system-ui'; g.textBaseline = 'middle'; g.lineWidth = 3; g.strokeStyle = 'rgba(0,0,0,.85)'; g.strokeText(mk.label, s[0] + 9, s[1]); g.fillStyle = '#ffe600'; g.fillText(mk.label, s[0] + 9, s[1]); }
    }
    if (picked) { const s = pv.toScreen(picked); g.strokeStyle = '#fff'; g.lineWidth = 2;
      g.beginPath(); g.moveTo(s[0] - 9, s[1]); g.lineTo(s[0] + 9, s[1]); g.moveTo(s[0], s[1] - 9); g.lineTo(s[0], s[1] + 9); g.stroke(); }
  }
  // ---- Per_1: the centers of period n (poly_info_n{n}.json from poly_info_n.py; n = 1 is c = 0)
  let PERN = [], PERNn = 0;
  async function setPerN(n) {
    PERNn = n; PERN = [];
    if (n === 1) PERN = [{ p: [0, 0], ks: [0, 0], label: '0', n: 1, title: 'center of period 1: c = 0', lines: ['z^2, angle 0'] }];
    else if (n >= 2) { const D = await fetch(`poly/poly_info_n${n}.json`).then(r => r.ok ? r.json() : null).catch(() => null); if (!D) return; const Q = 2**n - 1;
      PERN = D.map(P => ({ p: [Number(P.c_str[0]), Number(P.c_str[1])], ks: P.ks, label: String(P.label), n, kind: 'centre',
        title: `center ${P.label} of period ${n} (angles ${P.ks[0]}/${Q}, ${P.ks[1]}/${Q})`,
        lines: [`z^2 + c,  c = ${P.c_str[0]} ${P.c_str[1].startsWith('-') ? '-' : '+'} ${P.c_str[1].replace('-', '')}i`, `angles ${P.ks[0]}/${Q}, ${P.ks[1]}/${Q}`,
                `kneading ${P.kneading}`, `internal address ${P.address.join(' -> ')}`] })); }
    pv.redraw();
  }
  let PERSEL = null;                                   // the selected period-n center (its labels and rays highlighted)
  function perRays(m) {                                // parameter rays of the two characteristic angles (traced once, cached)
    if (!m.rays && PERNn >= 2) { const Q = BigInt(2**PERNn - 1); m.rays = m.ks.map(k => PM.traceRay('param', { p: BigInt(k), q: Q, text: `${k}/${2**PERNn - 1}` }, null, 8, 30*8).pts); }
    return m.rays || [];
  }
  function drawPerN(g, W, H) {
    if (!PERN.length) return; const Q = 2**PERNn - 1; g.font = '600 11px system-ui'; g.textBaseline = 'middle';
    const label = (m, s, dy) => { const t = PERNn === 1 ? '0' : `${m.ks[0]}/${Q}, ${m.ks[1]}/${Q}`; g.lineWidth = 3; g.strokeStyle = 'rgba(0,0,0,.85)'; g.strokeText(t, s[0] + 6, s[1] + dy); g.fillStyle = '#fff'; g.fillText(t, s[0] + 6, s[1] + dy); };
    if ($('perrays').checked) for (const m of PERN) { const sel = m === PERSEL;
      for (const r of perRays(m)) { g.beginPath(); r.forEach((q, i) => { const s = pv.toScreen(q); i ? g.lineTo(...s) : g.moveTo(...s); });
        g.strokeStyle = sel ? '#ffe650' : 'rgba(255,255,255,.55)'; g.lineWidth = sel ? 2.2 : 1; g.stroke(); } }
    const placed = [];                                   // labels: no two within 14 px vertically and 90 px horizontally (staggered)
    for (const m of PERN) { const s = pv.toScreen(m.p); if (s[0] < -10 || s[1] < -10 || s[0] > W + 10 || s[1] > H + 10) continue;
      g.beginPath(); g.arc(s[0], s[1], m === PERSEL ? 5.5 : 4, 0, 2*Math.PI); g.fillStyle = m === PERSEL ? '#ffe650' : '#fff'; g.fill(); g.strokeStyle = '#141414'; g.lineWidth = 1.2; g.stroke();
      if (m === PERSEL || $('perang').checked) { let dy = 0;
        for (let tries = 0; tries < 6 && placed.some(q => Math.abs(q[0] - s[0]) < 90 && Math.abs(q[1] - (s[1] + dy)) < 14); tries++) dy += (tries % 2 ? -1 : 1)*14*(tries + 1);
        if (m === PERSEL || !placed.some(q => Math.abs(q[0] - s[0]) < 90 && Math.abs(q[1] - (s[1] + dy)) < 14)) { label(m, s, dy); placed.push([s[0], s[1] + dy]); } } }
  }
  function hitPerRay(s) {                              // a period-n center whose ray passes within 6 px of screen point s
    if (!$('perrays').checked) return null; let best = null, bd = 6;
    for (const m of PERN) for (const r of m.rays || []) for (let i = 1; i < r.length; i++) {
      const a = pv.toScreen(r[i - 1]), b = pv.toScreen(r[i]), dx = b[0] - a[0], dy = b[1] - a[1], L2 = dx*dx + dy*dy || 1;
      const t = Math.max(0, Math.min(1, ((s[0] - a[0])*dx + (s[1] - a[1])*dy)/L2)), d = Math.hypot(a[0] + t*dx - s[0], a[1] + t*dy - s[1]);
      if (d < bd) { bd = d; best = m; } }
    return best;
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
    let x = null, y = null;
    if (torus) { const P = PM.torusXY(p); x = P.x; y = P.y; }
    // Caroline 2026-11-10: no card for plain parameters (the status line below shows the map and the free critical point's fate);
    // a marker shows its census data only
    if (marker) showCard(marker.title, marker.lines); else $('card').style.display = 'none';
    parent.postMessage({ e96: 'param', model: id, p, x, y, label: marker ? marker.title : null,
      marker: marker ? { kind: marker.kind, label: marker.label, n: marker.n, title: marker.title, lines: marker.lines,
                         ks: marker.ks || null, thetas: marker.thetas || null, pc: marker.pc || null } : null }, '*');
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
      if (E.rays) for (const m of PERN) { const q = pv.toScreen(m.p), d = Math.hypot(q[0] - s[0], q[1] - s[1]); if (d < bd) { bd = d; best = { m, o: [0, 0] }; } }
      if (E.rays && !best) { const r = hitPerRay(s); if (r) best = { m: r, o: [0, 0] }; }
      PERSEL = best && PERN.includes(best.m) ? best.m : null;
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
    const all = C.centers.map(c => ({ ...c, kind: 'centre' })).concat(C.punctures.map(p => ({ ...p, kind: 'puncture' })));
    for (const it of all) {
      const text = `${it.kind === 'centre' ? 'center' : it.kind} ${it.label}: ${it.summary || ''}`; if (q && !text.toLowerCase().includes(q)) continue;
      const d = document.createElement('div'); d.textContent = text; d.title = it.lines.join('\n');
      d.onclick = () => { if (it.p) { pv.ctr = it.p.slice(); pv.redraw(); pick(it.p, it); } else showCard(it.title, it.lines); };
      box.appendChild(d);
    }
  }
  $('mfilter').oninput = buildList; buildList();
  window.addEventListener('message', ev => { const d = ev.data;
    if (d && d.e96 === 'mark' && d.p) { MARKS.push({ p: d.p, label: d.label }); pv.redraw(); return; }
    if (d && d.e96 === 'goto' && d.label !== undefined) { const m = markers.find(q => q.kind === 'centre' && q.label === String(d.label)); if (m) { pv.ctr = m.p.slice(); pv.redraw(); pick(m.p, m); } return; }
    if (d && d.e96 === 'goto' && d.p) { pv.ctr = d.p.slice(); pv.redraw(); pick(d.p, null); } });
  if (E.rays) { $('perbox').style.display = ''; $('pern').onchange = e => setPerN(+e.target.value); $('perang').onchange = () => pv.redraw(); $('perrays').onchange = () => pv.redraw();
    $('raybox').style.display = ''; $('raygo').onclick = () => setRays($('rayin').value);
    $('rayin').addEventListener('keydown', e => { if (e.key === 'Enter') setRays(e.target.value); });
    $('rayclr').onclick = () => { RAYS = []; $('rayout').textContent = ''; parent.postMessage({ e96: 'rays', angles: [] }, '*'); pv.redraw(); }; }
  $('minb').onclick();                      // panel starts collapsed (Caroline 2026-11-10); + expands
  if (E.n < 2) $('matebox').style.display = 'none';
  $('mgo').onclick = () => parent.postMessage({ e96: 'mateReq', label: $('mf').value, angle: $('mg').value }, '*');
  $('mg').addEventListener('keydown', e => { if (e.key === 'Enter') $('mgo').click(); });
  window.addEventListener('keydown', e => { if (e.key === 'm' && e.target.tagName !== 'INPUT' && e.target.tagName !== 'SELECT') parent.postMessage({ e96: 'key', key: 'm' }, '*'); });
  window.E96READY = true;
})().catch(e => { document.body.innerHTML = '<pre style="color:#fff;padding:20px">' + e.stack + '</pre>'; console.error(e); });
