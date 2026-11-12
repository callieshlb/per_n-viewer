// E96: parameter pane for kind = plane (Per_1 .. Per_4) and kind = torus (Per_5). URL: plane.html?m=<model id>
'use strict';
(async function () {
  const id = new URLSearchParams(location.search).get('m') || 'per1';
  const man = await fetch('models/manifest.json?v=20261112225441').then(r => r.json());
  const E = man.models.find(m => m.id === id);
  const C = await fetch(E.census + '?v=20261112225441').then(r => r.json());
  const torus = E.kind === 'torus', T5 = PM.T5;
  const $ = s => document.getElementById(s);
  $('title').innerHTML = E.title; $('nf').innerHTML = E.normal; $('about').innerHTML = E.about || '';
  document.title = E.titleText || id;
  $('legend').innerHTML = GLSL.LEGEND(E.n, E.palette).map(([c, t]) => `<span class="sw" style="background:${c}"></span><span>${t}</span>`).join('');
  if (torus) ['gnl', 'gn', 'vil', 'vi', 'sml', 'sm'].forEach(k => { $(k).style.display = ''; });
  const opt = { mk: true, ml: true, gn: false, vi: false, sm: true, rc: false };      // real circles opt-in (Caroline 2026-11-10)
  if (E.real) ['rcl', 'rc'].forEach(k => { $(k).style.display = ''; });
  ['mk', 'ml', 'gn', 'vi', 'sm', 'rc'].forEach(k => { $(k).onchange = e => { opt[k] = e.target.checked; pv.redraw(); }; });
  $('pcl').innerHTML = GLSL.PERPAL.map((c, i) => `<span style="display:inline-block;min-width:13px;text-align:center;font-size:10px;background:rgb(${c.join(',')})">${i + 2}</span>`).join('');
  $('pc').onchange = e => { pv.perCol = e.target.checked; pv.redraw(); };
  // hotkeys (Caroline 2026-11-10): r real circles, g G_n edges, e violet arrows, p period colors, l labels, k markers, b bands, s seams, +/- iterates; t, i go to the dynamical pane
  const PKEY = { r: 'rc', g: 'gn', e: 'vi', p: 'pc', l: 'ml', k: 'mk', o: 'bands', s: 'sm', b: 'bt', j: 'bt' };   // b = bitransitive centers (Caroline 2026-11-12; j kept), o = escape-time bands (was b)
  function paneKey(key) {
    if (PKEY[key]) { const el = $(PKEY[key]); if (!el || el.style.display === 'none') return true; el.checked = !el.checked; el.dispatchEvent(new Event('change')); return true; }
    if (key === '+' || key === '=' || key === '-') { const el = $('pit'); el.value = Math.max(+el.min, Math.min(+el.max, +el.value + (key === '-' ? -1 : 1)*100)); el.dispatchEvent(new Event('input')); return true; }
    if (key === 'a') { $('side').classList.remove('min'); $('minb').textContent = '–'; window.focus(); $('bta').focus(); $('bta').select(); return true; }
    if (key === 't' || key === 'i' || key === 'u' || key === 'c' || key === 'h') { parent.postMessage({ e96: 'key', key }, '*'); return true; }
    return false;
  }
  window.addEventListener('keydown', e => { if (e.target.tagName === 'INPUT' && e.target.type !== 'checkbox' && e.target.type !== 'range' || e.target.tagName === 'SELECT' || e.metaKey || e.ctrlKey || e.altKey) return; if (paneKey(e.key)) e.preventDefault(); });
  window.addEventListener('message', ev => { const d = ev.data; if (d && d.e96 === 'key' && 'rgepjab'.includes(d.key)) paneKey(d.key); });
  window.E96TEST_bt = pk => { if (!BT || !BT.length) return null; let b = BT[0]; if (pk) b = BT.reduce((u, v) => Math.hypot(v.p[0] - pk[0], v.p[1] - pk[1]) < Math.hypot(u.p[0] - pk[0], u.p[1] - pk[1]) ? v : u); pv.ctr = b.p.slice(); pv.redraw(); return pv.toScreen(b.p); };
  window.E96TEST_view = () => ({ ctr: pv.ctr, upp: pv.upp });
  window.E96TEST_marks = () => MARKS.filter(m => m.vein).length;
  window.E96TEST_webs = () => WEBL.map(w => ({ id: w.id, paths: w.paths.length, dots: w.dots.length }));
  window.E96TEST_gp = () => ({ path: GPATH, edges: (C.edges || []).filter(e => GPATH.some((L, i) => i && ((GPATH[i - 1] === e.a && L === e.b) || (GPATH[i - 1] === e.b && L === e.a)))).length });
  window.E96TEST_pern = () => ({ n: PERN.length, tuned: PERN.filter(m => m.tuned).map(m => ({ label: m.label, p: m.p, line: m.lines.find(l => /TUNING/.test(l)).slice(0, 70) })) });
  window.E96TEST_pickPuncture = i => { const m = markers.filter(q => q.kind === 'puncture')[i]; if (m) pick(m.p, m); return !!m; };
  window.E96TEST_compr = () => COMPR && COMPR.rays.length === 2 ? { apart: +Math.hypot(COMPR.rays[0].at(-1)[0] - COMPR.rays[1].at(-1)[0], COMPR.rays[0].at(-1)[1] - COMPR.rays[1].at(-1)[1]).toExponential(1), fromCenter: +Math.hypot(COMPR.rays[0].at(-1)[0] - COMPR.c[0], COMPR.rays[0].at(-1)[1] - COMPR.c[1]).toExponential(1) } : null;
  $('minb').onclick = () => { $('side').classList.toggle('min'); $('minb').textContent = $('side').classList.contains('min') ? '+' : '–'; };
  let picked = null; const MARKS = [];
  // bitransitive centers from periodic captures (REGISTRATION addendum 2026-11-10 (h)); key j
  let BT = null, BTON = false;
  const negA = a => { const m = String(a).split('/'); return `−${+m[1] - +m[0]}/${m[1]}`; };
  const btText = b => `bitransitive center (ω = a${b.slot}): periodic captures ` + b.src.map(s => `${s.f} ^ ${s.angle} = ${s.f} ^ ${negA(s.angle)} [x${s.i}, ${s.side}]`).join(';  ');   // f ^ theta (Caroline 2026-11-12)
  // the three triangles of a bitransitive component (REGISTRATION addendum 2026-11-11 (j)): overlay image for the clicked center
  let TRI = null;
  function triangles(b) {
    const mapOf = p => { try { return PM.mapFrom(E.family, E.n, { p }); } catch (e) { return null; } }, R = BTRI.fill(mapOf, b.p, b.slot, { N: 281 }); if (R.fail) { TRI = null; return 'triangles: ' + R.fail; }
    const cvs = document.createElement('canvas'); cvs.width = cvs.height = R.N; const g = cvs.getContext('2d'), base = ((E.palette || GLSL.DEFPAL).concat(GLSL.DEFPAL))[(b.slot - 1) % 8], col = [base.map(v => Math.round(v + 0.45*(255 - v))), base, base.map(v => Math.round(0.55*v))];   // three shades of the component's own phase color (Caroline 2026-11-11)
    // paint(hl): hl(argW in degrees [0, 360)) -> true for the highlighted part (addendum (v)); null = no highlight
    const paint = hl => { const im = g.createImageData(R.N, R.N);
      for (let iy = 0; iy < R.N; iy++) for (let ix = 0; ix < R.N; ix++) { const i = iy*R.N + ix, o = 4*((R.N - 1 - iy)*R.N + ix); if (R.st[i] !== 1) continue; const wx = R.W[2*i], wy = R.W[2*i + 1], m = Math.hypot(wx, wy), [k, t] = BTRI.tri(wx, wy);
        let c = col[k], a = 215; const d0 = Math.min(t, 1 - t), dh = Math.abs(t - 0.5), ph = m*m*m;
        if (d0*m < 0.006) { c = [255, 255, 255]; a = 255; }                              // internal rays of angle 0 (triangle sides)
        else if (dh*m < 0.004) { c = [20, 20, 20]; a = 220; }                            // internal rays of angle 1/2
        else if (Math.abs(ph - 0.5) < 0.012 || Math.abs(ph - 0.25) < 0.008 || Math.abs(ph - 0.75) < 0.012) { c = c.map(v => Math.round(0.8*v)); a = 235; }   // |Phi| = 1/4, 1/2, 3/4
        if (hl) { let ar = Math.atan2(wy, wx)*180/Math.PI; if (ar < 0) ar += 360; if (hl(ar)) { c = c.map(v => Math.min(255, Math.round(v*0.55 + 255*0.45))); a = 255; } else a = Math.round(a*0.3); }
        im.data[o] = c[0]; im.data[o + 1] = c[1]; im.data[o + 2] = c[2]; im.data[o + 3] = a; }
      g.putImageData(im, 0, 0); };
    paint(null); TRI = { cvs, x0: R.x0, y0: R.y0, size: (R.N - 1)*R.dx, b, kappa: R.kappa, dx: R.dx, mapOf, paint, ends: null }; IRAYS = [];
    return `triangles (W = Φ^(1/3)): light / medium / dark shade of the component's phase color, white lines = internal rays of angle 0, dark = angle 1/2, rings |Φ| = 1/4, 1/2, 3/4; ${R.count} grid points, sup |W| = ${R.sup.toFixed(3)}`;
  }
  // (d) boundary point of the selected bitransitive component from an internal angle and a triangle index (REGISTRATION addendum 2026-11-11 (k)); key a
  let IRAYS = [];
  function internalRay(text) {
    const out = $('btao');
    if (!TRI) { out.textContent = 'select a bitransitive center first (toggle j, click a diamond)'; return; }
    const m = String(text).trim().match(/^(\S+)(?:\s+([012]))?$/), ang = m ? PM.parseAngle(m[1] === '0' ? '0/1' : m[1] === '1' ? '1/1' : m[1]) : null; if (!ang) { out.textContent = 'enter: angle [triangle 0, 1 or 2], e.g. 1/4 2'; return; }
    const ks = m[2] !== undefined ? [+m[2]] : [0, 1, 2], b = TRI.b, lines = []; let first = null;
    for (const k of ks) { const r = BTRI.boundaryPoint(TRI.mapOf, b.p, b.slot, k, ang, { kappa: TRI.kappa, dx: TRI.dx });
      if (r.fail) { lines.push(`triangle ${k}: ${r.fail}`); continue; } r.text = m[1]; IRAYS.push(r); if (!first) first = r;
      lines.push(`triangle ${k}, angle ${m[1]}: ` + (r.exact ? `boundary point ${PM.fmt(r.point, 10)} (critical orbit lands on a repelling cycle after ${r.pp[0]} return${r.pp[0] > 1 ? 's' : ''}; |multiplier| ${r.mult.toFixed(3)}, residual ${r.res.toExponential(1)}, ${r.gap.toExponential(1)} from the ray's end)`
        : (1 - r.rho > 1e-3 ? `NO boundary point: the ray could not be continued beyond |W| = ${r.rho.toFixed(3)}; it stops at ` : `ray followed to 1 − |W| = ${(1 - r.rho).toExponential(1)}: `) + `${PM.fmt(r.end, 8)}; still moving (${r.drift === null ? '?' : r.drift.toExponential(1)} over the last decade)` + (r.pp && r.pp[0] === 0 ? '; periodic angle: exact landing point not solved' : ''))); }
    out.textContent = lines.join('\n'); pv.redraw();
    if (first) { picked = first.point; let x = null, y = null; if (torus) { const P = PM.torusXY(first.point); x = P.x; y = P.y; }
      parent.postMessage({ e96: 'param', model: id, p: first.point, x, y, label: `internal ray of the bitransitive component ω = a${b.slot} at ${PM.fmt(b.p, 5)}\n` + lines.join('\n') }, '*'); }
  }
  function drawIRays(g, W, H, offs) { for (const o of offs) for (const r of IRAYS) { g.beginPath(); r.pts.forEach((q, i) => { const s = pv.toScreen([q[1] + o[0], q[2] + o[1]]); if (i) g.lineTo(s[0], s[1]); else g.moveTo(s[0], s[1]); });
      const e = pv.toScreen([r.point[0] + o[0], r.point[1] + o[1]]); g.lineTo(e[0], e[1]); g.lineWidth = 3.5; g.strokeStyle = 'rgba(0,0,0,.8)'; g.stroke(); g.lineWidth = 1.6; g.strokeStyle = '#ffe14d'; g.stroke();
      g.beginPath(); g.arc(e[0], e[1], 5, 0, 2*Math.PI); g.fillStyle = r.exact ? '#ffe14d' : 'rgba(0,0,0,0)'; g.fill(); g.lineWidth = 2; g.strokeStyle = r.exact ? '#141414' : '#ffe14d'; g.stroke();
      g.font = '500 10px system-ui'; g.textBaseline = 'middle'; g.lineWidth = 2.5; g.strokeStyle = 'rgba(0,0,0,.85)'; const t = `△${r.k} ${r.text}`; g.strokeText(t, e[0] + 8, e[1]); g.fillStyle = '#ffe14d'; g.fillText(t, e[0] + 8, e[1]); } }
  window.E96TEST_iray = t => { internalRay(t); return IRAYS.map(r => ({ k: r.k, t: r.t, exact: r.exact, point: r.point, rho: r.rho, n: r.pts.length })); };
  $('bta').addEventListener('keydown', e => { if (e.key === 'Enter') { internalRay($('bta').value); e.preventDefault(); } else if (e.key === 'Escape') $('bta').blur(); e.stopPropagation(); });
  $('btac').onclick = () => { IRAYS = []; $('btao').textContent = ''; pv.redraw(); };
  $('bt').onchange = async e => { BTON = e.target.checked; if (!BTON) { TRI = null; IRAYS = []; }
    if (BTON && !BT) { const D = await fetch(`models/bitrans_n${E.n}.json`).then(r => r.ok ? r.json() : []).catch(() => []); { const cnt = {}; D.forEach(b => { const k = cnt[b.slot] = (cnt[b.slot] || 0) + 1; b.name = 'b' + String(b.slot).replace(/\d/g, d => '⁰¹²³⁴⁵⁶⁷⁸⁹'[d]) + String(k).replace(/\d/g, d => '₀₁₂₃₄₅₆₇₈₉'[d]); }); }   /* b^j_k: j = slot (omega = a_j), k = running index within the slot, no meaning (Caroline 2026-11-12) */ BT = D.filter(b => b.p);   // b_1 .. b_N: the index is the position in the table and means nothing yet (Caroline 2026-11-12)
      $('btn').textContent = `${BT.length} shown` + (D.length > BT.length ? `, ${D.length - BT.length} at ∞ of this chart` : ''); }
    pv.redraw(); };
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
    if (torus && GPATH.length > 1 && C.edges) {   // shortest G_n path chosen in a center's card (lime)
      const want = new Set(); for (let i = 1; i < GPATH.length; i++) { want.add(GPATH[i - 1] + '|' + GPATH[i]); want.add(GPATH[i] + '|' + GPATH[i - 1]); }
      for (const o of offs) for (const e of C.edges) { if (!want.has(e.a + '|' + e.b)) continue; const a = pv.toScreen([e.P[0] + o[0], e.P[1] + o[1]]), b = pv.toScreen([e.Q[0] + o[0], e.Q[1] + o[1]]);
        if (Math.max(a[0], b[0]) < 0 || Math.min(a[0], b[0]) > W || Math.max(a[1], b[1]) < 0 || Math.min(a[1], b[1]) > H) continue;
        g.beginPath(); g.moveTo(...a); g.lineTo(...b); g.strokeStyle = 'rgba(0,0,0,.8)'; g.lineWidth = 6.5; g.stroke(); g.strokeStyle = '#8dff3a'; g.lineWidth = 3.5; g.stroke(); } }
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
        if (opt.ml && m.kind !== 'puncture') { g.lineWidth = 3; g.strokeStyle = 'rgba(0,0,0,.85)'; g.strokeText(m.label, s[0] + 8, s[1]);   // punctures: no number in the picture (Caroline 2026-11-12); the card has the data
          g.fillStyle = m.kind === 'centre' ? '#fff' : '#ffd0d8'; g.fillText(m.label, s[0] + 8, s[1]); }
      }
    }
    if (BTON && TRI) { for (const o of offs) { const a = pv.toScreen([TRI.x0 + o[0], TRI.y0 + TRI.size + o[1]]), q = pv.toScreen([TRI.x0 + TRI.size + o[0], TRI.y0 + o[1]]); if (q[0] < 0 || a[0] > W || q[1] < 0 || a[1] > H) continue; g.imageSmoothingEnabled = true; g.drawImage(TRI.cvs, a[0], a[1], q[0] - a[0], q[1] - a[1]); } }
    if (BTON && BT) { for (const o of offs) for (const b of BT) { const s = pv.toScreen([b.p[0] + o[0], b.p[1] + o[1]]); if (s[0] < -10 || s[1] < -10 || s[0] > W + 10 || s[1] > H + 10) continue;
        g.beginPath(); g.moveTo(s[0], s[1] - 5); g.lineTo(s[0] + 5, s[1]); g.lineTo(s[0], s[1] + 5); g.lineTo(s[0] - 5, s[1]); g.closePath(); g.fillStyle = '#00e5ff'; g.fill(); g.strokeStyle = '#141414'; g.lineWidth = 1; g.stroke();
        if (opt.ml) { const t = b.name; g.font = '500 10px system-ui'; g.textBaseline = 'middle'; g.lineWidth = 2.5; g.strokeStyle = 'rgba(0,0,0,.85)'; g.strokeText(t, s[0] + 7, s[1]); g.fillStyle = '#9ff3ff'; g.fillText(t, s[0] + 7, s[1]); } } }
    if (BTON && IRAYS.length) drawIRays(g, W, H, offs);
    if (COMPR) { const Q = 2**COMPR.n - 1; COMPR.rays.forEach((r, i) => { g.beginPath(); r.forEach((q, j) => { const s = pv.toScreen(q); j ? g.lineTo(...s) : g.moveTo(...s); }); g.strokeStyle = 'rgba(0,0,0,.7)'; g.lineWidth = 3.5; g.stroke(); g.strokeStyle = '#6ff'; g.lineWidth = 1.6; g.stroke();
        const e = pv.toScreen(r[Math.min(r.length - 1, Math.floor(r.length*0.35))]), t = `${COMPR.ks[i]}/${Q}`; g.font = '600 11px system-ui'; g.textBaseline = 'middle'; g.lineWidth = 3; g.strokeStyle = 'rgba(0,0,0,.85)'; g.strokeText(t, e[0] + 6, e[1]); g.fillStyle = '#9ff'; g.fillText(t, e[0] + 6, e[1]); });
      const s0 = pv.toScreen(COMPR.c); g.beginPath(); g.arc(s0[0], s0[1], 4.5, 0, 2*Math.PI); g.fillStyle = '#6ff'; g.fill(); g.strokeStyle = '#141414'; g.lineWidth = 1.2; g.stroke(); }
    for (const w of WEBL) for (const o of offs) {
      for (const pa of w.paths) { g.beginPath(); let prev = null; for (const q0 of pa.pts) { let q = q0;
          if (torus && prev) { const m = Math.round((prev[0] - q0[0])/T5.W1), l = Math.round((prev[1] - q0[1])/T5.W2); q = [q0[0] + m*T5.W1, q0[1] + l*T5.W2]; }      // nearest lattice translate: the chain is continuous on the torus
          const s = pv.toScreen([q[0] + o[0], q[1] + o[1]]); if (prev) g.lineTo(s[0], s[1]); else g.moveTo(s[0], s[1]); prev = q; }
        g.setLineDash(pa.dash || []); g.strokeStyle = 'rgba(0,0,0,.75)'; g.lineWidth = (pa.width || 2) + 2; g.stroke(); g.strokeStyle = pa.color; g.lineWidth = pa.width || 2; g.stroke(); g.setLineDash([]); }
      for (const d of w.dots) { const s = pv.toScreen([d.p[0] + o[0], d.p[1] + o[1]]); if (s[0] < -6 || s[1] < -6 || s[0] > W + 6 || s[1] > H + 6) continue;
        g.beginPath(); g.arc(s[0], s[1], d.r || 3, 0, 2*Math.PI); if (d.hollow) { g.fillStyle = 'rgba(0,0,0,.55)'; g.fill(); g.strokeStyle = d.color; g.lineWidth = 1.6; g.stroke(); } else { g.fillStyle = d.color; g.fill(); g.strokeStyle = '#141414'; g.lineWidth = 0.8; g.stroke(); } } }
    if (RAYS.length) drawRays(g, W, H);
    if (E.rays) drawPerN(g, W, H);
    const VL = [];                                        // vein labels already drawn (skip overlapping ones)
    for (const mk of MARKS) for (const o of torus ? lattice(W, H) : [[0, 0]]) {        // matings (yellow)
      const s = pv.toScreen([mk.p[0] + o[0], mk.p[1] + o[1]]); if (s[0] < -10 || s[1] < -10 || s[0] > W + 10 || s[1] > H + 10) continue;
      if (mk.vein) { g.beginPath(); g.arc(s[0], s[1], 3.2, 0, 2*Math.PI); g.fillStyle = '#ffe600'; g.fill(); g.strokeStyle = '#141414'; g.lineWidth = 1; g.stroke();
        if (opt.ml && !VL.some(q => Math.abs(q[0] - s[0]) < 46 && Math.abs(q[1] - s[1]) < 11)) { VL.push(s); g.font = '500 10px system-ui'; g.textBaseline = 'middle'; g.lineWidth = 2.5; g.strokeStyle = 'rgba(0,0,0,.85)'; g.strokeText(mk.label, s[0] + 6, s[1]); g.fillStyle = '#ffe600'; g.fillText(mk.label, s[0] + 6, s[1]); }
        continue; }
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
                `kneading ${P.kneading}`, `internal address ${P.address.join(' -> ')}`].concat((P.tuning || []).map(t => { const Qd = 2**t.d - 1, Qm = 2**t.inner_period - 1;
                  return `TUNING: (${t.base_ks[0]}/${Qd}, ${t.base_ks[1]}/${Qd})${t.base_name ? ' ' + t.base_name : ''} ⊥ (${t.inner_ks[0]}/${Qm}, ${t.inner_ks[1]}/${Qm})${t.inner_name ? ' ' + t.inner_name : ''}: in the small Mandelbrot copy of the period-${t.d} center ${t.base}`; })),
        tuned: (P.tuning || []).length > 0 })); }
    pv.redraw();
  }
  let COMPR = null;
  let GPATH = [];
  const WEBL = [];                                     // V_f / T'_f overlays (REGISTRATION addendum (u))
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
    for (const m of PERN) { if ($('pertun').checked && !m.tuned) continue; const s = pv.toScreen(m.p); if (s[0] < -10 || s[1] < -10 || s[0] > W + 10 || s[1] > H + 10) continue;
      g.beginPath(); g.arc(s[0], s[1], m === PERSEL ? 5.5 : (m.tuned ? 5 : 4), 0, 2*Math.PI); g.fillStyle = m === PERSEL ? '#ffe650' : (m.tuned ? '#ff7ad9' : '#fff'); g.fill();   // pink = a tuning (renormalizable) g.strokeStyle = '#141414'; g.lineWidth = 1.2; g.stroke();
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
      if (BTON && BT) { const [W0, H0] = pv.size(); for (const o of torus ? lattice(W0, H0) : [[0, 0]]) for (const b of BT) { const q = pv.toScreen([b.p[0] + o[0], b.p[1] + o[1]]);
          if (Math.hypot(q[0] - s[0], q[1] - s[1]) < 9) { picked = b.p; pv.redraw(); $('card').style.display = 'none'; let x = null, y = null; if (torus) { const P = PM.torusXY(b.p); x = P.x; y = P.y; }
            const tt = triangles(b); pv.redraw(); parent.postMessage({ e96: 'param', model: id, p: b.p, x, y, label: `bitransitive center (ω = a${b.slot})`, marker: { kind: 'bitrans', label: 'bt', slot: b.slot, src: b.src, title: `bitransitive center ${b.name}: ω = a${b.slot}`, lines: [tt] } }, '*'); return; } } }
      let best = null, bd = 12;
      if (opt.mk) { const [W, H] = pv.size(); for (const o of torus ? lattice(W, H) : [[0, 0]]) for (const m of markers) {
        const q = pv.toScreen([m.p[0] + o[0], m.p[1] + o[1]]), d = Math.hypot(q[0] - s[0], q[1] - s[1]); if (d < bd) { bd = d; best = { m, o }; } } }
      if (E.rays) for (const m of PERN) { if ($('pertun').checked && !m.tuned) continue; const q = pv.toScreen(m.p), d = Math.hypot(q[0] - s[0], q[1] - s[1]); if (d < bd) { bd = d; best = { m, o: [0, 0] }; } }
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
    if (d && d.e96 === 'clearMarks') { MARKS.length = 0; pv.redraw(); return; }
    if (d && d.e96 === 'gpath') { GPATH = (d.labels || []).map(String); pv.redraw(); return; }
    if (d && d.e96 === 'btHL') {                         // highlight the part of the selected bitransitive component attached to a capture (addendum (v), experimental rule)
      if (!TRI) return; if (!d.pts || !d.pts.length) { TRI.paint(null); pv.redraw(); return; }
      if (!TRI.ends) TRI.ends = [0, 1, 2].map(k => { const r = BTRI.ray(TRI.mapOf, TRI.b.p, TRI.b.slot, k, 0, { kappa: TRI.kappa, dx: TRI.dx }); return r.fail ? null : r.end; });
      const dist = (u, v) => { let dx = u[0] - v[0], dy = u[1] - v[1]; if (torus) { dx -= T5.W1*Math.round(dx/T5.W1); dy -= T5.W2*Math.round(dy/T5.W2); } return Math.hypot(dx, dy); };
      const info = d.pts.map(q => { const ds = TRI.ends.map(e => e ? dist(q, e) : Infinity), k = ds.indexOf(Math.min(...ds)), rest = ds.filter((_, i) => i !== k); return { k, d: ds[k], next: Math.min(...rest) }; });
      const ks = [...new Set(info.map(x => x.k))]; let hl;
      if (ks.length >= 2) { const a0 = 120*ks[0], a1 = 120*ks[1], lo = ((a1 - a0 + 360) % 360) === 120 ? a0 : a1; hl = ar => ((ar - lo + 360) % 360) <= 120; }      // the triangle between the two corners
      else { const c = 120*ks[0]; hl = ar => { const dd = Math.abs(((ar - c + 540) % 360) - 180); return dd <= 60; }; }                                         // the two hexants next to the corner
      TRI.paint(hl); pv.redraw(); parent.postMessage({ e96: 'btHLdone', key: d.key, info }, '*'); return; }
    if (d && d.e96 === 'webClear') { WEBL.length = 0; pv.redraw(); return; }
    if (d && d.e96 === 'web') { const k = WEBL.findIndex(w => w.id === d.id); if (k >= 0) WEBL.splice(k, 1); WEBL.push(d); pv.redraw(); return; }
    if (d && d.e96 === 'compRays') {                     // (s) the clicked component of Per_1: its center and the two rays landing at its root
      COMPR = { c: d.c, n: d.n, ks: d.ks, rays: [] }; if (d.ks && d.n >= 2) { const Q = BigInt(2**d.n - 1); COMPR.rays = d.ks.map(k => PM.traceRay('param', { p: BigInt(k), q: Q, text: `${k}/${2**d.n - 1}` }, null, 8, 30*8).pts); }
      if (d.go) { picked = d.c; } pv.redraw(); return; }
    if (d && d.e96 === 'mark' && d.p) { const k = MARKS.findIndex(q => q.label === d.label && !q.vein && !d.vein); if (k >= 0) MARKS.splice(k, 1);
      MARKS.push({ p: d.p, label: d.label, vein: !!d.vein }); if (d.go) pv.ctr = d.p.slice(); pv.redraw();
      if (d.go) parent.postMessage({ e96: 'markDone', label: d.label, ok: true, res: 0, tag: d.tag }, '*'); return; }
    if (d && d.e96 === 'goto' && d.label !== undefined) { const m = markers.find(q => q.kind === 'centre' && q.label === String(d.label)); if (m) { pv.ctr = m.p.slice(); pv.redraw(); pick(m.p, m); } return; }
    if (d && d.e96 === 'goto' && d.p) { pv.ctr = d.p.slice(); pv.redraw(); pick(d.p, null); } });
  if (E.rays) { $('perbox').style.display = ''; $('pern').onchange = e => setPerN(+e.target.value); $('perang').onchange = () => pv.redraw(); $('perrays').onchange = () => pv.redraw(); $('pertun').onchange = () => pv.redraw();
    $('raybox').style.display = ''; $('raygo').onclick = () => setRays($('rayin').value);
    $('rayin').addEventListener('keydown', e => { if (e.key === 'Enter') setRays(e.target.value); });
    $('rayclr').onclick = () => { RAYS = []; $('rayout').textContent = ''; parent.postMessage({ e96: 'rays', angles: [] }, '*'); pv.redraw(); }; }
  $('minb').onclick();                      // panel starts collapsed (Caroline 2026-11-10); + expands
  if (E.n < 2) $('matebox').style.display = 'none';
  $('mgo').onclick = () => parent.postMessage({ e96: 'mateReq', label: $('mf').value, angle: $('mg').value, K: $('mvn').checked ? +$('mvk').value || 0 : 0, mis: $('mvm').checked }, '*');
  $('mvk').addEventListener('keydown', e => { e.stopPropagation(); if (e.key === 'Enter') { $('mvn').checked = true; $('mgo').click(); } }); $('mvk').addEventListener('input', () => { $('mvn').checked = true; });
  $('mclr').onclick = () => { MARKS.length = 0; pv.redraw(); };
  if (E.n >= 3) { $('webbox').style.display = ''; const wgo = k => { const l = parseInt($('wf').value, 10); if (isFinite(l)) parent.postMessage({ e96: 'webReq', kind: k, label: l }, '*'); };
    $('wv').onclick = () => wgo('V'); $('wt').onclick = () => wgo('T'); $('wc').onclick = () => { WEBL.length = 0; pv.redraw(); }; $('wf').addEventListener('keydown', e => e.stopPropagation()); }
  $('mdlgb').onclick = () => parent.postMessage({ e96: 'key', key: 'm' }, '*'); $('cdlgb').onclick = () => parent.postMessage({ e96: 'key', key: 'c' }, '*');
  $('mg').addEventListener('keydown', e => { if (e.key === 'Enter') $('mgo').click(); });
  window.addEventListener('keydown', e => { if ((e.key === 'm' || e.key === 'v') && e.target.tagName !== 'INPUT' && e.target.tagName !== 'SELECT') parent.postMessage({ e96: 'key', key: e.key }, '*'); });
  window.E96READY = true;
})().catch(e => { document.body.innerHTML = '<pre style="color:#fff;padding:20px">' + e.stack + '</pre>'; console.error(e); });
