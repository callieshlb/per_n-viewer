// E96 (2026-11-17, request relayed by E97): Per_1 .. Per_5 drawn on a 3-D surface, colored live by the same fragment shader as the flat panes.
// Genus 0 (Per_1 .. Per_4): the pane's coordinate p on the Riemann sphere, p = c0 + s0 * conj(X + iY)/(1 - Z) (stereographic from the north pole,
// conjugated so that the picture is seen from outside with its own orientation); the home view of the flat pane fills the southern hemisphere.
// Genus 1 (Per_5): a torus of revolution conformal to the lattice (W1, i W2) of the flat torus: big-circle angle a <-> Re u, conformal small-circle
// coordinate phi(v) = int r dv/(R + r cos v) <-> Im u, R = 3, r = R q/sqrt(1 + q^2), q = W2/W1 (so phi(2 pi) = 2 pi q).
'use strict';
(async function () {
  const id = new URLSearchParams(location.search).get('m') || 'per3', $ = s => document.getElementById(s);
  const man = await fetch('models/manifest.json').then(r => r.json()), E = man.models.find(m => m.id === id), T5 = PM.T5, torus = E.kind === 'torus';
  $('title').innerHTML = (E.short || id) + '(0) in 3D';
  const c0 = E.view.ctr, s0 = 350*E.view.upp;
  $('about').innerHTML = torus ? `The torus Per<sub>5</sub>(0) as a torus of revolution conformal to its lattice (${T5.W1.toFixed(4)}, ${T5.W2.toFixed(4)}i); the long period goes around the big circle. Colors as in the flat pane.`
    : `The sphere: stereographic projection of the pane's coordinate; the flat pane's home view (center ${c0[0]}${c0[1] ? ' + ' + c0[1] + 'i' : ''}, radius ${s0.toFixed(2)}) is the lower hemisphere. Colors as in the flat pane.`;
  const gc = $('gc'), gl = gc.getContext('webgl2', { antialias: true, preserveDrawingBuffer: true }); if (!gl) { $('about').textContent = 'This browser has no WebGL2.'; return; }
  // fragment shader: the pane's shader, with the point p taken from the surface instead of the screen
  let fs = GLSL.fs(); const P0 = 'vec2 p = ctr + (gl_FragCoord.xy - .5*res)*upp;', O0 = 'o = vec4(fateColour(z0, z0i), 1.);';
  if (fs.indexOf(P0) < 0 || fs.indexOf(O0) < 0) { $('about').textContent = 'shader template changed: 3D view unavailable'; return; }
  fs = fs.replace('out vec4 o;', 'out vec4 o; in vec3 vS; in vec2 vU; in float vSh; uniform int geom; uniform vec2 sC; uniform float sR;')
         .replace(P0, 'vec2 p; if (geom == 1) p = vec2(vU.x*W1, vU.y*W2); else { vec3 s = normalize(vS); float d = max(1. - s.z, 1e-7); p = sC + sR*vec2(s.x, -s.y)/d; }')
         .replace(O0, 'o = vec4(fateColour(z0, z0i)*vSh, 1.);');
  const vs = `#version 300 es
in vec3 P; in vec2 UV; uniform mat4 MV; uniform mat4 PR; out vec3 vS; out vec2 vU; out float vSh;
void main(){ vS = P; vU = UV; vec4 q = MV*vec4(P, 1.); gl_Position = PR*q;
  vec3 n = normalize(mat3(MV)*NRM); vSh = .55 + .45*max(0., dot(n, normalize(-q.xyz))); }`;
  const sh = (t, src) => { const s = gl.createShader(t); gl.shaderSource(s, src); gl.compileShader(s); if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s)); return s; };
  const prog = gl.createProgram(); gl.attachShader(prog, sh(gl.VERTEX_SHADER, vs.replace('in vec2 UV;', 'in vec2 UV; in vec3 NRM;'))); gl.attachShader(prog, sh(gl.FRAGMENT_SHADER, fs)); gl.linkProgram(prog);
  if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(prog)); gl.useProgram(prog);
  const U = n => gl.getUniformLocation(prog, n);
  // mesh
  const pos = [], uv = [], nrm = [], idx = []; let NA, NB;
  if (torus) { NA = 240; NB = 120; const R = 3, q = T5.W2/T5.W1, r = R*q/Math.sqrt(1 + q*q), phi = [0];
    for (let j = 1; j <= NB*16; j++) { const v0 = 2*Math.PI*(j - 1)/(NB*16), v1 = 2*Math.PI*j/(NB*16), f = v => r/(R + r*Math.cos(v)); phi.push(phi[j - 1] + (v1 - v0)*(f(v0) + 4*f((v0 + v1)/2) + f(v1))/6); }
    const PH = phi[NB*16];
    for (let i = 0; i <= NA; i++) for (let j = 0; j <= NB; j++) { const a = 2*Math.PI*i/NA, v = 2*Math.PI*j/NB, w = R + r*Math.cos(v);
      pos.push(w*Math.cos(a)/R, w*Math.sin(a)/R, r*Math.sin(v)/R); nrm.push(Math.cos(v)*Math.cos(a), Math.cos(v)*Math.sin(a), Math.sin(v)); uv.push(i/NA, phi[j*16]/PH); }
    window.E96TORUS = { R, r, q, PH, phiRatio: PH/(2*Math.PI*q) };
  } else { NA = 256; NB = 160;
    for (let i = 0; i <= NA; i++) for (let j = 0; j <= NB; j++) { const a = 2*Math.PI*i/NA, t = Math.PI*j/NB, x = Math.sin(t)*Math.cos(a), y = Math.sin(t)*Math.sin(a), z = Math.cos(t);
      pos.push(x, y, z); nrm.push(x, y, z); uv.push(0, 0); } }
  for (let i = 0; i < NA; i++) for (let j = 0; j < NB; j++) { const k = i*(NB + 1) + j, k2 = (i + 1)*(NB + 1) + j; idx.push(k, k2, k + 1, k + 1, k2, k2 + 1); }
  const buf = (data, loc, n) => { const b = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, b); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(data), gl.STATIC_DRAW); const l = gl.getAttribLocation(prog, loc); gl.enableVertexAttribArray(l); gl.vertexAttribPointer(l, n, gl.FLOAT, false, 0, 0); };
  const vao = gl.createVertexArray(); gl.bindVertexArray(vao); buf(pos, 'P', 3); buf(uv, 'UV', 2); buf(nrm, 'NRM', 3);
  const ib = gl.createBuffer(); gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, ib); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(idx), gl.STATIC_DRAW);
  // uniforms of the pane
  gl.uniform1i(U('fam'), E.famId); GLSL.setPalette(gl, U, E.palette, E.flat); gl.uniform1f(U('W1'), T5.W1); gl.uniform1f(U('W2'), T5.W2); gl.uniform1f(U('qn'), Math.exp(-Math.PI*T5.W2/T5.W1)); gl.uniform1f(U('e1'), 0.75);
  gl.uniform1i(U('geom'), torus ? 1 : 0); gl.uniform2f(U('sC'), c0[0], c0[1]); gl.uniform1f(U('sR'), s0);
  gl.uniform1f(U('delta'), 1e-3); gl.uniform1f(U('tolP'), 1e-5); gl.uniform1i(U('plain'), 0); gl.uniform1i(U('perCol'), 0); gl.uniform2f(U('res'), 1, 1); gl.uniform2f(U('ctr'), 0, 0); gl.uniform1f(U('upp'), 1);
  // view: rotation (quaternion-free: yaw, pitch), distance
  const home = torus ? { yaw: 0, pitch: -0.95, dist: 5.6 } : { yaw: 0, pitch: 2.6, dist: 3.4 }; let V = { ...home }, maxIt = 600, bands = true, dirty = true, interacting = false, itimer = null;
  const mul = (A, B) => { const C = new Array(16).fill(0); for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) for (let k = 0; k < 4; k++) C[i + 4*j] += A[i + 4*k]*B[k + 4*j]; return C; };
  const rx = t => [1, 0, 0, 0, 0, Math.cos(t), Math.sin(t), 0, 0, -Math.sin(t), Math.cos(t), 0, 0, 0, 0, 1], rz = t => [Math.cos(t), Math.sin(t), 0, 0, -Math.sin(t), Math.cos(t), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
  function draw() { const dpr = Math.min(2, window.devicePixelRatio || 1), q = interacting ? 0.5 : 1, W = gc.clientWidth, H = gc.clientHeight, GW = Math.max(1, Math.round(W*dpr*q)), GH = Math.max(1, Math.round(H*dpr*q));
    if (gc.width !== GW || gc.height !== GH) { gc.width = GW; gc.height = GH; } gl.viewport(0, 0, GW, GH);
    const asp = GW/GH, f = 1/Math.tan(0.35), n0 = 0.05, f0 = 50, PR = [f/asp, 0, 0, 0, 0, f, 0, 0, 0, 0, (f0 + n0)/(n0 - f0), -1, 0, 0, 2*f0*n0/(n0 - f0), 0];
    const T = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, -V.dist, 1], MV = mul(T, mul(rx(V.pitch), rz(V.yaw)));
    gl.uniformMatrix4fv(U('MV'), false, MV); gl.uniformMatrix4fv(U('PR'), false, PR); gl.uniform1i(U('maxIt'), interacting ? Math.min(maxIt, 150) : maxIt); gl.uniform1i(U('bands'), bands ? 1 : 0);
    gl.enable(gl.DEPTH_TEST); gl.clearColor(0.1, 0.1, 0.1, 1); gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT); gl.drawElements(gl.TRIANGLES, idx.length, gl.UNSIGNED_INT, 0); }
  const poke = () => { dirty = true; }, act = () => { interacting = true; clearTimeout(itimer); itimer = setTimeout(() => { interacting = false; dirty = true; }, 250); };
  (function loop() { if ($('spin').checked) { V.yaw += 0.004; dirty = true; act(); } if (dirty) { dirty = false; draw(); } requestAnimationFrame(loop); })();
  let drag = null; gc.addEventListener('pointerdown', e => { drag = [e.clientX, e.clientY]; gc.setPointerCapture(e.pointerId); gc.style.cursor = 'grabbing'; });
  gc.addEventListener('pointermove', e => { if (!drag) return; V.yaw += (e.clientX - drag[0])*0.008; V.pitch = Math.max(-Math.PI, Math.min(Math.PI, V.pitch + (e.clientY - drag[1])*0.008)); drag = [e.clientX, e.clientY]; act(); poke(); });
  gc.addEventListener('pointerup', () => { drag = null; gc.style.cursor = 'grab'; });
  gc.addEventListener('wheel', e => { e.preventDefault(); V.dist = Math.min(20, Math.max(1.15, V.dist*Math.exp(e.deltaY*0.0012))); act(); poke(); }, { passive: false });
  window.addEventListener('resize', poke);
  $('bands').onchange = e => { bands = e.target.checked; poke(); }; $('it').oninput = e => { maxIt = +e.target.value; $('itv').textContent = maxIt; poke(); };
  $('reset').onclick = () => { V = { ...home }; poke(); }; $('minb').onclick = () => { $('panel').classList.toggle('min'); $('minb').textContent = $('panel').classList.contains('min') ? '+' : '–'; };
  window.addEventListener('keydown', e => { if (e.key === 'w' && e.target.tagName !== 'INPUT') { $('bands').checked = !$('bands').checked; $('bands').onchange({ target: $('bands') }); } });
  window.E96GLOBE = { V, set: v => { Object.assign(V, v); poke(); }, err: gl.getError(), verts: pos.length/3 };
})().catch(e => { document.getElementById('about').textContent = '3D view failed: ' + e.message; window.E96GLOBE_ERR = e.message; });
