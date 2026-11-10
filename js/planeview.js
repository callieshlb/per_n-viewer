// E96: a pannable / zoomable WebGL2 plane (parameter plane, torus cover, or dynamical plane) with a 2-D overlay canvas.
// View state in double: center ctr = [re, im], upp = units per CSS pixel. Drag pans, wheel / pinch zooms, a click
// (no drag) calls onPick(point). Low resolution while moving, full resolution when idle.
'use strict';
function PlaneView(host, opts) {
  var self = this;
  host.innerHTML = '';
  var gc = document.createElement('canvas'), oc = document.createElement('canvas');
  [gc, oc].forEach(function (c) { c.style.cssText = 'position:absolute;inset:0;width:100%;height:100%'; host.appendChild(c); });
  oc.style.pointerEvents = 'none'; gc.style.cursor = 'crosshair';
  var gl = gc.getContext('webgl2', { antialias: false, preserveDrawingBuffer: true });
  this.gl = gl; this.overlay = oc.getContext('2d'); this.host = host; this.canvas = gc;
  if (!gl) { host.innerHTML = '<p style="padding:20px">This page needs WebGL2.</p>'; return; }
  function sh(type, src) { var s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s)); return s; }
  var prog = gl.createProgram(); gl.attachShader(prog, sh(gl.VERTEX_SHADER, GLSL.vs)); gl.attachShader(prog, sh(gl.FRAGMENT_SHADER, GLSL.fs()));
  gl.linkProgram(prog); if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(prog));
  gl.useProgram(prog);
  var buf = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, buf); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
  var loc = gl.getAttribLocation(prog, 'p'); gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, 2, gl.FLOAT, false, 0, 0);
  var U = function (n) { return gl.getUniformLocation(prog, n); };
  this.U = U;
  this.ctr = (opts.ctr || [0, 0]).slice(); this.upp = opts.upp || 0.01;
  this.maxIt = opts.maxIt || 400; this.bands = true; this.plain = false; this.delta = opts.delta || 1e-3; this.tolP = 1e-5;
  this.setUniforms = opts.setUniforms || function () {};
  this.drawOverlay = opts.drawOverlay || function () {};
  this.onPick = opts.onPick || function () {};
  this.onView = opts.onView || function () {};
  var dpr = 1, W = 0, H = 0, dirty = true, interacting = false, idleT = null;
  function resize() { dpr = window.devicePixelRatio || 1; W = host.clientWidth; H = host.clientHeight;
    oc.width = Math.round(W*dpr); oc.height = Math.round(H*dpr); dirty = true; }
  this.resize = resize;
  if (window.ResizeObserver) new ResizeObserver(resize).observe(host); else window.addEventListener('resize', resize);
  resize();
  function poke() { interacting = true; clearTimeout(idleT); idleT = setTimeout(function () { interacting = false; dirty = true; }, 220); }
  this.redraw = function () { dirty = true; };
  this.toScreen = function (p) { return [W/2 + (p[0] - self.ctr[0])/self.upp, H/2 - (p[1] - self.ctr[1])/self.upp]; };
  this.toPlane = function (sx, sy) { return [self.ctr[0] + (sx - W/2)*self.upp, self.ctr[1] - (sy - H/2)*self.upp]; };
  this.size = function () { return [W, H, dpr]; };
  function render() {
    var q = interacting ? 0.4 : 1, GW = Math.max(1, Math.round(W*dpr*q)), GH = Math.max(1, Math.round(H*dpr*q));
    if (gc.width !== GW || gc.height !== GH) { gc.width = GW; gc.height = GH; }
    gl.viewport(0, 0, GW, GH);
    gl.uniform2f(U('res'), GW, GH); gl.uniform2f(U('ctr'), self.ctr[0], self.ctr[1]); gl.uniform1f(U('upp'), self.upp*W/Math.max(1, GW)*1);
    gl.uniform1i(U('maxIt'), interacting ? Math.min(self.maxIt, 150) : self.maxIt); gl.uniform1f(U('delta'), self.delta); gl.uniform1f(U('tolP'), self.tolP);
    gl.uniform1i(U('bands'), self.bands ? 1 : 0); gl.uniform1i(U('plain'), self.plain ? 1 : 0);
    self.setUniforms(gl, U);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
    var g = self.overlay; g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, W, H);
    self.drawOverlay(g, W, H);
    self.onView();
  }
  function loop() { if (dirty && W > 0 && H > 0) { dirty = false; render(); } requestAnimationFrame(loop); }
  requestAnimationFrame(loop);
  var down = null, last = null, moved = 0, pts = new Map(), pinch0 = null;
  gc.addEventListener('pointerdown', function (e) { pts.set(e.pointerId, [e.clientX, e.clientY]); gc.setPointerCapture(e.pointerId);
    if (pts.size === 1) { down = [e.clientX, e.clientY]; last = down; moved = 0; } else { pinch0 = null; moved = 99; } });
  gc.addEventListener('pointermove', function (e) {
    if (!pts.has(e.pointerId)) return; pts.set(e.pointerId, [e.clientX, e.clientY]);
    if (pts.size === 2) { var a = Array.from(pts.values()), d = Math.hypot(a[0][0] - a[1][0], a[0][1] - a[1][1]);
      if (pinch0) { self.upp *= pinch0/d; dirty = true; poke(); } pinch0 = d; return; }
    if (!down) return;
    moved = Math.max(moved, Math.hypot(e.clientX - down[0], e.clientY - down[1]));
    self.ctr[0] -= (e.clientX - last[0])*self.upp; self.ctr[1] += (e.clientY - last[1])*self.upp; last = [e.clientX, e.clientY];
    dirty = true; poke(); });
  gc.addEventListener('pointerup', function (e) { pts.delete(e.pointerId); if (pts.size) return;
    var wasDown = down; down = null; if (!wasDown || moved > 4) return;
    var r = gc.getBoundingClientRect(); self.onPick(self.toPlane(e.clientX - r.left, e.clientY - r.top), [e.clientX - r.left, e.clientY - r.top]); });
  // hotkeys (Caroline 2026-11-10): z = zoom in, x = zoom out, at the cursor; self.mouse = plane point under the cursor
  self.mouse = null; var lastXY = null;
  gc.addEventListener('pointermove', function (e) { var r = gc.getBoundingClientRect(); lastXY = [e.clientX - r.left, e.clientY - r.top]; self.mouse = self.toPlane(lastXY[0], lastXY[1]); });
  gc.addEventListener('pointerleave', function () { lastXY = null; self.mouse = null; });
  window.addEventListener('keydown', function (e) {
    if (!lastXY || e.target.tagName === 'INPUT' || e.target.tagName === 'SELECT' || (e.key !== 'z' && e.key !== 'x')) return;
    var before = self.toPlane(lastXY[0], lastXY[1]); self.upp *= e.key === 'z' ? 1/1.6 : 1.6; self.upp = Math.min(10, Math.max(1e-7, self.upp));
    var after = self.toPlane(lastXY[0], lastXY[1]); self.ctr[0] += before[0] - after[0]; self.ctr[1] += before[1] - after[1]; dirty = true; poke(); });
  gc.addEventListener('wheel', function (e) { e.preventDefault();
    var r = gc.getBoundingClientRect(), sx = e.clientX - r.left, sy = e.clientY - r.top, before = self.toPlane(sx, sy);
    self.upp = Math.min(10, Math.max(1e-7, self.upp*Math.exp(e.deltaY*0.0015)));
    var after = self.toPlane(sx, sy); self.ctr[0] += before[0] - after[0]; self.ctr[1] += before[1] - after[1];
    dirty = true; poke(); }, { passive: false });
}
