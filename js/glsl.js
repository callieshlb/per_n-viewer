// E96: fragment shaders (WebGL2, single precision) for the plane / torus parameter panes and the dynamical plane.
// Same colouring as the E92-E94 live viewers: capture into the marked cycle coloured by phase (j - t) mod n and shaded by
// time t; attracting fixed point dark grey; other attracting cycle light grey; undecided (bifurcation locus / Julia set)
// near-black; degenerate / overflow dark red. Port of PM.fz / PM.fate (pern_math.js, double) with infinity flags.
'use strict';
var GLSL = (function () {
  var HEAD = `#version 300 es
precision highp float;
uniform vec2 res; uniform vec2 ctr; uniform float upp;           // view: centre, units per pixel
uniform int maxIt; uniform float delta; uniform float tolP; uniform int fam; uniform int nM;
uniform vec2 uA, uB, uC, uD, uE, uF; uniform vec2 uMk[12]; uniform int uMkI[12];   // dynamical-plane map (mode 0)
uniform float W1, W2, qn, e1; uniform int invz; uniform vec2 mA, mB, mC, mD;   // invz: pixel w -> z = (mA w + mB)/(mC w + mD)                                    // Per_5 torus
out vec4 o;
vec2 cm(vec2 a, vec2 b){ return vec2(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x); }
vec2 cd(vec2 a, vec2 b){ return vec2(a.x*b.x + a.y*b.y, a.y*b.x - a.x*b.y)/dot(b, b); }
vec2 cinv(vec2 z){ return vec2(z.x, -z.y)/dot(z, z); }
float chF(vec2 a, vec2 b){ float A = dot(a, a), B = dot(b, b);
  if (A > 1. && B > 1.) { vec2 ia = cinv(a), ib = cinv(b); return length(ia - ib)/sqrt((1. + dot(ia, ia))*(1. + dot(ib, ib))); }
  return length(a - b)/sqrt((1. + A)*(1. + B)); }
float chI(vec2 a){ return 1./sqrt(1. + dot(a, a)); }
float chN(vec2 a, bool ai, vec2 b, bool bi){ if (ai && bi) return 0.; if (ai) return chI(b); if (bi) return chI(a); return chF(a, b); }
vec2 A_, B_, C_, D_, E_, F_; vec2 mk[12]; bool mi[12]; int nm;
// z' = f(z); infinity as a flag
void fz(inout vec2 z, inout bool zi){
  vec2 num, den;
  if (zi || dot(z, z) > 1.) { vec2 w = zi ? vec2(0.) : cinv(z), w2 = cm(w, w);
    num = A_ + cm(B_, w) + cm(C_, w2); den = D_ + cm(E_, w) + cm(F_, w2); }
  else { vec2 z2 = cm(z, z); num = cm(A_, z2) + cm(B_, z) + C_; den = cm(D_, z2) + cm(E_, z) + F_; }
  float ln = length(num), ld = length(den);
  if (ld <= 1e-30*ln || ld == 0.) { zi = true; return; }
  if (ld < ln) { vec2 r = cd(den, num); if (dot(r, r) < 1e-28) { zi = true; return; } z = cinv(r); zi = false; return; }
  z = cd(num, den); zi = false;
}
void corMap(vec2 x, vec2 y, int n, out vec2 om){
  vec2 one = vec2(1., 0.), w = cd(x + y - one, x - one), s = x + w, p = cm(x, w);
  A_ = one; B_ = -s; C_ = p; D_ = one; E_ = vec2(0.); F_ = vec2(0.); om = cd(2.*p, s);
  nm = n; mk[0] = vec2(0.); mi[0] = false; mk[1] = vec2(0.); mi[1] = true; mk[2] = one; mi[2] = false;
  if (n >= 4) { mk[3] = y; mi[3] = false; }
  for (int k = 4; k < 11; k++) { if (k >= n - 1) break; mk[k] = mk[k - 1]; mi[k] = mi[k - 1]; fz(mk[k], mi[k]); }
  if (n >= 5) { mk[n - 1] = x; mi[n - 1] = false; }
}
vec2 csin(vec2 z){ return vec2(sin(z.x)*cosh(z.y), cos(z.x)*sinh(z.y)); }
vec2 ccos(vec2 z){ return vec2(cos(z.x)*cosh(z.y), -sin(z.x)*sinh(z.y)); }
vec2 wpX(vec2 u, out bool pole){                 // Per_5 torus: wp via theta functions (pern_math.js wpX)
  float a = u.x - W1*floor(u.x/W1 + .5), b = u.y - W2*floor(u.y/W2 + .5);
  pole = (a*a + b*b) < 1e-14; vec2 v = 3.14159265358979*vec2(a, b)/W1, t1 = vec2(0.), t2 = vec2(0.); float t3 = 1., t4 = 1.;
  for (int n = 0; n < 5; n++) { float fn = float(n), qq = pow(qn, (fn + .5)*(fn + .5)); vec2 arg = (2.*fn + 1.)*v;
    t1 += csin(arg)*2.*qq*((n & 1) == 1 ? -1. : 1.); t2 += ccos(arg)*2.*qq;
    if (n > 0) { float q2 = pow(qn, fn*fn); t3 += 2.*q2; t4 += 2.*q2*((n & 1) == 1 ? -1. : 1.); } }
  vec2 R = cd(t2, t1)*t3*t4; float k = 3.14159265358979/W1;
  return vec2(e1, 0.) + cm(R, R)*k*k;
}
vec2 csq(vec2 z){ float r = length(z), t = atan(z.y, z.x)*.5; return sqrt(r)*vec2(cos(t), sin(t)); }
bool torusXY(vec2 u, out vec2 x, out vec2 y){
  bool pole, p1, p2; vec2 X = wpX(u, pole); if (pole) { x = vec2(1., 0.); y = vec2(1., 0.); return true; }
  vec2 r = csq(4.*cm(X, cm(X, X)) - 2.75*X + vec2(.375, 0.));
  float h = 2e-3; vec2 fd = (wpX(u + vec2(h, 0.), p1) - wpX(u - vec2(h, 0.), p2))/(2.*h);
  if (length(fd - r) > length(fd + r)) r = -r;
  vec2 one = vec2(1., 0.), s = X + vec2(1.25, 0.); if (dot(s, s) < 1e-20) return false;
  y = one - cinv(s); vec2 y1 = y - one, vv = cm(r, cm(y1, y1)), y2 = cm(y, y);
  vec2 A = y2 + y - one, B = 2.*y2 + y - one, C = 2.*y2 - y;
  if (length(B + vv) >= length(B - vv)) { if (dot(A, A) < 1e-24) return false; x = cd(B + vv, 2.*A); } else x = cd(2.*C, B - vv);
  return true;
}
uniform vec3 CAP[8]; uniform int flatPal;                       // phase palette (0..255), flat = no time shading
uniform int bands; uniform int plain;                                               // escape-time bands (Caroline 2026-11-10)
vec3 shadeT(vec3 c, float t){
  float s = flatPal == 1 ? 1. : max(.3, min(1., 1. - .6*log(1. + t)/log(81.)));
  if (bands == 1) s *= .72 + .28*cos(.55*t);
  return (c*s + 20.*(1. - s))/255.; }
vec3 greyT(float base, float t){ float v = base; if (bands == 1) v = base*(.78 + .22*cos(.5*t)) + 12.*(1. - exp(-.02*t)); return vec3(v)/255.; }
vec3 fateColour(vec2 z, bool zi){
  float th[12];
  for (int j = 0; j < 12; j++) { if (j >= nm) break; float m = 1.; for (int i = 0; i < 12; i++) { if (i >= nm) break; if (i != j) m = min(m, chN(mk[j], mi[j], mk[i], mi[i])); } th[j] = delta*m; }
  vec2 ring[32]; bool ringI[32]; float tconv = -1.;
  for (int it = 0; it < 20000; it++) {
    if (it >= maxIt) break;
    fz(z, zi);
    for (int j = 0; j < 12; j++) { if (j >= nm) break;
      float dj = chN(z, zi, mk[j], mi[j]);
      if (dj < th[j]) { int ph = ((j - 1 - (it + 1)) % nm + nm) % nm;   // critical-value convention     // smooth capture time (superattracting cycle: distances square)
        float ts = float(it + 1) - log2(max(1., log(max(dj, 1e-30))/log(th[j]))); return shadeT(CAP[ph % 8], max(ts, 0.)); } }
    if (tconv < 0. && it > 0 && chN(z, zi, ring[(it - 1) & 31], ringI[(it - 1) & 31]) < 1e-3) tconv = float(it);
    if (tconv >= 0. && it > 0 && chN(z, zi, ring[(it - 1) & 31], ringI[(it - 1) & 31]) > 1e-2) tconv = -1.;
    ring[it & 31] = z; ringI[it & 31] = zi;
  }
  for (int per = 1; per <= 24; per++) {
    int a = (maxIt - 1 - per + 3200) & 31, b = (maxIt - 1) & 31;
    if (chN(ring[b], ringI[b], ring[a], ringI[a]) < tolP) return per == 1 ? greyT(85., max(tconv, 0.)) : greyT(163., float(maxIt));
  }
  return vec3(18.)/255.;
}
`;
  var MAIN = `
void main(){
  vec2 p = ctr + (gl_FragCoord.xy - .5*res)*upp;
  if (plain == 1) { o = vec4(.94, .93, .90, 1.); return; }        // plain surface: the curve alone, no dynamics
  vec2 om; bool ok = true; vec2 z0; bool z0i = false;
  if (fam == 0) {                          // dynamical plane: the pixel is the starting point
    A_ = uA; B_ = uB; C_ = uC; D_ = uD; E_ = uE; F_ = uF; nm = nM;
    for (int j = 0; j < 12; j++) { mk[j] = uMk[j]; mi[j] = uMkI[j] == 1; }
    z0 = p; if (invz == 1) { vec2 nu = cm(mA, p) + mB, de = cm(mC, p) + mD; if (dot(de, de) < 1e-30*max(dot(nu, nu), 1e-30)) z0i = true; else z0 = cd(nu, de); }
  } else {
    if (fam == 1) { A_ = vec2(1., 0.); B_ = vec2(0.); C_ = p; D_ = vec2(0.); E_ = vec2(0.); F_ = vec2(1., 0.); nm = 1; mk[0] = vec2(0.); mi[0] = true; om = vec2(0.); }
    else if (fam == 2) { A_ = vec2(1., 0.); B_ = vec2(0.); C_ = -p; D_ = vec2(1., 0.); E_ = vec2(0.); F_ = vec2(-1., 0.); nm = 2;   // (z^2 - c)/(z^2 - 1)
      mk[0] = vec2(0.); mi[0] = true; mk[1] = vec2(1., 0.); mi[1] = false; om = vec2(0.); }
    else if (fam == 3) corMap(cd(p, vec2(2., 0.) - p), vec2(0.), 3, om);                 // x = c/(2 - c); om = c
    else if (fam == 4) { vec2 y4 = cd(vec2(.1, 0.), p) + vec2(.6180339887498949, 0.); corMap(y4, y4, 4, om); }
    else if (fam == 5) { vec2 x, y; ok = torusXY(p, x, y); if (ok) corMap(x, y, 5, om); }
    z0 = om;
    if (!ok || dot(om, om) > 1e30 || dot(B_, B_) > 1e30 || dot(C_, C_) > 1e30 || !(dot(om, om) >= 0.)) { o = vec4(85., 20., 20., 255.)/255.; return; }
  }
  o = vec4(fateColour(z0, z0i), 1.);
}`;
  return { fs: function () { return HEAD + MAIN; },
           vs: '#version 300 es\nin vec2 p; void main(){ gl_Position = vec4(p, 0., 1.); }',
           DEFPAL: [[42, 120, 214], [235, 104, 52], [27, 175, 122], [237, 161, 0], [232, 123, 164], [127, 63, 191], [160, 82, 45], [0, 150, 160]],
           setPalette: function (gl, U, pal, flat) { var p = (pal || this.DEFPAL).concat(this.DEFPAL).slice(0, 8), a = new Float32Array(24);
             p.forEach(function (c, i) { a.set(c, 3*i); }); gl.uniform3fv(U('CAP'), a); gl.uniform1i(U('flatPal'), flat ? 1 : 0); },
           LEGEND: function (n, pal) {
             var C = (pal || this.DEFPAL).map(function (c) { return 'rgb(' + c.join(',') + ')'; }), L = [];
             for (var j = 0; j < Math.min(n, 8); j++) L.push([C[j], 'captured, phase ' + j]);
             L.push(['#555555', 'attracting fixed point'], ['#a3a3a3', 'other attracting cycle'], ['#121212', 'undecided'], ['#551414', 'degenerate / overflow']);
             return L; } };
})();
if (typeof module !== 'undefined') module.exports = GLSL;
