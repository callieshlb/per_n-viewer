<h1>The Per<sub>n</sub>(0) atlas</h1>
<p>Pick a curve in the bar above: Per<sub>1</sub> … Per<sub>9</sub>, the PrePer<sub>k,n</sub> grid, or the ℝ² slice. Click a point of a curve to see its rational map in the dynamical plane on the right. This page is the default view; the title <b>Per<sub>n</sub>(0) atlas</b> brings it back.</p>
<h2>Where these curves live</h2>
<p><b>rat<sub>2</sub></b> is the moduli space of quadratic rational maps up to Möbius conjugacy; it is isomorphic to ℂ² (<a href="https://arxiv.org/abs/math/9209221" target="_blank" rel="noopener">Milnor, <i>Geometry and dynamics of quadratic rational maps</i></a>). A quadratic rational map has two critical points. Marking them, i.e. remembering which critical point is which, gives the branched double cover <b>rat<sub>2</sub><sup>cm</sup></b>, the critically marked moduli space.</p>
<p><i>Why "double cover"?</i> A puzzle: I have two US coins adding up to 30 cents, and one of them isn't a nickel. How? The other one is a nickel. "One of them" does not say <i>which</i> one, and the two critical points of a map are like the two coins: rat<sub>2</sub> only sees the unordered pair; rat<sub>2</sub><sup>cm</sup> also records which is which.</p>
<p>Inside rat<sub>2</sub><sup>cm</sup> we look at <b>critical orbit relation (COR) curves</b>: one-dimensional algebraic curves on which the marked critical point ω has a prescribed partial dynamics, while the other, <i>free</i>, critical point ω′ moves.</p>
<ul>
<li><b>Per<sub>n</sub>(0)</b>: ω is periodic of exact period n.</li>
<li><b>PrePer<sub>k,n</sub></b>: ω is strictly preperiodic; its orbit first lands on a cycle of exact period n after k steps.</li>
</ul>
<h2>What is known</h2>
<ul>
<li><b>Per<sub>1</sub>(0)</b> is the classical moduli space of quadratic polynomials {z² + c}: a rational map with a fully ramified fixed point (a fixed critical point) is conjugate to a polynomial. It contains the Mandelbrot set.</li>
<li><b>Per<sub>2</sub>(0)</b> has been studied by Wittner (thesis, Cornell 1988), Luo (thesis, Cornell 1995), <a href="https://arxiv.org/abs/1112.4780" target="_blank" rel="noopener">Dudko, <i>Matings with laminations</i></a>, and <a href="https://arxiv.org/abs/math/0610343" target="_blank" rel="noopener">Aspenberg–Yampolsky, <i>Mating non-renormalizable quadratic polynomials</i></a>.</li>
<li><b>Per<sub>3</sub>(0)</b>: Wittner, and Rees (<i>A partial description of the parameter space of rational maps of degree two</i>, Acta Math. 168 (1992); <i>Views of parameter space: Topographer and Resident</i>, Astérisque 288 (2003)).</li>
<li><b>Per<sub>n</sub>(0), n ≥ 4</b>: not much discussed. Epstein showed Per<sub>n</sub>(0) is smooth, so irreducibility, a longstanding open question, is the same as connectivity. Milnor showed every connected component contains a polynomial.</li>
<li><b>PrePer<sub>k,n</sub></b> curves are also mysterious. PrePer<sub>k,1</sub> was studied by <a href="https://arxiv.org/abs/1806.11221" target="_blank" rel="noopener">Buff–Epstein–Koch</a> (published as <i>Prefixed curves in moduli space</i>, Amer. J. Math. 144 (2022)), who proved irreducibility by algebraic means; a dynamical explanation is less clear. On the (k,1) pages, "Rays in the landing region" draws the landing region around its one puncture.</li>
</ul>
<h2>A fun fact</h2>
<p>The <b>symmetry locus</b> (maps commuting with a Möbius involution, e.g. f(z) = k(z + 1/z), which commutes with z ↦ −z) is a 2:1 cover of PrePer<sub>2,1</sub>: dividing out by the symmetry (w = z²) gives g(w) = k²(w + 2 + 1/w), whose critical point −1 maps to 0, then to the fixed point ∞.</p>
<h2>How the pages are made</h2>
<p>Genus 0 curves are drawn on their parameter plane, genus 1 on a flat torus, and genus ≥ 2 as conformal hyperbolic models (the Poincaré disk tiled by a Dirichlet polygon). Each pixel is coloured by the fate of the free critical point ω′, computed afresh. Everything here is numerical; the checks behind each page are recorded with the computations.</p>
