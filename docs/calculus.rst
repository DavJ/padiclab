Analytic series and local calculus
==================================

``AnalyticSeries`` represents a local function, or a family of functions
consistent with finite coefficient data and an explicit bound for all
omitted coefficients. It supports Qp and the two certified quadratic field
engines. It uses only Python's standard library.

Use ``from_polynomial`` for a finite polynomial and ``from_function`` for
exp, sin, cos, log1p or log. These factories establish their own tail bounds.
A custom prefix requires an explicit caller-supplied coefficient bound.

Exact polynomials and primitives
--------------------------------

Coefficients are in ascending powers of h=x-center. ``from_polynomial``
accepts a polynomial in the original variable x and shifts it exactly.
Exact rational coefficients and quadratic coordinate pairs remain exact;
finite coefficient balls retain their uncertainty. Centers require exact
source data, rather than an arbitrarily chosen representative of a ball.

.. doctest::

   >>> K = Qp(3, 3)
   >>> f = AnalyticSeries.from_polynomial(K, Polynomial([1, 2, 3]), center=2, precision=18)
   >>> f.coefficients
   (Fraction(17, 1), Fraction(14, 1), Fraction(3, 1))
   >>> f.derivative()(5).contains(32)
   True
   >>> F = f.antiderivative(constant=7)
   >>> F(2).contains(7)
   True
   >>> F.derivative().coefficients == f.coefficients
   True
   >>> f(5, precision=30).precision
   30

``antiderivative`` constructs the local primitive with F(center)=constant;
``integral`` is its alias. These methods do not evaluate an integral against
a measure. The separate :doc:`measures`, :doc:`transforms` and
:doc:`harmonic` APIs supply those operations.

A finite polynomial has ``tail_bound=None`` and may use ``radius=None``
for the whole field. Its evaluation can exceed the parent's default
precision when its coefficients and argument come from exact sources.

Analytic function factories
---------------------------

The factory's exact centers at zero (or at one for log) produce exact
rational Taylor coefficients. Other exponential/trigonometric centers use
local Taylor coefficient balls with guard digits; their uncertainty is
retained during calculus. Logarithmic coefficients of positive degree are
exact reciprocal powers of the exact scale, divided by the degree. Only
the logarithmic constant value needs a ball. Differentiation can therefore
remove that uncertainty. The logarithm uses the default log(p)=0 branch.

.. doctest::

   >>> L = Q2i(3)
   >>> s = AnalyticSeries.from_function(exp, L, order=48, precision=21)
   >>> s.radius
   3
   >>> x = L(4, 4, absprec=25)
   >>> s(x, precision=21).is_congruent(exp(x), 21)
   True
   >>> s.derivative(2)(x, precision=21).is_congruent(exp(x), 21)
   True
   >>> s.antiderivative()(x, precision=21).is_congruent(exp(x)-1, 21)
   True
   >>> logarithm = AnalyticSeries.from_function(log1p, L, order=96, precision=21)
   >>> logarithm.radius
   1
   >>> pi = L(0, 1, absprec=27)
   >>> logarithm.derivative()(pi, precision=21).is_congruent(1/(1+pi), 21)
   True
   >>> expected = (1+pi)*log1p(pi)-pi
   >>> logarithm.antiderivative()(pi, precision=21).is_congruent(expected, 21)
   True
   >>> local_log = AnalyticSeries.from_function(log, K, center=3, order=20, precision=9)
   >>> local_log.derivative()(3, precision=60).contains(Fraction(1, 3))
   True

``radius`` is a uniformizer exponent: r means v_pi(x-center)>=r.
It is not a real-valued radius. The whole input ball must belong to this
closed disc. By default, exponential factories use r>e/(p-1). Logarithmic
factories use r>v_pi(scale), with scale=center for log and scale=1+center
for log1p. Thus log1p supports the full principal displacement disc, even
in wild dyadic extensions. Read :doc:`extension_analytic` for conventions.

Known prefix and the infinite tail
----------------------------------

For order m and ``tail_bound=(A,B)``, the contract is

.. math::

   f(c+h)=\sum_{n=0}^{m-1}a_nh^n+\sum_{n\geq m}a_nh^n,
   \qquad v_\pi(a_n)\geq A+Bn\quad\text{for every }n\geq m.

A and B are exact rationals. The declared domain must have r+B>0. For
s=v_pi(h)>=r the entire omitted tail belongs to

.. math::

   O\left(\pi^{\lceil A+(s+B)m\rceil}\right).

The strict positive slope proves convergence and bounds every later term.
A bound supplied by a caller is a mathematical premise about the intended
function. It cannot be proved by examining a finite prefix. The returned
enclosure covers all coefficient sequences satisfying that premise.

For example, a prefix of ones and integral omitted coefficients includes
the geometric series 1/(1-x), as well as many other functions:

.. doctest::

   >>> family = AnalyticSeries(K, [1]*24, radius=1, tail_bound=(0, 0), precision=18)
   >>> family(3, precision=18).contains(Fraction(-1, 2))
   True
   >>> short = AnalyticSeries(K, [1], radius=1, tail_bound=(0, 0), precision=18)
   >>> short(3).precision
   1
   >>> short(3, precision=2)
   Traceback (most recent call last):
   ...
   padic_tools.core.PrecisionError: requested precision exceeds coefficient, input or tail information
   >>> short[7].precision
   0
   >>> short.tail_precision(0) is None
   True
   >>> short.tail_precision(K.O(2))
   2

``coefficient(n)`` (or indexing) returns a known coefficient, an exact zero
beyond a finite polynomial, or an unknown ball from the coefficient bound.
Unknown coefficients are never silently filled with exact zeros.
At an exact center the positive-order tail vanishes. A zero-containing
argument ball still represents nonzero displacements and retains its tail.
An empty prefix is supported: its constant coefficient is also unknown.

``precision`` is the default absolute output target, independent of order m.
Evaluation propagates coefficient and input balls through Horner arithmetic
and adds the full-tail enclosure. The default target is capped by available
information. An explicitly requested target above it raises PrecisionError.
This bound can be conservative because arithmetic does not track correlated
errors in coefficients or repeated occurrences of the input.

Derivatives
------------

On a convergent local power series, differentiation is termwise, as in
`Conrad, section 5 <https://kconrad.math.uconn.edu/blurbs/gradnumthy/infseriespadic.pdf>`_.
For derivative order d, each coefficient is multiplied by a falling
factorial, which is integral. Therefore the omitted coefficients satisfy

.. math::

   v_\pi\left(\frac{(n+d)!}{n!}a_{n+d}\right)
      \geq A+Bd+Bn.

The returned series has at most max(m-d,0) known coefficients and tail
bound (A+Bd,B). Known falling factorials propagate their actual valuation
gain; the generic tail bound conservatively ignores that gain. Even when
d>=m, the derivative's unknown tail remains represented.

Primitives and denominator spikes
----------------------------------

Formal coefficients transform to b_0=constant and b_k=a_(k-1)/k. If
v_pi(p)=e, division by k loses e*v_p(k) absolute coefficient digits.
This is visible with finite input coefficients:

.. doctest::

   >>> uncertain = AnalyticSeries.from_polynomial(L, [0, L(1, absprec=9)], precision=20)
   >>> primitive = uncertain.antiderivative()
   >>> primitive[2].precision
   7
   >>> primitive.derivative()[1].precision
   9

For the infinite tail we derive a new linear coefficient bound rather
than differentiate or integrate a single numerical error term. Put
delta=(r+B)/2>0 and define the full lower envelope

.. math::

   T(K,\delta)=\min_{k\geq K}
     \left(k\delta-e\lfloor\log_p k\rfloor\right).

Between powers of p the expression increases. Its values at successive
p-powers eventually increase with difference (p-1)*p**j*delta-e. Thus this
minimum is evaluated with finitely many checks, including every later
denominator spike. Since v_p(k)<=floor(log_p(k)), all k>=m+1 satisfy

.. math::

   \begin{aligned}
      v_\pi(b_k)&\geq A-B+Bk-e\lfloor\log_p k\rfloor\\
                &\geq A'+B'k,\\
      A'&=A-B+T(m+1,\delta),\\
      B'&=B-\delta.
   \end{aligned}

The new gap r+B'=(r+B)/2 is positive, so the same declared disc remains
valid. This is a conservative certificate; repeated integration can weaken
it more than an optimal function-specific bound. A finite polynomial needs
no infinite-tail bound, and its exact coefficient divisions remain exact.

Arithmetic, domains and scope
-----------------------------

Addition, subtraction, products and scalar multiplication preserve
enclosures. Series must have the same exact center and field presentation;
their result uses the intersection of the declared discs. A common slope
and global coefficient lower bounds control the omitted convolution terms.
The ultrametric inequality introduces no real-valued term-count factor.
``restrict(r)`` can shrink a disc but cannot enlarge it.

A structural zero polynomial kills a product exactly. A stored zero ball
contains unknown nonzero coefficients and keeps its uncertainty. Arithmetic
can weaken a certificate, so an algebraically equivalent representation
can return a coarser enclosure. ``from_function`` should be used when its
specific certificate is available.

Formal ``PowerSeries`` retains its separate rational/x-adic semantics;
compare `Sage's power-series methods
<https://doc.sagemath.org/html/en/reference/power_series/sage/rings/power_series_ring_element.html>`_.
Its prefix alone does not imply a convergent AnalyticSeries. Supply a proven
tail bound when transferring a custom prefix.

This layer is univariate and uses commutative local fields. General series
inversion/composition, multivariate differential calculus and ODE solvers
remain subsequent work. The :doc:`volkenborn` integrator now accepts
polynomials and these bounded analytic families. Local
primitives and measure integrals have different contracts. General p-adic
functions can have derivative zero without being constant; their global
integration theory cannot be inferred from this analytic-disc construction.

Independent tests use SymPy rational functions and multiplication matrices,
whole-input/coefficient perturbations and admissible omitted monomials at
later p-power denominators. PARI checks values, first/second derivatives and
normalized primitives at translated centers for p=2,3,5, plus both Gaussian
orientations over Q2(i) and Q3(i). The quadratic checks use independent
scalar identities, not native ramified POLMOD exp/log.

Run ``python examples/local_calculus.py`` for a dependency-free example.
