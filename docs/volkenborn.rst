Volkenborn integration
======================

``volkenborn_integral`` evaluates the normalized functional

.. math::

   I_p(f)=\lim_{N\to\infty}\frac{1}{p^N}
             \sum_{j=0}^{p^N-1}f(j),

for exact polynomials or the coefficient-bounded analytic families of
:doc:`calculus`. Its parameter belongs to Zp; values may belong to Qp or
a certified quadratic coefficient field. For an extension field the
parameter is still the embedded Zp, rather than its whole ring of integers.

Polynomial moments
-------------------

With B_1=-1/2, the moments are

.. math::

   I_p(x^n)=B_n,\qquad I_p\left(\sum c_nx^n\right)=\sum c_nB_n.

The result is independent of p for rational polynomials. The finite
power-sum identity in `DLMF 24.4.7 <https://dlmf.nist.gov/24.4.E7>`_
expresses the normalized sum as a polynomial in M=p**N whose constant
coefficient is B_n. Since M tends to zero p-adically, its limit is B_n.
This also fixes the sign convention at degree one. Related Volkenborn
identities are discussed by `Simsek <https://arxiv.org/abs/1910.09296>`_.

.. doctest::

   >>> x = Polynomial([0, 1])
   >>> volkenborn_integral(x)
   Fraction(-1, 2)
   >>> volkenborn_integral(x*x)
   Fraction(1, 6)
   >>> K = Qp(3, 3)
   >>> volkenborn_integral(x*x, K, precision=21).contains(Fraction(1, 6))
   True

Without an output field, a Polynomial with exact rational affine data
returns an exact Fraction. Supplying a field returns an absolute-precision
ball. AnalyticSeries infers its coefficient field and preserves finite
coefficient uncertainty. Formal PowerSeries prefixes and arbitrary
callables are rejected: they do not establish the necessary tail bound.

Affine parametrization and translation
--------------------------------------

``shift=a, scale=b`` means I_p(t -> f(a+b*t)). No Jacobian or volume
factor is inserted. The whole affine chart must belong to the declared
series disc: both v_pi(a-center) and v_pi(b) must be at least its radius
exponent. Exact affine data stays exact; compatible finite balls retain
their uncertainty in field-valued evaluation.

This functional is not translation invariant. For example, both t and
t+1 range over Zp but their Volkenborn integrals differ:

.. doctest::

   >>> volkenborn_integral(x, shift=1)
   Fraction(1, 2)
   >>> volkenborn_integral(x, shift=1)-volkenborn_integral(x)
   Fraction(1, 1)
   >>> volkenborn_integral(Polynomial([1]), scale=5)
   Fraction(1, 1)

For polynomials, I_p(f(t+1))-I_p(f(t))=f'(0). Thus changing an affine
representative can change the result even when its image is the same ball.
Keep this API distinct from real Haar integration, bounded p-adic measures
and local primitives. Their definitions are in :doc:`harmonic`,
:doc:`measures` and :doc:`calculus`.

Analytic integrands
-------------------

An AnalyticSeries with m known coefficients and tail_bound=(A,B) controls
every omitted coefficient, not just observed terms. Set
s=min(v_pi(shift-center),v_pi(scale)). Domain checking gives s>=radius
and therefore s+B>0. Every omitted integrated term has valuation at least

.. math::

   v_\pi\bigl(a_n I_p((d+bt)^n)\bigr)\geq A+(B+s)n-e,
   \qquad d=\mathrm{shift}-\mathrm{center},\quad v_\pi(p)=e.

The entire integrated tail is consequently contained in

.. math::

   O\left(\pi^{\lceil A+(B+s)m-e\rceil}\right).

For a custom series this guarantee is conditional on its explicit
coefficient-bound premise. Function factories prove their own bounds.
Default output precision is capped by coefficient, affine-parameter and
tail uncertainty; an unsupported explicit target raises PrecisionError.

.. doctest::

   >>> L = Q2i(3)
   >>> s = AnalyticSeries.from_function(exp, L, order=100, precision=21)
   >>> value = volkenborn_integral(s, scale=(4, 4), precision=21)
   >>> z = L(4, 4, absprec=40)
   >>> value.is_congruent(z/(exp(z)-1), 21)
   True
   >>> volkenborn_integral(s)
   Traceback (most recent call last):
   ...
   padic_tools.core.DomainError: the whole affine Zp chart must lie in the declared series disc

The example uses the identity I_p(exp(z*t))=z/(exp(z)-1) on its
convergence disc. At z=0 the constant function has integral one; no
division by a zero-containing exponential difference is required by the
integrator itself. A translation multiplies this example by exp(shift).

Why the bound covers the limit
------------------------------

The integer sum of j**k over j=0,...,p**N-1 is divisible by p**(N-1).
This follows by induction: split j=u+p**(N-1)*v into p blocks and expand
binomially. Thus normalized monomial averages have valuation at least -1
in Qp, uniformly in N and k. Expanding (d+b*t)**n and using the
ultrametric inequality gives

.. math::

   v_\pi\left(\frac{1}{p^N}\sum_{j=0}^{p^N-1}(d+bj)^n\right)
       \geq ns-e.

Multiplication by the omitted coefficient bound yields the preceding
positive-slope estimate, uniformly in N. The integrated tail tends to
zero uniformly, so finite moment limits can be interchanged with the
convergent series. This proves both existence for the supported family
and the returned error enclosure. A first zero Bernoulli moment does not
annihilate subsequent omitted terms.

Finite precision and constant charts
-------------------------------------

Bernoulli denominators can cost coefficient digits. Exact zeros among the
moments can instead remove a coefficient's uncertainty entirely. An
exactly zero scale is a constant chart and uses the sharper function-value
tail bound. A zero-containing scale ball remains uncertain.

.. doctest::

   >>> f = AnalyticSeries.from_polynomial(L, [0, L(1, absprec=9)], precision=20)
   >>> volkenborn_integral(f).precision
   7
   >>> family = AnalyticSeries(L, [0], radius=1, tail_bound=(0, 0), precision=30)
   >>> volkenborn_integral(family, scale=0, precision=60).precision
   60
   >>> volkenborn_integral(family, scale=L.O(7)).precision
   5

Independent tests compare exact moments with symbolic Faulhaber limits,
actual finite Riemann sums, translation identities and rational matrix
models. Affine/coefficient perturbations, omitted monomials and dyadic
precision losses exercise whole families. PARI checks polynomial moments
and exponential integral identities, including both Gaussian orientations
over Q2(i) and Q3(i), through independent scalar p-adic identities.

The native moment/substitution algorithms target modest degrees. General
integrable callables, unbounded distribution APIs and multivariate
integration remain subsequent work. Run ``python examples/volkenborn.py``.
