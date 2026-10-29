Reciprocals and centered composition
====================================

AnalyticSeries supports multiplicative reciprocals and centered composition
over Qp and both certified quadratic fields. The output is an enclosure of
all allowed input functions, including finite coefficient balls and every
omitted coefficient under its stated tail premise. It carries a sufficient
closed disc and a proved whole-tail bound, as in :doc:`calculus`.

Multiplicative reciprocals
--------------------------

``f.reciprocal()`` represents 1/f. Its constant coefficient must be resolved
and nonzero. Exact zero raises ZeroDivisionError; a ball containing zero
raises PrecisionError. The automatic disc ensures constant-term dominance
uniformly over the whole input family, so all represented denominators are
nonzero there.

.. doctest::

   >>> K = Qp(3, 3)
   >>> f = AnalyticSeries.from_polynomial(K, [1, -1])
   >>> inverse = f.reciprocal(order=48, precision=21)
   >>> inverse.radius
   1
   >>> inverse(9, precision=21).contains(Fraction(-1, 8))
   True
   >>> inverse(0, precision=60).contains(1)
   True
   >>> inverse(1)
   Traceback (most recent call last):
   ...
   padic_tools.core.DomainError: the whole argument ball must lie in the declared series disc

``order`` defaults to the positive length of the input prefix, or one for
an empty prefix. ``precision`` is a default evaluation target. Available
input and tail information can cap the result; an unsupported explicit
evaluation target raises PrecisionError. Exact source data can support
targets above the field's default, while existing balls retain uncertainty.
Absolute precision and radius exponents always use uniformizer units.

The reciprocal also preserves exact quadratic source coordinates:

.. doctest::

   >>> L = Q2i(3)
   >>> f = AnalyticSeries(L, [1, (-1, 1), 1], tail_bound=None)
   >>> inverse = f.reciprocal(order=32)
   >>> z = L(4, 4, absprec=45)
   >>> expected = 1/(1+L(-1, 1, absprec=45)*z+z*z)
   >>> inverse((4, 4), precision=21).is_congruent(expected, 21)
   True

Here (-1,1) is i=pi-1 in the ramified basis. A constant polynomial's
reciprocal has an exact zero tail, even when the nonzero constant is a
finite ball. A nonconstant denominator generally has an infinite inverse.

Reciprocal certificate
----------------------

Write f(c+h)=sum a_n*h**n and t=v_pi(a_0). Choose a rational lambda such
that every positive coefficient satisfies

.. math::

   v_\pi(a_n)\geq t+\lambda n,\qquad n\geq1.

The implementation takes the minimum of (v_pi(a_n)-t)/n over the known
positive coefficients. For an infinite tail (A,B), also include B and
B+(A-t)/m_plus, where m_plus=max(1,order). These bound the entire tail
because B+(A-t)/n is monotone. Exact zero coefficients impose no constraint.

Take r at least the input radius and floor(-lambda)+1. Then
q=(f-a_0)/a_0 has positive valuation throughout v_pi(h)>=r, uniformly in
the coefficient family, and its geometric reciprocal converges. Coefficients
can be computed from the finite recurrence

.. math::

   \begin{aligned}
   b_0&=a_0^{-1},\\
   b_n&=-a_0^{-1}\sum_{j=1}^n a_jb_{n-j}.
   \end{aligned}

Induction gives the global bound below. For M stored coefficients and
s=v_pi(h)>=r, the whole omitted reciprocal tail consequently obeys

.. math::

   \begin{aligned}
   v_\pi(b_n)&\geq-t+\lambda n,\\
   R_M(h)&\in O\left(\pi^{\lceil-t+(s+\lambda)M\rceil}\right).
   \end{aligned}

Since r+lambda>0, this proves analytic convergence and the enclosure for
every allowed denominator. The sufficient disc can be smaller than the
actual pole-free domain. If all positive coefficients are exact structural
zeros and the input is polynomial, no extra radius restriction is needed.

Centered composition
---------------------

``outer.compose(inner)`` represents outer(inner(x)). The fields must match.
Centers can differ, but the exact inner constant must equal the exact outer
center. A ball with that center does not establish exact equality and is
rejected, even if it contains the required value. The output center is the
inner series' center.

.. doctest::

   >>> f = AnalyticSeries(K, [1, -1], tail_bound=None).reciprocal(order=48)
   >>> g = AnalyticSeries(K, [0, 1, 1], tail_bound=None)
   >>> composed = f.compose(g)
   >>> composed(9, precision=21).contains(Fraction(-1, 89))
   True
   >>> exponential = AnalyticSeries.from_function(exp, K, order=48)
   >>> logarithm = AnalyticSeries.from_function(log, K, center=1, order=48)
   >>> logarithm.compose(exponential)(9, precision=21).contains(9)
   True
   >>> f.compose(AnalyticSeries(K, [K(0), 1], tail_bound=None))
   Traceback (most recent call last):
   ...
   TypeError: composition requires an exact inner constant coefficient

The default order is the full possible degree plus one for finite polynomial
composition; it is the infinite input's prefix length when the other input
is polynomial, or the smaller prefix length when both have infinite tails.
A constant outer polynomial needs only one coefficient. At least one
coefficient is stored. Explicit ``order`` can request a different positive
length, retaining unknown input coefficients as bounded balls.

.. doctest::

   >>> f = AnalyticSeries(K, [2, 3, 1], center=5, tail_bound=None)
   >>> g = AnalyticSeries(K, [5, 1, 1], center=-1, tail_bound=None)
   >>> result = f.compose(g)
   >>> result.order, result.center, result.tail_bound, result.radius
   (5, Fraction(-1, 1), None, None)
   >>> result(0, precision=60).contains(12)
   True

Fully stored polynomial results retain zero tails. A deliberately shorter
prefix instead gets a whole-tail bound on a declared disc. An exact constant
inner function gives a constant result, including the outer constant's
uncertainty when that coefficient is unknown. Observed coefficient
cancellation does not determine later coefficients.

Disc image and composition majorant
-----------------------------------

Set h=x-c_inner and d(h)=inner(x)-c_outer. Exact alignment makes d_0=0.
For every positive known inner coefficient record (k,v_pi(d_k)). For an
infinite inner tail (A_g,B_g), also record its first positive tail index
m_plus=max(1,inner.order) with bound A_g+B_g*m_plus. Define C(w) as the
minimum of v+k*w over these pairs. When the inner tail is infinite,
w+B_g>0 ensures its whole weighted tail is bounded by that endpoint.

The actual integer radius R must respect the inner declared domain and
the entire image must belong to the outer disc: C(R)>=r_outer when that
disc is restricted. These conditions apply to every represented inner
function and argument, including zero-containing argument balls.

Separately choose an auxiliary rational weight rho<R. If the inner tail is
infinite, require rho+B_g>0. If the outer tail is infinite with slope B_f,
require C(rho)+B_f>0. These are strict rational linear inequalities in rho;
the implementation picks an admissible integer R and a rational midpoint
when strict bounds are present. For a truncated polynomial prefix, any
auxiliary rho<R suffices when there are no strict constraints.
The auxiliary weight need not satisfy the outer declared domain: only the
actual argument disc is evaluated. Set C=C(rho) and

.. math::

   H=\inf_{n\geq0}\{v_\pi(a_n)+Cn\},

using the full coefficient bound for the outer tail. Weighted multiplication
gives a lower bound j*C-rho*n for [h**n]d(h)**j. Since d_0=0, only j<=n
can affect a fixed coefficient. The finite sum then proves

.. math::

   \begin{aligned}
   v_\pi([h^n]f(g(x)))&\geq H-\rho n,\\
   R_M(h)&\in O\left(\pi^{\lceil H+(s-\rho)M\rceil}\right),\\
   s&=v_\pi(h)\geq R>\rho.
   \end{aligned}

This covers every omitted output coefficient, including those contributed
by unknown coefficients of either input. The positive gap R-rho proves
convergence and a whole-tail enclosure. Exact coefficient arithmetic or
ball arithmetic computes the prefix by truncated Horner substitution.
The native algorithm is cubic in the requested order; no fast composition
algorithm is claimed.

For a fully stored polynomial or constant result, only the declared input
and image-domain restrictions are needed; the infinite-tail weight is
unnecessary. A zero polynomial still respects its declared input domain.
Finite uncertainty is retained throughout. Independent coefficient balls
can include functions beyond the correlated reciprocal/composition family,
so exact functional identities do not in general imply an exact zero tail.

Scope and sources
-----------------

This step supplies multiplicative inversion and exactly aligned composition.
General recentering, uncertain inner constants, compositional reversion,
nonlinear ODEs and optimal correlated precision remain subsequent work.
For formal coefficient counterparts and related non-Archimedean norms see
Sage's `lazy series composition
<https://doc.sagemath.org/html/en/reference/power_series/sage/rings/lazy_series.html>`_
and `Tate algebras
<https://doc.sagemath.org/html/en/reference/power_series/sage/rings/tate_algebra.html>`_.
The certificates used here are derived above. Tests execute independent
exact matrix/binomial oracles and PARI/GP rational/transcendental identities;
they do not execute Sage. Run ``python examples/analytic_operations.py``.
