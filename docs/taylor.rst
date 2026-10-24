Taylor expansions
==================

The center, the formal order and the numerical precision are distinct:

.. math::

   f(a+h)=\sum_{n=0}^{N-1}\frac{f^{(n)}(a)}{n!}h^n+R_N(h).

Exact symbolic coefficients
---------------------------

For a rational polynomial and a rational center, ``taylor`` returns a
``PowerSeries`` in h=x-center. ``Polynomial.shift`` gives the whole exact
polynomial f(center+h). Both use exact rational arithmetic.

.. doctest::

   >>> f = Polynomial([1, 2, 3])
   >>> f.shift(2).coefficients
   (Fraction(17, 1), Fraction(14, 1), Fraction(3, 1))
   >>> f.taylor(2, order=4).coefficients
   (Fraction(17, 1), Fraction(14, 1), Fraction(3, 1), Fraction(0, 1))
   >>> taylor(exp, order=8) == exp_series(8)
   True
   >>> taylor(log, 1, order=8) == log1p_series(8)
   True

Exact formal jets are also provided for sin, cos and log1p at zero.
Their coefficients are rational. A general nonzero analytic center requires
an explicit p-adic ball because its function value need not be rational.
An arbitrary numerical callable is rejected; finite differences do not
provide a certified Taylor expansion.

Local numerical expansions
---------------------------

At a Qp or certified quadratic center, ``taylor`` returns a ``TaylorExpansion``. Supported
functions are polynomials, exp, log, log1p, sin and cos. The logarithm uses
the same default Iwasawa branch log(p)=0 as the analytic API. Absolute
precision counts the target uniformizer pi, with v_pi(p)=e.

.. doctest::

   >>> K = Qp(5, 18)
   >>> jet = taylor(exp, K(5), order=24)
   >>> jet.order
   24
   >>> jet(K(10)).is_congruent(exp(K(10)), precision=18)
   True
   >>> short = taylor(exp, K(5), order=3)
   >>> short(K(10)).precision < jet(K(10)).precision
   True
   >>> logarithm = taylor(log, K(5), order=24)
   >>> logarithm(K(30)).is_congruent(log(K(30)))
   True
   >>> f.taylor(K(2), order=2)(K(7)).contains(f(Fraction(7)))
   True
   >>> L = Q2i(3)
   >>> center = L(4, 4, absprec=21)
   >>> point = center+L(0, 16, absprec=23)
   >>> taylor(exp, center, order=24)(point).is_congruent(exp(point), 21)
   True

``coefficients[n]`` encloses f^(n)(center)/n!. Evaluating the object
propagates their uncertainty and the argument's uncertainty, then adds the
entire omitted tail. ``tail_precision(x)`` reports only the tail bound;
it cannot restore digits already lost in coefficients or input. A short
expansion therefore returns a coarse but valid ball.

Remainder bounds and convergence
--------------------------------

For exp, sin and cos, both the center and displacement must lie in the
supported exponential disc. Set h=x-center and v=v_pi(h)>e/(p-1).
Every derivative value at the center is p-integral, and Legendre's bound
gives, for all k>=N,

.. math::

   v_\pi\left(\frac{f^{(k)}(a)h^k}{k!}\right)
      \geq k\left(v-\frac{e}{p-1}\right),\qquad
   R_N(h)\in O\left(\pi^{\left\lceil N(v-e/(p-1))\right\rceil}\right).

For log at a, put u=h/a. For log1p at a, put u=h/(1+a).
Require v_pi(u)>=1. The omitted log(1+u) terms satisfy

.. math::

   v_\pi(u^k/k)\geq k\,v_\pi(u)-e\lfloor\log_p k\rfloor,
   \qquad
   R_N(h)\in O\left(\pi^{T(N,v_\pi(u))}\right).

Here T(N,s) is the minimum of ``k*s-e*floor(log_p(k))`` over every k>=N,
computed by inspecting later powers of p as described in
:doc:`extension_analytic`. For e=1 the expression is already nondecreasing;
for ramified fields it can decrease later. The envelope bounds the whole
tail, including those denominator spikes.
For polynomials, the finite omitted coefficients supply the bound directly;
there is no convergence restriction. If none are omitted,
``tail_precision`` returns ``None``: only coefficient/input uncertainty remains.

The analytic domains follow the sources in :doc:`references`; the explicit
inequalities above state the numerical contracts used by this implementation.
