Analysis over quadratic local fields
====================================

``exp``, ``log1p``, ``log``, ``sin``, ``cos``, ``power`` and local
``taylor`` accept the certified :doc:`unramified <extensions>` and
:doc:`Eisenstein <ramified>` quadratic fields, including Q2(i). The
implementation needs only Python's standard library. A formal quadratic
coefficient algebra is not a supported analytic input: its field and
valuation would need to be certified first.

Valuation and domains
---------------------

Write L for the target field, pi for its chosen uniformizer, e for its
ramification index, and f for its residue degree. We always use
v_pi(pi)=1 and v_pi(p)=e. Thus ``absprec=N`` describes the whole ball
``c+O(pi**N)``. These are target digits, rather than base-p digits.

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Function
     - Supported domain
     - Default output precision
   * - exp, sin, cos
     - v_pi(x)>e/(p-1), for the whole ball
     - N
   * - log1p
     - v_pi(x)>=1
     - B(N), defined below
   * - log
     - Ball excluding zero
     - Relative input precision and branch uncertainty, as below
   * - power(base, exponent)
     - v_pi(base-1)>e/(p-1), exponent integral in L
     - Propagated through log1p, multiplication and exp

The trigonometric domain is a sufficient supported domain. The construction
follows the local-field exponential and logarithm conventions documented in
`Sage's p-adic element reference
<https://doc.sagemath.org/html/en/reference/padics/sage/rings/padics/padic_generic_element.html>`_.
The explicit tail and input-error bounds below specify this implementation's
precision contract.

.. doctest::

   >>> L = Q2i(precision=3)
   >>> x = L(4, 4, absprec=21)
   >>> (x.valuation(), exp(x).precision)
   (4, 21)
   >>> log(exp(x)).is_congruent(x, precision=21)
   True
   >>> (sin(x)**2+cos(x)**2).is_congruent(1, precision=21)
   True
   >>> exp(L(0, 2, absprec=21)).precision
   21
   >>> exp(L(2, absprec=21))
   Traceback (most recent call last):
   ...
   padic_tools.core.DomainError: requires v_pi(x) > e/(p-1), including the whole input ball

Here 2*pi has valuation three and is supported; 2 has valuation two and is
on the excluded exponential boundary. The requested 21 pi digits are
independent of the parent's default of three base digits. Existing balls
cannot be re-evaluated at higher precision without exact source data.

Certified series and input uncertainty
--------------------------------------

The series evaluator multiplies exact rational coordinate pairs modulo the
defining polynomial and rounds its sum only once. It does not repeatedly
round individual terms. For s=v_pi(x)>e/(p-1), Legendre's inequality gives

.. math::

   v_\pi(x^k/k!)\geq k\left(s-\frac{e}{p-1}\right).

The positive slope bounds every omitted exponential, sine or cosine term.
On this disc exp is a unit and exp(delta)-1 has valuation v_pi(delta);
the addition formulas give the same input-error bound for sine and cosine.
An input ideal O(pi**N) therefore remains within O(pi**N).

For log1p, all terms converge when s>=1, but the denominator bound can
decrease at a later power of p. Define a lower envelope for the entire tail:

.. math::

   T(K,s)=\min_{k\geq K}\bigl(ks-e\lfloor\log_p k\rfloor\bigr),
   \qquad B(N)=T(1,N).

This minimum is computed finitely. Between consecutive powers of p the
expression increases. Check K and the following powers p**j until the
successive minima increase; their difference is
``(p-1)*p**j*s-e``, which is increasing. Stop the series only when T(K,s)
reaches the requested precision. In the wild dyadic case e=2, s=1, the
k=3 term has valuation three while the k=4 term has valuation zero.

The same envelope controls input uncertainty. For a whole input ball
x=c+delta, with v_pi(delta)>=N and v_pi(c)>=1,

.. math::

   \log(1+c+\delta)-\log(1+c)
     =\log\left(1+\frac{\delta}{1+c}\right)\in O(\pi^{B(N)}).

Hence B(N) is a conservative absolute output precision. In our quadratic
fields it equals N, except that p=e=2 and N=1 gives B(1)=0. This is an
enclosure, not a claim of optimal precision after torsion cancellation.
An explicit ``precision`` above this bound raises ``PrecisionError``.

.. doctest::

   >>> log1p(L.O(1)).precision
   0
   >>> log1p(L.O(2)).precision
   2
   >>> log1p(L.O(1), precision=1)
   Traceback (most recent call last):
   ...
   padic_tools.core.PrecisionError: function evaluation cannot invent input digits

Logarithm branches and torsion
------------------------------

For a nonzero ball write x=pi**v*u with u a unit. Its known relative
precision is R=N-v. Put m=p**f-1; this rational integer is a unit, and
u**m-1 belongs to the maximal ideal. Define

.. math::

   \log_U(u)=\frac{\log(1+(u^m-1))}{m},\qquad
   \log(x)=\log_U(u)+v\log(\pi).

The unit part has precision B(R). This construction kills all roots of
unity, including p-power torsion, without selecting a public Teichmueller
lift in the extension. The identity follows from the convergent principal
logarithm's multiplicativity; torsion is killed because its finite order
times its logarithm is zero in characteristic zero.

The default branch is **log(p)=0**, also in a ramified field. In general
this does not mean log(pi)=0. The defining data determine the exact unit
u0=pi**e/p, so

.. math::

   \log(\pi)=\frac{\operatorname{p\_branch}+\log_U(u_0)}{e}.

Alternatively set ``pi_branch`` directly. Supply at most one branch.
Branches can be exact rational scalars, exact coordinate pairs, or finite
balls in the same field or its Qp base. Foreign primes and polynomial
presentations are rejected. Units do not depend on the branch.

.. doctest::

   >>> pi = L(0, 1, absprec=19)
   >>> log(pi).is_congruent(0)
   True
   >>> log(pi, p_branch=1).contains(Fraction(1, 2))
   True
   >>> log(L(2, absprec=19), pi_branch=1).is_congruent(2)
   True
   >>> E = EisensteinQuadraticField(Qp(3, 3), [-3, 3, 1])
   >>> rho = E(0, 1, absprec=21)
   >>> log(rho).is_zero()
   False
   >>> (2*log(rho)).is_congruent(log(rho*rho/3))
   True
   >>> log(E(3, absprec=21)).is_congruent(0)
   True

In Q2(i), pi**2=2*i and log(i)=0, so the default also has log(pi)=0.
The second example shows why that conclusion cannot be used for every
Eisenstein polynomial. As usual, log and exp are inverse on their common
deep local disc; exp(log(x)) need not equal x outside it.

Division by e costs ``e*v_p(e)`` pi digits. The exact defining unit u0
permits guard digits; a supplied finite branch does not. For example:

.. doctest::

   >>> branch = L(1, 1, absprec=7)
   >>> log(pi, p_branch=branch).precision
   5
   >>> log(pi, pi_branch=branch).precision
   7
   >>> log(L(Fraction(1, 4), absprec=13)).precision
   17

The last value gains absolute precision by removing the known negative
valuation: R=N-v. It does not gain significant input digits. Numerical
branch addition can reduce the final precision further, and the returned
ball always records it. ``precision=...`` only truncates known digits.

Powers, Taylor expansions and the bridge
----------------------------------------

``power`` uses exp(exponent*log1p(base-1)) on the stated deep principal
disc. Exponents may be integral extension elements, Qp elements of the
same prime, or exact rational scalars. Nonintegral exponents are rejected.
Integer ``x**n`` remains ordinary algebraic powering on all nonzero balls
for negative n, and on all balls for nonnegative n.

.. doctest::

   >>> (power(1+x, Fraction(1, 3))**3).is_congruent(1+x)
   True
   >>> jet = taylor(exp, x, order=24)
   >>> point = x+L(0, 16, absprec=23)
   >>> jet(point).is_congruent(exp(point), precision=21)
   True

``taylor`` accepts the same quadratic centers for polynomials and all five
elementary functions. ``tail_precision`` counts pi digits. Factorial tails
use e in Legendre's inequality; logarithmic tails use T(order,v_pi(ratio)).
The coefficients and argument retain their own uncertainty. Its logarithm
uses the default log(p)=0. See :doc:`taylor` for the coefficient contract.

An exact Gaussian source can first be mapped using :doc:`local_bridges`,
then evaluated in that selected local field:

.. doctest::

   >>> place = GaussianLocalPlace(Qp(2, 3), root_sign=1)
   >>> local = place(Gaussian(4, 4), absprec=21)
   >>> log(exp(local)).is_congruent(local, precision=21)
   True

The complex image of that source is a separate evaluation. No conversion
of a finite local result into a complex transcendental value is implied.
Quaternion and Clifford adapters continue to map coefficients algebraically;
noncommutative analytic functions are a separate future API. Public sqrt
and teichmuller, theta and KL-zeta evaluation still support Qp only.

Independent checks
-------------------

Tests use fixed, exact SymPy polynomials in rational 2x2 multiplication
matrices, whole-ideal perturbations, boundary balls, torsion, both branch
normalizations, negative valuations and explicit precision above the parent
default. They check norm/trace identities and Taylor tails as well as local
inverse identities.

Optional PARI/GP tests evaluate Gaussian exponentials using its scalar
p-adic exp, sin and cos, and Gaussian logarithms using scalar log of the
norm plus a separately truncated rational arctangent series. They compare
both coordinates and both Gaussian orientations over Q3(i) and Q2(i).
This is not a claim that PARI's generic polynomial-modulus exp/log provides
a native ramified extension oracle. CI installs GP and runs these checks.
