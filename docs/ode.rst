Local linear differential equations
===================================

``solve_linear_ode`` encloses the analytic power-series solution of

.. math::

   y'(x)=a(x)y(x)+b(x),\qquad y(c)=y_0.

The multiplier ``a`` is an AnalyticSeries over Qp or either certified
quadratic field. ``forcing=b`` is a compatible series at the same exact
center, or a scalar constant. ``initial`` defaults to one and ``forcing``
to zero. Finite coefficient and initial balls retain their uncertainty.
The returned AnalyticSeries supplies the coefficient prefix, a closed disc
and a proved whole-tail bound. Custom input tail bounds remain mathematical
premises, as in :doc:`calculus`.

Solving and evaluating
----------------------

The example y'=-y+2, y(c)=3 has solution 2+exp(-(x-c)). Its exact
source coefficients permit evaluation above the field's default precision:

.. doctest::

   >>> K = Qp(3, 3)
   >>> c = Fraction(1, 3)
   >>> a = AnalyticSeries.from_polynomial(K, [-1], center=c)
   >>> y = solve_linear_ode(a, 3, forcing=2, order=64, precision=21)
   >>> y.radius
   1
   >>> y(c+9, precision=21).is_congruent(2+exp(K(-9, absprec=45)), 21)
   True
   >>> y(c, precision=60).contains(3)
   True

``order`` is the positive number of stored solution coefficients, not a
number of p-adic digits. ``precision`` selects the usual default output
target; evaluation caps it by the available coefficient, input and tail
information. An unsupported explicit evaluation target raises
PrecisionError. Increasing ``order`` cannot recover missing input digits.

Variable coefficients and quadratic values
------------------------------------------

For a=1+h, b=-1-3h-h**2 and y(0)=3, the solution is
2+h+exp(h+h**2/2). Coefficients use h=x-c:

.. doctest::

   >>> a = AnalyticSeries(K, [1, 1], tail_bound=None)
   >>> b = AnalyticSeries(K, [-1, -3, -1], tail_bound=None)
   >>> y = solve_linear_ode(a, 3, forcing=b, order=64)
   >>> y[3]
   Fraction(2, 3)
   >>> h = K(9, absprec=45)
   >>> y(9, precision=21).is_congruent(2+h+exp(h+h*h/2), 21)
   True
   >>> L = Q2i(3)
   >>> a = AnalyticSeries(L, [(-1, 1)], tail_bound=None)
   >>> y = solve_linear_ode(a, order=64)
   >>> z = L(-1, 1, absprec=45)*L(4, 4, absprec=45)
   >>> y((4, 4), precision=21).is_congruent(exp(z), 21)
   True

In the last example (-1,1) is i=pi-1 in the ramified basis pi=1+i.
The same implementation accepts genuinely quadratic variable coefficients,
unramified bases, translated centers and finite coefficient balls.

Coefficient recurrence and convergence certificate
--------------------------------------------------

Write a, b and y in powers of h=x-c. Comparing coefficients gives

.. math::

   \begin{aligned}
   (n+1)y_{n+1}&=b_n+\sum_{j=0}^n a_jy_{n-j},\\
   y_0&=\mathrm{initial}.
   \end{aligned}

Known coefficients use exact source arithmetic or ball arithmetic. Input
coefficients beyond a stored prefix are balls from its full tail bound;
they are never set to exact zero. Division by n+1 loses e*v_p(n+1)
coefficient digits, where v_pi(p)=e.

The following factorial majorant controls every later coefficient. Choose
a rational L such that v_pi(a_j)>=L*(j+1) for all j. The implementation
takes the minimum of each known coefficient's valuation divided by j+1,
and, for an infinite multiplier tail (A_a,B_a) starting at m_a, the two
quantities B_a and (A_a+B_a*m_a)/(m_a+1). The latter work because
(A_a+B_a*j)/(j+1) is monotone. Exact zero coefficients impose no constraint.
Also include the forcing tail's slope B_b when present. If this candidate
set is empty, L=0 suffices.

Let C_b be a global intercept with v_pi(b_j)>=C_b+L*j. It is the
minimum of the known bounds v_pi(b_j)-L*j and the forcing tail bound
A_b+(B_b-L)*m_b, when present. Our choice L<=B_b makes this valid at
every later index. Set

.. math::

   H=\min\{v_\pi(y_0),C_b-L\},\qquad E=\frac{e}{p-1}.

With z_n=n!*y_n, the recurrence becomes

.. math::

   \begin{aligned}
   z_{n+1}&=n!b_n\\
           &\quad+\sum_{j=0}^n\frac{n!}{(n-j)!}a_jz_{n-j}.
   \end{aligned}

The falling factorials are integers. Induction gives
v_pi(z_n)>=H+L*n: each multiplier term has valuation at least
H+L*(n+1), and the forcing term does too because H<=C_b-L.
The initial case follows from H<=v_pi(y_0). Finally,
v_p(n!)<=n/(p-1) gives the global certificate

.. math::

   v_\pi(y_n)\geq H+(L-E)n.

The automatic radius is the largest input radius exponent or
floor(E-L)+1, whichever is larger. An explicit radius can only shrink
this disc. On it, r+L-E>0, so the series and its derivative converge
and satisfy the equation term by term. For M stored coefficients and
s=v_pi(x-c)>=r, the entire omitted solution tail lies in

.. math::

   O\left(\pi^{\lceil H+(s+L-E)M\rceil}\right).

This proves existence and an error enclosure; the coefficient recurrence
proves uniqueness within analytic power series at the chosen center for
each fixed input. The convergence disc is sufficient and can be smaller
than the solution's actual disc. General locally analytic functions can
have different constants on disjoint subdiscs, so that wider function
class is not the uniqueness claim here.

Families, zero solutions and limitations
----------------------------------------

For uncertain inputs the output encloses the corresponding solutions.
Its independent coefficient balls and tail also represent other functions:
correlations between solution coefficients are not retained. A residual
ball containing zero is a consistency check, rather than an independent
existence or uniqueness proof.

.. doctest::

   >>> a = AnalyticSeries(K, [1], radius=1, tail_bound=(0, 0), precision=30)
   >>> enclosure = solve_linear_ode(a, order=30)
   >>> enclosure[2].contains(1)
   True
   >>> enclosure(9).precision < 21
   True
   >>> enclosure(9, precision=21)
   Traceback (most recent call last):
   ...
   padic_tools.core.PrecisionError: requested precision exceeds coefficient, input or tail information

Here one allowed multiplier is 1/(1-h), whose solution with initial one
is 1/(1-h). Other omitted coefficients allow other solutions; this is why
knowing more computed solution terms alone does not imply 21 digits.

If the multiplier is structurally zero and polynomial forcing is fully
covered by ``order``, the solver returns a polynomial with an exact zero
tail. A structural zero initial condition and forcing also give an exact
zero solution for every allowed multiplier. Zero-containing balls do not
trigger either exact-zero argument.

.. doctest::

   >>> zero = AnalyticSeries.from_polynomial(K, [])
   >>> b = AnalyticSeries.from_polynomial(K, [1, 1])
   >>> y = solve_linear_ode(zero, 2, forcing=b, order=3)
   >>> y.tail_bound is None and y.radius is None
   True
   >>> y(1, precision=60).contains(Fraction(7, 2))
   True

This API handles scalar first-order linear equations at ordinary points.
The :doc:`ode_system` API handles coupled linear systems and matrix initial
data. Nonlinear equations, singular points, analytic continuation and
optimal differential/lattice precision remain subsequent work. Complexity
is quadratic in the prefix length; no fast Newton solver is claimed.

Related primary sources are Sage's
`formal solve_linear_de documentation
<https://doc.sagemath.org/html/en/reference/power_series/sage/rings/power_series_ring_element.html#sage.rings.power_series_ring_element.PowerSeries.solve_linear_de>`_
and `Lairez and Vaccon, On p-adic differential equations with separation
of variables <https://arxiv.org/abs/1602.00244>`_ for sharper precision methods
in that equation class. The conservative majorant used here is proved
above; tests execute SymPy and PARI/GP, not Sage. Run
``python examples/linear_ode.py`` for a dependency-free example.
