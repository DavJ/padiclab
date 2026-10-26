Linear differential systems
===========================

``solve_linear_system`` encloses local solutions of

.. math::

   Y'(x)=A(x)Y(x)+B(x),\qquad Y(c)=Y_0,

over Qp and both certified quadratic fields. A is square, while Y and B
can have several columns. A vector is represented by one column. This
extends :doc:`ode` to coupled equations, fundamental matrices and several
initial conditions, with the same exact/ball arithmetic and explicit
coefficient-tail premises as :doc:`calculus`.

Vectors and fundamental matrices
--------------------------------

Entries of A may be scalars or AnalyticSeries at the same exact center.
``initial`` is a scalar vector or a Matrix. A homogeneous example is the
coupled oscillator y_1'=y_2, y_2'=-y_1:

.. doctest::

   >>> K = Qp(3, 3)
   >>> y = solve_linear_system([[0, 1], [-1, 0]], [1, 0], field=K, order=64)
   >>> y.shape
   (2, 1)
   >>> x = K(9, absprec=45)
   >>> y(9, precision=21).is_congruent(Matrix(K, [[cos(x)], [-sin(x)]]), 21)
   True
   >>> y[1, 0].derivative()(9, precision=21).is_congruent(-cos(x), 21)
   True

``y[i,j]`` is one component AnalyticSeries. Calling the solution returns
a Matrix of local balls; ``precisions`` on that Matrix reports the available
digits separately. The result's ``field``, ``center``, ``radius``, ``order``,
``precision`` and ``tail_bound`` describe its shared context and certificate.
``entries`` and row indexing expose its immutable rows of component series.

Use a Matrix for multiple initial columns. This avoids ambiguity with exact
quadratic coordinate pairs in a vector. A rational identity initial matrix
gives the normalized fundamental solution:

.. doctest::

   >>> phi = solve_linear_system([[0, 1], [-1, 0]], Matrix.identity(exact, 2), field=K, order=64)
   >>> expected = Matrix(K, [[cos(x), sin(x)], [-sin(x), cos(x)]])
   >>> phi(9, precision=21).is_congruent(expected, 21)
   True
   >>> phi(0, precision=60).contains([[1, 0], [0, 1]])
   True

Rectangular or singular initial matrices are allowed. The solver does not
invert Y_0. Matrix(exact, ...) preserves rational source data; a local-field
Matrix retains its existing balls, including finite uncertainty at zero
centers. Default output precision comes from the supplied series, or from
the field if there are none. An explicit solver precision sets a default
target, and an unsupported explicit evaluation target raises PrecisionError.

Field and exact center are inferred from series entries. A local Matrix
base can also supply the field; with only exact scalars, specify ``field``.
The center defaults to zero when no series supplies it. ``field`` and
``center`` can be explicit, and every series entry must match them.

Variable coefficients and forcing
---------------------------------

``forcing=None`` denotes zero. For a vector solution, supply a vector of
scalar/series entries or a matching Matrix. For multiple columns, supply
rectangular rows of scalar/series entries or a matching Matrix. Forcing and
initial shapes must agree.

The noncommuting example

.. math::

   A(h)=\begin{pmatrix}0&1\\h&0\end{pmatrix},\qquad h=x-c,

reduces to the normalized Airy equation u''=h*u. Its rational normalized
solutions are the two coefficient series in `DLMF 9.4.1
<https://dlmf.nist.gov/9.4.E1>`_; the equation itself is in
`DLMF 9.2.1 <https://dlmf.nist.gov/9.2.E1>`_. Here these rational series
are evaluated on the certified p-adic disc. No complex Ai(0) or Ai'(0)
normalization constant is imported into the local field.

.. doctest::

   >>> h = AnalyticSeries(K, [0, 1], tail_bound=None, precision=30)
   >>> airy = solve_linear_system([[0, 1], [h, 0]], Matrix.identity(exact, 2), order=64)
   >>> airy[0, 0][3]
   Fraction(1, 6)
   >>> airy(9, precision=21).determinant().contains(1)
   True
   >>> b = [[1, -h*h], [-h*h, 2*h]]
   >>> forced = solve_linear_system([[0, 1], [h, 0]], Matrix.identity(exact, 2), forcing=b, order=64)
   >>> (forced(9, precision=21)-airy(9, precision=21)).contains([[9, 0], [0, 81]])
   True

The forced solution is the fundamental matrix plus diag(h,h**2).
For variable A, the recurrence preserves the order A_j*Y_(n-j).
In this example, exp(integral(A)) would give the incorrect cubic coefficient
1/4 in the first entry; the correct value is 1/6. No commutation of
coefficient matrices is assumed.

Quadratic coefficients use the same interface:

.. doctest::

   >>> L = Q2i(3)
   >>> y = solve_linear_system([[0, (-1, 1)], [(1, -1), 0]], [1, 0], field=L, order=64)
   >>> z = L(-1, 1, absprec=45)*L(4, 4, absprec=45)
   >>> y((4, 4), precision=21)[0, 0].is_congruent(cos(z), 21)
   True

Here (-1,1) is i=pi-1 in the ramified basis. Exact coordinate pairs in
initial vectors remain scalars; multiple initial columns require Matrix.

Matrix majorant and proof
-------------------------

Set v(M) to the minimum v_pi of its entries, with v(0)=infinity. For
compatible finite matrices, the ultrametric inequality gives

.. math::

   \begin{aligned}
   v(UV)&\geq v(U)+v(V),\\
   v(U+V)&\geq\min\{v(U),v(V)\}.
   \end{aligned}

There is no real-norm dimension factor in this valuation estimate. The
coefficient recurrence and factorial-normalized form are

.. math::

   \begin{aligned}
   (n+1)Y_{n+1}&=B_n\\
                &\quad+\sum_{j=0}^n A_jY_{n-j},\\
   Z_{n+1}&=n!B_n\\
           &\quad+\sum_{j=0}^n\frac{n!}{(n-j)!}A_jZ_{n-j},\\
   Z_n&=n!Y_n.
   \end{aligned}

Collect the scalar multiplier constraints of :doc:`ode` over every entry
of A, and all infinite forcing-tail slopes. Their minimum L bounds every
matrix coefficient. Exact structural zeros impose no constraint; if the
set is empty, choose L=0. The global forcing intercept is the minimum
of the scalar global intercepts of all B entries. Thus

.. math::

   \begin{aligned}
   v(A_j)&\geq L(j+1),\\
   v(B_j)&\geq C_B+Lj,\\
   H&=\min\{v(Y_0),C_B-L\},\\
   E&=e/(p-1),\qquad v_\pi(p)=e.
   \end{aligned}

The falling factorials are integral, and matrix multiplication satisfies
the preceding valuation inequality without changing product order.
Induction therefore gives v(Z_n)>=H+L*n. Since v_p(n!)<=n/(p-1),
all entries of every omitted coefficient obey the common certificate

.. math::

   \begin{aligned}
   v(Y_n)&\geq H+(L-E)n,\\
   R_M(h)&\in O\left(\pi^{\lceil H+(s+L-E)M\rceil}\right),\\
   s&=v_\pi(h).
   \end{aligned}

The automatic radius exponent is the maximum input radius and
floor(E-L)+1. On this disc, s+L-E>0, so the whole omitted tail tends
to zero uniformly. This proves analytic existence; the coefficient
recurrence gives uniqueness for each fixed set of input functions and
initial values within analytic power series at c. An explicit radius can
shrink the sufficient disc. The bound need not be optimal.

Uncertain coefficients, initial entries and argument balls are propagated
entrywise. Unknown input coefficients are bounded balls from their stated
tail premises. Returned component balls discard correlations and can
represent additional functions beyond the actual solution family. The
common majorant can also hide stronger bounds of individual components.
An exact zero initial matrix and forcing give an exact zero solution.
An exact zero multiplier and fully covered polynomial forcing give a
polynomial solution with zero tail. Zero-containing balls retain uncertainty.

Independent tests use rational block-matrix exponentials, normalized Airy
hypergeometric coefficients, whole input perturbations and later p-power
denominator spikes. PARI checks Airy Maclaurin sums for p=2,3,5 and coupled
Gaussian oscillators in both local orientations using scalar analytic
identities. These checks do not call generic ramified POLMOD analytics.

Both matrix dimensions are limited to 64. Arithmetic takes O(M**2*d**2*m)
scalar operations for d equations, m columns and M stored coefficients.
Nonlinear or singular systems, analytic continuation, optimal lattice
precision and accelerated solvers remain subsequent work. Run
``python examples/linear_system.py`` for a dependency-free example.
