# padiclab

A Python library for finite precision **p-adic analysis**, exact **symbolic
algebra**, and explicit bridges to **complex numbers, quaternions,
biquaternions and Clifford algebras**.

This is a usable **0.1 research library**, not a claim to implement every part
of p-adic mathematics. Domains, precision bounds and remaining scope are
documented. The numerical core uses only Python's standard library.

Repository: [DavJ/padiclab](https://github.com/DavJ/padiclab).
The current Python distribution is `padic-tools`, imported as `padic_tools`.

```bash
python -m pip install .
# Optional symbolic expressions, complex evaluation and HTML documentation:
python -m pip install '.[symbolic,complex,docs]'
```

```python
from fractions import Fraction
from padic_tools import Qp, exp, log, theta3, zeta_p

K = Qp(5, precision=20)
x = K(Fraction(1, 3))
assert (x * 3).is_congruent(1)
assert log(exp(K(5))).is_congruent(5)

t = theta3(K(5))
z = zeta_p(-3, K, branch=0)  # Euler-corrected zeta(-3), on this branch
assert z.is_congruent(Fraction(-31, 30))
```

| Layer | Implemented in 0.1 |
|---|---|
| Numerical | Qp balls, valuations, digits, inverses, precision propagation, simple-root Hensel lifting |
| Local fields | Certified unramified and Eisenstein quadratic extensions, Q3(i) and Q2(i), ideal-ball precision, exact residue fields and norm/trace; lifted Frobenius for unramified fields |
| Analytic | exp, log1p, sin, cos and principal-unit powers over Qp and certified quadratic fields; explicit log(p)/log(pi) branches; Qp square roots including p=2 |
| Special functions | Multiplicative theta3, theta4 constants, exact Bernoulli/zeta values, numerical even KL zeta branches on Zp away from s=1 |
| Symbolic | Rational polynomials, Euclidean division/gcd, formal series, composition/calculus, optional SymPy function nodes/evaluation |
| Transforms | Exact atomic and bounded rational Amice measures, moments/polynomial integration, unit restriction, certified branched Mellin, local Laplace, Amice evaluation and additive convolution |
| Taylor | Exact polynomial shifts and rational jets; local p-adic exp/log/log1p/sin/cos and polynomial expansions with whole-tail enclosures |
| Calculus | Analytic series over Qp and certified quadratic fields; exact/ball coefficients, full-tail bounds, derivatives/primitives, multiplicative reciprocals and centered composition |
| Volkenborn integration | Exact polynomial/Bernoulli moments and certified analytic-series integrals, explicit affine charts and extension-valued coefficients |
| Local ODEs | Scalar and coupled linear equations with noncommuting analytic coefficient matrices over Qp/quadratic fields, vector/matrix initial data and whole-solution tail certificates |
| Harmonic analysis | Exact additive characters, real Haar integration, finite Fourier transforms, compactly supported test functions, inversion, Parseval and convolution |
| Wavelets and operators | Exact finite Kozyrev bases and reconstruction, whole-Qp integer-order Vladimirov images with analytic tails, conservative ball generators and numerical diffusion |
| Algebraic bridges | Q and Q[i] to independent completions; chosen split, unramified and ramified local-field images; rational quaternions and biquaternions |
| Geometric algebra | Diagonal Clifford algebras, grades, geometric/wedge products, involutions, explicit Hamilton-even-Clifford isomorphism |
| Linear algebra | Exact Q/Q[i] and local-field matrices; division-free determinants, certified square solves/inverses and condition exponents |

## Unramified quadratic fields

```python
from padic_tools import UnramifiedQuadraticField

L = UnramifiedQuadraticField(Qp(3, 12), [1, 0, 1], name="i")
x = L(1, 2)
assert (L.generator**2).is_congruent(-1)
assert (x*x.inverse()).is_congruent(1)
assert x.norm().contains(5) and x.trace().contains(2)
assert L.residue_field().order == 9
assert x.frobenius().residue() == x.residue()**3
```

The constructor checks monicity, p-integral exact coefficients and
irreducibility modulo p. `[1, 0, 1]` defines Q3(i); it is rejected over Q5
(split) and Q2 (ramified). `[1, 1, 1]` instead defines the unramified
quadratic extension of Q2. Each element is a whole ideal ball with common
coordinate precision, including negative valuations. `norm()` is algebraic;
`absolute_value()` is the p-adic magnitude. See [the extension chapter](docs/extensions.rst)
for error bounds, exact residues and scope. Elementary analytic functions
and dedicated field adapters are included; higher degrees remain future work.

## Ramified Q2(i) and Eisenstein fields

```python
from padic_tools import Q2i, EisensteinQuadraticField

L2i = Q2i(precision=12)
pi = L2i.uniformizer()  # pi = 1+i, with pi^2-2*pi+2 = 0
i = pi-1
assert (i*i).is_congruent(-1)
assert pi.valuation() == 1 and L2i(2).valuation() == 2
assert pi.base_valuation() == Fraction(1, 2)
assert pi.norm().contains(2) and pi.trace().contains(2)
assert L2i.residue_field().order == 2

ball = L2i(1, 1, absprec=5)  # Precision is in powers of pi
assert tuple(c.precision for c in ball.coefficients) == (3, 2)
R = EisensteinQuadraticField(Qp(3, 12), [-3, 0, 1])
assert (R.generator**2).is_congruent(3)
```

Coordinates use the basis `(1, pi)`. Odd and negative uniformizer precisions
give staggered base-coordinate precisions. Norm and trace return Qp balls
in base p units, with bounds that account for the different, including wild
dyadic ramification. `absolute_value_squared()` gives the exact rational
square of the base-compatible magnitude. See [the ramified chapter](docs/ramified.rst)
for the derivation and explicit Gaussian source conversion.

## Chosen Gaussian and quaternion local places

```python
from padic_tools import Gaussian, GaussianBridge, GaussianLocalPlace, Biquaternions
from padic_tools import QuaternionBridge, CliffordAlgebra

place = GaussianLocalPlace(Qp(2, 12), root_sign=1)  # i -> pi-1 in Q2(i)
source = GaussianBridge(Fraction(1, 3), 2)
image = source.local(place, absprec=19)  # Absolute powers of pi
assert image.precision == 19 and image.contains(Fraction(-5, 3), 2)
assert image.norm().contains(source.value.norm())

q = Biquaternions(Gaussian(1, 2), Gaussian(3, -1), 0, 1)
local_q = QuaternionBridge(q).local(place)
C = CliffordAlgebra((1, 1, 1), place.field)
assert q.to_clifford().map_coefficients(C, place).is_congruent(local_q.to_clifford())
```

The same API selects Qp at split primes and an unramified quadratic field
at primes congruent to 3 modulo 4. `root_sign` explicitly chooses the
generator image. See [chosen local places](docs/local_bridges.rst) for
guard digits, conjugations, algebraic norms and the Clifford commuting square.
The field images support the elementary functions described below.

## Analysis over quadratic local fields

```python
from padic_tools import exp, log, log1p, power, taylor

L = Q2i(precision=3)
x = L(4, 4, absprec=21)
assert log(exp(x)).is_congruent(x, precision=21)
assert (power(1+x, Fraction(1, 3))**3).is_congruent(1+x)
pi = L(0, 1, absprec=19)
assert log(pi, p_branch=1).contains(Fraction(1, 2))
point = x+L(0, 16, absprec=23)
assert taylor(exp, x, order=24)(point).is_congruent(exp(point), 21)
```

Domains use `v_pi(pi)=1`, `v_pi(p)=e`: exp/sin/cos require
`v_pi(x)>e/(p-1)`, log1p requires `v_pi(x)>=1`, and log accepts balls
excluding zero. Logarithm defaults to `log(p)=0`; choose at most one of
`p_branch` and `pi_branch`. Finite branch uncertainty and all later
denominator spikes are included in the output ball. See
[quadratic analysis](docs/extension_analytic.rst) for the full proof and
precision contracts. Extension sqrt/Teichmueller lifts and theta/zeta
evaluation in extension fields remain future work.

## Differential and integral calculus

```python
from padic_tools import AnalyticSeries, Q2i, exp

L = Q2i(precision=3)
s = AnalyticSeries.from_function(exp, L, order=48, precision=21)
x = L(4, 4, absprec=25)
assert s.derivative(2)(x, precision=21).is_congruent(exp(x), 21)
assert s.antiderivative()(x, precision=21).is_congruent(exp(x)-1, 21)
```

The series retains bounds for every omitted coefficient. Derivatives and
primitives transform those bounds, including all later p-power denominator
spikes. Exact coefficients remain exact; finite balls preserve their
uncertainty. A custom prefix needs an explicit coefficient-tail premise.
See [local calculus](docs/calculus.rst) and run
`python examples/local_calculus.py`. Local primitives, Haar integration and
integration against p-adic measures have separate APIs and meanings.

## Analytic reciprocals and composition

```python
from fractions import Fraction
from padic_tools import AnalyticSeries, Qp

K = Qp(3, 3)
inverse = AnalyticSeries(K, [1, -1], tail_bound=None).reciprocal(order=48)
inner = AnalyticSeries(K, [0, 1, 1], tail_bound=None)
assert inverse.compose(inner)(9, precision=21).contains(Fraction(-1, 89))
```

Reciprocals certify a nonvanishing disc for every allowed denominator.
Composition checks exact center alignment and the whole disc image, with
a bound on all omitted output coefficients. See [analytic operations](docs/series_ops.rst)
and run `python examples/analytic_operations.py`.

## Volkenborn integration

```python
from fractions import Fraction
from padic_tools import Polynomial, volkenborn_integral

x = Polynomial([0, 1])
assert volkenborn_integral(x) == Fraction(-1, 2)
assert volkenborn_integral(x*x) == Fraction(1, 6)
```

The integrator also accepts `AnalyticSeries` with explicit affine parameters
and a proof-backed whole-tail bound. See [Volkenborn integration](docs/volkenborn.rst)
for domains, finite precision and the distinction from Haar integration.

## Local differential equations

```python
from padic_tools import AnalyticSeries, Qp, exp, solve_linear_ode

K = Qp(3, 3)
a = AnalyticSeries.from_polynomial(K, [-1])
y = solve_linear_ode(a, 3, forcing=2, order=64, precision=21)
assert y(0, precision=60).contains(3)
assert y(9, precision=21).is_congruent(2 + exp(K(-9, absprec=45)), 21)
```

This solves y'=-y+2, y(0)=3. Variable coefficients and uncertain inputs use
the same API. The returned `AnalyticSeries` includes a sufficient convergence
disc and a bound for every omitted solution coefficient. See
[local linear ODEs](docs/ode.rst) for the proof and scope, and run
`python examples/linear_ode.py`.

For coupled equations and fundamental matrices:

```python
from padic_tools import Matrix, exact, solve_linear_system

phi = solve_linear_system([[0, 1], [-1, 0]], Matrix.identity(exact, 2),
                          field=K, order=64)
assert phi(0, precision=60).contains([[1, 0], [0, 1]])
```

See [linear differential systems](docs/ode_system.rst) for variable
coefficients, forcing, the matrix majorant proof and examples.

## Matrices and certified linear systems

```python
from padic_tools import Matrix, exact

A = Matrix(exact, [[1, 2], [3, 7]])
assert A.solve([1, 0]) == (7, -3)
assert A @ A.inverse() == Matrix.identity(exact, 2)

D = Matrix(Qp(3, 12), [[9, 0], [0, 243]], absprec=21)
report = D.solve_with_info([1, 1])
assert report.pivot_valuations == (2, 5)
assert report.condition_exponent == 3
assert report.precisions == (14, 11)
assert D.inverse().precisions == ((17, 14), (14, 11))
```

Local matrices represent whole families of entry balls. Complete pivoting
checks valuation lower bounds, preserves uncertain zero entries and reports
insufficient precision explicitly. Exact sources can be reevaluated through
the Gaussian local places. See [linear algebra](docs/linalg.rst) for the
enclosure proof, condition normalization, scope and independent checks.
Run `python examples/linear_systems.py` for a standalone example.

## Transforms and Taylor expansions

```python
from padic_tools import AtomicMeasure, RationalAmiceMeasure, taylor

mu = AtomicMeasure(K, [(1, 1), (2, 2)])
assert mu.mellin(6, branch=2).is_congruent(mu.moment(6))
assert mu.amice(exp(K(-5))-1).is_congruent(mu.laplace(K(5)))

jet = taylor(exp, K(5), order=30)
assert jet(K(10)).is_congruent(exp(K(10)), precision=20)

rational = RationalAmiceMeasure(K, [1], [2, 1])  # Complete A(T)=1/(2+T)
assert rational.moment(3) == Fraction(1, 8)
assert rational.laplace(K(5)).is_congruent(1/(1+exp(K(-5))))

units = rational.restrict_units()  # Exact rational transform of the restriction
assert units.moment(3) == Fraction(-31, 2)
assert units.mellin(3, branch=3).contains(units.moment(3))
assert units.mellin(Fraction(1, 2), branch=3).precision == 20
```

The measure and kernel convention are explicit. Mellin support must lie in
Zp units; Laplace checks the local p-adic exponential domain. Atomic measures
are finite exact data. Rational Amice measures additionally describe complete
bounded moment sequences with a checked coefficient bound. Their Mellin
evaluator proves an infinite Mahler-tail bound, including the dyadic case;
`precision` requests an absolute output order and `max_terms` limits work.
General bounded-measure oracles and transform inverses are future extensions.
Formal series order and p-adic numerical precision are
tracked separately; a numerical Taylor expansion includes its omitted tail.

## Exact harmonic analysis

```python
from padic_tools import PadicTestFunction, additive_phase

K3 = Qp(3)
assert additive_phase(Fraction(1, 6), K3) == Fraction(2, 3)
f = PadicTestFunction(K3, [1, 2, 4])  # Supported on Z3, constant modulo 3
assert f.integral().as_fraction() == Fraction(7, 3)
spectrum = f.fourier()
assert spectrum.fourier(inverse=True) == f
assert spectrum.norm_squared() == f.norm_squared()
assert PadicTestFunction.ball(K3).fourier() == PadicTestFunction.ball(K3)
```

The coefficients are exact rational combinations of complex p-power roots of
unity; `.complex(dps=...)` provides an optional numerical view. Haar weights
are real, with volume(Zp)=1. The transform retains support/constancy metadata
for the dual Qp lattice. Raw `finite_fourier` arrays have an explicit 1/N
forward normalization. Native work limits are documented in the harmonic
chapter; accelerated backends and general complex samples remain future work.

## Wavelets and conservative diffusion

```python
from padic_tools import KozyrevBasis, kozyrev_wavelet, vladimirov, BallVladimirov

basis = KozyrevBasis(K3, constancy_exponent=2)
signal = PadicTestFunction(K3, range(9), constancy_exponent=2)
assert basis.reconstruct(basis.decompose(signal)) == signal
wave = kozyrev_wavelet(K3, scale=1)
assert vladimirov(wave, order=2).as_test_function() == wave * 81

operator = BallVladimirov(K3, constancy_exponent=2)
initial = (PadicTestFunction.ball(K3, exponent=2)*9).on_lattice(
    support_exponent=0, constancy_exponent=2)
state = operator.evolve(initial, Fraction(1, 10), dps=60)  # complex extra
assert abs(state.integral()-1) < 10**-50
```

Wavelets retain their exact squared norms. `vladimirov` includes the entire
exterior tail for nonzero-mass input on Qp. `BallVladimirov` instead uses the
regional conservative boundary on its declared ball. Its exact spectrum and
kernel matrix determine mass-preserving heat evolution. Numerical outputs
use isolated mpmath working precision. Run `python examples/padic_diffusion.py`
for a reproducible positivity, mass, energy and grid/boundary comparison.

## The bridge is included

```python
from padic_tools import Gaussian, GaussianBridge, Hamilton, QuaternionBridge

g = GaussianBridge(Gaussian(Fraction(1, 3), 2))
complex_image = g.complex(dps=60)
local_algebra_image = g.padic_algebra(Qp(3, 20))
chosen_split_image = g.split_padic(Qp(5, 20), root_sign=1)

q = Hamilton(1, 2, 3, 4)
q5 = QuaternionBridge(q).padic(Qp(5, 20))
assert (q*q).to_clifford() == q.to_clifford()*q.to_clifford()
```

The common source is an exact algebra over **Q**; its complex and p-adic
images are calculated separately. A finite p-adic residue is not a complex
approximation. At split primes the quadratic quotient has zero divisors.
The complex Riemann zeta and a KL branch are distinct functions with a
specified special-value interpolation relation. See the bridge chapter for
proofs, conventions and examples.

## Documentation and validation

```bash
python -m pip install '.[dev]'
python -m pytest
python -m sphinx -W --keep-going -b html docs docs/_build/html
python -m sphinx -W --keep-going -b doctest docs docs/_build/doctest
```

Open **docs/_build/html/index.html**. The guide includes a Czech introduction,
an API reference, mathematical derivations, executable examples, a feature
roadmap and links to primary sources. CI uploads the generated HTML manual as
the `html-documentation` artifact; it does not publish a website.

For sharing or opening a single file, build the portable manual instead:

```bash
python scripts/build_manual.py
```

Open **docs/_build/padic-tools-manual.html**. This file includes every chapter,
its styles and pre-rendered mathematical formulas. Chapter and language links
stay inside the file, so it works offline without companion files or
JavaScript. Building this portable edition additionally requires Node.js;
reading it does not. CI uploads it as `portable-html-manual`.

The regular multi-page HTML edition must be served as a complete directory
or opened from the full downloaded `html-documentation` archive. Sharing
only its `index.html` loses the other pages and assets. Both editions provide
a separate **English · Čeština** language navigation bar; the Czech guide
is kept out of the English chapter menu.

Start with [the Czech guide](docs/cesky.rst),
[the numerical tutorial](docs/tutorial.rst),
[the special-function definitions](docs/functions.rst), and
[the complex/quaternion/Clifford bridge](docs/bridges.rst).

The 2022 `main.py` helpers remain as an installed-source compatibility demo.
`padic()` now returns a ball; new code should import `padic_tools`.
The unfinished float-based arithmetic has been replaced with exact rational
conversion. The project does not currently grant an open-source license;
licensing remains the repository owner's decision.

