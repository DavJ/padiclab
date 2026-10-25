from fractions import Fraction as F
from math import ceil, factorial

import pytest
import sympy as sp

from padic_tools import (
    AnalyticSeries as S, DomainError, EisensteinQuadraticField, Polynomial,
    PowerSeries, PrecisionError, Q2i, Qp, UnramifiedQuadraticField, exp,
    solve_linear_ode as solve,
)
from padic_tools.core import PadicField


FIELDS = [Qp(2, 3), Qp(3, 3), Qp(5, 3), Q2i(3),
          UnramifiedQuadraticField(Qp(2, 3), [1, 1, 1]),
          UnramifiedQuadraticField(Qp(3, 3), [1, 0, 1]),
          EisensteinQuadraticField(Qp(3, 3), [-3, 0, 1]),
          EisensteinQuadraticField(Qp(3, 3), [-3, 3, 1]),
          EisensteinQuadraticField(Qp(2, 3), [2, 2, 1])]


def matrix(K, value):
    if isinstance(K, PadicField):
        return sp.Matrix([[sp.Rational(value)]])
    a, b = value if isinstance(value, (tuple, list)) else (value, 0)
    C, B, _ = K.polynomial.coefficients
    return sp.Matrix([[a, -C*b], [b, a-B*b]])


def image(M):
    return F(M[0, 0]) if M.rows == 1 else (F(M[0, 0]), F(M[1, 0]))


def pi_power(K, n):
    if isinstance(K, PadicField) or K.ramification_index == 1:
        return F(K.p)**n
    t = sp.Symbol('t')
    f = sum(sp.Rational(c)*t**j for j, c in enumerate(K.polynomial.coefficients))
    q = sp.rem(t**n, f, t) if n >= 0 else sp.invert(t**(-n), f, t)
    return tuple(F(q.coeff(t, j)) for j in range(2))


def contains_coefficient(K, coefficient, source):
    if hasattr(coefficient, 'contains'):
        assert coefficient.contains(source)
    else:
        assert matrix(K, coefficient) == matrix(K, source)


@pytest.mark.parametrize('K', FIELDS)
def test_constant_equation_against_exact_multiplication_matrix_exponential(K):
    A = matrix(K, pi_power(K, -1))
    G = matrix(K, 3 if isinstance(K, PadicField) else (1, 1))
    Y = matrix(K, 2 if isinstance(K, PadicField) else (2, -1))
    center = pi_power(K, -1)
    a = S(K, [image(A)], tail_bound=None, center=center)
    y = solve(a, image(Y), forcing=image(-A*G), order=64, precision=21)
    assert y.order == 64 and y.center == a.center
    for n in range(80):
        expected = Y if n == 0 else A**n*(Y-G)/factorial(n)
        contains_coefficient(K, y[n], image(expected))
    h = matrix(K, pi_power(K, y.radius+2))
    oracle = K(image(G), absprec=45)+K(image(Y-G), absprec=45)*exp(K(image(A*h), absprec=45))
    point = image(matrix(K, center)+h)
    assert y(point, precision=21).is_congruent(oracle, 21)
    assert y(center, precision=80).contains(image(Y))


@pytest.fixture(scope='module')
def integrating_factor_coefficients():
    t = sp.Symbol('t')
    primitive = 2*t-sp.Rational(3, 2)*t*t+sp.Rational(4, 3)*t**3
    inverse = sp.exp(-primitive).series(t, 0, 12).removeO()
    integral = sp.integrate(inverse*(3-4*t*t), t)
    factor = sp.exp(primitive).series(t, 0, 12).removeO()
    solution = sp.expand(factor*(sp.Rational(3, 5)+integral))
    return [F(solution.coeff(t, n)) for n in range(12)]


@pytest.mark.parametrize('K', FIELDS)
def test_variable_equation_against_symbolic_integrating_factor(K, integrating_factor_coefficients):
    center = pi_power(K, 1)
    a = S(K, [2, -3, 4], center=center, tail_bound=None)
    b = S(K, [3, 0, -4], center=center, tail_bound=None)
    y = solve(a, F(3, 5), forcing=b, order=12)
    for coefficient, oracle in zip(y.coefficients, integrating_factor_coefficients):
        contains_coefficient(K, coefficient, oracle)
    residual = y.derivative()-a*y-b
    for coefficient in residual.coefficients:
        contains_coefficient(K, coefficient, 0)


@pytest.mark.parametrize('K', FIELDS)
@pytest.mark.parametrize('mode', ['multiplier', 'forcing'])
@pytest.mark.parametrize('bound', [(4, F(1)), (-3, F(-1, 2))])
def test_whole_input_tail_families_include_later_denominator_spikes(K, mode, bound):
    A, B = bound
    r = 0 if B > 0 else 2
    family = S(K, [0, 0], radius=r, tail_bound=bound, precision=30)
    zero = S(K, [0], radius=r, tail_bound=None, precision=30)
    y = solve(family if mode == 'multiplier' else zero,
              1 if mode == 'multiplier' else 0,
              forcing=family if mode == 'forcing' else 0, order=12)
    h = matrix(K, pi_power(K, y.radius+7))
    value = y(image(h), precision=21)
    # Each omitted input monomial is a full admissible analytic function.
    # Its primitive is exact, including n+1=p^k denominator spikes; the
    # multiplier-only solution is its exponential, not the solver recurrence.
    for n in sorted({2, 3, K.p-1, K.p**2-1, 2*K.p**2-1}):
        if n < 2:
            continue
        primitive = matrix(K, pi_power(K, ceil(A+B*n)))/(n+1)
        argument = primitive*h**(n+1)
        expected = exp(K(image(argument), absprec=45)) if mode == 'multiplier' else image(argument)
        assert value.is_congruent(expected, 21)
        for j in range(1, 4):
            index = j*(n+1)
            coefficient = primitive**j/factorial(j) if mode == 'multiplier' else (
                primitive if j == 1 else sp.zeros(primitive.rows))
            contains_coefficient(K, y[index], image(coefficient))


@pytest.mark.parametrize('K', FIELDS)
def test_short_nonzero_input_prefix_is_not_padded_with_exact_zeros(K):
    a = S(K, [1], radius=1, tail_bound=(0, 0), precision=30)
    y = solve(a, order=30)
    for coefficient in y.coefficients:
        contains_coefficient(K, coefficient, 1)
    h = matrix(K, pi_power(K, y.radius+2))
    one = sp.eye(h.rows)
    expected = image((one-h).inv())  # a=1/(1-h), y=1/(1-h).
    value = y(image(h))
    assert value.contains(expected)
    with pytest.raises(PrecisionError):
        y(image(h), precision=value.precision+1)


@pytest.mark.parametrize('K', FIELDS)
def test_finite_coefficient_initial_and_input_ball_families(K):
    A, B, Y = (matrix(K, pi_power(K, -1)), matrix(K, pi_power(K, 1)),
               matrix(K, 2 if isinstance(K, PadicField) else (2, 1)))
    a = S(K, [K(image(A), absprec=11)], tail_bound=None, precision=30)
    y = solve(a, K(image(Y), absprec=12), forcing=K(image(B), absprec=13), order=40)
    X = matrix(K, pi_power(K, y.radius+1))
    value = y(K(image(X), absprec=14))
    assert value.precision <= 12
    with pytest.raises(PrecisionError):
        y(K(image(X), absprec=14), precision=13)
    with pytest.raises(PrecisionError):
        y(0, precision=60)
    for j in (-1, 0, 1):
        aa = A+j*matrix(K, pi_power(K, 11))
        bb = B+(j+1)*matrix(K, pi_power(K, 13))
        yy = Y+(j-1)*matrix(K, pi_power(K, 12))
        xx = X+j*matrix(K, pi_power(K, 14))
        P, oracle, power = aa*xx, yy, sp.eye(A.rows)
        for n in range(1, 40):
            oracle += (power*P*yy+power*bb*xx)/factorial(n)
            power *= P
        # This fixed rational matrix exponential has omitted terms of
        # valuation >70 on these charts; the compared balls have <=12 digits.
        assert value.contains(image(oracle))


@pytest.mark.parametrize('K', FIELDS)
def test_polynomial_and_structural_zero_solutions_keep_zero_tail(K):
    center = pi_power(K, -1)
    zero = S(K, [], center=center, tail_bound=None)
    b = S.from_polynomial(K, Polynomial([1, -3, 2]), center=center)
    y = solve(zero, 3, forcing=b, order=6)
    primitive = b.antiderivative(3)
    assert y.tail_bound is None and y.radius is None
    point = pi_power(K, -3)
    assert y(point, precision=60).is_congruent(primitive(point, precision=60), 60)
    uncertain = S(K, [K.O(-2)], center=center, radius=2, tail_bound=(0, 0))
    zero_solution = solve(uncertain, 0, order=1)
    assert zero_solution.tail_bound is None and zero_solution.radius == 2
    assert zero_solution(center, precision=80).contains(0)
    # A zero-containing initial ball is not structural zero.
    family = solve(uncertain, K.O(5), order=1)
    assert family.tail_bound is not None
    with pytest.raises(PrecisionError):
        family(center, precision=6)
    short = solve(zero, 3, forcing=b, order=1)
    assert short.tail_bound is not None  # The uncomputed primitive is retained.


@pytest.mark.parametrize('K', FIELDS)
def test_larger_than_unit_discs_and_precision_limits(K):
    a = S.from_polynomial(K, [K.p**4])
    y = solve(a, order=80)
    assert y.radius < 0
    point = pi_power(K, y.radius+1)
    oracle = exp(K(image(matrix(K, K.p**4)*matrix(K, point)), absprec=45))
    assert y(point, precision=21).is_congruent(oracle, 21)
    assert y(K.O(y.radius+1)).contains(1)
    with pytest.raises(PrecisionError):
        y(K.O(y.radius+1), precision=100)
    with pytest.raises(DomainError):
        y(K.O(y.radius-1))
    with pytest.raises(DomainError):
        solve(a, radius=y.radius-1)
    restricted = solve(a.restrict(2), forcing=S(K, [0], radius=5, tail_bound=None))
    assert restricted.radius == 5
    with pytest.raises(DomainError):
        solve(a.restrict(2), radius=1)


def test_argument_validation_and_one_term_tail():
    K = Qp(3, 3)
    a = S.from_polynomial(K, [1])
    for bad in (1, Polynomial([1]), PowerSeries([1], 20), lambda x:x):
        with pytest.raises(TypeError):
            solve(bad)
    for kwargs in ({'order': True}, {'radius': 1.5}, {'precision': True},
                   {'initial': 0.5}, {'forcing': 0.5},
                   {'forcing': S.from_polynomial(Qp(5), [1])},
                   {'forcing': S.from_polynomial(K, [1], center=1)}):
        with pytest.raises(TypeError):
            solve(a, **kwargs)
    for order in (-1, 0):
        with pytest.raises(ValueError):
            solve(a, order=order)
    y = solve(a, order=1)
    assert y.tail_bound is not None
    assert y(0, precision=80).contains(1)
    with pytest.raises(PrecisionError):
        y(3, precision=21)
    with pytest.raises(DomainError):
        solve(a, radius=0)
