from fractions import Fraction as F
from math import ceil, comb

import pytest
import sympy as sp

from padic_tools import (
    AnalyticSeries as S, DomainError, EisensteinQuadraticField, PrecisionError,
    Q2i, Qp, UnramifiedQuadraticField, exp, log, log1p,
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


def contains(K, coefficient, source):
    if hasattr(coefficient, 'contains'):
        assert coefficient.contains(source)
    else:
        assert matrix(K, coefficient) == matrix(K, source)


@pytest.mark.parametrize('K', FIELDS)
def test_reciprocal_geometric_coefficients_and_entire_tail(K):
    f = S(K, [pi_power(K, 2), image(-matrix(K, pi_power(K, 1)))], tail_bound=None)
    r = f.reciprocal(order=24, precision=21)
    assert r.radius == 2 and r.tail_bound == (-2, -1)
    for n in range(40):
        contains(K, r[n], pi_power(K, -2-n))
    x = pi_power(K, r.radius+2)
    oracle = (matrix(K, pi_power(K, 2))-matrix(K, pi_power(K, 1))*matrix(K, x)).inv()
    assert r(x, precision=21).contains(image(oracle))
    assert r(0, precision=80).contains(pi_power(K, -2))


@pytest.mark.parametrize('K', FIELDS)
def test_quadratic_denominator_against_closed_binomial_matrix_formula(K):
    G = matrix(K, pi_power(K, -1))
    one = sp.eye(G.rows)
    f = S(K, [1, image(G), 1], tail_bound=None)
    r = f.reciprocal(order=24)
    for n in range(36):
        # [h^n] 1/(1+G*h+h^2), from a geometric/binomial expansion.
        oracle = sum(((-1)**(n-k)*comb(n-k, k)*G**(n-2*k)
                      for k in range(n//2+1)), sp.zeros(G.rows))
        contains(K, r[n], image(oracle))
    X = matrix(K, pi_power(K, r.radius+3))
    assert r(image(X), precision=21).contains(image((one+G*X+X*X).inv()))


@pytest.mark.parametrize('K', FIELDS)
@pytest.mark.parametrize('bound,radius', [((2, F(-1, 2)), 1), ((-3, F(1, 2)), 0)])
def test_reciprocal_whole_tail_families_do_not_pad_unknowns_with_zero(K, bound, radius):
    A, B = bound
    f = S(K, [1], radius=radius, tail_bound=bound, precision=30)
    r = f.reciprocal(order=20)
    X = matrix(K, pi_power(K, r.radius+3))
    one = sp.eye(X.rows)
    value = r(image(X))
    assert value.precision < 21
    with pytest.raises(PrecisionError):
        r(image(X), precision=21)
    for n in sorted({1, 2, K.p, K.p**2}):
        C = matrix(K, pi_power(K, ceil(A+B*n)))
        assert value.contains(image((one+C*X**n).inv()))
        for j in (1, 2, 3):
            contains(K, r[n*j], image((-C)**j))


@pytest.mark.parametrize('K', [FIELDS[1], FIELDS[3], FIELDS[5]])
def test_reciprocal_constant_coefficient_and_argument_uncertainty(K):
    source = [pi_power(K, 1), 1, -1]
    f = S(K, [K(v, absprec=12+n) for n, v in enumerate(source)], tail_bound=None, precision=30)
    r = f.reciprocal(order=32)
    X = matrix(K, pi_power(K, r.radius+3))
    value = r(K(image(X), absprec=15))
    assert K.precision < value.precision <= 10  # Inverting a_0 costs 2*v(a_0).
    with pytest.raises(PrecisionError):
        r(0, precision=11)
    for j in (-1, 0, 1):
        cs = [matrix(K, v)+(j+n)*matrix(K, pi_power(K, 12+n)) for n, v in enumerate(source)]
        x = X+j*matrix(K, pi_power(K, 15))
        assert value.contains(image((cs[0]+cs[1]*x+cs[2]*x*x).inv()))


@pytest.mark.parametrize('K', FIELDS)
def test_full_polynomial_composition_and_truncation_on_large_disc(K):
    C, D = pi_power(K, -1), pi_power(K, -2)
    G = matrix(K, pi_power(K, 1))
    one = sp.eye(G.rows)
    f = S(K, [1, image(G), -1, 1], center=C, tail_bound=None)
    g = S(K, [C, image(G), 1], center=D, tail_bound=None)
    result = f.compose(g)
    assert result.order == 7 and result.tail_bound is None and result.radius is None
    assert result.center == g.center
    t = sp.Symbol('t')
    T = G*t+one*t*t
    oracle = one+G*T-T*T+T**3
    for n in range(7):
        contains(K, result[n], image(oracle.applyfunc(lambda v:sp.expand(v).coeff(t, n))))
    X = matrix(K, pi_power(K, -3))
    point = image(matrix(K, D)+X)
    exact_value = oracle.subs(t, X[0, 0]) if isinstance(K, PadicField) else (
        one+G*(G*X+X*X)-(G*X+X*X)**2+(G*X+X*X)**3)
    assert result(point, precision=60).contains(image(exact_value))
    short = f.compose(g, order=3, radius=-3)
    assert short.tail_bound is not None and short.radius == -3
    assert short(point).contains(image(exact_value))


@pytest.mark.parametrize('K', FIELDS)
def test_composition_of_infinite_families_against_closed_rational_coefficients(K):
    alpha, beta, gamma = (matrix(K, pi_power(K, n)) for n in (-1, 2, -1))
    C, D = pi_power(K, -1), pi_power(K, -2)
    one = sp.eye(alpha.rows)
    f = S(K, [image(alpha**n) for n in range(6)], center=C, radius=2, tail_bound=(0, -1), precision=30)
    g = S(K, [C]+[image(beta*gamma**(n-1)) for n in range(1, 4)], center=D,
          radius=2, tail_bound=(3, -1), precision=30)
    result = f.compose(g, order=24)
    assert result.radius == 2
    for n in range(40):
        oracle = one if n == 0 else alpha*beta*(gamma+alpha*beta)**(n-1)
        contains(K, result[n], image(oracle))
    X = matrix(K, pi_power(K, result.radius+4))
    T = beta*X*(one-gamma*X).inv()
    point = image(matrix(K, D)+X)
    assert result(point, precision=21).contains(image((one-alpha*T).inv()))


@pytest.mark.parametrize('K', FIELDS)
def test_functional_identities_and_center_alignment(K):
    exponential = S.from_function(exp, K, order=32, precision=21)
    logarithm = S.from_function(log1p, K, order=32, precision=21)
    f = exponential.compose(logarithm)
    assert f.tail_bound is not None  # Bounds still allow omitted perturbations.
    for n in range(10):
        contains(K, f[n], int(n in (0, 1)))
    point = pi_power(K, f.radius+3)
    assert f(point, precision=21).contains(image(matrix(K, 1)+matrix(K, point)))
    outer_log = S.from_function(log, K, center=1, order=32, precision=21)
    g = outer_log.compose(exponential)
    assert g.center == exponential.center
    for n in range(10):
        contains(K, g[n], int(n == 1))
    point = pi_power(K, g.radius+3)
    assert g(point, precision=21).contains(point)


@pytest.mark.parametrize('K', FIELDS)
def test_missing_outer_constant_and_constant_inner_are_families(K):
    f = S(K, [], radius=0, tail_bound=(-2, 1), precision=30)
    g = S(K, [0, pi_power(K, 2)], tail_bound=None)
    result = f.compose(g, order=20)
    assert result.radius == -2  # Actual image lies on the outer disc boundary.
    assert result(pi_power(K, -2)).precision <= -2
    with pytest.raises(PrecisionError):
        result(pi_power(K, -2), precision=-1)
    constant = S(K, [0], radius=1, tail_bound=None)
    collapsed = f.compose(constant)
    assert collapsed.tail_bound is None and collapsed.radius == 1
    assert collapsed(0).contains(pi_power(K, -2))
    with pytest.raises(PrecisionError):
        collapsed(0, precision=-1)


@pytest.mark.parametrize('K', [FIELDS[1], FIELDS[3], FIELDS[5]])
def test_polynomial_composition_preserves_all_finite_ball_uncertainties(K):
    C, D = pi_power(K, 1), pi_power(K, -1)
    fs, gs = [1, pi_power(K, -1), 1], [C, pi_power(K, 1), 1]
    f = S(K, [K(v, absprec=11+n) for n, v in enumerate(fs)], center=C, tail_bound=None, precision=30)
    g = S(K, [C]+[K(v, absprec=12+n) for n, v in enumerate(gs[1:], 1)], center=D, tail_bound=None, precision=30)
    result = f.compose(g)
    assert result.tail_bound is None
    X = matrix(K, pi_power(K, 3))
    point = image(matrix(K, D)+X)
    value = result(K(point, absprec=14))
    assert K.precision < value.precision <= 11
    for j in (-1, 0, 1):
        a = [matrix(K, v)+(j+n)*matrix(K, pi_power(K, 11+n)) for n, v in enumerate(fs)]
        b = [matrix(K, v)+(j+n)*matrix(K, pi_power(K, 12+n)) for n, v in enumerate(gs) if n]
        x = X+j*matrix(K, pi_power(K, 14))
        T = b[0]*x+b[1]*x*x
        assert value.contains(image(a[0]+a[1]*T+a[2]*T*T))


@pytest.mark.parametrize('K', FIELDS)
def test_constant_reciprocal_and_whole_disc_validation(K):
    a = pi_power(K, 1)
    f = S(K, [a, 0, 0], tail_bound=None)
    r = f.reciprocal()
    assert r.tail_bound is None and r.radius is None
    assert r(pi_power(K, -10), precision=80).contains(pi_power(K, -1))
    family = S(K, [K(a, absprec=11)], tail_bound=None).reciprocal()
    assert family(K.O(-10)).precision == K.precision  # Default targets use uniformizer units.
    with pytest.raises(PrecisionError):
        family(0, precision=10)
    geometric = S(K, [1, -1], tail_bound=None).reciprocal(order=32)
    with pytest.raises(DomainError):
        geometric(K.O(0))
    assert geometric(K.O(1)).contains(1)
    with pytest.raises(PrecisionError):
        geometric(K.O(1), precision=2)


def test_invalid_arguments_and_unsupported_alignment():
    K = Qp(3)
    f = S(K, [1, -1], tail_bound=None)
    g = S(K, [0, 1], tail_bound=None)
    for method in (f.reciprocal, lambda **kw:f.compose(g, **kw)):
        for kwargs in ({'order': True}, {'radius': 0.5}, {'precision': True}):
            with pytest.raises(TypeError):
                method(**kwargs)
        for order in (0, -1):
            with pytest.raises(ValueError):
                method(order=order)
    for source in (S(K, [], tail_bound=None), S(K, [0, 1], tail_bound=None)):
        with pytest.raises(ZeroDivisionError):
            source.reciprocal()
    for source in (S(K, [K.O(3)], tail_bound=None), S(K, [], radius=1, tail_bound=(0, 0))):
        with pytest.raises(PrecisionError):
            source.reciprocal()
    with pytest.raises(DomainError):
        f.reciprocal(radius=0)
    for inner in (3, S(Qp(5), [0, 1], tail_bound=None), S(K, [K(0), 1], tail_bound=None)):
        with pytest.raises(TypeError):
            f.compose(inner)
    with pytest.raises(DomainError):
        f.compose(S(K, [1, 1], tail_bound=None))
    with pytest.raises(DomainError):
        f.restrict(2).compose(g, radius=1)
