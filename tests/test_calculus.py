"""Independent polynomial/matrix oracles and adversarial coefficient families."""
from fractions import Fraction as F
from math import ceil

import pytest
import sympy as sp

from padic_tools import (
    AnalyticSeries as S, DomainError, EisensteinQuadraticField, Polynomial,
    PrecisionError, Q2i, Qp, UnramifiedQuadraticField, cos, exp, log, log1p, sin,
    valuation,
)
from padic_tools.core import PadicField


FIELDS = [Qp(2, 3), Qp(3, 3), Qp(5, 3),
          UnramifiedQuadraticField(Qp(2, 3), [1, 1, 1]),
          UnramifiedQuadraticField(Qp(3, 3), [1, 0, 1]),
          Q2i(3), EisensteinQuadraticField(Qp(2, 3), [-2, 0, 1]),
          EisensteinQuadraticField(Qp(3, 3), [-3, 3, 1]),
          EisensteinQuadraticField(Qp(5, 3), [5, 5, 1])]


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
    f = sum(sp.Rational(a)*t**j for j, a in enumerate(K.polynomial.coefficients))
    value = sp.rem(t**n, f, t) if n >= 0 else sp.invert(t**(-n), f, t)
    return tuple(F(value.coeff(t, j)) for j in range(2))


def polynomial_value(expression, point):
    t = sp.Symbol('t')
    M = sp.zeros(point.rows)
    for c in sp.Poly(expression, t).all_coeffs():
        M = M*point+c*sp.eye(point.rows)
    return M


@pytest.mark.parametrize('K', FIELDS)
def test_polynomial_calculus_against_sympy_at_nonintegral_points_and_centers(K):
    t = sp.Symbol('t')
    f = F(1, K.p**2)-7*t+F(2, K.p)*t**3+3*t**5
    cs = list(reversed(sp.Poly(f, t).all_coeffs()))
    cs = [F(c) for c in cs]
    # Nonrational center in quadratic fields, with no rounding at parent defaults.
    c = F(1, K.p) if isinstance(K, PadicField) else (F(1, K.p), F(1, K.p))
    series = S.from_polynomial(K, cs, center=c, precision=25)
    assert series.radius is None and series.tail_bound is None
    for x in (F(2), F(-2, K.p**2), F(1, K.p)):
        X = matrix(K, x)
        for d in (0, 1, 2, 5, 7):
            value = series.derivative(d)(x)
            assert value.precision == 25
            assert value.contains(image(polynomial_value(sp.diff(f, t, d), X)))
        primitive = sp.integrate(f, t)
        C = matrix(K, c)
        expected = polynomial_value(primitive, X)-polynomial_value(primitive, C)+7*sp.eye(X.rows)
        assert series.antiderivative(7)(x).contains(image(expected))
        assert series.tail_precision(x) is None
    assert S.from_polynomial(K, Polynomial(), center=c)(2, precision=70).contains(0)


@pytest.mark.parametrize('K', FIELDS)
def test_geometric_series_derivatives_against_independent_rational_matrices(K):
    r = S(K, [1]*35, radius=1, tail_bound=(0, 0), precision=22)
    x = pi_power(K, 2)
    X = matrix(K, x)
    I = sp.eye(X.rows)
    for d in (0, 1, 2, 4):
        # SymPy differentiates the rational function, rather than our prefix.
        t = sp.Symbol('t')
        expression = sp.diff(1/(1-t), t, d)
        numerator, denominator = sp.fraction(expression)
        expected = polynomial_value(numerator, X)*polynomial_value(denominator, X).inv()
        assert r.derivative(d)(x).contains(image(expected))
    # A separately fixed 100-term log polynomial is amply beyond 22 digits.
    expected = sum((X**n/sp.Integer(n) for n in range(1, 100)), sp.zeros(X.rows))
    assert r.antiderivative()(x).contains(image(expected))


@pytest.mark.parametrize('K', FIELDS)
@pytest.mark.parametrize('function', [exp, sin, cos, log1p, log])
def test_function_factories_at_standard_and_translated_centers(K, function):
    for translated in (False, True):
        c = (1+K.p**2 if function is log else K.p**2) if translated else (1 if function is log else 0)
        series = S.from_function(function, K, center=c, order=96, precision=21)
        step = K.p**3 if isinstance(K, PadicField) else (0, K.p**3)
        point = image(matrix(K, c)+matrix(K, step))
        x, a = K(point, absprec=28), K(c, absprec=28)
        assert series(point, precision=21).is_congruent(function(x), 21)
        assert series(K(point, absprec=22)).is_congruent(function(x))
        derivative = (exp(x) if function is exp else cos(x) if function is sin else
                      -sin(x) if function is cos else 1/(1+x) if function is log1p else 1/x)
        assert series.derivative()(point, precision=21).is_congruent(derivative, 21)
        primitive = (exp(x)-exp(a) if function is exp else cos(a)-cos(x) if function is sin else
                     sin(x)-sin(a) if function is cos else
                     (1+x)*log1p(x)-x-((1+a)*log1p(a)-a) if function is log1p else
                     x*log(x)-x-(a*log(a)-a))
        assert series.antiderivative()(point).is_congruent(primitive)


@pytest.mark.parametrize('K', FIELDS)
@pytest.mark.parametrize('radius,slope', [(1, F(-1, 2)), (-1, F(2)), (0, F(1, 2))])
def test_later_denominator_spikes_in_all_admissible_omitted_monomials(K, radius, slope):
    p = K.p
    r = S(K, [0]*3, radius=radius, tail_bound=(-3, slope), precision=30)
    x = pi_power(K, radius)
    X = matrix(K, x)
    # Each monomial is a different complete function allowed by the same bound.
    indices = sorted(set(range(3, 11)) | {p**j-1 for j in (1, 2, 3) if p**j-1 >= 3})
    for n in indices:
        exponent = ceil(-3+slope*n)
        coefficient = matrix(K, pi_power(K, exponent))
        assert r(x).contains(image(coefficient*X**n))
        for d in (1, 2, 5):
            t = sp.Symbol('t')
            expected = coefficient*polynomial_value(sp.diff(t**n, t, d), X)
            assert r.derivative(d)(x).contains(image(expected))
        expected = coefficient*X**(n+1)/sp.Integer(n+1)
        assert r.antiderivative()(x).contains(image(expected))
        expected = coefficient*X**(n+2)/(sp.Integer(n+1)*(n+2))
        assert r.antiderivative().antiderivative()(x).contains(image(expected))
    # The primitive stays on the same disc, despite its weaker tail slope.
    assert r.antiderivative().radius == radius


@pytest.mark.parametrize('K', FIELDS)
def test_finite_coefficient_precision_division_and_exact_center_values(K):
    e = 1 if isinstance(K, PadicField) else K.ramification_index
    for n in (K.p, K.p**2):
        cs = [0]*(n-1)+[K(1, absprec=11)]
        r = S(K, cs, tail_bound=None, precision=30)
        primitive = r.antiderivative(K(3, absprec=7))
        assert primitive[n].precision == 11-e*valuation(n, K.p)
        assert primitive.derivative()[n-1].precision == 11
        assert primitive(0).precision == 7
        with pytest.raises(PrecisionError):
            primitive(0, precision=8)
    source = S.from_polynomial(K, [0, K(1, absprec=11), 0, 3], precision=30)
    # Every nonconstant term is exactly zero at the exact center, even with
    # uncertain coefficients. A zero-containing argument ball is different.
    assert source(0, precision=60).precision == 60
    assert source(K.O(4)).precision == 4
    square = S.from_polynomial(K, [0, 0, 1], precision=30)
    x = K(pi_power(K, 1), absprec=7)
    assert square(x, precision=8).precision == 8
    with pytest.raises(PrecisionError):
        square(x, precision=9)


@pytest.mark.parametrize('K', FIELDS)
def test_whole_input_and_coefficient_ball_perturbations_against_exact_matrices(K):
    sources, precisions = [1, F(1, K.p), -3], [12, 14, 17]
    cs = [K(c, absprec=N) for c, N in zip(sources, precisions)]
    point, N = F(1, K.p**2), 16
    r = S.from_polynomial(K, cs, center=pi_power(K, -1), precision=40)
    answer = r(K(point, absprec=N))
    for j in (-2, -1, 0, 1, 2):
        X = matrix(K, point)+j*matrix(K, pi_power(K, N))
        coefficients = [matrix(K, c)+(j+i)*matrix(K, pi_power(K, n))
                        for i, (c, n) in enumerate(zip(sources, precisions))]
        expected = sum((a*X**i for i, a in enumerate(coefficients)), sp.zeros(X.rows))
        assert answer.contains(image(expected))


@pytest.mark.parametrize('K', FIELDS)
def test_arithmetic_domain_intersections_and_unknown_zero_balls(K):
    f = S(K, [1]*30, radius=1, tail_bound=(0, 0), precision=20)
    g = S.from_polynomial(K, [1, -2, 3], precision=20)
    x = pi_power(K, 2)
    X = matrix(K, x)
    I = sp.eye(X.rows)
    Fm, Gm = (I-X).inv(), I-2*X+3*X**2
    for model, expected in ((f+g, Fm+Gm), (f*g, Fm*Gm), (f*f, Fm*Fm),
                            (2*f-3, 2*Fm-3*I), (3-f, 3*I-Fm), (-f, -Fm)):
        assert model(x).contains(image(expected))
        assert model.radius == 1
    assert (f+f.restrict(3)).radius == 3
    with pytest.raises(DomainError):
        (f+f.restrict(3))(x)
    zero = S.from_polynomial(K, [], precision=30)
    uncertain_zero = S.from_polynomial(K, [K.O(3)], precision=30)
    assert (f*zero).tail_bound is None
    assert (f*zero)(x, precision=90).precision == 90
    assert (f*uncertain_zero).tail_bound is not None
    assert (f*uncertain_zero)(x).precision == 3
    with pytest.raises(PrecisionError):
        (f*uncertain_zero)(x, precision=4)


def test_empty_prefix_unknown_tail_and_no_manufactured_digits():
    K = Qp(3, 3)
    with pytest.raises(TypeError):
        S(K, [1])
    unknown_constant = S(K, [], radius=1, tail_bound=(0, 0), precision=20)
    assert unknown_constant(0).precision == 0
    assert unknown_constant[0].precision == 0
    assert unknown_constant.antiderivative(7)(0, precision=30).contains(7)
    r = S(K, [1], radius=1, tail_bound=(0, 0), precision=20)
    assert r.tail_precision(0) is None
    assert r.tail_precision(K.O(2)) == 2
    assert r(K.O(2)).precision == 2
    assert r[7].is_zero() and r[7].precision == 0
    assert r.derivative(3).order == 0
    assert r.derivative(3).tail_bound is not None
    with pytest.raises(PrecisionError):
        r(3, precision=2)
    assert r(3).precision == 1
    assert r(0, precision=25).precision == 25
    primitive = S(Q2i(3), [0]*3, radius=1, tail_bound=(0, 0), precision=30).antiderivative()
    assert primitive.tail_precision((0, 1)) == 0
    # The allowed omitted monomial pi_variable**3 integrates to pi_variable**4/4.
    assert primitive((0, 1)).contains((-1, 0))  # (1+i)**4/4 = -1


@pytest.mark.parametrize('K', [FIELDS[0], FIELDS[2], FIELDS[4], FIELDS[5]])
def test_arithmetic_with_distinct_slopes_and_admissible_nonzero_omitted_coefficients(K):
    f = S(K, [1, pi_power(K, -2)], center=pi_power(K, -1),
          radius=2, tail_bound=(-4, F(-1, 2)), precision=30)
    g = S(K, [K(1, absprec=8), 0, pi_power(K, -3), F(1, K.p)],
          center=pi_power(K, -1), radius=1, tail_bound=(1, F(1, 3)), precision=30)
    H = matrix(K, pi_power(K, 8))
    point = image(matrix(K, pi_power(K, -1))+H)
    for n in (8, 15, 26):
        actual_f = [matrix(K, 1), matrix(K, pi_power(K, -2))]
        actual_f += [sp.zeros(H.rows) for _ in range(n+1-len(actual_f))]
        actual_f[n] = matrix(K, pi_power(K, ceil(-4-F(n, 2))))
        actual_g = [matrix(K, 1)+matrix(K, pi_power(K, 8)), sp.zeros(H.rows),
                    matrix(K, pi_power(K, -3)), matrix(K, F(1, K.p))]
        actual_g += [sp.zeros(H.rows) for _ in range(n+2-len(actual_g))]
        actual_g[n+1] = matrix(K, pi_power(K, ceil(1+F(n+1, 3))))
        Fm = sum((a*H**j for j, a in enumerate(actual_f)), sp.zeros(H.rows))
        Gm = sum((a*H**j for j, a in enumerate(actual_g)), sp.zeros(H.rows))
        assert (f+g)(point).contains(image(Fm+Gm))
        assert (f*g)(point).contains(image(Fm*Gm))
        assert (f-g)(point).contains(image(Fm-Gm))
        primitive = sum((a*H**(j+1)/sp.Integer(j+1) for j, a in enumerate(actual_f)), sp.zeros(H.rows))
        primitive += sum((a*H**(j+1)/sp.Integer(j+1) for j, a in enumerate(actual_g)), sp.zeros(H.rows))
        assert (f+g).antiderivative()(point).contains(image(primitive))


def test_constructor_and_factory_domain_errors_are_explicit():
    K, L = Qp(2), Q2i(3)
    for field in (2, object(), Polynomial([1, 1])):
        with pytest.raises(TypeError):
            S(field, [1], tail_bound=None)
    for bound in ((0, 0), (0, -1)):
        with pytest.raises(DomainError):
            S(K, [1], radius=0, tail_bound=bound)
    with pytest.raises(DomainError):
        S(K, [1], tail_bound=(0, 0))
    for bound in (0, (0,), (0, 1, 2), ('x', 0), (0, 0.5)):
        with pytest.raises((TypeError, ValueError)):
            S(K, [1], radius=2, tail_bound=bound)
    for kwargs in ({'radius': True}, {'precision': 2.0}, {'center': K(0)}):
        with pytest.raises(TypeError):
            S(K, [1], tail_bound=None, **kwargs)
    for field, coefficients in ((K, [0.5]), (K, [Qp(3)(1)]), (L, [K(1), Qp(3)(2)])):
        with pytest.raises(TypeError):
            S(field, coefficients, tail_bound=None)
    r = S.from_function(exp, K, order=30)
    for operation in (lambda: r(Qp(3)(9)), lambda: r(K.O(1)),
                      lambda: r.restrict(1), lambda: r.derivative(-1),
                      lambda: r[-1], lambda: r(4, precision=True),
                      lambda: r+S.from_polynomial(K, [1], center=1),
                      lambda: r*S.from_polynomial(Qp(3), [1])):
        with pytest.raises((TypeError, ValueError)):
            operation()
    for function, field, kwargs in ((exp, L, {'center': 2}),
                                    (exp, L, {'radius': 2}),
                                    (log1p, L, {'center': 1}),
                                    (log, L, {'center': 0}),
                                    (log, L, {'center': 2, 'radius': 2}),
                                    (lambda x:x, L, {}),
                                    (exp, L, {'order': 0}), (exp, L, {'order': True})):
        with pytest.raises((TypeError, ValueError)):
            S.from_function(function, field, **kwargs)


@pytest.mark.parametrize('K', FIELDS)
def test_logarithm_factories_at_nonunit_centers_and_full_shallow_disc(K):
    for scale in (-3, -2, -1, 1, 2, 3):
        c = pi_power(K, scale)
        r = S.from_function(log, K, center=c, order=96, precision=18)
        point = image(matrix(K, c)*(sp.eye(matrix(K, c).rows)+matrix(K, pi_power(K, 1))))
        assert r.radius == scale+1
        assert r(point, precision=18).is_congruent(log(K(point, absprec=28)))
        assert r.derivative()(point).is_congruent(1/K(point, absprec=28))
        with pytest.raises(DomainError):
            r(0)
    r = S.from_function(log1p, K, order=96, precision=18)
    point = pi_power(K, 1)
    assert r(point, precision=18).is_congruent(log1p(K(point, absprec=28)))


@pytest.mark.parametrize('K', FIELDS[3:])
def test_genuinely_quadratic_centers_and_shallow_logarithms(K):
    for function in (exp, sin, cos, log1p, log):
        c = (1+K.p**2, K.p**2) if function is log else (K.p**2, K.p**2)
        point = (c[0]+K.p**3, c[1]+K.p**3)
        r = S.from_function(function, K, center=c, order=96, precision=21)
        assert r(point, precision=21).is_congruent(function(K(point, absprec=30)), 21)
    c = pi_power(K, 1)
    point = image(2*matrix(K, c))
    r = S.from_function(log1p, K, center=c, order=96, precision=21)
    assert r(point, precision=21).is_congruent(log1p(K(point, absprec=30)), 21)


@pytest.mark.parametrize('K', FIELDS)
@pytest.mark.parametrize('function', [log, log1p])
def test_logarithmic_derivatives_retain_exact_center_source_information(K, function):
    center = pi_power(K, 1)
    r = S.from_function(function, K, center=center, order=20, precision=9)
    C = matrix(K, center)
    scale = C if function is log else sp.eye(C.rows)+C
    # The value at the center is still a finite ball, but every positive
    # degree coefficient is exact source data, even in a quadratic field.
    assert hasattr(r[0], 'precision')
    assert all(not hasattr(a, 'precision') for a in r.coefficients[1:])
    for d in (1, 2, 5):
        t = sp.Symbol('t')
        derivative = sp.diff(sp.log(t), t, d)
        numerator, denominator = sp.fraction(derivative)
        expected = polynomial_value(numerator, scale)*polynomial_value(denominator, scale).inv()
        assert r.derivative(d)(center, precision=60).contains(image(expected))
    with pytest.raises(PrecisionError):
        r(center, precision=60)
