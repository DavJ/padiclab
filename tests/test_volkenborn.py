from fractions import Fraction as F
from math import ceil

import pytest
import sympy as sp

from padic_tools import (
    AnalyticSeries as S, DomainError, EisensteinQuadraticField, Polynomial as P,
    PowerSeries, PrecisionError, Q2i, Qp, UnramifiedQuadraticField, exp,
    volkenborn_integral as I,
)
from padic_tools.core import PadicField


FIELDS = [Qp(2, 3), Qp(3, 3), Qp(5, 3), Q2i(3),
          UnramifiedQuadraticField(Qp(3, 3), [1, 0, 1]),
          EisensteinQuadraticField(Qp(3, 3), [-3, 3, 1])]


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


def matrix_integral(M, t):
    # Independent CAS Bernoulli polynomial convention B_1(0)=-1/2.
    def entry(v):
        return sum(c*sp.bernoulli(power[0], 0) for power, c in sp.Poly(v, t).terms())
    return M.applyfunc(entry)


@pytest.mark.parametrize('degree', range(13))
def test_exact_moments_against_independent_symbolic_faulhaber_limits(degree):
    j, M = sp.symbols('j M', integer=True, positive=True)
    average = sp.summation(j**degree, (j, 0, M-1))/M
    oracle = F(sp.limit(average, M, 0))
    assert I(P([0]*degree+[1])) == oracle


@pytest.mark.parametrize('p,depth', [(2, 9), (3, 7), (5, 5)])
def test_against_actual_normalized_finite_riemann_sums(p, depth):
    f = P([1, 3, -2, 1])
    count = p**depth
    finite = sum((f(F(j)) for j in range(count)), F(0))/count
    result = I(f, Qp(p, 3), precision=depth-2)
    assert result.contains(finite)
    assert result.contains(I(f))


def test_exact_affine_parameters_and_translation_identity():
    x = P([0, 1])
    assert I(x) == F(-1, 2)
    assert I(x*x) == F(1, 6)
    assert I(x, shift=1)-I(x) == 1
    assert I(P([1]), scale=5) == 1  # No hidden chart-volume factor.
    for f in (P(), P([1, 2, -3, 4]), (x*x-2*x+F(1, 3))**5):
        assert I(f.shift(1))-I(f) == f.derivative()(F(0))
        for a, b in ((F(2, 3), F(3, 5)), (F(-2), F(0))):
            t = sp.Symbol('t')
            expression = sum(sp.Rational(c)*(a+b*t)**n for n, c in enumerate(f.coefficients))
            oracle = matrix_integral(sp.Matrix([[expression]]), t)[0, 0]
            assert I(f, shift=a, scale=b) == F(oracle)


@pytest.mark.parametrize('K', FIELDS)
def test_affine_quadratic_coefficients_and_finite_chart_uncertainty(K):
    sources = [1, pi_power(K, -1), pi_power(K, 2), -2]
    r = S.from_polynomial(K, [K(c, absprec=18+n) for n, c in enumerate(sources)],
                          center=pi_power(K, -1), precision=30)
    shift, scale = pi_power(K, 1), pi_power(K, 2)
    answer = I(r, shift=K(shift, absprec=17), scale=K(scale, absprec=16))
    t = sp.Symbol('t')
    for j in (-1, 0, 1):
        A = matrix(K, shift)+j*matrix(K, pi_power(K, 17))
        B = matrix(K, scale)+(j+1)*matrix(K, pi_power(K, 16))
        coefficients = [matrix(K, c)+(j+n)*matrix(K, pi_power(K, 18+n))
                        for n, c in enumerate(sources)]
        polynomial = sum((c*(A+B*t)**n for n, c in enumerate(coefficients)), sp.zeros(A.rows))
        assert answer.contains(image(matrix_integral(polynomial, t)))
    exact_model = S.from_polynomial(K, sources, center=pi_power(K, -1), precision=30)
    left = I(exact_model, shift=image(matrix(K, shift)+matrix(K, scale)), scale=scale)
    right = I(exact_model, shift=shift, scale=scale)
    expected = matrix(K, scale)*sum(
        (n*matrix(K, c)*matrix(K, shift)**(n-1) for n, c in enumerate(sources) if n),
        sp.zeros(matrix(K, shift).rows))
    assert (left-right).contains(image(expected))


@pytest.mark.parametrize('K', FIELDS)
def test_analytic_exponential_integral_and_explicit_domains(K):
    r = S.from_function(exp, K, order=100, precision=21)
    scale = pi_power(K, r.radius+2)
    shift = pi_power(K, r.radius+1)
    z, a = K(scale, absprec=45), K(shift, absprec=45)
    oracle = z/(exp(z)-1)
    assert I(r, scale=scale, precision=21).is_congruent(oracle, 21)
    assert I(r, shift=shift, scale=scale, precision=21).is_congruent(exp(a)*oracle, 21)
    with pytest.raises(DomainError):
        I(r)  # The unscaled Zp parameter contains points outside exp's disc.
    short = S.from_function(exp, K, order=2, precision=21)
    with pytest.raises(PrecisionError):
        I(short, scale=scale, precision=21)


@pytest.mark.parametrize('K', FIELDS)
def test_unknown_tail_includes_later_nonzero_bernoulli_moments(K):
    e = 1 if isinstance(K, PadicField) else K.ramification_index
    r = S(K, [0]*3, radius=0, tail_bound=(0, 1), precision=30)
    result = I(r)
    assert result.precision == 3-e  # B_3=0 does not annihilate the whole tail.
    for n in (3, 4, 6, 8, 12, 20):
        coefficient = matrix(K, pi_power(K, n))
        assert result.contains(image(coefficient*sp.bernoulli(n, 0)))
    with pytest.raises(PrecisionError):
        I(r, precision=4-e)
    empty = S(K, [], radius=0, tail_bound=(0, 1), precision=30)
    assert I(empty).contains(1)
    assert I(empty, scale=0).precision == 0


@pytest.mark.parametrize('K', FIELDS)
def test_coefficient_precision_bernoulli_losses_and_constant_charts(K):
    e = 1 if isinstance(K, PadicField) else K.ramification_index
    loss = e if K.p == 2 else 0
    r = S.from_polynomial(K, [0, K(1, absprec=11)], precision=30)
    result = I(r)
    assert result.precision == 11-loss
    with pytest.raises(PrecisionError):
        I(r, precision=result.precision+1)
    odd = S.from_polynomial(K, [0, 0, 0, K.O(3)])
    assert I(odd, precision=80).precision == 80  # Every allowed a_3 times B_3=0.
    family = S(K, [0], radius=1, tail_bound=(0, 0), precision=30)
    assert I(family, scale=0, precision=80).precision == 80
    uncertain = I(family, scale=K.O(7))
    assert uncertain.precision == 7-e
    assert I(family, shift=pi_power(K, 2), scale=0).is_congruent(family(pi_power(K, 2)))
    with pytest.raises(DomainError):
        I(family, shift=1, scale=0)


def test_argument_and_field_validation():
    f, K = P([1, 2]), Qp(3, 3)
    for bad in (lambda x:x, [1, 2], PowerSeries([1, 2], 10), 3):
        with pytest.raises(TypeError):
            I(bad)
    with pytest.raises(TypeError):
        I(f, precision=10)
    for kwargs in ({'shift': 0.5}, {'scale': 0.5}, {'shift': K(0)}):
        with pytest.raises(TypeError):
            I(f, **kwargs)
    for kwargs in ({'field': 3}, {'field': K, 'precision': True},
                   {'field': K, 'shift': Qp(5)(1)}):
        with pytest.raises(TypeError):
            I(f, **kwargs)
    r = S.from_polynomial(K, [1])
    for field in (Qp(5), Q2i(), UnramifiedQuadraticField(Qp(3), [1, 0, 1])):
        with pytest.raises(TypeError):
            I(r, field)
    assert I(f, K, precision=50).precision == 50
