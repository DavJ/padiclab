"""Quadratic analytic balls checked with exact rational multiplication matrices.

SymPy evaluates fixed polynomials, independently of the pair recurrence and
the stopping rule in the implementation. The fixed order 80 is ample here:
even the slowest tail k-2*log_2(k) is above our 21-digit targets for k>=80.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial

import pytest
import sympy as sp

from padic_tools import (
    DomainError, EisensteinQuadraticField, Polynomial, PrecisionError, Q2i, Qp,
    UnramifiedQuadraticField, cos, exp, log, log1p, power, sin, sqrt, taylor,
    teichmuller, valuation,
)
from padic_tools.analytic import _log_tail_bound


FIELDS = [
    UnramifiedQuadraticField(Qp(2, 3), [1, 1, 1]),
    UnramifiedQuadraticField(Qp(3, 3), [1, 0, 1]),
    Q2i(3),
    EisensteinQuadraticField(Qp(2, 3), [-2, 0, 1]),
    EisensteinQuadraticField(Qp(3, 3), [-3, 3, 1]),
    EisensteinQuadraticField(Qp(5, 3), [5, 5, 1]),
]


def exact_pi_power(L, n):
    if L.ramification_index == 1:
        return F(L.p)**n, F(0)
    t = sp.Symbol("t")
    f = sum(sp.Rational(c)*t**j for j, c in enumerate(L.polynomial.coefficients))
    value = sp.rem(t**n, f, t) if n >= 0 else sp.invert(t**(-n), f, t)
    return tuple(F(value.coeff(t, j)) for j in range(2))


@lru_cache(None)
def matrix_series(L, coordinates, function, terms=80):
    a, b = map(sp.Rational, coordinates)
    C, B, _ = map(sp.Rational, L.polynomial.coefficients)
    M = sp.Matrix([[a, -C*b], [b, a-B*b]])
    if function is log1p:
        cs = [sp.Rational(0)] + [sp.Rational((-1)**(j+1), j) for j in range(1, terms)]
    else:
        cs = [sp.Rational(
            1 if function is exp else
            (-1)**((j-1)//2) if function is sin and j % 2 else
            (-1)**(j//2) if function is cos and not j % 2 else 0,
            factorial(j)) for j in range(terms)]
    result = sp.zeros(2)
    for c in reversed(cs):
        result = result*M+c*sp.eye(2)
    return F(result[0, 0]), F(result[1, 0])


@pytest.mark.parametrize("L", FIELDS)
@pytest.mark.parametrize("function", [exp, sin, cos, log1p])
def test_exact_matrix_oracle_and_precision_above_parent(L, function):
    coordinates = (F(L.p**2), F(L.p**2))
    oracle = matrix_series(L, coordinates, function)
    for N in (11, 18):
        value = function(L(coordinates, absprec=N))
        assert value.precision == N
        assert value.contains(oracle)
        assert function(L(coordinates, absprec=26)).truncate(N).is_congruent(value, N)


@pytest.mark.parametrize("L", FIELDS)
def test_exponential_boundary_including_whole_zero_ball(L):
    s = L.ramification_index//(L.p-1)+1
    for function in (exp, sin, cos):
        source = exact_pi_power(L, s)
        value = function(L(source, absprec=17))
        assert value.precision == 17
        assert value.contains(matrix_series(L, source, function))
        assert function(L.O(s)).precision == s
        for x in (L(exact_pi_power(L, s-1), absprec=17), L.O(s-1)):
            with pytest.raises(DomainError):
                function(x)
    for x in (L(1, absprec=17), L.O(0)):
        with pytest.raises(DomainError):
            log1p(x)
    with pytest.raises(PrecisionError):
        exp(L.O(s), precision=s+1)
    source = exact_pi_power(L, 1)
    assert log1p(L(source, absprec=17)).contains(matrix_series(L, source, log1p))


@pytest.mark.parametrize("L", FIELDS)
def test_homomorphisms_conjugation_and_extension_valued_exponents(L):
    x, y = L(L.p**2, L.p**2, absprec=21), L(0, L.p**3, absprec=21)
    assert log(exp(x)).is_congruent(x, 21)
    assert exp(log(1+x)).is_congruent(1+x, 21)
    assert exp(x+y).is_congruent(exp(x)*exp(y), 21)
    assert log((1+x)*(1+y)).is_congruent(log(1+x)+log(1+y), 21)
    assert (sin(x)**2+cos(x)**2).is_congruent(1, 21)
    for function in (exp, log1p, sin, cos):
        assert function(x.conjugate()).is_congruent(function(x).conjugate(), 21)
    assert power(1+x, 3).is_congruent((1+x)**3, 21)
    assert (power(1+x, F(1, L.p+1))**(L.p+1)).is_congruent(1+x, 21)
    a, b = L(1, 1, absprec=19), L(2, 1, absprec=19)
    assert power(1+x, a+b).is_congruent(power(1+x, a)*power(1+x, b), 19)
    assert power(1+x, Qp(L.p, 18)(2)).is_congruent((1+x)**2)


@pytest.mark.parametrize("L", FIELDS)
@pytest.mark.parametrize("v", [-3, 1, 2, 4])
def test_global_log_two_branch_normalizations_and_norm_trace(L, v):
    x = L(exact_pi_power(L, v), absprec=23)*(1+L(L.p**2, L.p**2, absprec=25))
    pi = L(exact_pi_power(L, 1), absprec=30)
    assert log(x).is_congruent(log(x, pi_branch=log(pi)))
    assert log(x, pi_branch=(1, 1)).is_congruent(
        log(x, pi_branch=0)+v*L(1, 1, absprec=25))
    assert log(x, p_branch=2).is_congruent(log(x)+v*F(2, L.ramification_index))
    assert log(L(L.p, absprec=23)).is_congruent(0)
    assert log(L(L.p, absprec=23), p_branch=(1, 1)).is_congruent((1, 1))
    assert log(x).trace().is_congruent(log(x.norm()))


def test_shallow_wild_log_encloses_ideal_perturbations_and_later_spikes():
    L = EisensteinQuadraticField(Qp(2, 3), [-2, 0, 1])
    assert _log_tail_bound(3, 1, 2, 2) == 0  # k=4 is worse than k=3
    for N in range(1, 8):
        a, b = exact_pi_power(L, N)
        value = log1p(L(0, 1, absprec=N))
        assert value.precision == (0 if N == 1 else N)
        for r, s in ((0, 0), (1, 0), (-1, 0), (0, 1), (1, 1)):
            perturbed = (a*r+2*b*s, 1+b*r+a*s)
            assert value.contains(matrix_series(L, perturbed, log1p))
    with pytest.raises(PrecisionError):
        log1p(L.O(1), precision=1)
    assert log1p(L.O(1)).precision == 0
    assert log1p(L.O(1), precision=-2).precision == -2


@pytest.mark.parametrize("p,e", [(2, 1), (2, 2), (3, 1), (3, 2), (5, 2)])
def test_log_tail_bound_against_all_later_terms_and_power_spikes(p, e):
    for s in (1, 2, 4):
        for start in range(1, 30):
            bound = _log_tail_bound(start, s, p, e)
            assert all(bound <= k*s-e*valuation(k, p) for k in range(start, 250))
            assert all(bound <= p**j*s-e*j for j in range(1, 12) if p**j >= start)


def test_torsion_finite_branch_uncertainty_and_rejected_domains():
    L, N = Q2i(3), 19
    i, pi, p = L(-1, 1, absprec=N), L(0, 1, absprec=N), L(2, absprec=N)
    assert log(i).is_congruent(0, N)
    assert log(pi).is_congruent(0, N-1)
    assert log(p, pi_branch=1).is_congruent(2)
    assert log(pi, p_branch=1).contains(F(1, 2))
    finite = L(1, 1, absprec=7)
    assert log(pi, p_branch=finite).precision == 5
    assert log(pi, pi_branch=finite).precision == 7
    with pytest.raises(PrecisionError):
        log(pi, p_branch=finite, precision=6)
    assert log(i, p_branch=finite).precision == N
    for function in (exp, log1p, log, sin, cos):
        with pytest.raises(TypeError):
            function(F(1, 3))
        with pytest.raises(TypeError):
            function(L(4, 4, absprec=N), precision=True)
    for branch in (Qp(3)(1), FIELDS[0](1), 0.5):
        with pytest.raises(TypeError):
            log(pi, p_branch=branch)
        with pytest.raises(TypeError):
            log(i, pi_branch=branch)
    with pytest.raises(ValueError):
        log(pi, p_branch=0, pi_branch=0)
    with pytest.raises(PrecisionError):
        log(L.O(10))
    for function in (sqrt, teichmuller):
        with pytest.raises(TypeError):
            function(i)
    for exponent in (F(1, 2), L(0, F(1, 2)), Qp(3)(1), 0.5):
        with pytest.raises((TypeError, DomainError)):
            power(L(5, absprec=N), exponent)
    with pytest.raises(DomainError):
        power(L(1, 1, absprec=N), 2)


@pytest.mark.parametrize("L", FIELDS)
@pytest.mark.parametrize("function", [exp, sin, cos, log1p, log])
def test_quadratic_taylor_and_whole_remainder(L, function):
    center = L(L.p**2, L.p**2, absprec=21)
    point = center+L(0, L.p**4, absprec=23)
    jet = taylor(function, center, order=24)
    expected = center.precision-center.v if function is log else center.precision
    assert jet(point).precision == expected
    assert jet(point).is_congruent(function(point))
    assert jet.tail_precision(point) >= 21
    for order in (1, 2, 3, 6):
        assert taylor(function, center, order=order)(point).is_congruent(function(point))


def test_wild_taylor_tail_and_polynomial_rational_weights():
    L = EisensteinQuadraticField(Qp(2, 3), [-2, 0, 1])
    jet = taylor(log1p, L(0, absprec=21), order=3)
    point = L(0, 1, absprec=21)
    assert jet.tail_precision(point) == 0
    assert jet(point).contains(matrix_series(L, (F(0), F(1)), log1p))
    f = Polynomial([F(1, 4), 0, F(1, 2), 0, 3])
    center = L(0, F(1, 2), absprec=20)
    for order in (1, 2, 5, 7):
        assert taylor(f, center, order=order)(point).is_congruent(f(point))
    assert taylor(Polynomial(), center, order=3)(point).contains(0)
    assert taylor(f, center, order=5).tail_precision(point) is None
