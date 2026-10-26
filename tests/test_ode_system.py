from fractions import Fraction as F
from math import ceil, factorial

import pytest
import sympy as sp

from padic_tools import (
    AnalyticSeries as S, DomainError, EisensteinQuadraticField, LinearSystemSolution,
    Matrix, PrecisionError, Q2i, Qp, UnramifiedQuadraticField, exact,
    solve_linear_ode, solve_linear_system as solve,
)
from padic_tools.core import PadicField


FIELDS = [Qp(2, 3), Qp(3, 3), Qp(5, 3), Q2i(3),
          UnramifiedQuadraticField(Qp(2, 3), [1, 1, 1]),
          UnramifiedQuadraticField(Qp(3, 3), [1, 0, 1]),
          EisensteinQuadraticField(Qp(3, 3), [-3, 0, 1]),
          EisensteinQuadraticField(Qp(3, 3), [-3, 3, 1]),
          EisensteinQuadraticField(Qp(2, 3), [2, 2, 1])]


def scalar(K, value):
    if isinstance(K, PadicField):
        return sp.Matrix([[sp.Rational(value)]])
    a, b = value if isinstance(value, (tuple, list)) else (value, 0)
    C, B, _ = K.polynomial.coefficients
    return sp.Matrix([[a, -C*b], [b, a-B*b]])


def image(M):
    return F(M[0, 0]) if M.rows == 1 else (F(M[0, 0]), F(M[1, 0]))


def blocks(K, rows):
    return sp.BlockMatrix([[scalar(K, x) for x in row] for row in rows]).as_explicit()


def unblocks(K, M, shape):
    w = 1 if isinstance(K, PadicField) else 2
    return [[image(M[w*i:w*(i+1), w*j:w*(j+1)]) for j in range(shape[1])]
            for i in range(shape[0])]


def pi_power(K, n):
    if isinstance(K, PadicField) or K.ramification_index == 1:
        return F(K.p)**n
    t = sp.Symbol('t')
    f = sum(sp.Rational(c)*t**j for j, c in enumerate(K.polynomial.coefficients))
    q = sp.rem(t**n, f, t) if n >= 0 else sp.invert(t**(-n), f, t)
    return tuple(F(q.coeff(t, j)) for j in range(2))


def contains(K, coefficient, expected):
    if hasattr(coefficient, 'contains'):
        assert coefficient.contains(expected)
    else:
        assert scalar(K, coefficient) == scalar(K, expected)


@pytest.mark.parametrize('K', FIELDS)
def test_constant_three_by_two_system_against_rational_block_matrix_exponential(K):
    g = 2 if isinstance(K, PadicField) else (1, 1)
    raw = [[pi_power(K, -1), 1, 0], [g, -1, 1], [0, g, 1]]
    initial = Matrix(exact, [[1, 0], [0, 1], [1, 1]])
    force = [[1, g], [g, 0], [0, 1]]
    y = solve(raw, initial, forcing=force, field=K, order=32, precision=21)
    assert y.shape == (3, 2) and y.order == 32
    A, Y, B = blocks(K, raw), blocks(K, initial.rows), blocks(K, force)
    power = sp.eye(A.rows)
    for n in range(40):
        expected = Y if not n else (power*A*Y+power*B)/factorial(n)
        if n:
            power *= A
        rows = unblocks(K, expected, y.shape)
        for i in range(3):
            for j in range(2):
                contains(K, y[i, j][n], rows[i][j])
    x = pi_power(K, y.radius+2)
    X = sp.diag(*[scalar(K, x)]*3)
    Xout = sp.diag(*[scalar(K, x)]*2)
    P, Q, oracle, power = A*X, B*Xout, Y, sp.eye(A.rows)
    for n in range(1, 32):
        oracle += (power*P*Y+power*Q)/factorial(n)
        power *= P
    # The fixed block exponential's later terms have valuation >60 here.
    assert y(x, precision=21).contains(unblocks(K, oracle, y.shape))


def airy_coefficient(n, second=False):
    residue = 1 if second else 0
    if n < residue or (n-residue) % 3:
        return F(0)
    k = (n-residue)//3
    return F(1/(9**k*sp.factorial(k)*sp.rf(sp.Rational(4 if second else 2, 3), k)))


@pytest.mark.parametrize('K', FIELDS)
@pytest.mark.parametrize('inhomogeneous', [False, True])
def test_noncommuting_variable_system_against_closed_airy_hypergeometric_coefficients(K, inhomogeneous):
    c = pi_power(K, -1)
    h = S(K, [0, 1], center=c, tail_bound=None, precision=30)
    force = [[1, -h*h], [-h*h, 2*h]] if inhomogeneous else None
    y = solve([[0, 1], [h, 0]], Matrix.identity(exact, 2), forcing=force, order=24)
    assert y.precision == 30 and y.center == h.center
    contains(K, y[0, 0][3], F(1, 6))  # exp(integral(A)) would incorrectly give 1/4.
    for n in range(40):
        expected = [[airy_coefficient(n), airy_coefficient(n, True)],
                    [(n+1)*airy_coefficient(n+1), (n+1)*airy_coefficient(n+1, True)]]
        if inhomogeneous:
            expected[0][0] += int(n == 1)
            expected[1][1] += int(n == 2)
        for i in range(2):
            for j in range(2):
                contains(K, y[i, j][n], expected[i][j])
    if not inhomogeneous:
        point = image(scalar(K, c)+scalar(K, pi_power(K, y.radius+2)))
        assert y(point, precision=21).determinant().contains(1)


@pytest.mark.parametrize('K', FIELDS)
@pytest.mark.parametrize('mode', ['multiplier', 'forcing'])
def test_empty_input_prefix_families_and_later_denominator_spikes(K, mode):
    family = S(K, [], radius=2, tail_bound=(-3, F(-1, 2)), precision=30)
    a = [[0, family], [0, 0]] if mode == 'multiplier' else [[0, 0], [0, 0]]
    force = None if mode == 'multiplier' else [family, 0]
    initial = [0, 1] if mode == 'multiplier' else [0, 0]
    y = solve(a, initial, forcing=force, field=K, order=16)
    x = pi_power(K, y.radius+3)
    value = y(x)
    assert value.shape == (2, 1) and value[0, 0].precision < 21
    with pytest.raises(PrecisionError):
        y(x, precision=21)
    for n in sorted({0, 1, K.p-1, K.p**2-1, 2*K.p**2-1}):
        primitive = scalar(K, pi_power(K, ceil(-3-F(n, 2))))/(n+1)
        expected = image(primitive*scalar(K, x)**(n+1))
        assert value.contains([[expected], [initial[1]]])
        contains(K, y[0, 0][n+1], image(primitive))


@pytest.mark.parametrize('K', [FIELDS[1], FIELDS[3], FIELDS[5]])
def test_whole_matrix_initial_coefficient_forcing_and_argument_uncertainty(K):
    raw_a, raw_y, raw_b = [[1, 1], [0, -1]], [[1], [2]], [[0], [1]]
    a = Matrix(K, raw_a, absprec=11)
    initial = Matrix(K, raw_y, absprec=12)
    b = Matrix(K, raw_b, absprec=13)
    y = solve(a, initial, forcing=b, order=40, precision=30)
    x = pi_power(K, y.radius+2)
    value = y(K(x, absprec=14))
    assert all(p <= 12 for row in value.precisions for p in row)
    assert min(p for row in value.precisions for p in row) > K.precision
    for j in (-1, 0, 1):
        da, dy, db, dx = (scalar(K, pi_power(K, k)) for k in (11, 12, 13, 14))
        aa = blocks(K, [[image(scalar(K, v)+j*da) for v in row] for row in raw_a])
        yy = blocks(K, [[image(scalar(K, v)+(j-1)*dy) for v in row] for row in raw_y])
        bb = blocks(K, [[image(scalar(K, v)+(j+1)*db) for v in row] for row in raw_b])
        xx = scalar(K, x)+j*dx
        P = aa*sp.diag(xx, xx)
        Q = bb*xx
        oracle, power = yy, sp.eye(aa.rows)
        for n in range(1, 32):
            oracle += (power*P*yy+power*Q)/factorial(n)
            power *= P
        assert value.contains(unblocks(K, oracle, y.shape))
    with pytest.raises(PrecisionError):
        y(0, precision=13)


@pytest.mark.parametrize('K', FIELDS)
def test_scalar_adapter_and_quadratic_coordinate_vectors(K):
    a = S(K, [1, 2], tail_bound=None, precision=31)
    b = S(K, [3, 1], tail_bound=None, precision=30)
    native = solve_linear_ode(a, 2, forcing=b, order=20)
    y = solve([[a]], [2], forcing=[b], order=20)
    assert (y.radius, y.tail_bound, y.precision) == (native.radius, native.tail_bound, native.precision)
    for actual, expected in zip(y[0, 0].coefficients, native.coefficients):
        contains(K, actual, expected)
    if not isinstance(K, PadicField):
        z = solve([[0, 0], [0, 0]], [(1, 1), (2, -1)], field=K, order=1)
        assert z.shape == (2, 1) and z.tail_bound is None
        assert z(0, precision=60).contains([[(1, 1)], [(2, -1)]])


@pytest.mark.parametrize('K', FIELDS)
def test_zero_polynomial_and_larger_than_unit_discs(K):
    a = [[K.p**4, 0], [0, K.p**4]]
    y = solve(a, Matrix.identity(exact, 2), field=K, order=40)
    assert y.radius < 0  # Exact off-diagonal zeros do not impose rate zero.
    with pytest.raises(DomainError):
        solve(a, [1, 0], field=K, radius=y.radius-1)
    zero = [[0, 0], [0, 0]]
    h = S(K, [0, 1], tail_bound=None)
    p = solve(zero, [1, 2], forcing=[1+h, h], order=4)
    assert p.tail_bound is None and p.radius is None
    x = pi_power(K, -2)
    X, one = scalar(K, x), scalar(K, 1)
    assert p(x, precision=60).contains([[image(one+X+X*X/2)], [image(2*one+X*X/2)]])
    uncertain = Matrix(K, zero, absprec=3)
    assert solve(uncertain, [1, 0], order=1).tail_bound is not None
    exact_zero = solve(uncertain, [0, 0], order=1)
    assert exact_zero.tail_bound is None
    assert exact_zero(x, precision=60).contains([[0], [0]])


def test_context_shape_and_parameter_validation():
    K = Qp(3, 3)
    a = [[0, 1], [-1, 0]]
    for bad in ([], [[1, 2]], [[1, 2], [3]], [[1], [2]]):
        with pytest.raises(ValueError):
            solve(bad, [1, 0], field=K)
    with pytest.raises(TypeError):
        solve(a, [1, 0])
    for kwargs in ({'order': True}, {'radius': 1.5}, {'precision': True}, {'field': 3}):
        options = {'field': K, **kwargs}
        with pytest.raises(TypeError):
            solve(a, [1, 0], **options)
    for kwargs in ({'initial': [1]}, {'forcing': [1]}, {'order': 0}):
        options = {'initial': [1, 0], 'field': K, **kwargs}
        with pytest.raises(ValueError):
            solve(a, **options)
    local = S(K, [1], center=1, tail_bound=None, precision=31)
    inferred = solve([[0, local], [0, 0]], [1, 0], order=3)
    assert inferred.center == 1 and inferred.precision == 31
    with pytest.raises(TypeError):
        solve([[0, local], [0, 0]], [1, 0], center=0)
    with pytest.raises(TypeError):
        solve([[0, local], [0, 0]], [1, 0], forcing=[S.from_polynomial(Qp(5), [1]), 0])
    with pytest.raises(TypeError):
        solve(a, [Qp(5)(1), 0], field=K)
    with pytest.raises(DomainError):
        solve([[0]*65 for _ in range(65)], [0]*65, field=K)
    with pytest.raises(ValueError):
        LinearSystemSolution([[local, local.restrict(2)]])
    with pytest.raises(TypeError):
        LinearSystemSolution([[1]])
