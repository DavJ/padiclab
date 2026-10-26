"""Coupled oscillator and noncommuting Airy fundamental solution."""
from fractions import Fraction

from padic_tools import AnalyticSeries, Matrix, Q2i, Qp, cos, exact, sin, solve_linear_system


K = Qp(3, 3)
phi = solve_linear_system([[0, 1], [-1, 0]], Matrix.identity(exact, 2), field=K, order=64)
x = K(9, absprec=45)
expected = Matrix(K, [[cos(x), sin(x)], [-sin(x), cos(x)]])
assert phi(9, precision=21).is_congruent(expected, 21)
print("Oscillator fundamental matrix:", phi(9, precision=21))

h = AnalyticSeries(K, [0, 1], tail_bound=None, precision=21)
airy = solve_linear_system([[0, 1], [h, 0]], Matrix.identity(exact, 2), order=64)
assert airy[0, 0][3] == Fraction(1, 6)
assert airy(9, precision=21).determinant().contains(1)
print("Noncommuting Airy system, cubic coefficient:", airy[0, 0][3])

L = Q2i(3)
y = solve_linear_system([[0, (-1, 1)], [(1, -1), 0]], [1, 0], field=L, order=64)
z = L(-1, 1, absprec=45)*L(4, 4, absprec=45)
assert y((4, 4), precision=21)[0, 0].is_congruent(cos(z), 21)
print("Quadratic oscillator vector:", y((4, 4), precision=21))
