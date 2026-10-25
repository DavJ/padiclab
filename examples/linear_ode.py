"""Local y'=a(x)y+b(x), including a genuinely quadratic multiplier."""
from fractions import Fraction

from padic_tools import AnalyticSeries, Q2i, Qp, exp, solve_linear_ode


K = Qp(3, 3)
c = Fraction(1, 3)
a = AnalyticSeries.from_polynomial(K, [-1], center=c)
y = solve_linear_ode(a, 3, forcing=2, order=64, precision=21)
print("y'=-y+2, y(1/3)=3:", y(c+9))
print("Disc exponent and tail bound:", y.radius, y.tail_bound)
assert y(c+9, precision=21).is_congruent(2+exp(K(-9, absprec=45)), 21)

L = Q2i(3)
a = AnalyticSeries(L, [(-1, 1)], tail_bound=None)  # i=pi-1.
y = solve_linear_ode(a, order=64, precision=21)
z = L(-1, 1, absprec=45)*L(4, 4, absprec=45)
print("y'=i*y, y(0)=1:", y((4, 4)))
assert y((4, 4), precision=21).is_congruent(exp(z), 21)
