"""Reciprocals and exactly centered composition, without optional packages."""
from fractions import Fraction

from padic_tools import AnalyticSeries, Q2i, Qp, exp, log


K = Qp(3, 3)
inverse = AnalyticSeries(K, [1, -1], tail_bound=None).reciprocal(order=48, precision=21)
assert inverse(9, precision=21).contains(Fraction(-1, 8))
g = AnalyticSeries(K, [0, 1, 1], tail_bound=None)
composed = inverse.compose(g)
assert composed(9, precision=21).contains(Fraction(-1, 89))
print("1/(1-h-h^2):", composed(9, precision=21))

e = AnalyticSeries.from_function(exp, K, order=48)
l = AnalyticSeries.from_function(log, K, center=1, order=48)
assert l.compose(e)(9, precision=21).contains(9)
print("log(exp(h)):", l.compose(e)(9, precision=21))

L = Q2i(3)
f = AnalyticSeries(L, [1, (-1, 1), 1], tail_bound=None).reciprocal(order=32)
z = L(4, 4, absprec=45)
assert f((4, 4), precision=21).is_congruent(1/(1+L(-1, 1, absprec=45)*z+z*z), 21)
print("Quadratic reciprocal:", f((4, 4), precision=21))
