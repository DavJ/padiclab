"""Run with: python examples/local_calculus.py (no optional dependencies)."""
from fractions import Fraction

from padic_tools import AnalyticSeries, Polynomial, Q2i, Qp, exp, log1p


def main():
    for K in (Qp(3, 3), Q2i(3)):
        series = AnalyticSeries.from_function(exp, K, order=48, precision=21)
        x = K(K.p**2, absprec=25)
        derivative, primitive = series.derivative(2), series.antiderivative()
        assert derivative(x, precision=21).is_congruent(exp(x), 21)
        assert primitive(x, precision=21).is_congruent(exp(x)-1, 21)
        print(f"p={K.p}: second derivative = {derivative(x)}")
        print(f"p={K.p}: primitive with F(0)=0 = {primitive(x)}")

    K = Qp(3, 3)
    f = AnalyticSeries.from_polynomial(K, Polynomial([1, 2, 3]), center=2, precision=18)
    assert f.derivative()(5).contains(32)
    assert f.antiderivative(7)(2).contains(7)
    print("Shifted exact polynomial coefficients:", f.coefficients)

    family = AnalyticSeries(K, [1], radius=1, tail_bound=(0, 0), precision=18)
    print("One known coefficient, integral unknown tail at x=3:", family(3))
    assert family(3).precision == 1
    assert family(3).contains(Fraction(-1, 2))  # The geometric series is allowed.

    L = Q2i(3)
    logarithm = AnalyticSeries.from_function(log1p, L, order=96, precision=21)
    pi = L(0, 1, absprec=27)
    assert logarithm.derivative()(pi, precision=21).is_congruent(1/(1+pi), 21)
    assert logarithm.antiderivative()(pi, precision=21).is_congruent((1+pi)*log1p(pi)-pi, 21)
    uncertain = AnalyticSeries.from_polynomial(L, [0, L(1, absprec=9)])
    assert uncertain.antiderivative()[2].precision == 7
    print("Division by 2 loses two pi digits:", uncertain.antiderivative()[2])


if __name__ == '__main__':
    main()
