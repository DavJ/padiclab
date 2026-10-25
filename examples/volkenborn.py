"""Run with python examples/volkenborn.py; no optional dependencies needed."""
from fractions import Fraction

from padic_tools import AnalyticSeries, Polynomial, Q2i, Qp, exp, volkenborn_integral


def main():
    x = Polynomial([0, 1])
    assert volkenborn_integral(x) == Fraction(-1, 2)
    assert volkenborn_integral(x*x) == Fraction(1, 6)
    print('I(x), I(x^2):', volkenborn_integral(x), volkenborn_integral(x*x))
    print('Translation changes the functional:', volkenborn_integral(x, shift=1))
    for K in (Qp(3, 3), Q2i(3)):
        s = AnalyticSeries.from_function(exp, K, order=100, precision=21)
        z = K(K.p**2, absprec=45)
        value = volkenborn_integral(s, scale=K.p**2, precision=21)
        assert value.is_congruent(z/(exp(z)-1), 21)
        print(f'p={K.p}: I(exp(z*t)) =', value)


if __name__ == '__main__':
    main()
