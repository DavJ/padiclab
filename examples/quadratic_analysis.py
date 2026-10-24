"""Run with: python examples/quadratic_analysis.py (no optional dependencies)."""
from fractions import Fraction

from padic_tools import (
    EisensteinQuadraticField, Gaussian, GaussianLocalPlace, Q2i, Qp,
    UnramifiedQuadraticField, exp, log, log1p, power, taylor,
)


def main():
    for L in (UnramifiedQuadraticField(Qp(3, 3), [1, 0, 1]), Q2i(3)):
        x = L(L.p**2, L.p**2, absprec=21)
        point = x+L(0, L.p**4, absprec=23)
        value = exp(x)
        assert log(value).is_congruent(x, 21)
        assert (power(1+x, Fraction(1, L.p+1))**(L.p+1)).is_congruent(1+x, 21)
        assert taylor(exp, x, order=24)(point).is_congruent(exp(point), 21)
        print(f"p={L.p}, e={L.ramification_index}: exp(x) = {value}")

    L = Q2i(3)
    pi = L(0, 1, absprec=19)
    print("Q2(i), log(pi) with log(2)=1:", log(pi, p_branch=1))
    assert log(pi, p_branch=1).contains(Fraction(1, 2))
    print("Shallow log1p input:", L.O(1), "->", log1p(L.O(1)))

    E = EisensteinQuadraticField(Qp(3, 3), [-3, 3, 1])
    rho = E(0, 1, absprec=21)
    print("Another Eisenstein field, log(3)=0, log(rho):", log(rho))
    assert not log(rho).is_zero()
    assert (2*log(rho)).is_congruent(log(rho*rho/3))

    place = GaussianLocalPlace(Qp(2, 3), root_sign=-1)
    local = place(Gaussian(4, 4), absprec=21)
    assert log(exp(local)).is_congruent(local, 21)
    print("Chosen Gaussian local image:", local)


if __name__ == "__main__":
    main()
