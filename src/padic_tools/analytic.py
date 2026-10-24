"""p-adic analytic functions with checked domains and certified series tails.

The native algorithms prioritize transparent error bounds over asymptotic speed.
No series is stopped merely because one term happens to vanish modulo p**N.
"""
from __future__ import annotations

from fractions import Fraction

from .core import DomainError, Padic, PrecisionError, exact, p_power, valuation
from .extensions import _QuadraticElement, _multiply


def _is_ball(x):
    return isinstance(x, (Padic, _QuadraticElement))


def _check_ball(x):
    if not _is_ball(x):
        raise TypeError("first construct a ball in Qp or a certified quadratic local field")


def _ramification_index(x):
    return 1 if isinstance(x, Padic) else x.field.ramification_index


def _target(x, precision, known=None):
    known = x.precision if known is None else known
    n = known if precision is None else precision
    if type(n) is not int:
        raise TypeError("precision must be an integer")
    if n > known:
        raise PrecisionError("function evaluation cannot invent input digits")
    return n


def _exp_domain(x):
    _check_ball(x)
    if x.v * (x.p - 1) <= _ramification_index(x):
        raise DomainError("requires v_pi(x) > e/(p-1), including the whole input ball")


def _pair_product(field, a, b):
    return _multiply(a, b, field.polynomial[1], field.polynomial[0])


def _pair_power(field, a, exponent):
    if exponent < 0:
        x, y = a
        b, c = field.polynomial[1], field.polynomial[0]
        norm = x*x-b*x*y+c*y*y
        a, exponent = ((x-b*y)/norm, -y/norm), -exponent
    result = (Fraction(1), Fraction(0))
    while exponent:
        if exponent & 1:
            result = _pair_product(field, result, a)
        exponent >>= 1
        if exponent:
            a = _pair_product(field, a, a)
    return result


def _uniformizer_power(field, exponent):
    pi = ((Fraction(field.p), Fraction(0)) if field.ramification_index == 1
          else (Fraction(0), Fraction(1)))
    return _pair_power(field, pi, exponent)


def _factorial_series(x, kind, precision=None):
    _exp_domain(x)
    nout = _target(x, precision)
    if isinstance(x, _QuadraticElement):
        centre = x.representative()
        term, total, k = (Fraction(1), Fraction(0)), (Fraction(0), Fraction(0)), 0
        slope = x.v*(x.p-1)-_ramification_index(x)
        while k == 0 or k*slope < nout*(x.p-1):
            sign = (1 if kind == "exp" else
                    (-1)**((k-1)//2) if kind == "sin" and k % 2 else
                    (-1)**(k//2) if kind == "cos" and k % 2 == 0 else 0)
            total = tuple(a+sign*b for a, b in zip(total, term))
            k += 1
            term = tuple(c/k for c in _pair_product(x.field, term, centre))
        return x.field(total, absprec=nout)
    centre, term, total, k = x.representative(), Fraction(1), Fraction(0), 0
    # v(x^k/k!) >= k*(v(x)-1/(p-1)); a bound for the entire tail.
    slope = x.v * (x.p - 1) - 1
    while k == 0 or k * slope < nout * (x.p - 1):
        if kind == "exp":
            total += term
        elif kind == "sin" and k % 2:
            total += (-1)**((k-1)//2) * term
        elif kind == "cos" and k % 2 == 0:
            total += (-1)**(k//2) * term
        k += 1
        term = term * centre / k
    return x.field(total, absprec=nout)


def exp(x: Padic | _QuadraticElement, *, precision=None):
    """Exponential on v_pi(x)>e/(p-1), for Qp or certified quadratic fields.

    Output precision is absolute in the target uniformizer, at most the
    input precision. In Qp, e=1. The whole input ball must lie in the domain.
    """
    return _factorial_series(x, "exp", precision)


def sin(x: Padic | _QuadraticElement, *, precision=None):
    """Sine on the same sufficient convergence domain as exp."""
    return _factorial_series(x, "sin", precision)


def cos(x: Padic | _QuadraticElement, *, precision=None):
    """Cosine on the same sufficient convergence domain as exp."""
    return _factorial_series(x, "cos", precision)


def _floor_log(n, p):
    r = 0
    while n >= p:
        n //= p
        r += 1
    return r


def _log_tail_bound(start, v, p, e):
    """Lower bound for k*v-e*floor(log_p(k)) for EVERY k>=start.

    Between p-powers this expression increases. Inspect the first interval
    and the following p-powers until their successive minima increase.
    """
    j = _floor_log(start, p)
    bound = start*v-e*j
    power = p**(j+1)
    j += 1
    while True:
        bound = min(bound, power*v-e*j)
        if (p-1)*power*v >= e:
            return bound
        power *= p
        j += 1


def log1p(x: Padic | _QuadraticElement, *, precision=None):
    """log(1+x) on v_pi(x)>=1 in Qp or certified quadratic fields.

    Both input uncertainty and the full tail account for denominator spikes
    at p-powers. In ramified fields a shallow ball can yield fewer absolute
    output digits than input digits. Precision counts target uniformizers.
    """
    _check_ball(x)
    if x.v < 1:
        raise DomainError("log1p requires v_pi(x) >= 1")
    if isinstance(x, _QuadraticElement):
        e = x.field.ramification_index
        known = _log_tail_bound(1, x.precision, x.p, e)
        nout = _target(x, precision, known)
        c, power_x = x.representative(), x.representative()
        total, k = (Fraction(0), Fraction(0)), 1
        while _log_tail_bound(k, x.v, x.p, e) < nout:
            total = tuple(a+Fraction((-1)**(k+1), k)*b for a, b in zip(total, power_x))
            k += 1
            power_x = _pair_product(x.field, power_x, c)
        return x.field(total, absprec=nout)
    nout = _target(x, precision)
    c, power_x, total, k = x.representative(), x.representative(), Fraction(0), 1
    while k*x.v - _floor_log(k, x.p) < nout:
        total += (-1)**(k+1) * power_x / k
        k += 1
        power_x *= c
    return x.field(total, absprec=nout)


def teichmuller(x: Padic) -> Padic:
    """Root-of-unity part of a known unit.

    For odd p this lifts the residue to a (p-1)-st root of unity. At p=2
    returns the sign determined by the residue modulo 4; two digits are needed.
    """
    if not isinstance(x, Padic):
        raise TypeError("teichmuller currently supports Qp arguments")
    if x.is_zero() or x.v != 0:
        raise DomainError("teichmuller requires a known p-adic unit")
    n, p = x.precision, x.p
    if p == 2:
        if n < 2:
            raise PrecisionError("2-adic sign requires a residue modulo 4")
        return x.field(1 if x.residue(2) == 1 else -1, absprec=n)
    modulus, a = p**n, x.residue(1)
    for _ in range(n):
        a = pow(a, p, modulus)
    return x.field(a, absprec=n)


def _unit_log(x):
    m = x.field.residue_field().order-1
    return log1p(x**m-1)/m


def _branch_ball(field, value, precision):
    coordinates, known = field._data(value)
    return field(coordinates, absprec=precision if known is None else known)


def _quadratic_log(x, p_branch, pi_branch, precision):
    field, e, v = x.field, x.field.ramification_index, x.v
    unit = x/_uniformizer_power(field, v)
    answer = _unit_log(unit)
    if v:
        # Exact source data pi**e/p permits guard digits before division by e.
        work = max(e//(x.p-1)+1, answer.precision+e*valuation(e, x.p))
        if pi_branch is not None:
            branch = _branch_ball(field, pi_branch, work)
        else:
            branch_p = _branch_ball(field, 0 if p_branch is None else p_branch, work)
            u0 = tuple(c/x.p for c in _uniformizer_power(field, e))
            branch = (branch_p+_unit_log(field(u0, absprec=work)))/e
        answer = answer+v*branch
    else:
        # Branches do not affect units, but their scalar domain is checked.
        if p_branch is not None:
            field._data(p_branch)
        if pi_branch is not None:
            field._data(pi_branch)
    return answer if precision is None else answer.truncate(precision)


def log(x: Padic | _QuadraticElement, *, p_branch=None, pi_branch=None, precision=None):
    """Logarithm on nonzero balls in Qp and certified quadratic fields.

    The default is log(p)=0. Specify at most one of ``p_branch=log(p)``
    and ``pi_branch=log(pi)``. For quadratic fields these can be exact field
    coordinates or compatible balls. The relation pi**e=p*u determines
    log(pi)=(log(p)+log(u))/e. Roots of unity have logarithm zero; log and
    exp are inverse only on their common local domain.
    """
    _check_ball(x)
    if p_branch is not None and pi_branch is not None:
        raise ValueError("specify at most one of p_branch and pi_branch")
    if precision is not None and type(precision) is not int:
        raise TypeError("precision must be an integer")
    if x.is_zero():
        raise PrecisionError("logarithm of a ball containing zero is undefined")
    if isinstance(x, _QuadraticElement):
        return _quadratic_log(x, p_branch, pi_branch, precision)
    branch = pi_branch if pi_branch is not None else (0 if p_branch is None else p_branch)
    v = x.v
    unit = x / p_power(x.p, v)
    principal = unit / teichmuller(unit)
    answer = log1p(principal - 1) + v*exact(branch)
    return answer if precision is None else answer.truncate(precision)


def power(base: Padic | _QuadraticElement, exponent):
    """Principal-unit power exp(exponent*log(base)) for exponent in Z_p.

    Requires v_pi(base-1)>e/(p-1). Quadratic targets accept integral
    exponents in that field, in its Qp base, or exact rational exponents.
    No implicit fractional power of a root of unity is chosen.
    """
    _check_ball(base)
    delta = base - 1
    _exp_domain(delta)
    if isinstance(base, _QuadraticElement):
        if _is_ball(exponent):
            exponent = base.field(exponent)
            v = exponent.v
        else:
            exponent = exact(exponent)
            v = base.field.ramification_index*valuation(exponent, base.p)
        if v < 0:
            raise DomainError("exponent must be integral in the target field")
        return exp(log1p(delta)*exponent)
    if isinstance(exponent, Padic):
        base._check_prime(exponent)
        if exponent.v < 0:
            raise DomainError("exponent must lie in Z_p")
    else:
        if valuation(exact(exponent), base.p) < 0:
            raise DomainError("exponent must lie in Z_p")
    return exp(log1p(delta) * exponent)


def _tonelli_shanks(a, p):
    if a % p == 0:
        return 0
    if pow(a, (p-1)//2, p) != 1:
        raise DomainError("unit is not a square modulo p")
    if p % 4 == 3:
        r = pow(a, (p+1)//4, p)
        return min(r, p-r)
    q, s = p-1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p-1)//2, p) != p-1:
        z += 1
    c, r, t, m = pow(z, q, p), pow(a, (q+1)//2, p), pow(a, q, p), s
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = tt*tt % p
            i += 1
        b = pow(c, 2**(m-i-1), p)
        r, t, c, m = r*b % p, t*b*b % p, b*b % p, i
    return min(r, p-r)


def sqrt(x: Padic) -> Padic:
    """One square root in Q_p, with the derivative's precision loss.

    Odd p: choose the smaller nonzero initial residue. At p=2 choose the
    unit root congruent to 1 modulo 4; a square unit must be 1 modulo 8.
    The negative gives the other root. Zero balls need further information.
    """
    if not isinstance(x, Padic):
        raise TypeError("sqrt currently supports Qp arguments")
    if x.is_zero():
        raise PrecisionError("sqrt of a zero ball cannot select a root branch")
    if x.v % 2:
        raise DomainError("a square in Q_p has even valuation")
    k, rel, p, u = x.v//2, x.relative_precision, x.p, x.unit
    if p == 2:
        if rel < 3:
            raise PrecisionError("a 2-adic square unit needs three known digits")
        if u % 8 != 1:
            raise DomainError("a 2-adic square unit must be 1 modulo 8")
        r = 1
        for j in range(3, rel):
            r += (((u-r*r) >> j) & 1) * 2**(j-1)
        nout = x.precision - k - 1
    else:
        r, modulus = _tonelli_shanks(u % p, p), p
        for _ in range(1, rel):
            digit = ((u-r*r)//modulus * pow(2*r, -1, p)) % p
            r += digit * modulus
            modulus *= p
        nout = x.precision - k
    return x.field(r*p_power(p, k), absprec=nout)
