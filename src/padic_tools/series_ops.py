"""Reciprocal and centered composition enclosures for analytic series."""
from __future__ import annotations

from fractions import Fraction
from math import ceil, floor

from . import analytic
from .calculus import (
    AnalyticSeries, _add, _exact_zero, _integer, _multiply, _scale, _v, _zero,
)
from .core import DomainError, PadicField, PrecisionError


def _positive_terms(series):
    """(index, valuation bound) for positive coefficients, including the tail endpoint.

    At weight rho>-B, an infinite tail's minimum occurs at its first
    positive index. Its slope constraint is checked separately by callers.
    """
    terms = [(n, Fraction(_v(series.field, a))) for n, a in enumerate(series.coefficients)
             if n and not _exact_zero(series.field, a)]
    if series.tail_bound is not None:
        A, B = series.tail_bound
        n = max(1, series.order)
        terms.append((n, A+B*n))
    return terms


def _degree(series):
    return max((n for n, a in enumerate(series.coefficients)
                if not _exact_zero(series.field, a)), default=-1)


def _order(order, default):
    n = default if order is None else _integer(order, "series order")
    if n < 1:
        raise ValueError("series order must be positive")
    return n


def _radius(radius, minimum):
    if radius is not None:
        _integer(radius, "radius")
        if minimum is not None and radius < minimum:
            raise DomainError(f"radius must be at least {minimum} for the domain certificate")
    return minimum if radius is None else radius


def _reciprocal(series, *, order=None, radius=None, precision=None):
    K = series.field
    count = _order(order, max(1, series.order))
    target = series.precision if precision is None else _integer(precision, "precision")
    a0 = series[0]
    if _exact_zero(K, a0):
        raise ZeroDivisionError("the series constant is exactly zero")
    if analytic._is_ball(a0) and a0.v >= a0.precision:
        raise PrecisionError("reciprocal requires a resolved nonzero constant coefficient")
    t = _v(K, a0)
    inverse = (a0.inverse() if analytic._is_ball(a0) else
               1/a0 if isinstance(K, PadicField) else analytic._pair_power(K, a0, -1))
    terms = _positive_terms(series)
    rates = [(v-t)/n for n, v in terms]
    if series.tail_bound is not None:
        rates.append(series.tail_bound[1])
    bound = None
    minimum = series.radius
    if rates:
        slope = min(rates)
        required = floor(-slope)+1
        minimum = required if minimum is None else max(minimum, required)
        bound = (-t, slope)
    r = _radius(radius, minimum)
    cs = [inverse]
    for n in range(1, count):
        value = _zero(K)
        for j in range(1, n+1):
            value = _add(K, value, _multiply(K, series[j], cs[n-j]))
        cs.append(_scale(K, _multiply(K, inverse, value), -1))
    return AnalyticSeries(K, cs, center=series.center, radius=r,
                          tail_bound=bound, precision=target)


def _compose(outer, inner, *, order=None, radius=None, precision=None):
    if not isinstance(inner, AnalyticSeries):
        raise TypeError("composition requires an AnalyticSeries")
    key = lambda K: (PadicField, K.p) if isinstance(K, PadicField) else K._key()
    if key(outer.field) != key(inner.field):
        raise TypeError("composition requires the same coefficient field")
    K = outer.field
    constant = inner[0]
    if analytic._is_ball(constant):
        raise TypeError("composition requires an exact inner constant coefficient")
    if not _exact_zero(K, _add(K, constant, _scale(K, outer.center, -1))):
        raise DomainError("inner constant must equal the outer exact center")
    df, dg = _degree(outer), _degree(inner)
    if outer.tail_bound is None and df <= 0:
        default = 1
    elif outer.tail_bound is None and inner.tail_bound is None:
        default = max(1, df*max(0, dg)+1)
    elif inner.tail_bound is None:
        default = max(1, outer.order)
    elif outer.tail_bound is None:
        default = max(1, inner.order)
    else:
        default = max(1, min(outer.order, inner.order))
    count = _order(order, default)
    target = min(outer.precision, inner.precision) if precision is None else _integer(precision, "precision")
    terms = _positive_terms(inner)
    if not terms:  # An exact constant inner function; all higher terms vanish.
        r = _radius(radius, inner.radius)
        return AnalyticSeries(K, [outer[0]]+[_zero(K)]*(count-1), center=inner.center,
                              radius=r, tail_bound=None, precision=target)

    # The actual disc must map into the outer declared disc. For an
    # infinite inner tail its endpoint bound applies on the input disc.
    radii = [] if inner.radius is None else [inner.radius]
    if outer.radius is not None:
        radii.extend(ceil((outer.radius-v)/n) for n, v in terms)
    minimum = max(radii) if radii else None
    finite = (outer.tail_bound is None and
              (df <= 0 or (inner.tail_bound is None and count > df*max(0, dg))))
    bound = None
    if not finite:
        # Auxiliary weight rho controls coefficients, separately from the
        # actual image-domain test. Choose rho<R so the returned tail decays.
        strict = []
        if inner.tail_bound is not None:
            strict.append(-inner.tail_bound[1])
        if outer.tail_bound is not None:
            B = outer.tail_bound[1]
            strict.extend((-B-v)/n for n, v in terms)
        q = max(strict) if strict else None
        if q is not None:
            required = floor(q)+1
            minimum = required if minimum is None else max(minimum, required)
        r = _radius(radius, minimum)
        if r is None:
            r = 0
        rho = (q+r)/2 if q is not None else Fraction(r-1)
        C = min(v+rho*n for n, v in terms)
        H = outer._global_bound(-C)
        bound = (H, -rho)
    else:
        r = _radius(radius, minimum)

    # Centered substitution means outer terms of degree >=count cannot
    # affect the prefix. Truncated Horner arithmetic retains finite balls.
    cs = [_zero(K)]*count
    last = min(count, outer.order) if outer.tail_bound is None else count
    gs = [inner[n] for n in range(count)]
    for j in reversed(range(last)):
        out = [_zero(K)]*count
        for n in range(1, count):
            value = _zero(K)
            for k in range(1, n+1):
                value = _add(K, value, _multiply(K, gs[k], cs[n-k]))
            out[n] = value
        out[0] = outer[j]
        cs = out
    return AnalyticSeries(K, cs, center=inner.center, radius=r,
                          tail_bound=bound, precision=target)
