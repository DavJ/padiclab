"""Volkenborn integration of polynomials and bounded analytic-series families.

I(f)=lim_N p**(-N) sum_{j=0}^{p**N-1} f(j), with I(x**n)=B_n,
and B_1=-1/2. Affine parameters mean I(t -> f(shift+scale*t));
this functional is not translation-invariant Haar integration.
"""
from __future__ import annotations

from fractions import Fraction
from math import ceil, inf

from . import analytic
from .calculus import (
    AnalyticSeries, _add, _coefficient, _exact_zero, _field, _integer,
    _multiply, _scale, _v, _zero,
)
from .core import DomainError, PadicField, PrecisionError, exact
from .special import bernoulli
from .symbolic import Polynomial


def _key(field):
    return (PadicField, field.p) if isinstance(field, PadicField) else field._key()


def volkenborn_integral(function, field=None, *, shift=0, scale=1, precision=None):
    """Return I(t -> function(shift+scale*t)) on the standard Zp parameter.

    Polynomial with exact rational affine parameters returns an exact
    Fraction when field is omitted. Supplying Qp or a certified quadratic
    field returns a local ball. AnalyticSeries infers its coefficient field
    and retains finite coefficient and affine-parameter uncertainty.

    The whole affine chart must lie in the series' declared disc. For an
    infinite tail v_pi(a_n)>=A+B*n, m=order and
    s=min(v_pi(shift-center),v_pi(scale)), the omitted integral has
    valuation at least ceil(A+(B+s)*m-e), where v_pi(p)=e. This bound
    covers all represented functions; custom coefficient bounds remain
    caller premises. An exactly constant chart uses the sharper value bound.

    Default precision is capped by coefficient, parameter and tail data;
    explicit unsupported precision raises PrecisionError. A formal
    PowerSeries or arbitrary callable cannot establish an analytic tail.
    No factor 1/scale is inserted. Changing shift can change this functional
    even when two affine charts describe the same subset of Qp.
    """
    if isinstance(function, Polynomial):
        if field is None:
            if precision is not None:
                raise TypeError("a numerical precision requires an explicit output field")
            a, b = exact(shift), exact(scale)
            shifted = function.shift(a)
            return sum((c*b**n*bernoulli(n) for n, c in enumerate(shifted.coefficients)), Fraction(0))
        _field(field)
        function = AnalyticSeries.from_polynomial(field, function)
    elif isinstance(function, AnalyticSeries):
        if field is None:
            field = function.field
        else:
            _field(field)
            if _key(field) != _key(function.field):
                raise TypeError("the output field must match the analytic series field")
    else:
        raise TypeError("Volkenborn integration requires Polynomial or AnalyticSeries")

    target = function.precision if precision is None else _integer(precision, "precision")
    a, b = _coefficient(field, shift), _coefficient(field, scale)
    d = _add(field, a, _scale(field, _coefficient(field, function.center), -1))
    s = min(_v(field, d), _v(field, b))
    if function.radius is not None and s < function.radius:
        raise DomainError("the whole affine Zp chart must lie in the declared series disc")

    # Substitute h=d+b*t using exact/ball coefficient arithmetic. Finite
    # input uncertainty is retained; exact structural zero factors vanish.
    cs = []
    for coefficient in reversed(function.coefficients):
        out = [_zero(field) for _ in range(len(cs)+1)]
        for n, c in enumerate(cs):
            out[n] = _add(field, out[n], _multiply(field, c, d))
            out[n+1] = _add(field, out[n+1], _multiply(field, c, b))
        out[0] = _add(field, out[0], _coefficient(field, coefficient))
        cs = out
    value = _zero(field)
    for n, c in enumerate(cs):
        value = _add(field, value, _scale(field, c, bernoulli(n)))
    available = value.precision if analytic._is_ball(value) else inf
    if function.tail_bound is not None:
        if _exact_zero(field, b):
            tail = function._tail_precision(_v(field, d))
        else:
            A, B = function.tail_bound
            e = 1 if isinstance(field, PadicField) else field.ramification_index
            tail = ceil(A+(B+s)*function.order-e)
        if tail is not None:
            available = min(available, tail)
    if precision is not None and target > available:
        raise PrecisionError("Volkenborn precision exceeds coefficient, parameter or tail information")
    return field(value, absprec=min(target, available))
