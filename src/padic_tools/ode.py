"""Local linear differential equations with whole-solution tail bounds."""
from __future__ import annotations

from fractions import Fraction
from math import floor, inf

from .calculus import (
    AnalyticSeries, _add, _coefficient, _exact_zero, _integer, _multiply,
    _scale, _v,
)
from .core import DomainError, PadicField


def _rate(a, b):
    """Find L with v(a_j)>=L*(j+1), and L<=b's infinite-tail slope."""
    candidates = [Fraction(_v(a.field, c), j+1)
                  for j, c in enumerate(a.coefficients)
                  if not _exact_zero(a.field, c)]
    if a.tail_bound is not None:
        A, B = a.tail_bound
        candidates.extend((B, (A+B*a.order)/(a.order+1)))
    if b.tail_bound is not None:
        candidates.append(b.tail_bound[1])
    return min(candidates, default=Fraction(0))


def solve_linear_ode(coefficient, initial=1, *, forcing=0, order=20,
                     radius=None, precision=None):
    """Enclose the analytic solution of y'=coefficient*y+forcing, y(c)=initial.

    ``coefficient`` must be an AnalyticSeries over Qp or a certified
    quadratic field. ``forcing`` is a compatible series at the same exact
    center c, or a scalar constant. Exact initial data stays exact; finite
    coefficient/initial balls retain their uncertainty. ``order`` is the
    positive length of the computed coefficient prefix, not p-adic precision.
    Missing input coefficients are bounded balls, never exact zero padding.

    Returns AnalyticSeries with a proved whole-tail bound. The automatic
    disc intersects both input domains and a sufficient factorial-majorant
    convergence disc. An explicit radius may shrink this disc, not enlarge
    it. The disc and precision bounds are conservative. Evaluation uses the
    usual precision cap, or raises PrecisionError for an unsupported target.

    For fixed input data there is one analytic power-series solution. The
    returned coefficient/tail family encloses all allowed solutions and
    can contain other functions, since coefficient correlations are lost.
    Custom input tail certificates remain caller premises. This routine
    does not solve nonlinear equations, systems or singular initial points.
    """
    if not isinstance(coefficient, AnalyticSeries):
        raise TypeError("coefficient must be an AnalyticSeries")
    _integer(order, "ODE series order")
    if order < 1:
        raise ValueError("ODE series order must be positive")
    if radius is not None:
        _integer(radius, "radius")
    a, b = coefficient, coefficient._coerce(forcing)
    K = a.field
    y0 = _coefficient(K, initial)
    target = min(a.precision, b.precision) if precision is None else _integer(precision, "precision")
    radii = [r for r in (a.radius, b.radius) if r is not None]

    L = _rate(a, b)
    H = min(_v(K, y0), b._global_bound(L)-L)
    # A zero multiplier and a fully covered polynomial forcing have no
    # omitted solution coefficients. A structurally zero initial value
    # and forcing also give zero, even for an uncertain multiplier.
    finite = (a.tail_bound is None and all(_exact_zero(K, c) for c in a.coefficients)
              and b.tail_bound is None
              and all(_exact_zero(K, c) for c in b.coefficients[order-1:]))
    bound = None
    if not finite and H != inf:
        e = 1 if isinstance(K, PadicField) else K.ramification_index
        slope = L-Fraction(e, K.p-1)
        radii.append(floor(-slope)+1)
        bound = (H, slope)
    minimum = max(radii) if radii else None
    if radius is not None and minimum is not None and radius < minimum:
        raise DomainError(f"ODE radius must be at least {minimum} for input domains and the convergence bound")
    r = minimum if radius is None else radius

    cs = [y0]
    for n in range(order-1):
        value = b[n]
        for j in range(n+1):
            value = _add(K, value, _multiply(K, a[j], cs[n-j]))
        cs.append(_scale(K, value, Fraction(1, n+1)))
    return AnalyticSeries(K, cs, center=a.center, radius=r,
                          tail_bound=bound, precision=target)
