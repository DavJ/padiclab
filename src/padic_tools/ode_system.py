"""Matrix linear differential equations with a uniform solution majorant."""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import floor, inf

from .calculus import (
    AnalyticSeries, _FIELDS, _add, _center, _coefficient, _exact_zero,
    _field, _integer, _multiply, _scale, _v,
)
from .core import DomainError, PadicField
from .linalg import Matrix, _MAX_DIMENSION
from .ode import _multiplier_rates


def _rows(source, name):
    if isinstance(source, Matrix):
        return source.rows
    if isinstance(source, (str, bytes)):
        raise TypeError(f"{name} must be a matrix or an iterable of rows")
    rows = tuple(tuple(row) for row in source)
    if not rows or not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError(f"{name} must be a nonempty rectangular array")
    return rows


@dataclass(frozen=True, init=False)
class LinearSystemSolution:
    """Immutable matrix of AnalyticSeries returned by solve_linear_system.

    ``solution[i,j]`` is one component series; calling the solution returns
    a local Matrix, including for a vector (one column). Entries share field,
    exact center, order, disc, default precision and whole-tail certificate.
    Their coefficient balls discard correlations and can contain functions
    that do not solve the original equation. No matrix exponential or
    invertibility of the initial matrix is assumed.
    """
    entries: tuple

    def __init__(self, entries):
        rows = _rows(entries, "solution entries")
        if max(len(rows), len(rows[0])) > _MAX_DIMENSION:
            raise DomainError("native ODE systems support dimensions <= 64")
        if any(not isinstance(s, AnalyticSeries) for row in rows for s in row):
            raise TypeError("solution entries must be AnalyticSeries")
        first = rows[0][0]
        for row in rows:
            for s in row:
                first._coerce(s)
                if (s.order, s.radius, s.precision, s.tail_bound) != (
                        first.order, first.radius, first.precision, first.tail_bound):
                    raise ValueError("solution entries must share order, disc, precision and tail bound")
        object.__setattr__(self, "entries", rows)

    @property
    def shape(self):
        """The solution matrix shape; a vector has one column."""
        return len(self.entries), len(self.entries[0])

    @property
    def field(self):
        """Common certified coefficient field."""
        return self.entries[0][0].field

    @property
    def center(self):
        """Exact common expansion center."""
        return self.entries[0][0].center

    @property
    def radius(self):
        """Common uniformizer exponent of the declared closed disc."""
        return self.entries[0][0].radius

    @property
    def order(self):
        """Number of stored coefficients in each component."""
        return self.entries[0][0].order

    @property
    def precision(self):
        """Default target; evaluation caps each entry by its available digits."""
        return self.entries[0][0].precision

    @property
    def tail_bound(self):
        """Common (intercept, slope) solution-coefficient bound, or zero tail."""
        return self.entries[0][0].tail_bound

    def __getitem__(self, index):
        if isinstance(index, tuple):
            i, j = index
            return self.entries[i][j]
        return self.entries[index]

    def __iter__(self):
        return iter(self.entries)

    def __call__(self, x, *, precision=None):
        """Evaluate entrywise, preserving each component's available digits."""
        return Matrix(self.field, [[s(x, precision=precision) for s in row] for row in self.entries])


def solve_linear_system(coefficient, initial, *, forcing=None, field=None,
                        center=None, order=20, radius=None, precision=None):
    """Enclose Y'=A(x)Y+B(x), Y(c)=initial, over a certified local field.

    A is a nonempty square Matrix or rows of scalars/AnalyticSeries. The
    initial value is a vector of scalars or a Matrix with any positive
    number of columns. Use Matrix for multiple initial columns; scalar
    quadratic coordinate pairs in a vector remain unambiguous. Forcing is
    None (zero), a vector for a one-column solution, a matching Matrix, or
    rectangular rows of scalar/series entries for several columns.

    Field and exact center are inferred from series entries, then local
    Matrix bases; otherwise field is required and center defaults to zero.
    All series entries must have the same field and exact center. Exact
    source scalars stay exact; existing local Matrix entries remain balls.

    The returned LinearSystemSolution has order stored coefficients per
    entry, a common sufficient disc and a proved whole-solution tail bound.
    Product order is A_j*Y_(n-j), so variable coefficient matrices need not
    commute. The majorant uses the minimum entry valuation; it is conservative
    and discards correlations. Custom input tail bounds remain premises.
    Explicit evaluation targets above available information raise
    PrecisionError; the solver precision is a default evaluation target.
    Nonlinear/singular equations and optimal lattice precision are outside
    this API. Both matrix dimensions are limited to 64.
    """
    _integer(order, "ODE series order")
    if order < 1:
        raise ValueError("ODE series order must be positive")
    if radius is not None:
        _integer(radius, "radius")
    raw_a = _rows(coefficient, "coefficient")
    d = len(raw_a)
    if any(len(row) != d for row in raw_a):
        raise ValueError("the ODE coefficient matrix must be square")
    if d > _MAX_DIMENSION:
        raise DomainError("native ODE systems support dimensions <= 64")
    if isinstance(initial, (str, bytes)):
        raise TypeError("initial must be a vector or Matrix")
    raw_y = initial.rows if isinstance(initial, Matrix) else tuple((v,) for v in initial)
    if len(raw_y) != d or not raw_y or not raw_y[0]:
        raise ValueError("initial value must have one row per equation")
    m = len(raw_y[0])
    if m > _MAX_DIMENSION:
        raise DomainError("native ODE systems support dimensions <= 64")
    if forcing is None:
        raw_b = tuple((0,)*m for _ in range(d))
    elif isinstance(forcing, Matrix) or m > 1:
        raw_b = _rows(forcing, "forcing")
    else:
        if isinstance(forcing, (str, bytes)):
            raise TypeError("forcing must be a vector or Matrix")
        raw_b = tuple((v,) for v in forcing)
    if len(raw_b) != d or any(len(row) != m for row in raw_b):
        raise ValueError("forcing must have the same shape as the initial value")

    sources = [s for rows in (raw_a, raw_b) for row in rows for s in row
               if isinstance(s, AnalyticSeries)]
    if field is None:
        bases = [x.base for x in (coefficient, initial, forcing)
                 if isinstance(x, Matrix) and isinstance(x.base, _FIELDS)]
        field = sources[0].field if sources else (bases[0] if bases else None)
    _field(field)
    c = _center(field, sources[0].center if center is None and sources else (0 if center is None else center))
    target = (min((s.precision for s in sources), default=field.precision)
              if precision is None else _integer(precision, "precision"))
    prototype = AnalyticSeries(field, [], center=c, tail_bound=None, precision=target)
    a = tuple(tuple(prototype._coerce(s) for s in row) for row in raw_a)
    b = tuple(tuple(prototype._coerce(s) for s in row) for row in raw_b)
    aa, bb = [s for row in a for s in row], [s for row in b for s in row]
    y0 = [[_coefficient(field, v) for v in row] for row in raw_y]
    rates = [rate for s in aa for rate in _multiplier_rates(s)]
    rates.extend(s.tail_bound[1] for s in bb if s.tail_bound is not None)
    L = min(rates, default=Fraction(0))
    H = min(min(_v(field, v) for row in y0 for v in row),
            min(s._global_bound(L) for s in bb)-L)
    finite = (all(s.tail_bound is None and all(_exact_zero(field, v) for v in s.coefficients) for s in aa)
              and all(s.tail_bound is None and all(_exact_zero(field, v) for v in s.coefficients[order-1:]) for s in bb))
    radii = [s.radius for s in aa+bb if s.radius is not None]
    bound = None
    if not finite and H != inf:
        e = 1 if isinstance(field, PadicField) else field.ramification_index
        slope = L-Fraction(e, field.p-1)
        bound = (H, slope)
        radii.append(floor(-slope)+1)
    minimum = max(radii) if radii else None
    if radius is not None and minimum is not None and radius < minimum:
        raise DomainError(f"ODE radius must be at least {minimum} for input domains and the convergence bound")
    r = minimum if radius is None else radius
    coefficients = [y0]
    for n in range(order-1):
        out = [[s[n] for s in row] for row in b]
        for k in range(n+1):
            for i in range(d):
                for l in range(d):
                    v = a[i][l][k]
                    if _exact_zero(field, v):
                        continue
                    for j in range(m):
                        out[i][j] = _add(field, out[i][j], _multiply(field, v, coefficients[n-k][l][j]))
        coefficients.append([[_scale(field, v, Fraction(1, n+1)) for v in row] for row in out])
    return LinearSystemSolution([
        [AnalyticSeries(field, [cs[i][j] for cs in coefficients], center=c,
                        radius=r, tail_bound=bound, precision=target) for j in range(m)]
        for i in range(d)])
