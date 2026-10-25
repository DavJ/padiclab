"""Local analytic series with coefficient bounds, calculus and ball evaluation.

For n >= order, tail_bound=(A,B) asserts v_pi(a_n) >= A+B*n.
The declared closed disc v_pi(x-center) >= radius must have radius+B>0.
Unknown coefficients remain unknown after differentiation and integration.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import ceil, factorial, inf

from . import analytic
from .core import DomainError, Padic, PadicField, PrecisionError, exact, valuation
from .extensions import UnramifiedQuadraticField
from .ramified import EisensteinQuadraticField
from .symbolic import Polynomial


_FIELDS = (PadicField, UnramifiedQuadraticField, EisensteinQuadraticField)


def _field(field):
    if not isinstance(field, _FIELDS):
        raise TypeError("field must be Qp or a certified quadratic local field")


def _integer(n, name):
    if type(n) is not int:
        raise TypeError(f"{name} must be an integer")
    return n


def _coefficient(field, value):
    if isinstance(field, PadicField):
        return field(value) if isinstance(value, Padic) else exact(value)
    coordinates, known = field._data(value)
    return coordinates if known is None else field(value)


def _v(field, value):
    if analytic._is_ball(value):
        return value.v
    return valuation(value, field.p) if isinstance(field, PadicField) else field._valuation(value)


def _zero(field):
    return Fraction(0) if isinstance(field, PadicField) else (Fraction(0), Fraction(0))


def _exact_zero(field, value):
    return not analytic._is_ball(value) and _v(field, value) == inf


def _add(field, a, b):
    if analytic._is_ball(a):
        return a+b
    if analytic._is_ball(b):
        return b+a
    return a+b if isinstance(field, PadicField) else tuple(x+y for x, y in zip(a, b))


def _scale(field, a, q):
    if not q or _exact_zero(field, a):
        return _zero(field)
    if analytic._is_ball(a) or isinstance(field, PadicField):
        return a*q
    return tuple(x*q for x in a)


def _multiply(field, a, b):
    if _exact_zero(field, a) or _exact_zero(field, b):
        return _zero(field)
    if analytic._is_ball(a):
        return a*b
    if analytic._is_ball(b):
        return b*a
    return a*b if isinstance(field, PadicField) else analytic._pair_product(field, a, b)


def _center(field, center):
    value = _coefficient(field, center)
    if analytic._is_ball(value):
        raise TypeError("center requires exact source data, not a finite ball")
    return value


@dataclass(frozen=True, init=False)
class AnalyticSeries:
    """A coefficient-ball family sum a_n*(x-center)**n on a closed local disc.

    ``coefficients`` gives the known prefix, in ascending degree. Exact
    rational data and exact quadratic coordinate pairs remain exact;
    compatible local balls retain their known precision.

    Mandatory ``tail_bound=(A,B)`` asserts v_pi(a_n)>=A+B*n for EVERY
    n>=len(coefficients). A caller-supplied bound is a mathematical premise,
    not inferred from a prefix. Explicit ``tail_bound=None`` means a finite
    polynomial. Prefer ``from_function`` for internally proved bounds, and
    ``from_polynomial`` for a polynomial in the original variable x.

    ``radius`` is an integer uniformizer exponent, not a real radius.
    It declares v_pi(x-center)>=radius and must satisfy radius+B>0 for
    an infinite tail. None is allowed only for a polynomial on the whole
    field. ``center`` must be exact. ``precision`` is the default absolute
    output precision, in uniformizer units, capped by available information.
    """
    field: object
    coefficients: tuple
    center: object
    radius: int | None
    tail_bound: tuple | None
    precision: int

    def __init__(self, field, coefficients, *, tail_bound, center=0, radius=None,
                 precision=None):
        _field(field)
        c = _center(field, center)
        cs = tuple(_coefficient(field, a) for a in coefficients)
        target = field.precision if precision is None else _integer(precision, "precision")
        if radius is not None:
            _integer(radius, "radius")
        if tail_bound is not None:
            if not isinstance(tail_bound, (tuple, list)) or len(tail_bound) != 2:
                raise TypeError("tail_bound must be (intercept, slope), or explicit None")
            A, B = map(exact, tail_bound)
            if radius is None or radius+B <= 0:
                raise DomainError("an infinite tail requires radius + tail slope > 0")
            tail_bound = (A, B)
        for name, value in (("field", field), ("coefficients", cs), ("center", c),
                            ("radius", radius), ("tail_bound", tail_bound),
                            ("precision", target)):
            object.__setattr__(self, name, value)

    @property
    def order(self):
        """Length of the known coefficient prefix, distinct from output precision."""
        return len(self.coefficients)

    def coefficient(self, n):
        """Known coefficient, exact polynomial zero, or a bounded unknown ball."""
        _integer(n, "coefficient index")
        if n < 0:
            raise ValueError("coefficient index must be nonnegative")
        if n < self.order:
            return self.coefficients[n]
        if self.tail_bound is None:
            return _zero(self.field)
        A, B = self.tail_bound
        return self.field.O(ceil(A+B*n))

    __getitem__ = coefficient

    def _displacement(self, x):
        h = _add(self.field, _coefficient(self.field, x), _scale(self.field, self.center, -1))
        s = _v(self.field, h)
        if self.radius is not None and s < self.radius:
            raise DomainError("the whole argument ball must lie in the declared series disc")
        return h, s

    def _tail_precision(self, s):
        if self.tail_bound is None or (s == inf and self.order):
            return None
        A, B = self.tail_bound
        return ceil(A if not self.order else A+(s+B)*self.order)

    def tail_precision(self, x):
        """Entire-tail order in pi digits; None means no truncation error here.

        Coefficient and input uncertainty are separate. At an exact center,
        a positive-order omitted tail vanishes; a zero-containing input
        ball is not an exact center and retains its uncertainty.
        """
        _, s = self._displacement(x)
        return self._tail_precision(s)

    def __call__(self, x, *, precision=None):
        """Enclose all represented inputs and coefficients, including the tail.

        Default output is capped by the known information. An explicit
        precision above the coefficient/input/tail bound raises PrecisionError.
        Exact source data permits output above the parent's default precision.
        """
        target = self.precision if precision is None else _integer(precision, "precision")
        h, s = self._displacement(x)
        value = _zero(self.field)
        for a in reversed(self.coefficients):
            value = _add(self.field, _multiply(self.field, value, h), a)
        available = value.precision if analytic._is_ball(value) else inf
        tail = self._tail_precision(s)
        if tail is not None:
            available = min(available, tail)
        if precision is not None and target > available:
            raise PrecisionError("requested precision exceeds coefficient, input or tail information")
        return self.field(value, absprec=min(target, available))

    def derivative(self, order=1):
        """Local derivative of nonnegative integer order, with a transformed tail."""
        _integer(order, "derivative order")
        if order < 0:
            raise ValueError("derivative order must be nonnegative")
        if order == 0:
            return self
        cs = tuple(_scale(self.field, a, factorial(n)//factorial(n-order))
                   for n, a in enumerate(self.coefficients) if n >= order)
        bound = None
        if self.tail_bound is not None:
            A, B = self.tail_bound
            bound = (A+B*order, B)  # Falling factorials are integral.
        return type(self)(self.field, cs, center=self.center, radius=self.radius,
                          tail_bound=bound, precision=self.precision)

    def antiderivative(self, constant=0):
        """Local primitive with F(center)=constant; not a measure integral.

        Division by n costs e*v_p(n) coefficient digits. The transformed
        infinite-tail bound includes every later denominator spike. Its
        slope is weakened conservatively while keeping the declared disc.
        """
        cs = (_coefficient(self.field, constant),)+tuple(
            _scale(self.field, a, Fraction(1, n+1)) for n, a in enumerate(self.coefficients))
        bound = None
        if self.tail_bound is not None:
            A, B = self.tail_bound
            delta = (self.radius+B)/2
            e = 1 if isinstance(self.field, PadicField) else self.field.ramification_index
            minimum = analytic._log_tail_bound(self.order+1, delta, self.field.p, e)
            bound = (A-B+minimum, B-delta)
        return type(self)(self.field, cs, center=self.center, radius=self.radius,
                          tail_bound=bound, precision=self.precision)

    def integral(self, constant=0):
        """Alias for antiderivative; use a measure object for a measure integral."""
        return self.antiderivative(constant)

    def restrict(self, radius):
        """Restrict to a smaller concentric closed disc; never enlarge the domain."""
        _integer(radius, "radius")
        if self.radius is not None and radius < self.radius:
            raise DomainError("restrict cannot enlarge the declared series disc")
        return type(self)(self.field, self.coefficients, center=self.center, radius=radius,
                          tail_bound=self.tail_bound, precision=self.precision)

    @classmethod
    def from_polynomial(cls, field, polynomial, *, center=0, radius=None, precision=None):
        """Expand a finite polynomial f(x) in powers of x-center, with zero tail.

        Accepts Polynomial or a coefficient iterable. Exact source coordinates
        stay exact, and finite coefficient balls retain their uncertainty.
        """
        _field(field)
        c = _center(field, center)
        source = polynomial.coefficients if isinstance(polynomial, Polynomial) else polynomial
        cs = []
        for a in reversed(tuple(_coefficient(field, a) for a in source)):
            out = [_zero(field) for _ in range(len(cs)+1)]
            for n, b in enumerate(cs):
                out[n] = _add(field, out[n], _multiply(field, b, c))
                out[n+1] = _add(field, out[n+1], b)
            out[0] = _add(field, out[0], a)
            cs = out
        return cls(field, cs, center=c, radius=radius, tail_bound=None, precision=precision)

    @classmethod
    def from_function(cls, function, field, *, center=0, order=20, radius=None, precision=None):
        """Certified series for exp, sin, cos, log1p or log (default log(p)=0).

        Exact centers only. The standard rational centers retain exact
        coefficients; exponential/trigonometric centers use Taylor balls
        with guard digits. Logarithmic coefficients of positive degree stay
        exact. Logarithms require a principal displacement disc. A finite
        prefix can yield less than the desired default output precision.
        """
        from .taylor import taylor
        _field(field)
        _integer(order, "series order")
        if order < 1:
            raise ValueError("series order must be positive")
        c = _center(field, center)
        target = field.precision if precision is None else _integer(precision, "precision")
        e = 1 if isinstance(field, PadicField) else field.ramification_index
        if radius is not None:
            _integer(radius, "radius")
        if function in (analytic.exp, analytic.sin, analytic.cos):
            minimum_radius = e//(field.p-1)+1
            work = max(target, minimum_radius)
            analytic._exp_domain(field(c, absprec=work))
            r = minimum_radius if radius is None else radius
            bound = (0, -Fraction(e, field.p-1))
        elif function in (analytic.log, analytic.log1p):
            if function is analytic.log1p and _v(field, c) < 1:
                raise DomainError("log1p series requires v_pi(center)>=1")
            scale = c if function is analytic.log else _add(field, c, _coefficient(field, 1))
            t = _v(field, scale)
            if t == inf:
                raise DomainError("log series requires a nonzero exact center")
            minimum_radius = t+1
            r = minimum_radius if radius is None else radius
            if r < minimum_radius:
                raise DomainError("logarithmic displacement requires radius>v_pi(scale)")
            delta = Fraction(r-t, 2)
            bound = (analytic._log_tail_bound(order, delta, field.p, e), -t-delta)
            work = max(target+max(t, 0), t+2, 3)
        else:
            raise TypeError("supported functions are exp, sin, cos, log1p and log")
        if r < minimum_radius:
            raise DomainError("series disc exceeds the supported analytic displacement domain")
        if _exact_zero(field, c) and function is not analytic.log:
            coefficients = taylor(function, 0, order=order).coefficients
        elif function is analytic.log and _exact_zero(
                field, _add(field, c, _coefficient(field, -1))):
            coefficients = taylor(function, 1, order=order).coefficients
        elif function in (analytic.log, analytic.log1p):
            # The exact center makes 1/(n*scale**n) exact. Only log(center)
            # needs a numerical ball; derivatives must not inherit an
            # artificial rounding of the center into all these coefficients.
            inverse = (1/scale if isinstance(field, PadicField) else
                       analytic._pair_power(field, scale, -1))
            power = inverse
            coefficients = [function(field(c, absprec=work))]
            for n in range(1, order):
                coefficients.append(_scale(field, power, Fraction((-1)**(n+1), n)))
                power = _multiply(field, power, inverse)
        else:
            coefficients = taylor(function, field(c, absprec=work), order=order).coefficients
        return cls(field, coefficients, center=c, radius=r, tail_bound=bound, precision=target)

    def _coerce(self, other):
        if not isinstance(other, AnalyticSeries):
            return type(self)(self.field, [other], center=self.center, radius=self.radius,
                              tail_bound=None, precision=self.precision)
        key = lambda f: (PadicField, f.p) if isinstance(f, PadicField) else f._key()
        if key(self.field) != key(other.field) or self.center != other.center:
            raise TypeError("series arithmetic requires the same field and exact center")
        return other

    def _context(self, other):
        radii = [r for r in (self.radius, other.radius) if r is not None]
        return max(radii) if radii else None, min(self.precision, other.precision)

    def _global_bound(self, slope):
        bounds = [_v(self.field, a)-slope*n for n, a in enumerate(self.coefficients)]
        if self.tail_bound is not None:
            A, B = self.tail_bound
            bounds.append(A+(B-slope)*self.order)
        return min(bounds, default=inf)

    def _combine(self, other, multiply):
        other = self._coerce(other)
        field = self.field
        radius, precision = self._context(other)
        tails = [s.tail_bound for s in (self, other) if s.tail_bound is not None]
        if not tails:
            count = max(0, self.order+other.order-1) if multiply else max(self.order, other.order)
            if multiply and (not self.order or not other.order):
                count = 0
            bound = None
        else:
            count = min(s.order for s in (self, other) if s.tail_bound is not None)
            B = min(t[1] for t in tails)
            A1, A2 = self._global_bound(B), other._global_bound(B)
            A = A1+A2 if multiply else min(A1, A2)
            if A == inf:  # A structural zero polynomial, never a zero ball.
                return type(self)(field, [], center=self.center, radius=radius,
                                  tail_bound=None, precision=precision)
            bound = (A, B)
        cs = []
        for n in range(count):
            if multiply:
                value = _zero(field)
                for j in range(n+1):
                    value = _add(field, value, _multiply(field, self[j], other[n-j]))
            else:
                value = _add(field, self[n], other[n])
            cs.append(value)
        return type(self)(field, cs, center=self.center, radius=radius,
                          tail_bound=bound, precision=precision)

    def __add__(self, other):
        return self._combine(other, False)

    __radd__ = __add__

    def __mul__(self, other):
        return self._combine(other, True)

    __rmul__ = __mul__

    def __neg__(self):
        return self*-1

    def __sub__(self, other):
        return self+-self._coerce(other)

    def __rsub__(self, other):
        return self._coerce(other)+-self
