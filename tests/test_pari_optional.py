"""Independent optional PARI oracles for zeta, transforms and Taylor values.

Run where PARI/GP's gp command is installed. The normal exact-oracle tests
do not depend on it; skipping this file must never be reported as a pass
of the PARI cross-check.
"""
from fractions import Fraction
import shutil
import subprocess
import pytest
from padic_tools import AtomicMeasure, BallVladimirov, EisensteinQuadraticField, Gaussian, GaussianLocalPlace, PadicTestFunction, Qp, RationalAmiceMeasure, UnramifiedQuadraticField, exp, finite_fourier, log, taylor, zeta_p
from padic_tools import Matrix
from padic_tools import log1p, power
from padic_tools import AnalyticSeries, sin, cos
from padic_tools import solve_linear_ode


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p", [2, 3, 5])
def test_linear_ode_against_pari_closed_variable_coefficient_solution(p):
    script = (f"p={p};h=p^2+O(p^45);u=exp(h+h^2/2);\n"
              "print(lift(2+h+u));print(lift(1+(1+h)*u));quit\n")
    result = subprocess.run(["gp", "-q", "-f"], input=script, text=True,
                            capture_output=True, check=True, timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 2
    K, c = Qp(p, 3), Fraction(1, p)
    a = AnalyticSeries(K, [1, 1], center=c, tail_bound=None)
    b = AnalyticSeries(K, [-1, -3, -1], center=c, tail_bound=None)
    y = solve_linear_ode(a, 3, forcing=b, order=80)
    for model, expected in zip((y, y.derivative()), oracle):
        assert model(c+p*p, precision=21).contains(expected)


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,sign", [(2, -1), (2, 1), (3, -1), (3, 1)])
def test_quadratic_linear_ode_against_pari_scalar_exponential_identities(p, sign):
    # y=2+h+exp(i*(h+h^2)), with genuine Gaussian ODE coefficients.
    # Scalar exp/sin/cos provide the oracle, not native POLMOD analytics.
    script = (f"p={p};a=p^2+O(p^45);b=2*p^2+O(p^45);\n"
              "u=-b-2*a*b;v=a+a^2-b^2;E=exp(u)*cos(v);F=exp(u)*sin(v);\n"
              "print(lift(2+a+E));print(lift(b+F));\n"
              "print(lift(1-2*b*E-(1+2*a)*F));print(lift(-2*b*F+(1+2*a)*E));quit\n")
    result = subprocess.run(["gp", "-q", "-f"], input=script, text=True,
                            capture_output=True, check=True, timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 4
    def coordinates(a, b):
        return (a-sign*b, sign*b) if p == 2 else (a, sign*b)
    K = GaussianLocalPlace(Qp(p, 3), root_sign=sign).field
    a = AnalyticSeries(K, [coordinates(0, 1), coordinates(0, 2)], tail_bound=None)
    b = AnalyticSeries(K, [coordinates(1, -2), coordinates(0, -5), coordinates(0, -2)], tail_bound=None)
    y = solve_linear_ode(a, 3, forcing=b, order=80)
    point = coordinates(p*p, 2*p*p)
    for model, expected in zip((y, y.derivative()), zip(oracle[::2], oracle[1::2])):
        assert model(point, precision=21).contains(coordinates(*expected))
from padic_tools import Polynomial, volkenborn_integral


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p", [2, 3, 5])
def test_volkenborn_polynomial_and_exponential_against_pari(p):
    script = (
        f"p={p};f=1/p-2*(1+p*t)+3*(1+p*t)^3+(1+p*t)^4;\n"
        "print(sum(n=0,4,polcoef(f,n)*bernfrac(n)));\n"
        "z=p^2+O(p^45);q=z/(exp(z)-1);print(lift(q));print(lift(exp(z)*q));quit\n"
    )
    result = subprocess.run(["gp", "-q", "-f"], input=script, text=True,
                            capture_output=True, check=True, timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 3
    f = Polynomial([Fraction(1, p), -2, 0, 3, 1])
    assert volkenborn_integral(f, shift=1, scale=p) == oracle[0]
    s = AnalyticSeries.from_function(exp, Qp(p, 3), order=100, precision=21)
    assert volkenborn_integral(s, scale=p*p, precision=21).contains(oracle[1])
    assert volkenborn_integral(s, shift=p*p, scale=p*p, precision=21).contains(oracle[2])


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,sign", [(2, -1), (2, 1), (3, -1), (3, 1)])
def test_quadratic_volkenborn_exponential_against_pari_scalar_identities(p, sign):
    script = (
        f"p={p};X=p^2+O(p^45);Y=2*p^2+O(p^45);\n"
        "U=exp(X)*cos(Y)-1;V=exp(X)*sin(Y);D=U^2+V^2;\n"
        "A=(X*U+Y*V)/D;B=(Y*U-X*V)/D;print(lift(A));print(lift(B));\n"
        "C=exp(2*p^2+O(p^45))*cos(p^2+O(p^45));"
        "E=exp(2*p^2+O(p^45))*sin(p^2+O(p^45));\n"
        "print(lift(C*A-E*B));print(lift(C*B+E*A));quit\n"
    )
    result = subprocess.run(["gp", "-q", "-f"], input=script, text=True,
                            capture_output=True, check=True, timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 4
    def coordinates(a, b):
        return (a-sign*b, sign*b) if p == 2 else (a, sign*b)
    K = GaussianLocalPlace(Qp(p, 3), root_sign=sign).field
    s = AnalyticSeries.from_function(exp, K, order=100, precision=21)
    z, a = coordinates(p*p, 2*p*p), coordinates(2*p*p, p*p)
    assert volkenborn_integral(s, scale=z, precision=21).contains(coordinates(*oracle[:2]))
    assert volkenborn_integral(s, shift=a, scale=z, precision=21).contains(coordinates(*oracle[2:]))


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p", [2, 3, 5])
@pytest.mark.parametrize("function", [exp, sin, cos, log1p, log])
def test_analytic_series_calculus_against_pari_at_translated_centers(p, function):
    formulas = {
        exp: ("exp(x)", "exp(x)", "exp(x)", "exp(x)-exp(a)"),
        sin: ("sin(x)", "cos(x)", "-sin(x)", "cos(a)-cos(x)"),
        cos: ("cos(x)", "-sin(x)", "-cos(x)", "sin(x)-sin(a)"),
        log1p: ("log(1+x)", "1/(1+x)", "-1/(1+x)^2",
                "(1+x)*log(1+x)-x-((1+a)*log(1+a)-a)"),
        log: ("log(x)", "1/x", "-1/x^2", "x*log(x)-x-(a*log(a)-a)"),
    }
    script = f"p={p};a=p^2+O(p^45);x=a+p^3;\n"
    script += "".join(f"print(lift({f}));\n" for f in formulas[function])+"quit\n"
    result = subprocess.run(["gp", "-q", "-f"], input=script, text=True,
                            capture_output=True, check=True, timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 4
    series = AnalyticSeries.from_function(function, Qp(p, 3), center=p*p,
                                          order=100, precision=24)
    values = [series, series.derivative(), series.derivative(2), series.antiderivative()]
    for model, reference in zip(values, oracle):
        assert model(p*p+p**3, precision=18).contains(reference)


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,sign", [(2, -1), (2, 1), (3, -1), (3, 1)])
def test_quadratic_series_calculus_against_pari_gaussian_identities(p, sign):
    # Independent scalar Qp identities, not generic exp/log on a POLMOD.
    script = (
        f"p={p};a=p^2+O(p^45);b=2*p^2+O(p^45);x=a+p^3;y=b+p^3;\n"
        "emit(a,b)={print(lift(a));print(lift(b))};\n"
        "lg(a,b)=[log((1+a)^2+b^2)/2,sum(k=0,99,(-1)^k*(b/(1+a))^(2*k+1)/(2*k+1))];\n"
        "E=exp(x)*cos(y);F=exp(x)*sin(y);emit(E,F);\n"
        "emit(E-exp(a)*cos(b),F-exp(a)*sin(b));\n"
        "v=lg(x,y);w=lg(a,b);emit(v[1],v[2]);\n"
        "D=(1+x)^2+y^2;emit((1+x)/D,-y/D);\n"
        "emit((y^2-(1+x)^2)/D^2,2*(1+x)*y/D^2);\n"
        "emit((1+x)*v[1]-y*v[2]-x-((1+a)*w[1]-b*w[2]-a),"
        "(1+x)*v[2]+y*v[1]-y-((1+a)*w[2]+b*w[1]-b));quit\n"
    )
    result = subprocess.run(["gp", "-q", "-f"], input=script, text=True,
                            capture_output=True, check=True, timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 12
    pairs = [(a-sign*b, sign*b) if p == 2 else (a, sign*b)
             for a, b in zip(oracle[::2], oracle[1::2])]
    K = GaussianLocalPlace(Qp(p, 3), root_sign=sign).field
    def coordinates(a, b):
        return (a-sign*b, sign*b) if p == 2 else (a, sign*b)
    center, point = coordinates(p*p, 2*p*p), coordinates(p*p+p**3, 2*p*p+p**3)
    exponential = AnalyticSeries.from_function(exp, K, center=center, order=100, precision=25)
    logarithm = AnalyticSeries.from_function(log1p, K, center=center, order=100, precision=25)
    for d in (0, 1, 2):
        assert exponential.derivative(d)(point, precision=21).contains(pairs[0])
    assert exponential.antiderivative()(point, precision=21).contains(pairs[1])
    for d in (0, 1, 2):
        assert logarithm.derivative(d)(point, precision=21).contains(pairs[2+d])
    assert logarithm.antiderivative()(point, precision=21).contains(pairs[5])


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,sign", [(2, -1), (2, 1), (3, -1), (3, 1)])
def test_gaussian_analytic_functions_against_pari_scalar_identities(p, sign):
    # Do not call exp/log on a GP POLMOD: those are not a native ramified
    # analytic oracle. Use Qp exp/sin/cos, log(norm), and a separately fixed
    # rational arctangent polynomial. The omitted atan terms have v_p>300.
    script = (
        f"p={p};a=p^2+O(p^40);b=2*p^2+O(p^40);d=p+1;\n"
        "emit(a,b)={print(lift(a));print(lift(b))};\n"
        "emit(exp(a)*cos(b),exp(a)*sin(b));\n"
        "A=1+a;B=b;u=log(A^2+B^2)/2;"
        "v=sum(k=0,79,(-1)^k*(B/A)^(2*k+1)/(2*k+1));\n"
        "emit(u,v);emit(exp(u/d)*cos(v/d),exp(u/d)*sin(v/d));\n"
        "emit(exp(a)*cos(b+p^4),exp(a)*sin(b+p^4));\n"
        "B=b+p^4;u=log(A^2+B^2)/2;"
        "v=sum(k=0,79,(-1)^k*(B/A)^(2*k+1)/(2*k+1));emit(u,v);quit\n"
    )
    result = subprocess.run(["gp", "-q", "-f"], input=script, text=True,
                            capture_output=True, check=True, timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 10
    pairs = [(a-sign*b, sign*b) if p == 2 else (a, sign*b)
             for a, b in zip(oracle[::2], oracle[1::2])]
    place = GaussianLocalPlace(Qp(p, 3), root_sign=sign)
    x = place(Gaussian(p*p, 2*p*p), absprec=21)
    h = place(Gaussian(0, p**4), absprec=23)
    values = [exp(x), log1p(x), power(1+x, Fraction(1, p+1)),
              taylor(exp, x, order=30)(x+h),
              taylor(log, 1+x, order=30)(1+x+h)]
    assert log(1+x).contains(pairs[1])
    for value, pair in zip(values, pairs):
        assert value.precision == 21
        assert value.contains(pair)


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p", [2, 3, 5, 7])
def test_qp_matrix_inverse_solve_and_determinant_against_pari(p):
    script = (
        f"p={p};A=[p^2,1,1/p;0,p,2;1,0,1+p]+matrix(3,3,i,j,O(p^40));\n"
        "b=[1/p,2,3]~+vector(3,i,O(p^40))~;J=A^-1;x=J*b;\n"
        "print(lift(matdet(A)));for(i=1,3,for(j=1,3,print(lift(J[i,j]))));\n"
        "for(i=1,3,print(lift(x[i])));quit\n"
    )
    result = subprocess.run(["gp", "-q", "-f"], input=script, text=True,
                            capture_output=True, check=True, timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 13
    A = Matrix(Qp(p, 4), [[p*p, 1, Fraction(1, p)], [0, p, 2], [1, 0, 1+p]], absprec=17)
    assert A.determinant().contains(oracle[0])
    assert A.inverse().contains([oracle[1+i*3:4+i*3] for i in range(3)])
    assert all(x.contains(y) for x, y in zip(A.solve([Fraction(1, p), 2, 3]), oracle[10:]))


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,sign", [(2, -1), (2, 1), (3, -1), (3, 1)])
def test_quadratic_matrix_inverse_and_solve_against_pari(p, sign):
    polynomial, root = ("t^2-2*t+2", "(t-1)") if p == 2 else ("t^2+1", "t")
    script = (
        f"p={p};s={sign};r={root};f={polynomial};\n"
        "g(a,b)=Mod((a+O(p^40))+s*(b+O(p^40))*r,f);\n"
        "A=[g(1,1),g(2,-1);g(3,1),g(1,2)];b=[g(1,-1),g(2,1)]~;J=A^-1;x=J*b;d=matdet(A);\n"
        "print(lift(polcoef(lift(d),0)));print(lift(polcoef(lift(d),1)));\n"
        "for(i=1,2,for(j=1,2,print(lift(polcoef(lift(J[i,j]),0)));print(lift(polcoef(lift(J[i,j]),1)))));\n"
        "for(i=1,2,print(lift(polcoef(lift(x[i]),0)));print(lift(polcoef(lift(x[i]),1))));quit\n"
    )
    result = subprocess.run(["gp", "-q", "-f"], input=script, text=True,
                            capture_output=True, check=True, timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 14
    pairs = list(zip(oracle[::2], oracle[1::2]))
    place = GaussianLocalPlace(Qp(p, 4), root_sign=sign)
    source = Matrix(Gaussian, [[Gaussian(1, 1), Gaussian(2, -1)], [Gaussian(3, 1), Gaussian(1, 2)]])
    A = source.map_coefficients(place.field, lambda c: place(c, absprec=17))
    assert A.determinant().contains(pairs[0])
    assert A.inverse().contains([pairs[1:3], pairs[3:5]])
    rhs = [place(Gaussian(1, -1), absprec=17), place(Gaussian(2, 1), absprec=17)]
    assert all(x.contains(y) for x, y in zip(A.solve(rhs), pairs[5:]))


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,sign", [(5, 1), (5, -1), (13, 1), (13, -1)])
def test_split_gaussian_place_and_guard_digits_against_pari(p, sign):
    script = (
        f"p={p};r=sqrt(-1+O(p^60));if(lift(r)%p>p/2,r=-r);\n"
        f"z=1/p^6-2*{sign}*r/p^5;\n"
        f"print(lift(z));print(lift(1/z));print(lift({sign}*r));quit\n"
    )
    result = subprocess.run(["gp", "-q", "-f"], input=script, text=True,
                            capture_output=True, check=True, timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 3
    place = GaussianLocalPlace(Qp(p, 4), root_sign=sign)
    image = place(Gaussian(Fraction(1, p**6), Fraction(-2, p**5)), absprec=11)
    assert image.precision == 11 and image.contains(oracle[0])
    assert image.inverse().contains(oracle[1])
    assert place.generator_image(absprec=13).contains(oracle[2])


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,sign", [(2, 1), (2, -1), (3, 1), (3, -1)])
def test_nonsplit_gaussian_places_against_pari_quotients(p, sign):
    polynomial, root = ("t^2-2*t+2", "(t-1)") if p == 2 else ("t^2+1", "t")
    script = (
        f"p={p};z=Mod((1/p^3+O(p^40))+{sign}*(2/p^2+O(p^40))*{root},{polynomial});\n"
        "print(lift(polcoef(lift(z),0)));print(lift(polcoef(lift(z),1)));\n"
        "w=1/z;print(lift(polcoef(lift(w),0)));print(lift(polcoef(lift(w),1)));\n"
        "print(lift(norm(z)));print(lift(trace(z)));quit\n"
    )
    result = subprocess.run(["gp", "-q", "-f"], input=script, text=True,
                            capture_output=True, check=True, timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 6
    place = GaussianLocalPlace(Qp(p, 4), root_sign=sign)
    image = place(Gaussian(Fraction(1, p**3), Fraction(2, p**2)), absprec=19)
    assert image.precision == 19 and image.contains(tuple(oracle[:2]))
    assert image.inverse().contains(tuple(oracle[2:4]))
    assert image.norm().contains(oracle[4]) and image.trace().contains(oracle[5])


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,s", [(3,-1),(5,-3),(7,-5),(3,2),(5,3),(7,4),(2,3)])
def test_principal_zeta_against_pari(p,s):
    script=f"print(lift(zeta({s}+O({p}^24))))\nquit\n"
    result=subprocess.run(["gp","-q","-f"],input=script,text=True,capture_output=True,check=True,timeout=20)
    oracle=Fraction(result.stdout.strip())
    assert zeta_p(s,Qp(p,12),branch=0).contains(oracle)


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p", [2,3,5])
def test_transforms_and_taylor_against_pari(p):
    K = Qp(p,12)
    point = 3 if p == 2 else 2
    denominator = 3 if p == 2 else 2
    omega = "-1" if p == 2 else "teichmuller(x)"
    script = (
        f"x={point}+O({p}^32); w={omega}; s=1/{denominator};\n"
        f"print(lift(exp(-{p*p}+O({p}^32))/{p}+2*exp(-2*{p*p}+O({p}^32))));\n"
        f"print(lift(w*exp(s*log(x/w))/{p}));\n"
        f"print(lift(exp({p*p+p**3}+O({p}^32))));\n"
        f"print(lift(log({p*p+p**3}+O({p}^32))));\nquit\n"
    )
    result = subprocess.run(["gp","-q","-f"],input=script,text=True,
                            capture_output=True,check=True,timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 4
    values = [
        AtomicMeasure(K,[(1,Fraction(1,p)),(2,2)]).laplace(p*p),
        AtomicMeasure(K,[(point,Fraction(1,p))]).mellin(Fraction(1,denominator),branch=1),
        taylor(exp,K(p*p),order=30)(K(p*p+p**3)),
        taylor(log,K(p*p),order=30)(K(p*p+p**3)),
    ]
    for value,reference in zip(values,oracle):
        assert value.precision >= 11
        assert value.contains(reference)


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,c", [(2,3),(3,2),(5,2)])
def test_rational_measure_laplace_and_moments_against_pari(p,c):
    mu = RationalAmiceMeasure(Qp(p,14),[Fraction(1,p)],[c,1])
    script = (
        f"print(lift(1/({p}*({c}+exp(-{p*p}+O({p}^32))-1))));\n"
        f"f=1/({p}*({c}+exp(-x+O(x^10))-1));\n"
        "for(n=0,8,print(polcoef(f,n)))\nquit\n"
    )
    result = subprocess.run(["gp","-q","-f"],input=script,text=True,
                            capture_output=True,check=True,timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 10
    value = mu.laplace(p*p)
    assert value.precision >= 14
    assert value.contains(oracle[0])
    assert mu.laplace_series(9).coefficients == tuple(oracle[1:])
    from math import factorial
    for n,c in enumerate(oracle[1:]):
        assert mu.moment(n) == (-1)**n*factorial(n)*c


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p", [3,5,7])
def test_nonatomic_mellin_against_pari_kl_identity(p):
    # A(T)=1/(2+T), restricted to units, has the following analytic Mellin
    # identity on branch p-2. PARI evaluates its independent KL algorithm.
    mu = RationalAmiceMeasure(Qp(p,14),1,[2,1]).restrict_units()
    script = (
        f"s=1/2+O({p}^32);\n"
        f"print(lift((exp((s+1)*log(2+O({p}^32)))-1)*zeta(-s)));\nquit\n"
    )
    result = subprocess.run(["gp","-q","-f"],input=script,text=True,
                            capture_output=True,check=True,timeout=20)
    assert not result.stderr, result.stderr
    value = mu.mellin(Fraction(1,2),branch=p-2)
    assert value.precision == 14
    assert value.contains(Fraction(result.stdout.strip()))


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,n", [(2,4),(3,9),(5,5)])
def test_exact_fourier_complex_place_against_pari(p,n):
    mp = pytest.importorskip("mpmath")
    script = (
        "default(realprecision,90);\n"
        f"N={n};\n"
        "for(k=0,N-1,z=sum(j=0,N-1,(j^2-2*j+3)*exp(-2*Pi*I*j*k/N))/N;"
        "print(real(z));print(imag(z)))\nquit\n"
    )
    result = subprocess.run(["gp","-q","-f"],input=script,text=True,
                            capture_output=True,check=True,timeout=20)
    assert not result.stderr, result.stderr
    oracle = result.stdout.splitlines()
    assert len(oracle) == 2*n
    ctx = mp.mp.clone(); ctx.dps = 85
    values = finite_fourier([j*j-2*j+3 for j in range(n)],p)
    for k,value in enumerate(values):
        reference = ctx.mpc(ctx.mpf(oracle[2*k].replace(" ","")),
                            ctx.mpf(oracle[2*k+1].replace(" ","")))
        assert abs(value.complex(80)-reference) < ctx.mpf("1e-74")


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,a,b,order", [(2,0,2,1),(3,-1,1,2),(5,0,1,1)])
def test_regional_vladimirov_against_independent_pari_kernel(p,a,b,order):
    n = p**(b-a)
    script = (
        f"p={p};a={a};b={b};s={order};N={n};C=(p^s-1)/(1-p^(-1-s));\n"
        "M=matrix(N,N,j,k,if(j==k,0,-C*p^(-b)*p^((1+s)*(a+valuation(j-k,p)))));\n"
        "for(j=1,N,M[j,j]=-sum(k=1,N,M[j,k]));\n"
        "v=vector(N,j,(j-1)^2+1)~;w=M*v;for(j=1,N,print(w[j]));\nquit\n"
    )
    result = subprocess.run(["gp","-q","-f"],input=script,text=True,
                            capture_output=True,check=True,timeout=20)
    assert not result.stderr, result.stderr
    oracle = tuple(Fraction(line) for line in result.stdout.splitlines())
    assert len(oracle) == n
    op = BallVladimirov(p,order,support_exponent=a,constancy_exponent=b)
    f = PadicTestFunction(p,[j*j+1 for j in range(n)],support_exponent=a,constancy_exponent=b)
    assert tuple(value.as_fraction() for value in op.apply(f).values) == oracle


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,c,b", [(2,1,1),(3,1,0),(5,2,0),(7,1,0)])
def test_unramified_quadratic_against_pari_polymod(p,c,b):
    # PARI inverts a polynomial modulo f over its own Qp coefficient field;
    # it independently computes the norm and trace of that polymod.
    script = (
        f"f=t^2+{b}*t+{c};z=Mod((1/{p}+O({p}^32))+(3+O({p}^32))*t,f);\n"
        "w=lift(z^-1);print(lift(polcoef(w,0)));print(lift(polcoef(w,1)));\n"
        "print(lift(norm(z)));print(lift(trace(z)));\n"
        f"w=lift(z*Mod(1+{p}*t,f));print(lift(polcoef(w,0)));print(lift(polcoef(w,1)));\n"
        f"print(valuation(lift(z),{p}));quit\n"
    )
    result = subprocess.run(["gp","-q","-f"],input=script,text=True,
                            capture_output=True,check=True,timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 7
    L = UnramifiedQuadraticField(Qp(p,12),[c,b,1])
    x = L(Fraction(1,p),3)
    assert x.inverse().contains(oracle[:2])
    assert x.norm().contains(oracle[2])
    assert x.trace().contains(oracle[3])
    assert (x*L(1,p)).contains(oracle[4:6])
    assert x.valuation() == oracle[6]


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,c,b", [(2,2,-2),(2,-2,0),(3,-3,0),(5,5,5)])
def test_eisenstein_quadratic_against_pari_polymod_and_norm_valuation(p,c,b):
    points = [(0,1),(Fraction(1,p),3),(0,Fraction(1,p)),(p,p)]
    script = f"f=t^2+({b})*t+({c});\n"
    for a,d in points:
        script += (
            f"z=Mod(({a}+O({p}^40))+({d}+O({p}^40))*t,f);\n"
            "w=lift(z^-1);print(lift(polcoef(w,0)));print(lift(polcoef(w,1)));\n"
            "print(lift(norm(z)));print(lift(trace(z)));\n"
            f"w=lift(z*Mod(1+{p}*t,f));print(lift(polcoef(w,0)));print(lift(polcoef(w,1)));\n"
            f"print(valuation(norm(z),{p}));\n"
        )
    result = subprocess.run(["gp","-q","-f"],input=script+"quit\n",text=True,
                            capture_output=True,check=True,timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 7*len(points)
    L = EisensteinQuadraticField(Qp(p,12),[c,b,1])
    for j,(a,d) in enumerate(points):
        x, reference = L(a,d), oracle[7*j:7*j+7]
        assert x.inverse().contains(reference[:2])
        assert x.norm().contains(reference[2])
        assert x.trace().contains(reference[3])
        assert (x*L(1,p)).contains(reference[4:6])
        assert x.valuation() == reference[6]
