# Validation of the 0.1 implementation

Validated on 26 October 2026 with Python 3.12.14. This records the current
implemented scope; it is not a formal verification or performance benchmark.

- **1020 pytest cases passed**, with no skipped cases in this run. They
  include **90 independent PARI/GP cases**. GP 2.15.4 was run from a locally
  extracted Ubuntu package, whose SHA-256 matched its signed package index.
  CI installs GP explicitly; without it, these optional comparisons skip.
- Exact oracles cover rational/Qp arithmetic, ideal-ball perturbations,
  inverses, valuations, cancellation and negative precision; quadratic
  fields include unramified and tame/wild Eisenstein presentations.
  SymPy quotient-polynomial and rational multiplication-matrix models check
  field arithmetic independently of the native coordinate engine.
- **98 quadratic-analysis cases** cover exp/log/log1p/sin/cos, integral
  extension-valued powers, Taylor enclosures, convergence boundaries,
  whole-ideal perturbations, roots of unity, log(p)/log(pi) branches,
  finite branch uncertainty and precision above the parent's default.
  Fixed exact rational matrix polynomials provide separate series oracles.
- **156 local-calculus cases** check analytic coefficient-ball families,
  arbitrary-order derivatives, primitives, repeated primitives, exact
  centers, finite coefficient/input uncertainty, explicit precision,
  arithmetic with distinct tail slopes and disc intersections. Independent
  rational functions and multiplication matrices supply exact oracles;
  admissible omitted monomials probe later p-power denominator spikes on
  small, unit and larger-than-unit discs. A custom coefficient bound is an
  explicit mathematical premise, not inferred from a finite prefix.
  Logarithmic derivatives retain exact positive-degree coefficients from
  exact centers; tests request 60 digits even when the stored constant
  logarithm has only nine. Existing finite input balls cannot gain digits.
- **42 Volkenborn cases** compare polynomial moments with symbolic
  Faulhaber limits, actual normalized finite Riemann sums, translation
  identities and independent rational multiplication-matrix models.
  Coefficient/affine-parameter perturbations and admissible omitted
  monomials test the whole-family enclosure, including later nonzero
  moments, Bernoulli denominator losses and constant charts.
  Seven further PARI cases check polynomial moments and exponential
  integral identities for p=2,3,5 and both Gaussian local orientations.
- **91 linear-ODE cases** check variable coefficients, translated centers,
  genuinely quadratic coefficients, exact/ball initial values and whole
  coefficient/input families. Independent symbolic integrating factors,
  rational multiplication-matrix exponentials and admissible omitted input
  monomials exercise the solution prefix and whole tail, including later
  p-power denominators, negative disc exponents and structural zero cases.
  Seven additional PARI cases check variable-coefficient closed solutions
  and their derivatives for p=2,3,5 and both Gaussian local orientations.
  The conservative factorial-majorant proof supplies existence, uniqueness
  within analytic power series, and the full solution-tail certificate.
- **67 differential-system cases** check noncommuting variable coefficients,
  rectangular initial matrices, coupled fundamental solutions, quadratic
  coordinate vectors, finite matrix/initial/argument uncertainty and empty
  input prefixes. Independent rational block exponentials and normalized
  Airy hypergeometric coefficients check both stored and omitted terms.
  The matrix valuation proof extends the scalar factorial majorant without
  assuming commutativity of coefficient matrices. Seven new PARI cases
  compare Airy Maclaurin sums and coupled Gaussian oscillator identities.
- Four of the PARI cases check Gaussian analytic images in both orientations
  over Q3(i) and Q2(i). They use PARI's scalar p-adic exp/sin/cos, scalar
  log of the norm and a separately truncated arctangent polynomial. This
  does not claim a native ramified POLMOD exp/log oracle. The other PARI
  cases cover KL-zeta conventions, transforms, harmonic operators, chosen
  local places, quadratic arithmetic and matrix solves.
- Nineteen additional PARI cases independently check values, first/second
  derivatives and normalized primitives at translated centers. Scalar
  tests use p=2,3,5; Gaussian tests cover both orientations over Q2(i) and
  Q3(i), using scalar transcendental identities and a fixed arctangent
  polynomial rather than a native ramified POLMOD analytic oracle.
- Bridge, quaternion and Clifford tests check scalar maps, products,
  conjugations, norms and the explicit even-Clifford isomorphism. Matrix
  checks include exact and whole-ball solves, unresolved pivots,
  determinants, inverses and condition exponents. Measure and harmonic
  tests cover the documented transform, wavelet and diffusion contracts.
- **612 Sphinx documentation examples passed**. Both HTML builds completed
  with warnings treated as errors. The portable manual contains **660
  validated internal anchors and 77 pre-rendered native MathML formulas**.
- Chromium 156 checked both HTML editions offline at **1440x1000 and
  390x844**. The portable edition had JavaScript disabled; English/Czech
  switching, chapter links, visible formulas and absence of page-width
  overflow passed. The complete multi-page edition also passed language
  navigation, the new chapter link and local formula rendering.
- Wheel and source archives built successfully and **twine check passed**.
  Their analytic modules, documentation, tests and example were compared
  with the verified source. All CI wheel smoke checks passed in a fresh
  environment with **no optional dependencies**, including both quadratic
  analytic engines, branches, powers, Taylor evaluation, local calculus,
  public Volkenborn integration, scalar ODEs and coupled differential systems.
- GitHub Actions for the implementation commit passed all six configured
  jobs: Python **3.10, 3.12 and 3.13**, PARI, documentation and packages.

Reproduce with the commands in README and RELEASING. Examples are
`python examples/quadratic_analysis.py`, `python examples/local_calculus.py`,
`python examples/volkenborn.py`, `python examples/linear_ode.py` and
`python examples/linear_system.py`.
The precision/branch arguments are in `docs/extension_analytic.rst`; full
coefficient-tail, derivative and primitive proofs are in `docs/calculus.rst`.
The uniform normalized-sum argument and integration tail enclosure are in
`docs/volkenborn.rst`.
The coefficient recurrence, factorial majorant, sufficient ODE convergence
disc and analytic uniqueness argument are in `docs/ode.rst`.
The matrix majorant, preserved product order and system contracts are in
`docs/ode_system.rst`.

Public extension square roots/Teichmueller lifts, extension-valued theta
and zeta, general higher-degree fields, noncommutative analytic functions,
series inversion/composition, singular systems and nonlinear ODE solvers
remain outside this step. Agreement of algebraic images or special values
does not establish a complex/p-adic zero correspondence.
