# Validation of the 0.1 implementation

Validated on 24 October 2026 with Python 3.12.14. This records the current
implemented scope; it is not a formal verification or performance benchmark.

- **624 pytest cases passed**, with no skipped cases in this run. They
  include **50 independent PARI/GP cases**. GP 2.15.4 was run from a locally
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
- Four of the PARI cases check Gaussian analytic images in both orientations
  over Q3(i) and Q2(i). They use PARI's scalar p-adic exp/sin/cos, scalar
  log of the norm and a separately truncated arctangent polynomial. This
  does not claim a native ramified POLMOD exp/log oracle. The other PARI
  cases cover KL-zeta conventions, transforms, harmonic operators, chosen
  local places, quadratic arithmetic and matrix solves.
- Bridge, quaternion and Clifford tests check scalar maps, products,
  conjugations, norms and the explicit even-Clifford isomorphism. Matrix
  checks include exact and whole-ball solves, unresolved pivots,
  determinants, inverses and condition exponents. Measure and harmonic
  tests cover the documented transform, wavelet and diffusion contracts.
- **476 Sphinx documentation examples passed**. Both HTML builds completed
  with warnings treated as errors. The portable manual contains **583
  validated internal anchors and 55 pre-rendered native MathML formulas**.
- Chromium 150 checked both HTML editions offline at **1440x1000 and
  390x844**. The portable edition had JavaScript disabled; English/Czech
  switching, chapter links, visible formulas and absence of page-width
  overflow passed. The complete multi-page edition also passed language
  navigation, the new chapter link and local formula rendering.
- Wheel and source archives built successfully and **twine check passed**.
  Their analytic modules, documentation, tests and example were compared
  with the verified source. All CI wheel smoke checks passed in a fresh
  environment with **no optional dependencies**, including both quadratic
  analytic engines, branches, powers and Taylor evaluation.
- GitHub Actions for the implementation commit passed all six configured
  jobs: Python **3.10, 3.12 and 3.13**, PARI, documentation and packages.

Reproduce with the commands in README and RELEASING. The new example is
`python examples/quadratic_analysis.py`; the full precision and branch
arguments are in `docs/extension_analytic.rst`.

Public extension square roots/Teichmueller lifts, extension-valued theta
and zeta, general higher-degree fields and noncommutative analytic functions
remain outside this step. Agreement of algebraic images or special values
does not establish a complex/p-adic zero correspondence.
