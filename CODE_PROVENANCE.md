# Code provenance

This file records the origin and licensing of numerical implementation code.
Every imported, translated, or substantially adapted routine must be entered
before it is merged.

## Current implementation inventory

All files listed below are original SciFort implementations under the MIT
License. Version 0.1.0 contains no translated third-party source code.

- `src/math.f90`: stable elementary helpers implemented from standard
  power-series and algebraic identities.
- `src/stats/normal.f90`: standard normal `erfc` identity, a
  Mills-ratio asymptotic expansion for log-CDF tails, and safeguarded
  bisection for quantiles.
- `src/stats/uniform.f90`: direct implementation from the mathematical
  definition.
- `src/stats/exponential.f90`: direct implementation with stable
  `expm1` and `log1p` helpers.
- `src/stats/laplace.f90`: direct piecewise implementation from the
  mathematical definition.
- `src/stats/logistic.f90`: direct implementation using stable
  logistic and softplus identities.
- `src/stats/cauchy.f90`: direct implementation using reciprocal-angle
  identities and cancellation-resistant logarithmic tails.
- `src/c_api/api.f90`: interoperable wrapper layer over the native
  Fortran API.

No item above is a source translation from SciPy, Fortran `stdlib`, Netlib, or
another numerical package.

## Required entry format for future imports

Add a row or section containing:

- SciFort file and procedure names;
- upstream project, file path, and URL;
- exact version, tag, or commit hash;
- original authors and copyright holders;
- original license and location of retained license text;
- whether the code was copied, translated, adapted, or independently
  reimplemented from a publication;
- meaningful modernization changes;
- tests used to establish parity.
