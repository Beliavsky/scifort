# Language bindings

Language bindings are intentionally not implemented as separate numerical
libraries. They will wrap the stable C ABI declared in `include/scifort.h`.

Planned layout:

```text
bindings/python/
bindings/R/
bindings/matlab/
bindings/octave/
```

The current development version provides scalar C ABI entry points for PDF
or PMF, CDF, and PPF across every implemented distribution. Normal PDF, CDF,
and PPF additionally provide bulk array entry points with explicit lengths and
status codes. Generated normal-distribution scalar and bulk entry points now
also cover SF, log-PDF, log-CDF, log-SF, and ISF; see the
[generator configuration and coverage](cabi_config/README.md). Existing PDF,
CDF, and PPF symbols retain their original implementations and behavior.
For other families, survival, logarithmic, and inverse-survival methods remain
native-Fortran-only, as do special-function, random-state, random-variate, and
likelihood entry points.

The generator also provides an explicit-shape C adapter for the native
assumed-shape lower-triangular solver, `scifort_solve_lower_f64`. It uses ordinary
Fortran array declarations and documented dimensions; see the
[linear-algebra configuration](cabi_config/linalg.json) and
[array interface contract](cabi_config/README.md).

Wrapper implementations must preserve the numerical results of the Fortran
core and must document language-specific choices involving shapes, missing
values, exceptions, and memory ownership.
