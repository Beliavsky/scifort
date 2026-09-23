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
status codes. Survival, logarithmic, inverse-survival, special-function,
random-state, random-variate, and likelihood entry points remain
native-Fortran-only for now.

Wrapper implementations must preserve the numerical results of the Fortran
core and must document language-specific choices involving shapes, missing
values, exceptions, and memory ownership.
