# Language bindings

Language bindings are intentionally not implemented as separate numerical
libraries. They will wrap the stable C ABI declared in `include/scifort.h`.

Binding directories (Python implemented; others planned):

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
Outside the eleven generated families listed below, survival, logarithmic and
inverse-survival methods remain native-Fortran-only, as do special-function,
random-state, random-variate and likelihood entry points.

The generator also provides an explicit-shape C adapter for the native
assumed-shape lower-triangular solver, `scifort_solve_lower_f64`. It uses ordinary
Fortran array declarations and documented dimensions; see the
[linear-algebra configuration](cabi_config/linalg.json) and
[array interface contract](cabi_config/README.md).

Wrapper implementations must preserve the numerical results of the Fortran
core and must document language-specific choices involving shapes, missing
values, exceptions, and memory ownership.

The generated ABI also covers all eight PDF/log-PDF/CDF/SF/log-CDF/log-SF/PPF/ISF
methods for exponential, logistic, Laplace, Cauchy, Weibull, lognormal, gamma,
beta, Student t and chi-square. All eleven families also provide elementwise
C entry points for array-valued parameters; Python applies NumPy broadcasting.
Together with normal, these provide bulk entry points for all 56 distribution
array workloads in the SciPy timing comparison. See the
[family manifests and regeneration commands](cabi_config/README.md#benchmark-distribution-families).

A generated [NumPy/ctypes Python interface](python/README.md) now exposes these
eleven families, with scalar calls, same-shape arrays, explicit defaults and optional
caller-owned output memory. The focused development build and practical
Python-to-SciPy timing commands are documented there. R, MATLAB and Octave wrappers
remain planned.

The [discovery workflow](../tools/cabi_generator/README.md#api-inventory-and-draft-manifests)
reports public interfaces, existing C bindings and configured C/Python coverage.
The latest additions are available as `api.gamma`, `api.beta` and `api.t`, with
required shapes `a`, `a, b` and degrees of freedom `df`, respectively.

The [unified workflow](../tools/cabi_generator/README.md#unified-reviewed-interface-workflow)
now generates the dependency-ordered development build plan, aggregate C header
and [reviewed interface coverage](cabi_config/interface_coverage.md). Adding the
reviewed chi-square manifest did not require editing the development builder.
