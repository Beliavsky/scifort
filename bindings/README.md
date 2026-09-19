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

The current version provides only the shared C ABI vertical slice for normal
PDF, CDF, and PPF. It is sufficient to test ABI design before committing to
package-specific build systems.

Wrapper implementations must preserve the numerical results of the Fortran
core and must document language-specific choices involving shapes, missing
values, exceptions, and memory ownership.
