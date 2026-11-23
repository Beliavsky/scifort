# SciFort C ABI configuration

The reusable implementation is in [tools/cabi_generator](../../tools/cabi_generator/README.md).
This directory contains SciFort-specific interface and compatibility decisions.

`normal.json` exposes scalar and bulk `isf`, `logcdf`, `logpdf`, `logsf`, and
`sf` methods. Existing scalar and bulk `pdf`, `cdf`, and `ppf` wrappers remain
in `src/c_api/api.f90`, with their original symbols and behavior. They are
explicit exclusions in the generated coverage report, not missing methods.

The checked-in outputs are:

- [Fortran wrappers](../../src/c_api/generated/normal_generated.f90), compiled by FPM.
- [C header](../../include/scifort_generated_normal.h), included by `scifort.h`.
- [Coverage record](normal_coverage.json), accounting for all eight public exports
  and listing ten generated symbols.

Regenerate after changing the manifest or native declarations; review the diff
and run `--check` before committing. The manifest makes optional location and
scale arguments required at the C boundary. Bulk calls validate finite location
and finite positive scale before touching output. Their length/buffer/status
contract is documented in the generated header and generator guide.

`linalg.json` uses schema version 2 to wrap the existing assumed-shape
`solve_lower(l(:,:),b(:),x(:),status)` procedure. Its generated adapter declares
`l(n,n)`, `b(n)`, and `x(n)` and calls the native routine directly. It adds
`scifort_solve_lower_f64` without modifying the numerical implementation.

- [Explicit-shape Fortran adapter](../../src/c_api/generated/linalg_generated.f90).
- [C header](../../include/scifort_generated_linalg.h), included by `scifort.h`.
- [Coverage record](linalg_coverage.json), with explicit exclusions for other
  routines and native status constants.

The solver accepts a by-value C `int` order and contiguous column-major buffers;
all pointers must be valid and non-NULL. Orders below one return status 1 and
leave output unchanged. Otherwise, native status codes apply: 0 success and
5 singular/nonfinite diagonal, with the native routine's partial-output behavior.
No allocation, buffer copying, or pointer conversion occurs in this adapter.

Schema-version-1 elemental manifests retain their pointer-based default. New
version-2 elemental manifests default to explicit-shape bulk wrappers. The
normal manifest pins `array_style: "pointer"` to retain its existing NULL-buffer
contract; do not remove that policy when updating its schema version.

Future Python wrappers can provide defaults and array conventions over this
ABI. They are not implemented by this initial generator.

## Benchmark distribution families

[manifests.json](manifests.json) lists all eight configured modules. Regenerate
or check them together from the repository root:

```text
python tools/cabi_generator/generate.py bindings/cabi_config/manifests.json --batch
python tools/cabi_generator/generate.py bindings/cabi_config/manifests.json --batch --check
```

The [exponential](exponential.json), [logistic](logistic.json), [Laplace](laplace.json),
[Cauchy](cauchy.json), [Weibull](weibull.json), and [lognormal](lognormal.json)
manifests each account for all eight public probability methods: PDF, log-PDF,
CDF, SF, log-CDF, log-SF, PPF and ISF. Each emits five new scalar symbols and
eight bulk symbols, for 78 new symbols across these six families. PDF/CDF/PPF
scalar symbols remain in `src/c_api/api.f90`; `existing_scalar_symbol` records
them without emitting duplicate implementations or declarations.

Together with normal, all 56 distribution-array benchmark workloads now have
bulk entry points. New bulk adapters use explicit-shape Fortran arrays and require
valid non-NULL, nonoverlapping buffers even for zero-length calls. All use
`size_t` counts and return status 0 on success or 1 for invalid lengths/uniform
parameters, preserving output on call-level errors. Shapes (`c` for Weibull,
`s` for lognormal) and scale must be finite and positive; location must be finite.
Scalar calls retain native validation and NaN/endpoint behavior. The header for
each family documents the complete call contract and is included by `scifort.h`.

The focused family test compiles verbatim production scalar wrapper bodies
separately to avoid building unrelated distribution modules. It compares generated
bulk results to those scalar wrappers and direct native calls in checked and
optimized builds. This bounds validation to the changed ABI integration; it does
not replace full-project validation before release or merge.
