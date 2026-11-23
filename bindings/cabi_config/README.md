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
