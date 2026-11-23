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

Future Python wrappers can provide defaults and array conventions over this
ABI. They are not implemented by this initial generator.
