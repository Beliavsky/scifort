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

The generated [Python interface](../python/README.md) supplies explicit defaults
and array conventions over this ABI. [python.json](python.json) configures it;
its generation and read-only check are separate from the C manifest batch.

## Benchmark distribution families

[manifests.json](manifests.json) lists all eleven configured modules. Regenerate
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

## Gamma, beta and Student t

The [gamma](gamma.json), [beta](beta.json) and [Student t](student_t.json)
manifests were reviewed from discovery drafts. Each adds five scalar and eight
bulk symbols (39 new C symbols total); existing PDF/CDF/PPF scalar definitions
remain unchanged in `src/c_api/api.f90`. All three appear in the C batch and
[Python policy](python.json), bringing coverage to 80 scalar/bulk method pairs
and Python methods across ten families. This is configured family coverage,
not a claim of full-project parity. The current Python timing benchmark still
covers its original 56 workloads across seven families.

Gamma shape `a` and beta shapes `a, b` must be finite and positive. Student t
`df` must be positive but **may be positive infinity**, which delegates to the
normal distribution. Bulk validation intentionally omits `df` from the finite
list. All three require finite location and finite positive scale. Python
requires the shapes and defaults only location/scale to 0/1; it exposes Student
t as `api.t`, preserving the native `t_*` and existing `scifort_t_*` names.

A report for these three families can be regenerated without writing wrappers:

```text
python tools/cabi_generator/discover.py src/stats/gamma.f90 src/stats/beta.f90 src/stats/student_t.f90 --abi-source src/c_api --manifest-batch bindings/cabi_config/manifests.json --python-policy bindings/cabi_config/python.json --output build/cabi-discovery/after.json
```

It writes JSON and a sibling Markdown report. The reviewed scope has 24 public
methods, 24 observed scalar C wrappers, 24 observed procedural bulk wrappers,
24 configured scalar symbols, 24 configured bulk symbols and 24 Python methods.
Observed bindings come from syntactic imported references, not name matching;
configured coverage does not prove that generated artifacts or binaries are current.
Use the separate generation drift checks and focused integration tests.

A wider `src/stats` scan is available locally in
`build/cabi-discovery/distributions.json` and `.md`. It attempted 162 files,
parsed 156, recorded five parser failures and explicitly skipped
`qmc_sobol_data.f90` under a 1,000,000-byte limit. One parsed default-public module
also remains unresolved. It therefore marks itself incomplete and exits nonzero.
The resolved scope contains 1,273 public exports, with 909 supported scalar-function signatures
and 906 supported elemental signatures; these include methods beyond the eight
probability methods and are not counts of distributions or a SciPy parity score.
The report records 80 configured scalar/bulk/Python method triples. Regenerate:

```text
python tools/cabi_generator/discover.py src/stats --progress --max-source-bytes 1000000 --abi-source src/c_api --manifest-batch bindings/cabi_config/manifests.json --python-policy bindings/cabi_config/python.json --output build/cabi-discovery/distributions.json
```

The five parser failures are `foldcauchy.f90` and `skewcauchy.f90` (fparser's
intrinsic ATAN argument-count handling), plus `power_transform_engine.f90`,
`probability_plot.f90` and `sample_distance.f90` (Fortran 2018 `error stop ...,
quiet=` syntax outside the configured parser grammar). `scifort_power_transform`
has default-public access, which the current public-export resolver deliberately
rejects. These are discovery limitations; numerical source has not been rewritten
to satisfy the parser. The three reviewed families have no such inventory errors.

## Unified regeneration and chi-square

Use the reviewed [pipeline policy](pipeline.json) to regenerate/check all related
outputs together:

```text
python tools/cabi_generator/pipeline.py bindings/cabi_config/pipeline.json
python tools/cabi_generator/pipeline.py bindings/cabi_config/pipeline.json --check
```

The batch now includes twelve modules: eleven distribution families and the
existing lower-triangular solver. The [chi-square manifest](chi2.json) adds eight
bulk symbols and five scalar symbols while preserving its three manual scalar
symbols. Its degrees of freedom must be finite and positive, unlike Student t's
permitted infinite-DF limit. Python exposes it as `api.chi2` with required `df`
and explicit location/scale defaults 0/1.

The [combined coverage table](interface_coverage.md) reports 88 scalar/bulk/Python
method triples plus one array adapter. Constants/excluded solver exports remain
visible in the JSON. The generated development build plan discovers the
chi-square -> gamma -> incomplete-gamma/normal dependency chain and includes only
Python-selected families, not the full manual C API object. The old 56-workload
Python timing suite retains its existing scope and reports.
