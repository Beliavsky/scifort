# SciFort Python interface

This development interface exposes PDF, log-PDF, CDF, SF, log-CDF, log-SF, PPF,
and ISF for `normal`, `exponential`, `logistic`, `laplace`, `cauchy`, `weibull`,
`lognormal`, `gamma`, `beta`, Student `t`, and `chi2`. It uses NumPy arrays and the stable
binary64 C calling convention. Eight matrix methods are available under `api.linalg`,
and multivariate normal/Student-t PDF/log-PDF under `api.multivariate_normal`
and `api.multivariate_t`. `api.dirichlet` provides PDF, log-PDF, entropy, mean,
variance and covariance. Multinomial and multivariate hypergeometric methods are
available under `api.multinomial` and `api.multivariate_hypergeom`.
`api.dirichlet_multinomial` provides PMF, log-PMF, mean, variance and covariance.
Its numerical calculations run in Fortran. Python 3.10+ and NumPy 1.24+ are
required; SciPy is needed only for comparisons/tests. The generated package does
not depend on fparser at runtime.

## Try it from the repository

In CMD, with gfortran on PATH and NumPy installed:

```text
python bindings\python\example.py
```

The example builds a cached, optimized, eleven-family development shared library,
loads it, and demonstrates scalar calls, arrays, required shapes and reusable
output memory. Compilation is outside all timings.

For installation into your chosen Python environment:

```text
python -m pip install .\bindings\python
python bindings\python\build_library.py --directory build\python-bindings\release
```

Installation packages the Python layer; it does not compile or bundle a shared
library or Fortran runtime. Build and installation are separate steps. To load
the development library from the repository root:

```python
import json
from pathlib import Path
import numpy as np
from scifort import Scifort

metadata = json.loads(Path("build/python-bindings/release/build.json").read_text())
api = Scifort(metadata["library"], dll_directories=metadata["dll_directories"])
print(api.normal.cdf(np.array([-1., 0., 1.])))
print(api.normal.ppf(.975))
print(api.lognormal.sf([.1, 1., 10.], s=.8))
out = np.empty(3)
api.normal.logpdf([-1., 0., 1.], loc=.5, scale=2., out=out)
```

`Scifort` requires an explicit library path and validates all configured symbols
when loading it. On Windows, `dll_directories` supplies directories containing
Fortran runtime DLLs; the development builder records the compiler directory.
Keep that runtime installed. The Python layer can also load a full SciFort shared
library exporting these symbols; FPM's ordinary library build is not itself a
Python extension or a shared-library packaging step.

## Distribution contract

- All-scalar arguments return a Python float. Zero-dimensional NumPy arguments are scalar.
- Array-like real numeric input returns a native-endian float64 ndarray with the
  same shape. Lists, integer/float32 arrays, noncontiguous views, unaligned arrays
  and nonnative-endian arrays are converted when needed. Complex, string and
  object arrays are rejected. Empty arrays are supported and still validate
  uniform parameters.
- Location and scale default explicitly to 0 and 1 in the Python configuration.
  Weibull `c`, lognormal `s`, gamma `a`, beta `a, b`, and Student t/chi-square
  `df` are required. All eleven configured families accept array-valued parameters
  with NumPy broadcasting: `api.normal.cdf(0., loc=[-1., 0., 1.])` and
  `api.gamma.cdf([[.5], [2.]], a=[1., 2., 3.])` return shapes `(3,)` and `(2, 3)`.
  If any parameter is an array, every parameter is evaluated per element;
  invalid parameter elements return the native result (usually NaN) rather than
  raising a uniform-parameter error. Empty broadcasts evaluate no elements.
  Student t allows positive infinite `df` per element and returns the normal
  limit; chi-square requires finite positive `df`. Other shape parameters must
  also be finite and positive. Shapes, location and scale may broadcast together.
- Optional `out=` accepts a writable, aligned, C-contiguous, native float64 NumPy
  array of exactly the input or broadcast result shape. It is returned by identity. Wrong dtype,
  shape, layout, read-only storage and possible overlap with any input or parameter are
  rejected before the C call. All-scalar calls do not accept `out`.
- Scalar parameter errors retain native NaN results. Bulk call-level errors
  (invalid uniform parameters or rejected lengths) raise `ValueError`; output
  remains unchanged. These uniform-parameter checks apply when all parameters
  are scalar. Invalid individual points/probabilities retain native
  NaN/infinity/support behavior. An unknown ABI status raises `RuntimeError`.
- No numerical method is reimplemented in Python. Native behavior and intentional
  SciPy differences remain applicable. In particular SciPy 1.15.3 returns NaN
  for gamma/Weibull PDF/log-PDF at positive infinity, whereas SciFort returns the
  mathematical limits 0 and negative infinity.

## Matrix procedures

`api.linalg` exposes eight existing C adapters. Each returns an
`ArrayResult(status, outputs)`, exported from `scifort`. `status` is the native
integer status; `outputs` is a dictionary with the following keys:

| Method | Inputs | Output keys |
| --- | --- | --- |
| `cholesky_lower` | A(n,n) | `l` |
| `symmetric_eigen_jacobi` | A(n,n) | `eigenvalues`, `eigenvectors` |
| `psd_decompose` | A(n,n), `allow_singular=True` | `eigenvalues`, `eigenvectors`, `rank`, `log_pdet`, `cutoff` |
| `solve_lower` | L(n,n), b(n) | `x` |
| `spd_factor` | A(n,n) | `l`, `logdet` |
| `invert_lower` | L(n,n) | `linv` |
| `solve_lower_matrix` | L(n,n), B(n,k) | `x` |
| `solve_spd_matrix` | L(n,n), B(n,k) | `x` |

n >= 1; k >= 0. Eigenvectors and RHS vectors occupy columns. The SPD solve uses
a supplied lower Cholesky factor and solves `(L L^T) X = B`. Each call handles one
matrix, with exact shape relationships and no matrix-batch broadcasting.
See [linalg shape tables](../../docs/linalg_shapes.md) and
[matrix-helper shape tables](../../docs/matrix_helper_shapes.md).

```python
a = np.array([[4., 1.], [1., 3.]])
factorization = api.linalg.cholesky_lower(a)
if factorization.status:
    raise RuntimeError(f"Cholesky failed: {factorization.status}")
lower = factorization.outputs["l"]
solution = api.linalg.solve_spd_matrix(lower, np.array([[1.], [2.]]))
if solution.status:
    raise RuntimeError(f"Solve failed: {solution.status}")
print(solution.outputs["x"])
```

Inputs accept real numeric array-like data. Row-major, noncontiguous, unaligned,
integer/float32 and nonnative-endian arrays are converted to aligned,
Fortran-contiguous binary64 buffers as needed. Read-only inputs are accepted.
Complex/string/object inputs are rejected. Input buffers must not overlap each
other, matching the existing C adapter contracts; possible overlap is checked
against original arrays even when a conversion could make a copy.

By default array outputs are newly allocated Fortran-contiguous float64 arrays.
Pass `out={"l": buffer}` or another dictionary of array-output names to reuse
buffers; supplying only some outputs is permitted. Every supplied buffer must
have the exact shape, native-endian float64 dtype, Fortran-contiguous storage,
alignment and write permission. Outputs must not overlap inputs or other outputs.
The returned dictionary contains the supplied arrays by identity. Scalar outputs
are returned as Python float/int values and do not accept caller buffers.

Shape/type/layout/overlap errors raise before the C call and leave caller buffers
untouched. Native numerical failures return a nonzero `status` with the native
zeroed/partial/complete outputs, including changes to caller buffers. Always
check status before using a factorization or solution. All C statuses are returned
unchanged; the first four methods use `linalg_status_*`, while the last four use
`matrix_status_*`. See the [C status and failure contracts](../cabi_config/README.md).

`psd_decompose` explicitly defaults to allowing singular matrices, matching its
native procedure. `allow_singular` must be a Python bool or NumPy boolean scalar;
numeric flag values are rejected in Python. Allowed rank deficiency returns
status 0; disallowed rank deficiency returns status 5 with computed rank/cutoff,
complete eigenpairs and log_pdet=0. No Fortran logical representation crosses C.
Empty RHS batches return shape (n,0); helper routines skip diagonal evaluation
for those batches, even for otherwise-invalid diagonal values. Non-NULL scratch
buffers are supplied for all zero-size arrays.

The development library now includes both matrix modules/adapters and their
module dependencies through its generated build plan. Rebuild the library when
updating to this Python interface; older development DLLs may lack matrix symbols.

## Multivariate normal densities

```python
api.multivariate_normal.pdf([0., 0.])
api.multivariate_normal.logpdf([.5, -.2], mean=[0., 1.], cov=[[2., .3], [.3, 1.]])
api.multivariate_normal.pdf([1., 0.], cov=[[4., 0.], [0., 0.]], allow_singular=True)
```

Both methods accept one point `x(n)`, `mean(n)` and `cov(n,n)`, with n >= 1,
and return a Python float. `mean=None` explicitly supplies a zero vector;
`cov=None` explicitly supplies an identity matrix. `allow_singular=False` is
the default and accepts Python/NumPy boolean scalars only. These defaults are
reviewed Python policy; every input is required at the C boundary.

The matrix input conversion and nonoverlap rules above apply. Point and mean
must be vectors; covariance must be a square matrix. Scalar covariance, point
batches and NumPy broadcasting are not supported by these two methods. No `out`
argument is provided for scalar results. Shape/type/overlap errors raise before C.
See the [generated native shape tables](../../docs/multivariate_normal_shapes.md).

Numerical behavior follows the existing native procedures. Invalid covariance
or nonfinite mean returns NaN, as does singular covariance unless
`allow_singular=True`. A point outside a permitted singular covariance's support
has PDF 0 and log-PDF negative infinity. The lower covariance triangle is
authoritative, but **all entries must be finite**, including the upper triangle.
Native NaN/infinity handling and numerical limitations are retained; this wrapper
does not imply complete SciPy edge-case parity. Numerical NaN does not raise.

The C symbols `scifort_multivariate_normal_pdf_f64` and
`scifort_multivariate_normal_logpdf_f64` take dimension, point, mean, column-major
covariance, integer singularity flag, scalar output pointer and status pointer.
Status 1 rejects n < 1 or a flag other than 0/1, preserving the scalar output.
Status 0 means the call was valid, even when its numerical result is NaN.
All buffers/pointers must be valid, non-NULL, aligned and nonoverlapping.
Rebuild the development library to obtain the two new symbols.

## Multivariate Student-t densities

```python
api.multivariate_t.pdf([0., 0.])
api.multivariate_t.logpdf([.5, -.2], loc=[0., 1.], shape=[[2., .3], [.3, 1.]], df=5.)
api.multivariate_t.pdf([1., 0.], shape=[[4., 0.], [0., 0.]], df=5., allow_singular=True)
```

Both methods accept one point `x(n)`, `loc(n)` and `shape(n,n)`, with n >= 1,
and return a Python float. `loc=None` supplies zero location; `shape=None`
supplies identity shape; `df=1` is the default. Degrees of freedom must be a real
scalar greater than zero; positive infinity takes the native normal limit.
**Shape is not covariance**: for finite df > 2 the covariance is
`df/(df-2) * shape`. At infinite df the normal covariance is the supplied shape.

PDF defaults to `allow_singular=False`; a rejected singular matrix returns NaN.
Set `allow_singular=True` (Python/NumPy boolean scalar) to permit it. Log-PDF
always permits positive-semidefinite singular shape and has **no**
`allow_singular` argument, matching the existing native API. For finite df the
native density formula uses projected distances without rejecting points in
null-space directions. At infinite df it delegates to the native normal support
check, so points outside singular support return PDF 0/log-PDF negative infinity.
Both behaviors are retained by the wrappers.

All shape entries must be finite, including the unused upper triangle; the lower
triangle is authoritative. Invalid df, nonfinite location or invalid shape
returns native NaN without raising. Native numerical limitations are retained.
Input conversion, exact shape checking and nonoverlap rules match the multivariate
normal methods. No point batching, parameter broadcasting, scalar shape or `out`
argument is provided. See the [native shape tables](../../docs/multivariate_t_shapes.md).

The C symbols `scifort_multivariate_t_pdf_f64` and
`scifort_multivariate_t_logpdf_f64` require every native input. PDF takes the
integer 0/1 singularity flag; log-PDF has no flag. Status 1 rejects n < 1 or an
invalid PDF flag and preserves the scalar output. Otherwise status is 0, including
for numerical NaN. C buffer validity/storage rules match the normal adapters;
rebuild the development library to obtain these new symbols.

## Dirichlet distributions

```python
api.dirichlet.pdf([.2, .3, .5], alpha=[2., 3., 4.])
api.dirichlet.logpdf([.2, .3], alpha=[2., 3., 4.])  # infer final coordinate .5
api.dirichlet.entropy([2., 3., 4.])
api.dirichlet.mean([2., 3., 4.])
api.dirichlet.var([2., 3., 4.])
covariance = np.empty((3, 3), dtype=np.float64, order="F")
api.dirichlet.cov([2., 3., 4.], out=covariance)
```

`alpha(k)` contains finite positive concentrations, k >= 1. Density methods
accept a point `x(m)`, where m >= 1 and m is k or k-1. For m=k-1 the final
coordinate is `1-sum(x)`; this shortened representation requires k >= 2.
The full point must have finite coordinates in [0,1] and sum to 1 within the
unchanged native absolute tolerance 1e-9. Length/type/overlap errors raise before
C. Numerical parameter errors retain native NaN results without raising.
At a zero coordinate whose concentration is below 1, native PDF/log-PDF return
NaN; native boundary behavior is preserved. k=1 accepts only the full point [1].

PDF, log-PDF and entropy return Python floats. Mean and variance return arrays
of shape (k); covariance returns shape (k,k), with native-endian float64 and
Fortran-contiguous storage. The array methods accept `out=ndarray`, return that
same object, and use the matrix output rules above: exact shape, native float64,
alignment, write permission and no overlap with inputs. Scalar methods have no
`out` argument. Invalid concentrations fill native array outputs with NaN;
pre-call validation failures leave caller buffers untouched. There is no point
batching or parameter broadcasting. Inputs use the existing real-array conversion
rules. See the [native shape tables](../../docs/dirichlet_shapes.md).

The six `scifort_dirichlet_*_f64` C adapters require k, explicit density point
length m when applicable, native inputs, caller-owned result storage and status.
Array results need k or k*k doubles in column-major order. Status 1 rejects
invalid extents or point-length choices without changing result buffers;
status 0 retains native numerical values, including NaN. All C buffers must be
valid, aligned, contiguous, non-NULL and disjoint. These are call-level statuses,
not native `dirichlet_status_*` error values. Native numerical implementations
and their limitations are unchanged. Rebuild the development library for these
six symbols. RNG procedures remain native Fortran APIs; their array shapes and
empty-batch behavior are documented in the source.

## Multinomial and multivariate hypergeometric distributions

```python
api.multinomial.pmf([1., 3.], n=4, p=[.2, .8])
api.multinomial.logpmf([1., 3.], n=4, p=[.2, .8])
api.multinomial.entropy(n=4, p=[.2, .8])
api.multinomial.mean(n=4, p=[.2, .8])
api.multinomial.cov(n=4, p=[.2, .8])
api.multivariate_hypergeom.pmf([1, 3], m=[3, 5], n=4)
api.multivariate_hypergeom.logpmf([1, 3], m=[3, 5], n=4)
api.multivariate_hypergeom.mean(m=[3, 5], n=4)
api.multivariate_hypergeom.var(m=[3, 5], n=4)
api.multivariate_hypergeom.cov(m=[3, 5], n=4)
```

Both families use k >= 1 categories and return scalar floats or binary64 arrays.
Mean/variance return shape (k), covariance shape (k,k); array methods accept
`out=ndarray` with the same Fortran-contiguous storage, identity, alignment,
write-permission and nonoverlap rules as Dirichlet. There is no broadcasting or
point batching. See [multinomial shapes](../../docs/multinomial_shapes.md) and
[multivariate hypergeometric shapes](../../docs/multivariate_hypergeom_shapes.md).

Multinomial counts `x(k)` use the native real-array API, so floating-point input
is accepted. A valid count vector is finite, nonnegative, integer-valued and sums
to n; impossible counts return PMF 0/log-PMF negative infinity. Invalid n or
adjusted probabilities return native NaN. All input probabilities must be finite.
When their sum differs from 1 by more than `10*epsilon(binary64)`, native code
replaces the final probability with `1-sum(p[:-1])`, then checks the adjusted
vector lies in [0,1]. This can replace an originally negative/out-of-range final
entry. The wrapper retains this native behavior. Zero trials are supported;
entropy uses a native loop with work proportional to n times category count.

Trial/draw count n must be a Python/NumPy integer scalar. Hypergeometric `x(k)`
and population `m(k)` must contain integers; floating-point arrays (even integral
ones), booleans, strings and object arrays are rejected. Signed/unsigned and
nonnative-endian integer arrays are accepted after range checks and converted to
aligned native `np.intc` storage. Noncontiguous and read-only inputs are accepted.
No integer narrowing or float-to-integer truncation happens before validation.

The integer ABI uses C int, with Python values checked against signed 32-bit
bounds. These reviewed adapters further require n <= 2147483646, and the
hypergeometric arrays each satisfy `sum(abs(array)) <= 2147483646`, evaluated in
widened arithmetic. These call-level bounds prevent overflow in existing native
default-integer sums and count increments. Out-of-range values raise before C;
raw C guards return status 1 and preserve output. Smaller negative counts still
reach native numerical validation: invalid population/draw parameters yield NaN,
while impossible observations yield PMF 0/log-PMF negative infinity. Hypergeometric
zero population with n=0 retains native PMF 1 and zero moments.

Every C argument is required, with explicit category count k, native inputs,
caller-owned result buffers and status. Arrays must be valid, aligned, contiguous,
non-NULL and disjoint; real arrays/results use doubles, integer inputs use C int.
Status 1 denotes invalid extent or configured integer bounds; status 0 retains
native numerical results, including NaN. Numerical failures can write NaN into
caller buffers; Python validation failures leave them untouched. These call
statuses are separate from native sampling status constants. Rebuild the
development library for the ten new symbols. Native numerical limitations remain
unchanged; the source now documents sampling and private-helper shape contracts.

## Generation and shared-library boundary

[Python policy](../cabi_config/python.json) references the existing C manifests.
It provides explicit defaults and mappings for the three normal methods already
maintained outside the C generator. The reusable
[generator/runtime](../../tools/cabi_generator/README.md) emits method signatures,
ctypes metadata, a runtime copy and small supplementary Fortran adapters.

```text
python tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json
python tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json --check
```

Generation requires the pinned fparser dependency in the generator requirements.
All outputs are checked in, so building the development library does not require
that parser. The supplementary adapters under `native/` supply the scalar symbols
normally defined in `src/c_api/api.f90` and the three original normal bulk methods.
They call the same native procedures; they do not copy numerical algorithms.
They are compiled **only** by the distribution-only development builder, never together
with `src/c_api/api.f90`, which would duplicate C symbols. The normal supplementary
bulk calls validate lengths and permit NULL buffers for empty calls. These
supplements are a development build choice, not replacements for the established
full-library implementations. The normal methods generated in `src/` remain
unchanged. Supplementary headers and coverage files are internal build records.

The development library excludes other distributions, the generated linear solver,
callbacks, random states, allocatable results, R/MATLAB/Octave wrappers and wheel
bundling of native binaries. Its cache hashes compiler, flags and all Fortran
sources. Use `--checked` for bounds/runtime checks; `--compiler` chooses a compiler.
Only GNU Fortran 17 on Windows has been validated for this initial build path.

## Tests and practical Python timing

```text
python -m unittest discover -s bindings/python/tests -v
python -m unittest discover -s benchmarks -p test_compare_python_scipy.py -v
python benchmarks\compare_python_scipy.py --quick --check-only
python benchmarks\compare_python_scipy.py
```

The first command checks both runtime-checked and optimized native builds against
SciPy, plus scalar/bulk agreement, round trips, tails, parameters and memory layout.
The timing command compares Python calling the public `scipy.stats` distribution
methods with Python calling SciFort. Both allocate outputs; Python validation,
allocation, dispatch and FFI overhead are timed. Inputs, imports and builds are
excluded. It covers all 56 distribution-array workloads, validates both input
sets before timing, calibrates each backend, alternates sample order and consumes
checksums. It uses public `scipy.stats.norm.cdf/ppf`, rather than the direct special
kernels in the older native-throughput benchmark.

Results are saved to `build/python-scipy-speed/results.json` and `results.md`.
Markdown includes median, geometric mean, arithmetic mean, sample SD, minimum
and maximum of workload ratios. `--case NAME` selects workloads; `--list` lists
without compiling. `--quick` and short timing batches are smoke checks and should
not be published as representative speed measurements. Use a quiet machine for
comparison runs and preserve both JSON and Markdown when publishing results.

## Validation recorded for this change

On 2026-11-23, validation used Python 3.13.3, NumPy 2.2.6, SciPy 1.15.3,
fparser 0.2.5 and GNU Fortran/GCC/G++ 17.0.0 20260510 (experimental) on Windows.
The focused checks cover the new generator/runtime, package and development
library; numerical core implementations and the global FPM build are unchanged.
The full FPM suite was not rerun. Other compilers/platforms remain unverified.

Commands used (the isolated generator environment already contains fparser):

```text
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
python -m unittest discover -s bindings/python/tests -v
python -m unittest discover -s benchmarks -p test_compare_python_scipy.py -v
build\cabi-generator-venv\Scripts\python.exe tools/cabi_generator/generate_python.py bindings/cabi_config/python.json --check
build\cabi-generator-venv\Scripts\python.exe tools/cabi_generator/generate.py bindings/cabi_config/manifests.json --batch --check
python bindings/python/example.py
python bindings/python/build_library.py --directory build/python-bindings/release
python benchmarks/compare_python_scipy.py --quick --check-only
python benchmarks/compare_python_scipy.py --repeats 3 --target-seconds .01 --output build/python-scipy-speed/initial.json
python -m pip wheel --no-deps --no-build-isolation ./bindings/python -w build/python-wheel
```

The reusable suite passed 18 existing tests and four new Python-generator tests;
five Python integration tests and two benchmark tests also passed (29 total).
Compiled paths use both `-O0 -g -fcheck=all -fbacktrace` and `-O3`, with
`-std=f2018 -fimplicit-none -Werror=implicit-interface`. All 56 benchmark workloads
passed validation and the short trial timing run. The wheel was unpacked under
`build/` and imported in an isolated Python process, which successfully called
the native normal CDF. Generated-file drift, local links, Fortran ASCII/line
lengths and whitespace checks passed. Trial timing output is not a publishable
performance baseline.

## Additional reviewed families (2026-11-24)

The interface now has 80 methods across ten families. For example:

```python
api.gamma.cdf([0., 1., 4.], a=2.3)
api.beta.ppf([.1, .5, .9], a=2.3, b=3.1)
api.t.sf([-2., 0., 2.], df=5.)
api.t.cdf([-1., 0., 1.], df=float("inf"))  # exact native normal limit
```

Gamma/beta shapes require finite positive values. Student t accepts any positive
`df`, including infinity. Tests check this explicitly in scalar, array and empty
array calls. The Python timing harness still measures its original 56 workloads;
adding interface methods does not change those saved comparisons.

SciPy 1.15.3 has small Student t(5) quantile errors at probabilities .1 and .9.
Those two test references use an independent 80-digit incomplete-beta inversion
instead of relaxing the common tolerances. Existing numerical implementations
and public C scalar symbols are unchanged.

Focused validation commands for this extension:

```text
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s bindings/cabi_config/tests -p test_generated_shape_families.py -v
python -m unittest discover -s bindings/python/tests -v
python -m unittest discover -s benchmarks -p test_compare_python_scipy.py -v
```

The new C integration test compiles verbatim production PDF/CDF/PPF scalar bodies
and direct native oracles, then links real C clients and checks C++ headers. It
checks all 24 methods in runtime-checked and optimized profiles, invalid uniform
parameters, preservation of output, empty/oversized calls and infinite Student
degrees of freedom. The Python integration covers all 80 methods in both profiles.
Compiler, flags and environment remain those recorded above. The global FPM
suite is not needed for these isolated tooling/interface additions.

This extension passed 41 focused tests: 32 reusable generator/discovery tests,
one compiled C ABI integration test with 24 methods in both profiles, six Python
integration tests covering all 80 methods in both profiles, and two benchmark
harness tests. The original 56 benchmark workloads also passed validation with
the expanded development library. Generation drift checks covered 33 C artifacts
and 32 Python/supplement artifacts. The updated example, local links, generated
Fortran ASCII/line limits and whitespace checks passed. No full FPM run was made.

## Generated build plan and chi-square (2026-11-24)

The current interface has 88 methods across eleven families. `api.chi2.cdf(x, df)`
and its seven other methods require finite positive `df`. For example:

```python
api.chi2.cdf([0., 1., 4.], df=5.)
```

The development builder no longer contains distribution or source lists. It reads
[build_sources.json](build_sources.json), generated by the unified workflow above.
The plan records a deterministic module-dependency order and source hashes.
Changing any selected source (including a dependency body) requires running the
unified regeneration command before building this development library. A normal
build uses checked-in metadata and requires no parser. The FPM build is unchanged.

The unified `--check` also covers C wrappers, Python code, compatibility adapters,
aggregate headers and [combined coverage](../cabi_config/interface_coverage.md).
The original individual generators remain available, but alone do not refresh
the build plan. Supplemental adapters remain exclusive to the development build;
linking them with the full manual C API would duplicate symbols.

Validation commands for the unified workflow:

```text
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
python -m unittest discover -s bindings/python/tests -v
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s bindings/cabi_config/tests -p test_generated_shape_families.py -v
python -m unittest discover -s benchmarks -p test_compare_python_scipy.py -v
build\cabi-generator-venv\Scripts\python.exe tools/cabi_generator/pipeline.py bindings/cabi_config/pipeline.json --check
python benchmarks/compare_python_scipy.py --quick --check-only --output build/python-scipy-speed/pipeline-check.json
python bindings/python/example.py
```

These checks use the same Python 3.13.3, NumPy 2.2.6, SciPy 1.15.3, fparser
0.2.5 and GNU Fortran/GCC/G++ 17.0.0 20260510 (experimental) Windows environment
recorded above. The optional development build and generated ABI paths are the
changed scope; no numerical implementation or global FPM configuration changed.
No full FPM suite was run. Other compilers/platforms remain unverified.

This change passed 49 focused tests: 39 reusable generator/workflow tests,
seven Python integration tests, one compiled C integration test exercising
32 methods in two profiles, and two benchmark tests. The independent workflow
fixture compiles its generated source order under runtime checks and `-O3`.
The unified check covers 75 artifacts, and the original 56 benchmark workloads
passed validation with the new build plan. No timing baseline was refreshed.

## Array-valued parameters (2026-11-24)

Normal and gamma now expose 16 additional C elementwise entry points, used by
Python when any parameter is an array. Existing scalar and uniform-parameter
symbols and numerical implementations remain unchanged. This initial step covered two families; the expansion below covers the other nine. The Python contract above specifies broadcasting,
invalid-element behavior, output shape and overlap rules.

Validation used the first five commands in the preceding validation block:
41 generator/workflow tests, eight Python integration tests, one C integration
test and two benchmark tests (52 total), plus the 75-artifact read-only pipeline
check. Python and compiled C checks cover runtime-checked (`-O0 -g -fcheck=all
-fbacktrace`) and optimized (`-O3`) builds with GNU Fortran/GCC/G++ 17.0.0
20260510 (experimental) on Windows; Python/NumPy/SciPy/parser versions remain
those recorded above. New tests cover all eight normal/gamma methods with mixed
scalar/array parameters, multidimensional broadcasts, endpoints, invalid elements,
round trips, dtype/layout conversions, empty results and overlap with any input.
SciPy references explicitly use binary64 parameters to match this interface.
Changed local links, generated Fortran ASCII/line limits and whitespace passed.
No full FPM suite or timing-baseline refresh was run: changes are confined to
interface tooling/adapters and the optional development build's generated hashes.
Other compilers/platforms remain unverified.

## Broadcasting across all configured families (2026-11-24)

All 88 Python methods now accept array-valued parameters. The nine-family
expansion adds 72 elementwise C symbols through reviewed manifests, using the
existing generator/runtime without changes. Legacy scalar and uniform-bulk
wrapper bodies and parameter policies are unchanged. The current Python contract
above applies to every configured family.

Validation commands (CMD):

```text
python -m unittest discover -s bindings/python/tests -v
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s bindings/cabi_config/tests -p test_generated_shape_families.py -v
python -m unittest discover -s benchmarks -p test_compare_python_scipy.py -v
build\cabi-generator-venv\Scripts\python.exe tools/cabi_generator/pipeline.py bindings/cabi_config/pipeline.json --check
python bindings/python/example.py
```

The 14 focused tests passed: ten Python integration tests across all 88 methods,
two compiled C integration tests and two benchmark harness tests. Python and C
integration used checked (`-O0 -g -fcheck=all -fbacktrace`) and optimized (`-O3`)
builds with GNU Fortran/GCC/G++ 17.0.0 20260510 (experimental) on Windows. Python,
NumPy, SciPy and fparser versions remain those recorded above. C clients link the
elementwise gamma/beta/Student t/chi-square symbols and compile under C11 and
C++17. The updated beta-broadcast and varying-df examples ran successfully.

Coverage includes multidimensional and scalar-point broadcasts, both beta shapes,
invalid shapes/location/scale, support/probability endpoints, empty results,
round trips, output-buffer validation, and array-valued infinite Student t df.
SciPy's known inverse-t reference inaccuracies are handled by checking forward
probabilities and exact agreement with unchanged native scalar calls. No numerical
source, generator implementation or global build configuration changed, so the
full FPM suite and unchanged generator unit tests were not rerun. Other platforms
and compilers remain unverified; timing baselines were not refreshed.

The read-only pipeline check passed for all 75 artifacts. Changed local links,
Fortran ASCII/line limits, encoding and whitespace checks also passed.

## Matrix Python bindings and shape contracts (2026-11-24)

The matrix interface above is generated from two reviewed C manifests and eight
structured native shape contracts. The optional build plan includes both native
matrix modules, their dependencies and their existing C adapters. No numerical
implementation or existing C symbol changed.

Validation commands (CMD):

```text
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
python -m unittest discover -s bindings/python/tests -v
python -m unittest discover -s benchmarks -p test_compare_python_scipy.py -v
build\cabi-generator-venv\Scripts\python.exe tools/cabi_generator/pipeline.py bindings/cabi_config/pipeline.json --check
python bindings/python/example.py
```

All 69 focused tests passed: 52 generator/workflow tests, 15 Python integration
tests and two benchmark harness tests. Matrix integration tests cover all eight
methods in checked (`-O0 -g -fcheck=all -fbacktrace`) and optimized (`-O3`) builds,
including native statuses and partial outputs, singular matrices, empty RHS
batches, shape constraints, conversions, caller buffers and overlap rejection.
Validation used GNU Fortran/GCC/G++ 17.0.0 20260510 (experimental), Python 3.13.3,
NumPy 2.2.6, SciPy 1.15.3 and fparser 0.2.5 on Windows. The example ran successfully.
The full FPM suite was not run because changes affect optional interface tooling
and documentation, with unchanged numerical implementations and global build
configuration. Other platforms and compilers remain unverified.

The read-only pipeline check passed for all 80 artifacts. Changed local Markdown
links, UTF-8/BOM, Fortran ASCII/line limits and whitespace checks passed. The
native `linalg.f90` diff contains only comments. Timing baselines were not refreshed.

## Multivariate normal density interface validation (2026-11-24)

Validation commands (CMD):

```text
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
python -m unittest discover -s bindings/python/tests -v
python -m unittest discover -s benchmarks -p test_compare_python_scipy.py -v
build\cabi-generator-venv\Scripts\python.exe tools/cabi_generator/pipeline.py bindings/cabi_config/pipeline.json --check
python bindings/python/example.py
```

All 78 focused tests passed: 56 generator/workflow tests, 20 Python integration
tests and two benchmark harness tests. New independent fixtures cover function
result/status emission, required optional arguments, invalid policy rejection,
Python array defaults, method prefixes, build-plan/coverage integration, compiled
C11 clients and C++17 header checks. The five multivariate integration tests cover
reference values, default dimensions, singular support, invalid/nonfinite inputs,
layout conversions, shape/alias rejection, C output preservation and direct
native Fortran/C agreement. Numerical NaN/infinity values are retained.

Tests ran on Windows with GNU Fortran/GCC/G++ 17.0.0 20260510 (experimental),
Python 3.13.3, NumPy 2.2.6, SciPy 1.15.3 and fparser 0.2.5. The development library
was exercised in checked (`-O0 -g -fcheck=all -fbacktrace`) and optimized (`-O3`)
builds, also using `-std=f2018 -fimplicit-none -Werror=implicit-interface`.
C clients use `-std=c11 -Wall -Wextra -Werror`; headers also pass C++17 checks.
The example ran successfully. Other platforms/compilers remain unverified.

The development DLL now links with its final basename inside a temporary staging
directory before atomic publication. This preserves its embedded Windows import
name for C/Fortran clients; old cached builds are invalidated. Both client types
linked and ran successfully against the rebuilt DLLs.

Native numerical implementations and global FPM configuration are unchanged;
`multivariate_normal.f90` changes contain only documentation. Validation targets
the affected optional generator/build/runtime paths, so the full FPM suite was
not rerun. Timing baselines were not refreshed.

The read-only pipeline check passed for all 84 artifacts. Shape-policy agreement,
local Markdown links, encoding/BOM, Fortran ASCII/line limits and whitespace
checks passed. Reviewed coverage now includes 98 Python methods and ten C array
adapters; it remains a reviewed subset, not whole-project coverage.

## Multivariate Student t density interface validation (2026-11-24)

Validation commands (CMD):

```text
python -m unittest discover -s bindings/python/tests -v
build\cabi-generator-venv\Scripts\python.exe tools/cabi_generator/pipeline.py bindings/cabi_config/pipeline.json --check
python bindings/python/example.py
```

All 25 Python integration tests passed, including five new Student t tests.
Both checked (`-O0 -g -fcheck=all -fbacktrace`) and optimized (`-O3`) libraries
were exercised, using `-std=f2018 -fimplicit-none -Werror=implicit-interface`.
The new tests cover 1D/2D/3D reference values, explicit and default parameters,
finite df through 100, exact infinite-df agreement with the normal methods,
singular PDF/log-PDF policies, invalid/nonfinite inputs, layout conversions,
shape/alias rejection, and C call-status/output preservation. Compiled Fortran
clients compare native and C results, including NaN/infinity cases. C clients
compile/run under `-std=c11 -Wall -Wextra -Werror`; headers pass C++17 checks.
The example ran successfully with both Student t methods.

Validation used GNU Fortran/GCC/G++ 17.0.0 20260510 (experimental), Python 3.13.3,
NumPy 2.2.6, SciPy 1.15.3 and fparser 0.2.5 on Windows. Other platforms/compilers
remain unverified. No generator implementation, Python runtime, numerical code
or global FPM configuration changed. The native module diff contains comments
only; tests target the expanded optional build plan and its interfaces. The full
FPM suite and unchanged generator unit tests were not rerun. Timing baselines
were not refreshed.

The read-only pipeline check passed for all 88 artifacts. Native shape-policy
agreement, Markdown local links, encoding/BOM, Fortran ASCII/line limits and
whitespace checks passed. Reviewed coverage now includes 100 Python methods and
12 C array adapters; whole-project coverage remains broader than this subset.

## Dirichlet and array-result validation (2026-11-25)

Validation commands (CMD):

```text
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
python -m unittest discover -s bindings/python/tests -v
python -m unittest discover -s bindings/python/tests -p test_dirichlet.py -v
build\cabi-generator-venv\Scripts\python.exe tools/cabi_generator/pipeline.py bindings/cabi_config/pipeline.json --check
python bindings/python/example.py
```

All 60 generator/workflow tests passed, including independent compiled vector
and matrix result fixtures, C buffer preservation, length-choice guards, invalid
policies and source-contract drift. The 25 existing Python integration tests
passed, and all six new Dirichlet tests passed on the focused rerun after fixing
a ctypes scalar-output pointer type in the test. Together these cover 91 test
cases. The implementation required no correction from the numerical comparisons.
The example ran successfully with all six Dirichlet methods.

Dirichlet coverage includes full/shortened points, k=1, simplex boundaries and
sum tolerance, invalid/nonfinite concentrations, reference densities/entropy/
moments, covariance identities, dtype/layout conversions, shape/alias rejection,
caller-buffer identity and preservation, native NaN outputs and direct Fortran/C
agreement. C clients exercise all six methods under C11 with
`-Wall -Wextra -Werror`; headers also pass C++17 checks.

Tests used GNU Fortran/GCC/G++ 17.0.0 20260510 (experimental), Python 3.13.3,
NumPy 2.2.6, SciPy 1.15.3 and fparser 0.2.5 on Windows. The libraries were tested
in checked (`-O0 -g -fcheck=all -fbacktrace`) and optimized (`-O3`) builds, with
`-std=f2018 -fimplicit-none -Werror=implicit-interface`. Other platforms and
compilers remain unverified.

All 92 artifacts passed the read-only pipeline check. Native shape-policy
agreement, local Markdown links, encoding/BOM, Fortran ASCII/line limits and
whitespace checks passed. Stripping comments confirms unchanged executable
Dirichlet/Student t code. Sampling and private validation now have readable
shape documentation, including the native empty-batch exception; RNG C handles
are outside this batch. No new external dependency or numerical formula was
introduced. The full FPM numerical suite was not rerun: the changed generator,
optional ABI/runtime and dependency build plan are covered by the focused checks.
Timing baselines were not refreshed.

## Integer-input distribution interface validation (2026-11-25)

Validation commands (CMD):

```text
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
python -m unittest discover -s bindings/python/tests -v
python -m unittest discover -s bindings/python/tests -p test_discrete_multivariate.py -v
build\cabi-generator-venv\Scripts\python.exe tools/cabi_generator/pipeline.py bindings/cabi_config/pipeline.json --check
python bindings/python/example.py
```

All 64 generator/workflow tests passed. Independent integer fixtures compile and
run in checked and optimized profiles, including `-ftrapv` for the checked integer
guards, INT_MIN/INT_MAX arrays, widened sum limits and unchanged rejected outputs.
All 31 existing Python integration tests passed; the six new discrete-multivariate
cases passed on their focused rerun after correcting a C test's bitwise covariance
symmetry assumption to a numerical check. The native entries differed by one
rounding unit; numerical implementations were not changed. These checks cover
101 distinct test cases. The example ran successfully with all ten new methods.

The new cases cover SciPy reference values, mass/support behavior, zero trials
and populations, final-probability adjustment, entropy and covariance identities,
integer dtype/endianness/alignment conversions, rejected floats/booleans/overflow,
shape/alias validation, caller-buffer identity/preservation and direct native/C
agreement. Compiled C clients call all ten symbols under C11 with
`-Wall -Wextra -Werror`; headers pass C++17 checks.

Validation used GNU Fortran 17.0.0 20260510 (experimental), Python 3.13.3,
NumPy 2.2.6, SciPy 1.15.3 and fparser 0.2.5 on Windows, with GCC/G++ clients.
The shared libraries use checked (`-O0 -g -fcheck=all -fbacktrace`) and optimized
(`-O3`) profiles, plus `-std=f2018 -fimplicit-none -Werror=implicit-interface`.
Other platforms/compilers remain unverified. Native numerical limitations and
probability-adjustment behavior are retained; explicit ABI bounds protect native
integer arithmetic without extending its precision.

The read-only pipeline check passed for all 100 artifacts. Native source diffs
contain documentation only after stripping comments. Shape-policy agreement,
Markdown local links, encoding/BOM, Fortran ASCII/line limits and whitespace
checks passed. Coverage now records 116 Python methods and 28 array C adapters
within the reviewed subset. No new dependency or numerical algorithm was added.
The full FPM numerical suite was not rerun: focused checks cover the changed
optional generator/runtime, C integration and dependency build plan. Timing
baselines were not refreshed.

## Dirichlet-multinomial distribution

`api.dirichlet_multinomial.pmf(x, alpha, n)` and `logpmf(x, alpha, n)`
return a Python float. `mean(alpha, n)`, `var(alpha, n)` and `cov(alpha, n)`
return NumPy arrays and accept an optional `out=` array. `alpha` and `x` have
shape `(k,)`, with `k >= 1`; covariance has shape `(k,k)`. Matrices and caller
buffers follow the existing column-major float64 contract. These are single
observations without broadcasting or batching. See the generated
[shape tables](../../docs/dirichlet_multinomial_shapes.md).

`alpha` must be finite and positive. Counts are real-valued inputs constrained
to finite nonnegative integers. Native invalid counts return NaN; finite integer
counts whose sum differs from `n` return zero mass and negative-infinite log mass.
This retains SciFort behavior, including NaN for negative or fractional counts,
which can differ from SciPy support handling. Invalid concentrations or negative
`n` produce native NaNs. Zero trials are supported, whereas the SciPy 1.15.3
reference used here requires positive `n`. Concentration sums and
intermediate moment products must remain finite; this batch adds no extended
range algorithm for extreme concentrations.

Python validates shapes, types, aliases, output buffers and the signed C integer
range before calling Fortran. The reviewed ABI caps `n` at 2147483646, consistently
with the other discrete multivariate adapters. Call-level errors raise in Python;
direct C calls return status 1 and preserve the output buffer. Valid calls return
status 0 even when the numerical result is NaN or infinity. The native Fortran
API retains its default-integer range, including `huge(0)` after the count-increment
fix. All C buffers must be valid, nonnull, aligned and nonoverlapping.

```python
mass = api.dirichlet_multinomial.pmf([1, 3], [2., 3.], 4)
mean = api.dirichlet_multinomial.mean([2., 3.], 4)
cov = api.dirichlet_multinomial.cov([2., 3.], 4)
```

### Validation: Dirichlet-multinomial batch (2026-11-25)

Windows validation used Python 3.13.3, NumPy 2.2.6, SciPy 1.15.3 and GNU
Fortran/GCC/G++ 17.0.0 20260510 (experimental). Commands from the repository root:

```text
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json --check
python -m unittest discover -s bindings/python/tests -v
python bindings\python\example.py
```

All 42 Python binding integration tests passed on the final code/header state.
The five new cases cover reference masses and moments, zero trials, invalid
concentrations/counts, nonfinite/large counts, caller buffers, integer bounds,
output preservation and native/C agreement for all five methods. Each runs
against checked (`-O0 -g -fcheck=all -fbacktrace`) and optimized (`-O3`) libraries;
both use `-std=f2018 -fimplicit-none -Werror=implicit-interface -shared -fPIC`.
The numerical module and count regression program are also compiled directly
with `-ftrapv -ffpe-trap=invalid` in the checked case, and `-ffpe-trap=invalid`
in the optimized case. C clients use `-std=c11 -Wall -Wextra -Werror`; C++17
header syntax checks also passed. The example demonstrates all five methods. The read-only pipeline check passed
for all 104 artifacts. Encoding/BOM, Fortran ASCII/line lengths, Markdown local
links and whitespace checks passed.

An initial comparison incorrectly asked SciPy to accept zero trials; this is now
a native-behavior check. The broader integration run also caught and corrected
a duplicate header guard in the new manifest. No generator/runtime implementation
changed. All 19 reviewed header guards are distinct. Coverage now records 121
Python methods and 33 array adapters in the reviewed subset.

The full FPM suite was not rerun: executable changes are confined to this
numerical module and its optional interfaces, covered by direct native tests and
the binding integration suite. Other platforms and extreme concentration ranges
remain unverified; benchmarks were not refreshed.
