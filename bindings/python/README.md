# SciFort Python interface

This development interface exposes PDF, log-PDF, CDF, SF, log-CDF, log-SF, PPF,
and ISF for `normal`, `exponential`, `logistic`, `laplace`, `cauchy`, `weibull`,
`lognormal`, `gamma`, `beta`, Student `t`, and `chi2`. It uses NumPy arrays and the stable
binary64 C calling convention. Eight matrix methods are available under `api.linalg`,
and multivariate normal/Student-t PDF/log-PDF, entropy and marginal parameters under `api.multivariate_normal`
and `api.multivariate_t`. `api.dirichlet` provides PDF, log-PDF, entropy, mean,
variance and covariance. Multinomial and multivariate hypergeometric methods are
available under `api.multinomial` and `api.multivariate_hypergeom`.
`api.dirichlet_multinomial` provides PMF, log-PMF, mean, variance and covariance.
`api.wishart` and `api.invwishart` provide densities, entropy and matrix moments.
`api.matrix_normal` provides PDF, log-PDF and entropy; `api.matrix_t` provides
PDF and log-PDF for rectangular matrix variates.
`api.poisson_binom` provides eight mass, tail and quantile methods for a vector
of Bernoulli probabilities.
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
six symbols. The later [explicit RNG and Dirichlet sampling section](#explicit-rng-streams-and-dirichlet-sampling)
describes the added batch sampler and stream bridge.

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
unchanged; the source documents sampling and private-helper shape contracts.
The later [integer-count batch sampling section](#integer-count-batch-sampling)
describes the added samplers and their distinct status and integer-output contracts.

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

## Wishart and inverse-Wishart distributions

Both `api.wishart` and `api.invwishart` expose these methods:

| Python call | Result |
| --- | --- |
| `pdf(x, df, scale)` | Python float |
| `logpdf(x, df, scale)` | Python float |
| `entropy(df, scale)` | Python float |
| `mean(df, scale, out=None)` | float64 matrix of entrywise means |
| `mode(df, scale, out=None)` | float64 mode matrix |
| `var(df, scale, out=None)` | float64 matrix of entrywise variances, not a covariance tensor |

`p = scale.shape[0]`, with `p >= 1`; `scale`, `x` and matrix outputs have shape
`(p,p)`. Supply full, finite, symmetric positive-definite matrices. Degrees of
freedom are real and finite, with `df > p-1`. Inputs have no default scale or
batch/broadcast axis. Python accepts numeric matrix inputs in C or Fortran order
and converts them to aligned, column-major float64 storage when necessary.
`out=` must already be writable, aligned, Fortran-contiguous float64 of the exact
shape. It is returned by identity. Buffers may not overlap any input.

| Method | Additional existence condition | Native behavior when undefined |
| --- | --- | --- |
| Wishart mode | `df >= p+1` | matrix of NaNs; zero matrix at equality |
| Inverse-Wishart mean | `df > p+1` | matrix of NaNs |
| Inverse-Wishart variance | `df > p+3` | matrix of NaNs |

Invalid numerical parameters, singular/non-positive-definite matrices, and
nonfinite matrix entries produce native NaNs. In particular, singular evaluation
matrices produce NaN rather than a support-boundary density. This batch preserves
native behavior and does not add singular Wishart distributions. These choices
may differ from SciPy exceptions, zero densities or `None` for undefined moments.
A density may underflow to zero while `logpdf` remains finite.

Factorization uses lower triangles, but every matrix entry must be finite.
Densities and entropy use the scale's lower-triangle interpretation. The moment
formulas read the full scale entries; the Wishart density also reads the full
`x` matrix in its trace solve. Inverse-Wishart density factorization uses only
`x`'s lower triangle. No symmetry check or automatic symmetrization is performed:
supply symmetric, fully populated inputs for consistent results. These are
existing native conventions, retained by the adapters.

Python raises for invalid shapes, dtypes, aliasing and output buffers before
writing output. C uses explicit dimensions, column-major contiguous arrays and
caller-owned result buffers. All pointers must be valid, nonnull, aligned and
nonoverlapping. C status 1 means an invalid dimension and preserves the result
buffer. Status 0 means the call was valid, including numerical NaN/infinity
results. The later [batch sampling section](#wishart-and-inverse-wishart-batch-sampling)
describes the added samplers and their distinct status contracts.

See the generated [Wishart shapes](../../docs/wishart_shapes.md) and
[inverse-Wishart shapes](../../docs/invwishart_shapes.md). Native sampling and
private helper declarations also document dimensions and permitted empty batches.

```python
scale = [[2., .4], [.4, 1.]]
value = api.wishart.logpdf([[1.5, .2], [.2, 1.1]], 7., scale)
mean = api.invwishart.mean(7., scale)
variance = api.invwishart.var(7., scale)
```

### Validation: Wishart interface batch (2026-11-25)

Windows checks used GNU Fortran 17.0.0 20260510 (experimental), Python 3.13.3,
NumPy 2.2.6 and SciPy 1.15.3. Commands from the repository root:

```text
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json --check
python -m unittest discover -s bindings/python/tests -v
python bindings\python\example.py
```

All 48 binding integration tests passed, including six new cases exercising
all twelve methods against SciPy, degrees-of-freedom thresholds, undefined
moments, finite log densities when PDFs underflow, matrix-triangle conventions,
invalid matrices, layout conversion, aliases, caller buffers, status/output
preservation, and native/C agreement. Tests use both checked
`-O0 -g -fcheck=all -fbacktrace` and optimized `-O3` libraries, with
`-std=f2018 -fimplicit-none -Werror=implicit-interface -shared -fPIC` in both.
Native clients use the corresponding optimization/runtime-check flags; C clients
use `-std=c11 -Wall -Wextra -Werror`, with C++17 header syntax checks too.

The read-only consistency check passed for 112 generated artifacts, and the
example exercised all twelve new methods. All 21 reviewed header guards are
distinct. Encoding/BOM, ASCII Fortran, 132-column limits, Markdown local links
and whitespace checks passed. Reviewed coverage now includes 133 Python methods
and 45 array adapters.

After stripping comments, both native modules' executable code matches the
previous commit. The generator/runtime implementation is unchanged. Focused
native and complete binding integration checks cover the new manifests,
dependency build plan, aggregate header and output contracts; the full FPM
numerical suite was not rerun. Sampling remains native-only. Other platforms,
singular distribution extensions and extreme matrix/parameter ranges were not
newly validated; benchmark baselines were not refreshed.

## Matrix-normal and matrix-t distributions

| Python call | Result |
| --- | --- |
| `api.matrix_normal.pdf(x, mean=None, rowcov=None, colcov=None)` | Python float |
| `api.matrix_normal.logpdf(x, mean=None, rowcov=None, colcov=None)` | Python float |
| `api.matrix_normal.entropy(rowcov=None, colcov=None)` | Python float |
| `api.matrix_t.pdf(x, mean=None, row_spread=None, col_spread=None, df=1.0)` | Python float |
| `api.matrix_t.logpdf(x, mean=None, row_spread=None, col_spread=None, df=1.0)` | Python float |

For density calls, `m,n = x.shape`, with `m >= 1`, `n >= 1`. Evaluation and
mean matrices have shape `(m,n)`, the row covariance/spread `(m,m)`, and the column
covariance/spread `(n,n)`. Rectangular evaluation matrices are supported; there
is no batching or broadcasting. Missing mean becomes zero, and missing row or
column matrices become identities of the corresponding dimensions.

Entropy has no evaluation matrix: `m` comes from the supplied row covariance's
first extent, or is 1 when absent; `n` comes independently from the supplied
column covariance's first extent, or is 1 when absent. Thus `entropy()` is the
entropy of a 1-by-1 standard normal matrix. `entropy(colcov=np.eye(3))` uses a
1-by-3 matrix, while `entropy(np.eye(2), np.eye(3))` uses a 2-by-3 matrix.
These are explicit reviewed fallback dimensions in Python, matching native
Fortran defaults.

Every input entry must be finite. Covariance/spread matrices must be positive
definite; their lower triangles define the factorizations. Full upper triangles
must still contain finite values. Matrix-t `df` is finite and strictly positive;
its default is 1. Its spread matrices are not covariance matrices. For a 1-by-1
matrix with row spread `r` and column spread `c`, the ordinary Student-t scale is
`sqrt(r*c/df)`, preserving the native/SciPy matrix-t parameterization.

Numerical parameter errors retain native NaNs. Python rejects wrong ranks/shapes,
complex inputs and overlapping input buffers before the native call. Numeric
inputs are converted to aligned, column-major float64 storage when needed.
These scalar-result methods have no `out=` argument. Direct C adapters require
explicit dimensions and all arrays/scalars, including the promoted native
optional arguments. All pointers must be valid, nonnull, aligned and disjoint;
arrays use contiguous column-major storage. C status 1 rejects nonpositive
extents and preserves the scalar result. Status 0 retains numerical NaN/infinity
results. The later [batch sampling section](#matrix-normal-and-matrix-t-batch-sampling)
describes the added matrix samplers and their status contracts.

See the [matrix-normal shapes](../../docs/matrix_normal_shapes.md),
[matrix-t shapes](../../docs/matrix_t_shapes.md), and the generator's
[dimension-fallback policy](../../tools/cabi_generator/README.md#optional-array-dimension-defaults).
Shape comments also cover sampling and private helpers, including empty batches
and the covariance helper's permitted zero extent. The matrix-t sampler's existing
BSD-3-Clause notices and provenance are retained.

```python
x = [[.1, .2, .3], [.4, .5, .6]]
normal_logpdf = api.matrix_normal.logpdf(x)
t_logpdf = api.matrix_t.logpdf(x, df=3.)
entropy = api.matrix_normal.entropy(np.eye(2), np.eye(3))
```

### Validation: matrix-normal and matrix-t interface batch (2026-11-25)

Windows checks used GNU Fortran 17.0.0 20260510 (experimental), Python 3.13.3,
NumPy 2.2.6 and SciPy 1.15.3. Commands from the repository root:

```text
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json --check
python -m unittest discover -s bindings/python/tests -v
python bindings\python\example.py
```

All 66 generator tests and 54 binding integration tests passed. The two new
independent generator tests validate fallback policy errors and compile/run an
unrelated numerical fixture with omitted/supplied dimension sources. The six
new integration cases cover rectangular matrix-normal references, the repository's
saved SciPy 1.17 matrix-t references, Student-t reductions, transpose identities,
optional arrays, independently omitted entropy covariances, invalid inputs,
layout conversion, aliases, scalar-output preservation, and native/C agreement
for all five methods. The native client also compares absent-argument entropy
calls directly with the corresponding explicit C inputs. Installed SciPy 1.15.3
has no `matrix_t`; a live comparison with that method was not run.

Integration tests use both checked `-O0 -g -fcheck=all -fbacktrace` and optimized
`-O3` libraries; both use `-std=f2018 -fimplicit-none -Werror=implicit-interface
-shared -fPIC`. Native clients use the corresponding checked/optimized flags.
C clients use `-std=c11 -Wall -Wextra -Werror`, with C++17 header syntax checks.
The independent fallback fixture also runs at `-O0 -g -fcheck=all` and `-O3`.

The example exercised all five new methods and entropy with absent covariances.
The read-only check passed for all 120 generated artifacts. All 23 reviewed
header guards are distinct. Encoding/BOM, ASCII Fortran, 132-column limits,
Markdown local links and whitespace checks passed. Reviewed coverage now records
138 Python methods and 50 array adapters.

Both native numerical implementations match the previous commit after stripping
comments; the matrix-t sampler's license notices are retained. The changes to
shared Python interface tooling are covered by its complete generator and binding
integration suites. The full FPM numerical suite was not rerun because the native
numerical code and overall FPM build/test process are unchanged. Sampling remains
native-only; other platforms and extreme parameter/matrix ranges were not newly
validated. Benchmark baselines were not refreshed.

## Poisson-binomial distribution

`api.poisson_binom` provides `pmf(k, p, loc=0.0)`, `logpmf`, `cdf`, `sf`,
`logcdf`, `logsf`, `ppf(probability, p, loc=0.0)` and `isf`. Each returns one
Python float. `p` is a rank-one vector of length `n >= 1`, containing finite
Bernoulli success probabilities in `[0,1]`; zeros and ones are permitted.
`loc` is a real support shift, defaulting to zero. These eight methods have no batching,
broadcasting or `out=` buffer argument.

Mass functions require an integer-valued `k-loc`; fractional or out-of-support
counts give zero mass and negative-infinite log mass. Tail functions floor
`k-loc`. Large finite observations are handled without converting them to default
integers; infinite observations yield the corresponding tail endpoints. NaN
observations or invalid probability vectors produce native NaNs. Quantile inputs
must be finite and in `[0,1]`; invalid probabilities produce NaN. Quantile
endpoints use the nominal support: PPF(0) and ISF(1) return `loc-1`, while PPF(1)
and ISF(0) return `loc+n`, including deterministic trial vectors.

Use finite `loc` for meaningful distribution calculations. Existing handling of
nonfinite shifts is retained: mass/tail methods evaluate the resulting shifted
observation, while quantile methods return NaN. Log tails are computed in the
log domain and can remain finite when their probabilities underflow to zero.

Python validates ranks and extents before calling C, and converts numeric inputs
to aligned, contiguous float64 as needed. Empty vectors or wrong ranks raise;
invalid probability values return native NaN. C adapters return status 1 for a
nonpositive trial count, preserving the scalar result; a valid call returns
status 0 even for numerical NaN/infinity. All C pointers must be valid, nonnull,
aligned and disjoint. Dimensions use C int, with valid buffers and sufficient
memory for native work arrays required.

The new status-returning symbols are named
`scifort_poisson_binom_<method>_explicit_f64`. Existing published scalar-returning
`scifort_poisson_binom_pmf_f64`, `scifort_poisson_binom_cdf_f64` and
`scifort_poisson_binom_ppf_f64` retain their original calling conventions.
The development Python library uses the new explicit adapters. The six mass/tail
public generics select the real-count native specifics; native integer overloads
remain available to Fortran callers. PPF/ISF use their direct native procedures.

See the [shape tables](../../docs/poisson_binom_shapes.md) and
[generic selection policy](../../tools/cabi_generator/README.md#reviewed-public-generic-selections).

### Poisson-binomial score

**SciFort extension:** `api.poisson_binom.logpmf_derivative_p(count, p, *, out=None)`
returns `ArrayResult(status, outputs)`. `outputs['logf']` is a float and
`outputs['dpvec']` is the length-`n` gradient of log mass with respect to `p`.
The observation is unshifted; this method has no `loc` argument.

```python
gradient = np.empty(3, dtype=np.float64)
result = api.poisson_binom.logpmf_derivative_p(
    1, [.2, .5, .8], out={'dpvec': gradient})
print(result.status, result.outputs['logf'], gradient)
```

`p` and the gradient have shape `(n,)`, `n >= 1`. Caller-owned outputs must be
writable, aligned, native float64 and contiguous, with no overlap with inputs.
C adds a final call-status pointer without changing the native procedure.
Nonpositive extents return status 1 with log mass and gradient unchanged;
valid extents return status 0, including numerical NaNs/infinities.

Invalid probabilities, nonfinite counts and fractional counts give NaN log
mass and gradients. Finite integer counts outside support and zero-mass
observations give negative-infinite log mass and NaN gradients. At probability
0 or 1, gradients are inward one-sided derivatives when the observed mass is
positive. Endpoint scores use closed forms; other scores use signed log-domain
leave-one-out masses, so finite gradients survive probability underflow.
A mathematically finite score beyond binary64 range can still overflow; small
scores near cancellation are limited by log-mass rounding. The general method
uses cubic work and linear workspace in the trial count; this change does not
add a faster score algorithm. Private allocation bounds are documented too.

Native probability tables use quadratic work in the number of trials; quantile
searches rebuild tables for their tail evaluations. This batch adds no caching
or bulk evaluation algorithm.

```python
p = [.2, .6, .8]
mass = api.poisson_binom.pmf(2, p)
log_tail = api.poisson_binom.logsf(2, p)
quantile = api.poisson_binom.ppf(.95, p)
```

### Validation: Poisson-binomial interface batch (2026-11-25)

Windows checks used GNU Fortran 17.0.0 20260510 (experimental), Python 3.13.3,
NumPy 2.2.6, SciPy 1.15.3 and FPM 0.12.0 alpha. Commands:

```text
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -p test_generic_specifics.py -v
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json --check
python -m unittest discover -s bindings/python/tests -v
python bindings\python\example.py
```

All 69 generator tests passed; the three generic-selection cases were rerun after
adding the shape-table assertion. All 60 binding integration tests passed, with
six new cases covering eight SciPy methods, an equal-probability binomial
identity, independent outcome enumeration, quantile inequalities, endpoints,
invalid/nonfinite parameters, large observations, log-tail underflow, layout
conversion, shape guards, output preservation and native/C agreement.

Integration uses checked `-O0 -g -fcheck=all -fbacktrace` and optimized `-O3`
libraries, both with `-std=f2018 -fimplicit-none -Werror=implicit-interface
-shared -fPIC`. The changed numerical module and native client were also compiled
directly with `-ftrapv -ffpe-trap=invalid` in both configurations for the large
finite-count regressions. Native tests check integer/real overload agreement,
ordinary two-trial derivatives, a one-trial derivative and out-of-support large
positive/negative derivative counts. C11 clients use `-Wall -Wextra -Werror`;
C++17 header syntax checks also passed. The example exercised all eight methods.

Dependency inspection identified likelihood, score, random-variates and the
existing C API as native callers. The existing
`test_wallenius_poisson_binom` program passed after recompiling the changed
Poisson-binomial module and its test helpers and linking the cached native
library. This exercised the existing C ABI, likelihood scores, fitting and RNG
integration. The local cache-specific compile/run commands were:

```bat
gfortran -O0 -g -fcheck=all -fbacktrace ^
  -I build\gfortran_2806A042DAA8BBF7 -J build\poisson-binom-regression ^
  src\stats\poisson_binom.f90 test\test_checks.f90 ^
  test\wallenius_poisson_binom_reference.f90 ^
  test\test_wallenius_poisson_binom.f90 ^
  build\gfortran_3C159EB802C22FB9\scifort\libscifort.a ^
  -o build\poisson-binom-regression\dependent.exe
build\poisson-binom-regression\dependent.exe
```

`fpm test --target test_wallenius_poisson_binom --profile debug` was started but
stopped when it triggered a broad rebuild of unrelated library/test targets.
The direct native check above passed; the FPM invocation did not finish and is
not counted as a passing run. The full FPM suite was not run. Cached native
dependencies were reused; the new-module direct tests additionally cover checked
and optimized builds from the current source.

The read-only consistency check passed for all 124 generated artifacts. All 24
reviewed header guards are distinct. Encoding/BOM, ASCII Fortran, 132-column
limits, Markdown links and whitespace checks passed. Coverage now records 146
Python methods and 58 array adapters in the reviewed subset. Existing published
C prototypes and symbols are unchanged. Native changes are limited to the two
conversion-safety fixes and documentation. Other platforms, very large trial
vectors and extended-range derivatives were not newly validated; benchmark
baselines were not refreshed.

## Poisson-binomial score validation (2026-11-25)

Validation used GNU Fortran/GCC/G++ 17.0.0 20260510 (experimental), Python 3.13.3,
NumPy 2.2.6, SciPy 1.15.3 and fparser 0.2.5. Commands from the repository root:

```bat
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
python -m unittest discover -s bindings/python/tests -v
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json --check
python bindings\python\example.py
```

All 72 generator tests and all 64 binding tests passed. The three independent
synthetic-status generator tests also passed separately after the final emitter
formatting change. Integration compiled checked `-O0 -g -fcheck=all -fbacktrace`
and optimized `-O3` libraries, with `-std=f2018 -fimplicit-none
-Werror=implicit-interface -shared -fPIC`. Native/C clients passed in both
configurations; the Poisson-binomial native module/client additionally passed
with `-ftrapv -ffpe-trap=invalid`. C11 and C++17 clients used
`-Wall -Wextra -Werror`.

New cases cover derivative enumeration, centered finite differences, binomial
score identities, an independent 80-digit Decimal reference for heterogeneous
underflow, deterministic probabilities, zero-mass observations, nonfinite and
fractional counts, caller-owned buffers, alias rejection and unchanged outputs
on invalid extents. Near-zero gradients use an absolute rounding tolerance.
The example exercised the score with a caller-owned gradient.

The existing `test_wallenius_poisson_binom` native dependency test passed after
recompiling the final numerical module with the exact cached-library compile/run
commands in the preceding validation section. This covers likelihood scores,
fitting, random sampling and the existing manual C ABI. The full FPM suite was
not run: the numerical change is confined to this module, while all generator
and configured binding integrations cover the shared adapter change. Cached
native dependencies were reused. A final declaration edit made an intermediate
binding run stale; the source-hash guard stopped it. The 64-test final run above
used regenerated metadata and completed successfully.

All 124 generated artifacts passed the read-only consistency check. Encoding,
BOM, ASCII Fortran, 132-column limits, local Markdown links and whitespace checks
passed. All previous Poisson-binomial C prototypes remain unchanged. Coverage is
147 Python methods and 59 array adapters in the reviewed subset. General scores
still require cubic work; cancellation accuracy and binary64 overflow limits
are described above. Other platforms and large-vector performance were not
newly validated.

## Multivariate entropy and marginal parameters

Both `api.multivariate_normal` and `api.multivariate_t` now expose `entropy`
and `marginal`. Native procedures retain their names and numerical behavior.

| Python call | Result |
| --- | --- |
| `multivariate_normal.entropy(mean=None, cov=None)` | Python float |
| `multivariate_t.entropy(loc=None, shape=None, df=1.0)` | Python float |
| `multivariate_normal.marginal(dimensions, mean, cov, *, out=None)` | `ArrayResult(status, outputs)` |
| `multivariate_t.marginal(dimensions, loc, shape, *, out=None)` | `ArrayResult(status, outputs)` |

Entropy dimension is inferred first from a supplied mean/location vector, then
from the supplied matrix's first axis, otherwise defaults to one. Supplied
vectors must have shape `(n,)` and matrices `(n,n)`, `n >= 1`; mismatches and
empty dimensions raise. Absent vectors become zeros and absent matrices become
identity. Provided location values must be finite even though entropy does not
depend on their values. This preserves native SciFort dimension inference.
Local SciPy 1.15.3 Student-t entropy uses a one-dimensional default when shape
is omitted, even with a longer loc; cross-library references therefore supply
an explicit identity shape for that case.
C always receives explicit arrays, n, and Student-t df.
Invalid numerical parameters return NaN with call status 0. Nonpositive C
extents return status 1 with the scalar output unchanged.

Both entropy functions permit PSD matrices and use the lower triangle, while
requiring every entry to be finite. Normal entropy uses numerical rank and the
pseudo-determinant, with rank-zero entropy zero. Finite-df Student-t entropy uses
the ambient dimension n and pseudo-determinant, including singular matrices;
this is a documented computational convention rather than a full-dimensional
density entropy for a degenerate law. Infinite df delegates to normal entropy,
including its rank convention. Invalid df gives NaN; the Python default is 1.

Marginal `dimensions` is a nonempty, rank-one **zero-based** integer vector in
Python. C and Fortran keep **one-based** indices. Let m be its length; n is the
input vector length. Vectors have shape `(n,)`, matrices `(n,n)`, selected
vectors `(m,)` and selected matrices `(m,m)`. Indices must be unique and in
range, implying `1 <= m <= n`; requested order is preserved. Floats and booleans
are rejected, and out-of-range Python indices raise before C conversion.
Duplicate indices return native status 2 and preserve caller-owned outputs.
Locations and all matrix entries must be finite. Selection symmetrizes from
the lower triangle; it does not require the input matrix to be positive definite.
Student-t marginal selection returns location/shape parameters; df is unchanged
and is not an argument or output of this native routine.

The normal output keys are `marginal_mean` and `marginal_cov`; Student-t uses
`marginal_loc` and `marginal_shape`. Pass buffers with these names in `out`.
They must be aligned, writable native float64, column-major, correctly shaped
and disjoint from each other and all inputs. Other buffers are converted as
needed; caller index arrays are never modified. Native statuses are 0 success,
1 invalid shape and 2 invalid parameters. Invalid C extents are rejected before
native execution; all failure paths preserve caller-owned outputs. All C array
and scalar output pointers must be valid, aligned, nonnull and have sufficient
storage; these adapters cannot verify addresses or buffer capacity.

```python
print(api.multivariate_normal.entropy(cov=np.eye(3)))
print(api.multivariate_t.entropy(shape=np.eye(3), df=5.0))
selected_mean = np.empty(2)
selected_cov = np.empty((2, 2), order='F')
result = api.multivariate_normal.marginal(
    [2, 0], [1., 2., 3.], np.eye(3),
    out={'marginal_mean': selected_mean, 'marginal_cov': selected_cov})
print(result.status, selected_mean, selected_cov)
```

See the [normal shape tables](../../docs/multivariate_normal_shapes.md) and
[Student-t shape tables](../../docs/multivariate_t_shapes.md). C symbols follow
`scifort_multivariate_normal_entropy_f64`, `scifort_multivariate_t_entropy_f64`,
`scifort_multivariate_normal_marginal_f64` and `scifort_multivariate_t_marginal_f64`.

## Multivariate parameter-interface validation (2026-11-25)

Validated on Windows with GNU Fortran/GCC/G++ 17.0.0 20260510 (experimental),
Python 3.13.3, NumPy 2.2.6, SciPy 1.15.3 and fparser 0.2.5. Commands:

```bat
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
python -m unittest discover -s bindings/python/tests -v
python -m unittest discover -s bindings/python/tests -p test_multivariate_parameters.py -v
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json --check
python bindings\python\example.py
```

All 76 generator tests passed, including four new independent policy/compiled
fixture tests. The complete 70-test binding run passed all 64 existing tests;
two new assertions required correction: one expected SciPy's omitted-shape
Student-t dimension to match SciFort's, and one reused mutated C index test
storage across families. After those test-only fixes, all six new integration
tests passed in the focused rerun. No library behavior was changed to satisfy
those assertions.

Integration used checked `-O0 -g -fcheck=all -fbacktrace` and optimized `-O3`
libraries, both with `-std=f2018 -fimplicit-none -Werror=implicit-interface
-shared -fPIC`. Native clients compare all four adapters directly with native
procedures; C11 and C++17 clients use `-Wall -Wextra -Werror`. Tests cover entropy
references, independent defaults and dimension precedence, singular and normal
limits, invalid/nonfinite parameters, marginal ordering and lower triangles,
zero-/one-based conversion, read-only/endian/strided integer inputs, unsigned
range checks, unchanged inputs, output ownership/aliasing and failure preservation.
The updated example ran successfully. General Python-runtime changes justify
checking every existing configured interface; native numerical algorithms did
not change.

The existing native `test_multivariate` program also passed after compiling the
final native sources and linking cached dependencies:

```bat
gfortran -O0 -g -fcheck=all -fbacktrace ^
  -I build\gfortran_2806A042DAA8BBF7 -J build\multivariate-parameters-regression ^
  src\stats\multivariate_normal.f90 src\stats\multivariate_t.f90 ^
  test\test_checks.f90 test\multivariate_reference.f90 test\test_multivariate.f90 ^
  build\gfortran_3C159EB802C22FB9\scifort\libscifort.a ^
  -o build\multivariate-parameters-regression\dependent.exe
build\multivariate-parameters-regression\dependent.exe
```

The cache paths above are local validation artifacts, not portable build paths.
This checks existing multivariate density/entropy/marginal, fitting and random
sampling paths without rebuilding unrelated modules. The full FPM suite was not
run. It remains required before release or merge under repository policy.

All 124 generated artifacts passed the read-only check. Previous published C
prototypes are unchanged. Encoding/BOM, ASCII Fortran, 132-column source limits,
Markdown links and whitespace were checked. Reviewed coverage is now 151 Python
methods and 63 array adapters. Other platforms and timing baselines were not
newly validated in that batch. Normal fitting is described below; see the later
[explicit RNG streams](#explicit-rng-streams-and-dirichlet-sampling) and
[CDF interfaces](#multivariate-box-cdf-and-normal-log-cdf) sections for their contracts.

## Multivariate normal fitting

`api.multivariate_normal.fit(data, fixed_mean=None, fixed_cov=None, *, out=None)`
returns `ArrayResult(status, outputs)` with `outputs['mean']` and `outputs['cov']`.
The input has shape `(n,k)`: n variables and k **observations by column**, both
positive. Outputs have shapes `(n,)` and `(n,n)`. This preserves native Fortran
orientation; SciPy fitting takes observations by row, so compare with `data.T`.
The native names `fixed_mean`/`fixed_cov` are retained in Python; SciPy calls
these parameters `fix_mean`/`fix_cov`.

With no constraints, the result is the sample mean and maximum-likelihood
covariance, divided by k rather than k-1. A supplied fixed mean is returned and
used for centering if covariance is estimated. A supplied fixed covariance is
returned after symmetrizing its lower triangle; an absent mean is still estimated
from the data. Either, both or neither parameter can be fixed. A single
observation is supported, yielding zero unconstrained covariance. Fixed
covariance may be PSD/singular; its lower triangle is authoritative, but every
entry must be finite. None remains genuinely absent in the native call, including
when the other fixed parameter is supplied.

Optional arrays must have the matching shapes when supplied. Python converts
numeric arrays to aligned column-major float64 as needed; complex inputs, shape
mismatches and empty extents raise. `out={'mean': mean_buffer, 'cov': cov_buffer}`
uses caller-owned, writable, aligned native float64 buffers with the exact shapes.
They must not overlap each other or any provided input, including fixed arrays.
The input data/fixed arrays remain unchanged. The wrapper allocates full-sized
ignored placeholders for absent arrays; present optional inputs are copied inside
the C adapter. Numerical fitting formulas and their arithmetic limits are unchanged.

```python
data = np.array([[1., 3., 5.], [2., 6., 4.]])
mean = np.empty(2)
cov = np.empty((2, 2), order='F')
result = api.multivariate_normal.fit(
    data, fixed_mean=[2., 3.], out={'mean': mean, 'cov': cov})
print(result.status, mean, cov)
```

C calls `scifort_multivariate_normal_fit_f64` with arguments in this order:
`(n, k, data, mean, cov, status, has_fixed_mean, fixed_mean, has_fixed_cov, fixed_cov)`.
Lengths and presence flags are by-value C int; status is an int pointer. Arrays
are double pointers to contiguous column-major storage; n,k must be positive.
Both fixed-array buffers must be valid, aligned, nonnull, disjoint and fully
sized even when their flags are 0. Flags are 0 absent or 1 present; other values
return status 1 with mean/cov unchanged. Invalid C extents likewise preserve
outputs. All flags are checked before any optional-array copy/native call.

Native numerical failures retain native outputs rather than promising
preservation: nonfinite data returns status 1 and NaNs in both outputs;
invalid fixed mean returns status 2 and NaNs in both; invalid fixed covariance
returns status 2, retains the computed/fixed mean, and writes NaN covariance.
Success returns status 0. Python exposes these statuses and partial outputs.
See the [native shape table](../../docs/multivariate_normal_shapes.md),
[manifest](../cabi_config/multivariate_normal.json) and
[optional-array generator policy](../../tools/cabi_generator/README.md#presence-controlled-optional-arrays).

## Normal-fitting interface validation (2026-11-25)

Validated on Windows with GNU Fortran/GCC/G++ 17.0.0 20260510 (experimental),
Python 3.13.3, NumPy 2.2.6, SciPy 1.15.3 and fparser 0.2.5. Commands:

```bat
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
python -m unittest discover -s bindings/python/tests -v
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json --check
python bindings\python\example.py
```

All 79 generator tests and all 76 binding integration tests passed. Three new
independent generator tests validate policies, C/Python signatures and true
native presence/absence for all four combinations, including ignored NaN
placeholders. Six new fit integration cases compare SciPy and independently
computed centered cross-products, plus single observations, PSD/singular fixed
covariance, lower-triangle conventions, invalid/nonfinite inputs, optional shape
checks, layout conversion, input/output alias rejection, output ownership,
invalid flags/extents and native partial-output behavior.

Integration used checked `-O0 -g -fcheck=all -fbacktrace` and optimized `-O3`
libraries, both with `-std=f2018 -fimplicit-none -Werror=implicit-interface
-shared -fPIC`. Native and C11 clients exercise all four fitting combinations;
C11/C++17 checks use `-Wall -Wextra -Werror`. The updated example passed.
All existing binding tests ran because the validator, emitter and Python array
runtime are shared by configured interfaces. No native numerical algorithm was
changed.

The existing `test_multivariate` program passed with the final native source
and cached dependent modules/library, using:

```bat
gfortran -O0 -g -fcheck=all -fbacktrace ^
  -I build\gfortran_2806A042DAA8BBF7 -J build\multivariate-parameters-regression ^
  src\stats\multivariate_normal.f90 src\stats\multivariate_t.f90 ^
  test\test_checks.f90 test\multivariate_reference.f90 test\test_multivariate.f90 ^
  build\gfortran_3C159EB802C22FB9\scifort\libscifort.a ^
  -o build\multivariate-parameters-regression\dependent.exe
build\multivariate-parameters-regression\dependent.exe
```

Those cache paths are local validation artifacts. The full FPM suite was not
run; the native edits are documentation and explicit implicit-none only, while
focused tests cover fitting and the generated interface/runtime changes. The
full suite remains required before release or merge under repository policy.

All 124 generated artifacts passed the read-only consistency check. All previous
multivariate-normal C prototypes are preserved. Encoding/BOM, Fortran ASCII and
132-column limits, local Markdown links and whitespace checks passed. Coverage
is 152 Python methods and 64 array adapters in the reviewed subset. Optional
input copies and ignored full-sized placeholders add memory overhead. Other
platforms, compiler-specific behavior and timing baselines were not newly
validated.

## Explicit RNG streams and Dirichlet sampling

**SciFort extension:** caller-owned MRG32k3a streams and snapshot/restore API.

```python
stream = api.rng(seed=12345)
snapshot = stream.state              # independent int64[6] copy
uniforms = stream.uniform(10)
stream.set_state(snapshot)           # validate and restore, preserving on error
samples = api.dirichlet.rvs_array(stream, [2., 3., 4.], n_samples=100)
print(samples.status, samples.outputs['samples'].shape)  # 0, (3, 100)
```

`RandomState(api, seed)` is also exported. Every stream owns its six-component
native MRG32k3a state and retains its library owner; no RNG is hidden in module
state. Seeds span the full signed int64 range; floats and booleans are rejected.
The serialized layout is the existing native `rng_get_state` order: first three
components in `[0,4294967087)`, last three in `[0,4294944443)`, with neither
triple all zero. Snapshots can be restored to reproduce the native stream.
`set_state` checks shape, integer range and recurrence validity before replacing
the live state. Snapshot mutations do not affect the stream. Independent stream
objects can be interleaved; do not mutate the same stream concurrently. These
are SciFort native streams, not NumPy/SciPy bit generators; equal seeds do not
promise equal sequences between libraries.

`stream.uniform(n, out=None)` returns a length-n float64 array in `(0,1)`.
A caller-owned output must be native float64, aligned, writable, contiguous,
correctly shaped and disjoint from state. n must be a nonnegative C-int extent.
An empty call validates state and leaves a valid stream unchanged. Invalid
states/counts raise without changing outputs or the live stream.

`api.dirichlet.rvs_array(state, alpha, *, n_samples=1, out=None)` returns
`ArrayResult(status, {'samples': array})`. For k concentrations the output shape
is `(k,n_samples)`, with independent simplex draws **by column**. Supply
`out={'samples': buffer}` for an aligned, writable, column-major native float64
buffer, disjoint from inputs/state. Integer n_samples is keyword-only; booleans
are rejected. It must fit nonnegative C int. A state from a different `Scifort`
instance is rejected, even when that instance loads the same library path.

Empty batches succeed without validating concentration values, including an
empty alpha vector, and do not advance a valid stream. Nonempty batches require
k >= 1 and finite positive concentrations. Invalid parameters return native
status 1, NaN output values where the shape has elements, and unchanged state.
A native numerical generation failure can advance state: completed columns are
retained and the failed/remaining columns are NaN. For example, very small
positive concentrations can underflow all inverse-Gamma draws to zero. This
behavior is preserved, not replaced by a new sampler. Typical successful draws
consume k native uniforms per column; numerical failures stop at the native
failure point. Both checked and optimized builds use the same native stream.

C includes [scifort_rng.h](../../include/scifort_rng.h) and calls:

```c
void scifort_rng_seed_i64(int64_t seed, int64_t values[6], int *status);
void scifort_rng_uniform_f64(int n, int64_t values[6], double *output, int *status);
void scifort_dirichlet_rvs_array_f64(int k, int n_samples, int64_t *state,
                                   double *samples, const double *alpha, int *status);
```

All C buffers and scalar output pointers must be valid, aligned, contiguous,
fully sized and disjoint, including for empty calls (use nonnull sentinels for
empty arrays). State is always validated, even when n/n_samples is zero. Seed
returns 0. Uniform statuses are 0 success, 1 invalid extent, 3 invalid state.
Dirichlet statuses are 0 success, 1 native parameter/numerical failure, 2 invalid
extent/shape, 3 invalid state. Extent/state errors preserve buffers and state;
native numerical errors retain native partial outputs and advanced state.
No Fortran descriptors, derived-type layout, opaque handles or mutable global
state cross this ABI. The six-component MRG32k3a format is an ABI commitment;
a future different RNG needs a separately documented format/interface.

See the [generated sampler shapes](../../docs/dirichlet_shapes.md) and
[reviewed codec policy](../../tools/cabi_generator/README.md#reviewed-serialized-state-codecs).
Multivariate CDFs are described in the later box-integration section; additional
distribution sampler bindings remain future work.

## RNG and sampling interface validation (2026-11-26)

Validated on Windows with GNU Fortran/GCC/G++ 17.0.0 20260510 (experimental),
Python 3.13.3, NumPy 2.2.6 and fparser 0.2.5. Commands:

```bat
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
python -m unittest discover -s bindings/python/tests -p test_random.py -v
python -m unittest discover -s bindings/python/tests -v
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json --check
python bindings\python\example.py
```

All 82 generator tests passed, including three independent two-component state
codec/output-dimension tests. All seven new RNG/sampler integration tests passed
in a focused run, then all 83 binding tests passed together. Integration builds
use checked `-O0 -g -fcheck=all -fbacktrace` and optimized `-O3`, both with
`-std=f2018 -fimplicit-none -Werror=implicit-interface -shared -fPIC`.
C11 clients and C++17 header checks use `-Wall -Wextra -Werror`.

Native/C clients compare exact uniform draws, final recurrence components and
Dirichlet draws against the actual native procedures, including zero, negative,
large and signed-int64 endpoint seeds. Python checks stream splitting, snapshots,
restoration, independent ownership, open-interval uniforms, state/buffer validation,
empty batches and exact native uniform consumption. Dirichlet moments use analytic
means/variances with fixed-seed, six-standard-error mean bounds; simplex columns
sum to one. Tiny positive concentrations exercise native numerical underflow:
state is advanced through the failing draw and exported, while remaining output
columns become NaN. Invalid parameters do not advance state; invalid serialized
state preserves outputs. All configured bindings ran because the Python runtime
and generator are shared. The updated example also passed.

Existing native RNG reference-sequence/variates and multivariate sampling tests
passed using the updated Dirichlet source and cached native dependencies:

```bat
gfortran -O0 -g -fcheck=all -fbacktrace ^
  -I build\gfortran_2806A042DAA8BBF7 -J build\rng-c-api-regression ^
  src\stats\dirichlet.f90 test\test_random.f90 ^
  build\gfortran_3C159EB802C22FB9\scifort\libscifort.a ^
  -o build\rng-c-api-regression\random.exe
build\rng-c-api-regression\random.exe
gfortran -O0 -g -fcheck=all -fbacktrace ^
  -I build\gfortran_2806A042DAA8BBF7 -J build\multivariate-parameters-regression ^
  src\stats\dirichlet.f90 test\test_checks.f90 test\multivariate_reference.f90 ^
  test\test_multivariate.f90 ^
  build\gfortran_3C159EB802C22FB9\scifort\libscifort.a ^
  -o build\multivariate-parameters-regression\dependent.exe
build\multivariate-parameters-regression\dependent.exe
```

These cache paths are local validation artifacts, not portable build paths.
The full FPM suite was not run: native RNG/sampling algorithms are unchanged,
while the bounded interface change is covered by native stream comparisons and
all configured binding integrations. The full suite remains required before
release or merge. Other platforms, RNG algorithms, concurrency of a shared
stream and timing baselines were not newly validated.

Final regeneration updated only generated argument comments; the sampler's
executable-source hash (excluding comments/blank lines) matches the tested
version. No numerical tests were repeated for those documentation-only edits.
Previous Dirichlet C prototypes remain unchanged. Coverage is now 153 generated
Python methods and 65 array adapters; two manual RNG C entry points and the
Python stream factory are additional interfaces outside that generated count.

The final read-only check passed for all 124 generated artifacts. Encoding/BOM,
Fortran ASCII and 132-column limits, local Markdown links and whitespace checks
also passed. The new manual RNG header is included by the public umbrella header;
C/C++ integration clients exercised all three new C entry points.

## Integer-count batch sampling

```python
state = api.rng(123)
counts = api.multinomial.rvs_array(state, 12, [0.2, 0.3, 0.5], n_samples=100)
draws = api.multivariate_hypergeom.rvs_array(state, [10, 20, 30], 12, n_samples=100)
assert counts.status == draws.status == 0
assert counts.outputs['samples'].shape == (3, 100)
```

Both methods return `ArrayResult(status, {'samples': array})`. Each column is
one count vector and sums to n on success. Counts are nonnegative; hypergeometric
counts also cannot exceed the corresponding population entry. Output arrays
have native-endian `np.intc` dtype (C int), shape `(k,n_samples)` and Fortran
column-major layout. Supply `out={'samples': buffer}` to reuse aligned, writable,
exactly shaped `np.intc` storage disjoint from all inputs and state. Float64 and
int64 output buffers are rejected when they do not match native C int.
`n_samples` defaults to 1, is keyword-only, and must be a nonnegative C integer;
booleans are rejected. Integer population inputs accept signed/unsigned integer
arrays after range checks; floats and booleans are rejected.

Both use the [explicit RNG stream](#explicit-rng-streams-and-dirichlet-sampling)
and preserve the native algorithm and sample order. Splitting a batch preserves
the stream. Uniform consumption depends on remaining counts and parameters;
do not assume a fixed number of uniforms per column. NumPy/SciPy seed streams
are not promised to match. Snapshot/restore, library ownership and alias rules
are the same as for Dirichlet. No state is hidden or shared globally.

For multinomial, k=size(p) may be zero for storage, but nonempty batches require
k>=1, n>=0 and valid adjusted probabilities. The native final-entry probability
adjustment is preserved. Empty batches succeed without validating n or p values.
Hypergeometric requires k=size(m)>=1, m>=0, and 0<=n<=sum(m), even when the batch
is empty. Valid empty batches and rejected parameters do not advance state.
A single category or zero draws also consumes no uniforms. Parameter failures
return status 1 and zero all sample entries. Empty multinomial p with a nonempty
batch returns status 1.

Additional ABI safety limits always apply, including empty batches: n must be
in [-2147483648,2147483646]; hypergeometric requires sum(abs(m))<=2147483646,
computed in widened arithmetic before native default-integer sums. Python rejects
violations before calling C. C returns status 2, preserving sample buffers and
state. Invalid extents also return 2; invalid state returns 3 unchanged. State is
validated even for empty calls, after extent/integer guards. Successful calls
return 0. Every native call exports its resulting state, including on failure.

```c
void scifort_multinomial_rvs_array_i32(int k, int n_samples, int64_t *state,
                                     int *samples, int n, const double *p, int *status);
void scifort_multivariate_hypergeom_rvs_array_i32(int k, int n_samples, int64_t *state,
                                               const int *m, int n, int *samples, int *status);
```

The `_i32` names identify the reviewed integer-count ABI. Buffers use C int and
require the native default integer kind to match `c_int`. All C pointers must be
valid, aligned, disjoint, nonnull and fully sized; empty arrays need nonnull
sentinels. C callers supply column-major storage. See the generated
[multinomial shapes](../../docs/multinomial_shapes.md#multinomial_rvs_array) and
[hypergeometric shapes](../../docs/multivariate_hypergeom_shapes.md#multivariate_hypergeom_rvs_array).

## Integer-output sampler validation (2026-11-26)

Validation uses GNU Fortran/GCC/G++ 17.0.0 20260510 (experimental), Python 3.13.3,
NumPy 2.2.6 and fparser 0.2.5 on Windows. Commands:

```bat
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
python -m unittest discover -s bindings/python/tests -p test_random.py -v
python -m unittest discover -s bindings/python/tests -v
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json --check
```

The independent integer-output fixture compiles rank-one and rank-two outputs
for C11/Python clients, including empty arrays, preserved guard outputs and
invalid NumPy storage. Production integration checks exact native/C/Python
samples and final state, split-batch equivalence, signed seed endpoints,
nonnegative counts and totals, population bounds, analytic moments, invalid
parameters and output/state preservation. Checked builds use
`-O0 -g -fcheck=all -fbacktrace`; optimized builds use `-O3`. Both use
`-std=f2018 -fimplicit-none -Werror=implicit-interface -shared -fPIC`.
C11 and C++17 checks use `-Wall -Wextra -Werror`.

The dependent native multivariate regression was compiled against the existing
local FPM archive, with the two edited native modules compiled directly:

```bat
gfortran -O0 -g -fcheck=all -fbacktrace ^
  -I build\gfortran_2806A042DAA8BBF7 -J build\multivariate-parameters-regression ^
  src\stats\multinomial.f90 src\stats\multivariate_hypergeom.f90 ^
  test\test_checks.f90 test\multivariate_reference.f90 test\test_multivariate.f90 ^
  build\gfortran_3C159EB802C22FB9\scifort\libscifort.a ^
  -o build\multivariate-parameters-regression\dependent.exe
build\multivariate-parameters-regression\dependent.exe
```

These cache paths are local artifacts. The full FPM suite was not run: numerical
algorithms are unchanged, while the bounded generator/runtime change is checked
against every configured binding and independent compiled fixtures. The full
suite remains required before merge or release. Other platforms were not newly
validated. Existing C symbols are unchanged; generated coverage now includes
155 Python methods and 67 array adapters.

All 84 generator tests and all 86 binding integration tests passed, including
ten RNG/sampler checks and the updated density/moment C-guard regression. The
native multivariate regression passed. The final read-only pipeline check passed
for all 124 artifacts; encoding/BOM, Fortran ASCII/line limits, added local links
and whitespace checks passed. No numerical tests were repeated for the final
documentation-only edits.

## Multivariate vector batch sampling

```python
state = api.rng(123)
normal = api.multivariate_normal.rvs_array(
    state, mean=[1., 2.], cov=[[4., 1.], [1., 3.]], n=2, n_samples=100)
student = api.multivariate_t.rvs_array(
    state, loc=[1., 2.], shape=[[4., 1.], [1., 3.]], df=5, n=2, n_samples=100)
assert normal.status == student.status == 0
```

The signatures are:

```python
api.multivariate_normal.rvs_array(state, mean=None, cov=None, *, n=1, n_samples=1, out=None)
api.multivariate_t.rvs_array(state, loc=None, shape=None, df=1.0, *, n=1, n_samples=1, out=None)
```

Both return `ArrayResult(status, {'samples': array})`. Outputs are native-endian
float64, shape `(n,n_samples)`, with one vector per column. Explicit n defines
the dimension even when parameters are absent; it is not inferred from supplied
parameters or out buffers. Pass `n=2` for two-dimensional parameters. Both extent
keywords require nonnegative C integers and reject booleans; normal requires
positive n even for empty batches. Nonempty Student t batches require positive
n; empty Student t batches permit n=0. `out={'samples': buffer}` requires an
exactly shaped, aligned, writable, Fortran-contiguous float64 buffer disjoint
from inputs and state. Invalid Python storage raises before native execution.

Absent mean/loc selects native zeros, and absent cov/shape selects native
identity. None preserves native optional-array absence through presence flags;
it does not create an explicit default parameter array. Present vectors require
shape `(n)` and matrices `(n,n)`. All present array entries must be finite,
including the unused upper matrix triangle. The lower triangle is authoritative;
no symmetry check is performed. Singular positive-semidefinite matrices,
including all-zero matrices, are allowed using the existing native PSD cutoff.
Student t shape is not covariance: covariance is `df/(df-2)*shape` when df>2;
the mean exists for df>1. df must be positive or positive infinity. Python uses
the native default 1; C always supplies df explicitly. df=None is rejected.

These use the [explicit stream API](#explicit-rng-streams-and-dirichlet-sampling).
Successful normal sampling consumes n uniforms per column even for rank-zero
covariance. Student t consumes n+1 for finite df, n for infinite df. The
infinite-df samples agree with normal sampling up to rounding, and final state
matches exactly with the same native parameters and initial state. Native normal
and Student t add location at different stages. Splitting successful batches preserves stream
order. NumPy/SciPy seed sequences are not promised to match.

| Status | Meaning | Output and state behavior |
| ---: | --- | --- |
| 0 | Success | Sample columns and advanced native state |
| 1 | Invalid extent, shape or presence flag | Transport rejection preserves caller outputs and state |
| 2 | Invalid native parameters or Student t scale generation failure | Entire batch NaN; state retained at the native failure point |
| 3 | Native PSD decomposition failure | Entire batch NaN; state unchanged |
| 4 | Invalid serialized RNG state | Output and state unchanged |

Normal validates parameters and PSD covariance even for zero columns. Student t
returns success without distribution validation for zero columns, including
invalid df or matrix values. Both always validate state, extents, supplied array
shapes and C presence flags. Valid empty calls do not advance state. Native
parameter rejection occurs before sampling. Student t can fail after consuming
uniforms, for example when a tiny positive df underflows the chi-square scale;
then it fills the **entire** batch with NaN and exports the advanced state.
Thus completed columns are also discarded on a later native generation failure.
Student t n=0 with nonzero columns returns native status 1 with no output elements.
These existing numerical behaviors are preserved.

```c
void scifort_multivariate_normal_rvs_array_f64(
    int n, int n_samples, int64_t *state, double *samples,
    int has_mean, const double *mean, int has_cov, const double *cov, int *status);
void scifort_multivariate_t_rvs_array_f64(
    int n, int n_samples, int64_t *state, double *samples,
    int has_loc, const double *loc, int has_shape, const double *shape, double df, int *status);
```

Flags must be 0 absent or 1 present. All pointers remain required, nonnull,
aligned, disjoint and fully sized, even when a flag is zero or a batch is empty.
Absent parameters' buffer values are ignored. Use nonnull scratch buffers for
zero-size arrays. Arrays use column-major storage; no Fortran descriptors or
RNG object layout crosses C. State belongs to the loaded library instance and
must not be mutated concurrently. See the generated
[normal shapes](../../docs/multivariate_normal_shapes.md#multivariate_normal_rvs_array)
and [Student t shapes](../../docs/multivariate_t_shapes.md#multivariate_t_rvs_array).

## Vector sampler validation (2026-11-26)

Validation uses GNU Fortran/GCC/G++ 17.0.0 20260510 (experimental), Python 3.13.3,
NumPy 2.2.6 and fparser 0.2.5 on Windows. Commands:

```bat
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
python -m unittest discover -s bindings/python/tests -p test_vector_sampling.py -v
python -m unittest discover -s bindings/python/tests -v
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json --check
```

The independent counter-state fixture combines state codecs, optional-array
presence and promoted scalar inputs. Production clients compare exact native/C/
Python sample columns and final state for all four optional-array combinations,
finite/infinite df and native omitted df. Integration also checks split batches,
normal limits, analytic moments, singular/zero matrices, lower-triangle use,
empty-batch differences, tiny-df underflow, guard preservation and output buffers.
Checked builds use `-O0 -g -fcheck=all -fbacktrace`; optimized builds use `-O3`.
Both use `-std=f2018 -fimplicit-none -Werror=implicit-interface -shared -fPIC`.
Strict C11 clients and C++17 header checks use `-Wall -Wextra -Werror`.

The dependent native multivariate regression compiled both edited native modules
directly and reused the existing local FPM dependency archive:

```bat
gfortran -O0 -g -fcheck=all -fbacktrace ^
  -I build\gfortran_2806A042DAA8BBF7 -J build\multivariate-parameters-regression ^
  src\stats\multivariate_normal.f90 src\stats\multivariate_t.f90 ^
  test\test_checks.f90 test\multivariate_reference.f90 test\test_multivariate.f90 ^
  build\gfortran_3C159EB802C22FB9\scifort\libscifort.a ^
  -o build\multivariate-parameters-regression\dependent.exe
build\multivariate-parameters-regression\dependent.exe
```

These cache paths are local artifacts. The full FPM suite was not run: numerical
algorithms are unchanged, and the generator/interface effects are bounded by
compiled fixtures, all configured binding integrations and the native dependent
test. The full suite remains required before merge or release. Other platforms
were not newly validated. Existing C symbols are unchanged; generated coverage
now includes 157 Python methods and 69 array adapters. Optional present arrays
are copied by the existing bridge, and native Student t still factors the shape
for each column; this change does not optimize either implementation.

All 86 generator tests and all 92 binding tests passed, including the six new
sampler tests in checked and optimized builds and strict native/C/C++ clients.
The first focused run exposed an overly strict cross-family infinite-df sample
comparison (maximum difference 8.9e-16); the corrected test allows rounding while
requiring exact final state. Each method's native/C/Python stream comparison
remains exact. The native dependent multivariate regression also passed.
The final read-only check passed for all 124 generated artifacts. Encoding/BOM,
Fortran ASCII/line limits, added local links and whitespace checks passed.
No numerical tests were repeated for the final documentation-only edits.

## Wishart and inverse-Wishart batch sampling

```python
state = api.rng(123)
scale = [[2., .4], [.4, 1.]]
w = api.wishart.rvs_array(state, 7., scale, n_samples=100)
iw = api.invwishart.rvs_array(state, 7., scale, n_samples=100)
assert w.status == iw.status == 0
assert w.outputs['samples'].shape == (2, 2, 100)
```

Both signatures are `rvs_array(state, df, scale, *, n_samples=1, out=None)`.
They return `ArrayResult(status, {'samples': array})`. p is inferred from square
scale; sample shape is `(p,p,n_samples)`. Each matrix occupies one slice along
the third axis, with native-endian float64 Fortran layout. C storage places all
p*p entries of one matrix before the next. NumPy/SciPy commonly place the sample
axis first; use `np.moveaxis(samples,2,0)` when that orientation is needed.
`out={'samples': buffer}` requires exact shape/dtype, alignment, writeability,
Fortran contiguity and disjointness from scale/state. n_samples must be a
nonnegative C integer; floats and booleans are rejected. Inputs may be converted,
while caller outputs must already meet the storage contract.

Successful nonempty draws require p>=1, finite df>p-1 and finite SPD scale.
The lower scale triangle is authoritative; all entries, including the ignored
upper triangle, must be finite. Singular scale is rejected. The native Bartlett
Wishart and Wishart-inversion inverse-Wishart algorithms remain in use. They
share the [explicit RNG bridge](#explicit-rng-streams-and-dirichlet-sampling),
consume p*(p+1)/2 uniforms per successful matrix and preserve stream order when
splitting batches. Native/C/Python samples and final state agree exactly for
the same sampler, parameters and starting state. NumPy/SciPy seeded streams do
not match by contract. Snapshot/restore and library ownership rules are unchanged.

Empty batches skip native distribution validation, so invalid df or scale values
can return success without advancing a valid stream. Native code only checks
leading output extents on this path; C/Python additionally require square scale
and matching tensor shapes even when empty. Zero-dimensional `(0,0)` scale is
accepted for an empty batch; it returns native shape status 1 for a nonempty
batch, with no sample elements. All state/extent/buffer rules still apply.

| Status | Meaning | Output/state behavior |
| ---: | --- | --- |
| 0 | Native success | Generated matrices and advanced state |
| 1 | Invalid extent or shape | Transport rejection preserves caller output/state |
| 2 | Invalid native df or scale values | Failed matrix and all later matrices NaN; no draws consumed |
| 3 | Native SPD/solve failure | Completed matrices retained; failed and later matrices NaN; consumed state exported |
| 4 | Invalid serialized RNG state | Output/state unchanged |

Both native batch routines now explicitly fill the previously unwritten
remainder with NaN when a scalar draw fails. Successful sample generation is
unchanged. Native parameter validation failures occur before RNG consumption.
Inverse-Wishart may fail after drawing its internal Wishart matrix; state then
includes those uniforms. Existing extreme-df limitations remain: for p=1 and a
tiny positive df, chi-square inversion can produce a nonfinite draw. At df=1e-300
the checked and optimized tests observe NaN Wishart output with native status 0,
while inverse-Wishart reports status 3 when the subsequent factorization fails. These interfaces preserve
that numerical behavior; status 0 is not a universal finite/SPD guarantee for
extreme parameters. Ordinary-parameter tests check SPD, symmetry and moments.

```c
void scifort_wishart_rvs_array_f64(int p, int n_samples, int64_t *state,
                                 double df, const double *scale, double *samples, int *status);
void scifort_invwishart_rvs_array_f64(int p, int n_samples, int64_t *state,
                                    double df, const double *scale, double *samples, int *status);
```

All pointers must be valid, aligned, disjoint, nonnull and fully sized, including
empty calls; use nonnull scratch storage for zero-size arrays. State uses six
int64 components. Fortran descriptors stay inside the adapter. See the generated
[Wishart shapes](../../docs/wishart_shapes.md#wishart_rvs_array) and
[inverse-Wishart shapes](../../docs/invwishart_shapes.md#invwishart_rvs_array).

## Wishart sampler validation (2026-11-26)

Validation uses GNU Fortran/GCC/G++ 17.0.0 20260510 (experimental), Python 3.13.3,
NumPy 2.2.6 and fparser 0.2.5 on Windows. Commands:

```bat
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -p test_rank_three_outputs.py -v
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
python -m unittest discover -s bindings/python/tests -p test_wishart_sampling.py -v
python -m unittest discover -s bindings/python/tests -v
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json --check
```

The independent rank-three fixture checks third-axis extent parsing, tensor
layout, C/Python buffer ownership, zero extents and rejection of unsupported
rank/type/intent policies. Sampler integration compares exact native/C/Python
samples and final state, split and one-matrix batches, lower-triangle use,
symmetry/SPD, analytic mean/entrywise variance, empty and invalid parameters,
state/buffer guards and defined failure remainders. A deterministic later
inverse-Wishart failure checks retained completed matrices and advanced state:
seed 2**40+123, ten initial uniforms, df=0.0001, unit one-dimensional scale.
The first matrix succeeds and the second fails. A separate tiny-df check records
the existing Wishart NaN/status-0 limitation rather than assuming a zero draw.

Checked builds use `-O0 -g -fcheck=all -fbacktrace`; optimized builds use `-O3`.
Both use `-std=f2018 -fimplicit-none -Werror=implicit-interface -shared -fPIC`.
Native/C clients are compiled and run for both builds. C11 and C++17 header
checks use `-Wall -Wextra -Werror`. The native matrix-distribution regression
also compiled the edited modules directly, reusing local FPM dependencies:

```bat
gfortran -O0 -g -fcheck=all -fbacktrace ^
  -I build\gfortran_2806A042DAA8BBF7 -J build\multivariate-parameters-regression ^
  src\stats\wishart.f90 src\stats\invwishart.f90 test\test_checks.f90 ^
  test\test_matrix_distributions_reference.f90 test\test_matrix_distributions.f90 ^
  build\gfortran_3C159EB802C22FB9\scifort\libscifort.a ^
  -o build\multivariate-parameters-regression\matrix.exe
build\multivariate-parameters-regression\matrix.exe
```

These cache paths are local artifacts. The full FPM suite was not run: the native
change affects only previously unwritten batch-error remainders, while compiled
generator fixtures, all configured bindings and the native matrix-distribution
regression cover the affected interfaces and dependencies. The full suite
remains required before merge or release. Other platforms were not newly
validated. Existing C prototypes and successful sampler algorithms are unchanged;
coverage now includes 159 Python methods and 71 array adapters. Matrix-normal
and matrix-t batch samplers are described in the later sampling section.

All 88 generator tests and all 98 binding tests passed. A subsequent focused run
of all six sampler tests passed with the additional deterministic late-failure
native/C/Python checks, in checked and optimized builds. The native dependent
matrix-distribution regression passed. All 124 generated artifacts passed the
read-only consistency check. Encoding/BOM, Fortran ASCII/line limits, added local
links, preserved existing C prototypes and whitespace checks passed. No numerical
tests were repeated for final documentation-only edits.

## Matrix-normal and matrix-t batch sampling

```python
state = api.rng(123)
mean = [[1., 2., 3.], [4., 5., 6.]]
row = [[2., .3], [.3, 1.]]
col = [[1., .1, .2], [.1, 2., .3], [.2, .3, 3.]]
normal = api.matrix_normal.rvs_array(state, mean, row, col, m=2, n=3, n_samples=100)
student = api.matrix_t.rvs_array(state, mean, row, col, df=12, m=2, n=3, n_samples=100)
assert normal.status == student.status == 0
```

```python
api.matrix_normal.rvs_array(state, mean=None, rowcov=None, colcov=None,
                           *, m=1, n=1, n_samples=1, out=None)
api.matrix_t.rvs_array(state, mean=None, row_spread=None, col_spread=None, df=1.0,
                      *, m=1, n=1, n_samples=1, out=None)
```

Both return `ArrayResult(status, {'samples': array})`. Outputs have shape
`(m,n,n_samples)`, with one rectangular matrix per third-axis slice, native-endian
float64 and Fortran-contiguous layout. m/n are explicitly requested row/column
dimensions; neither is inferred from supplied parameters or out buffers. All
three extent keywords default to 1, accept nonnegative C integers and reject
booleans. Positive m and n are required for successful nonempty draws.
`out={'samples': buffer}` requires exact shape/dtype, alignment, writeability,
Fortran contiguity and disjointness from inputs/state. Empty array pointers still
use valid nonnull storage. `np.moveaxis(samples,2,0)` moves the sample axis first.

Optional mean has shape `(m,n)`, the row matrix `(m,m)` and the column matrix
`(n,n)`. None preserves native absence: zero mean and identity row/column matrices.
All present entries must be finite, including ignored upper triangles. Row and
column matrices must be positive definite; singular matrices are rejected.
Their lower triangles define symmetry; no explicit symmetry check is performed.
Python validates ranks/shapes/buffers before C. Numerical parameter failures
remain native statuses with native output values.

For column-major vectorization, matrix-normal covariance is
`kron(colcov,rowcov)`. Matrix-t uses the existing unscaled inverse-Wishart mixture:
row covariance is drawn from `InvWishart(df+m-1,row_spread)`, followed by a
matrix-normal draw with col_spread. Thus for df>2 the vectorized covariance is
`kron(col_spread,row_spread)/(df-2)`; the mean exists for df>1. This differs from
multivariate Student t's scale convention. Matrix-t df must be finite and positive;
positive infinity is invalid. Python uses native default 1; C always supplies df.

Both reuse the [explicit RNG stream](#explicit-rng-streams-and-dirichlet-sampling).
Normal sampling consumes m*n uniforms per matrix. Successful matrix-t sampling
consumes m*(m+1)/2+m*n per matrix. Splitting successful batches preserves samples
and final state. Native/C/Python streams agree exactly for the same method,
parameters and starting state; NumPy/SciPy seed streams are not promised to match.
Each stream belongs to its loaded library instance and must not be mutated
concurrently. Snapshot/restore and ownership rules are unchanged.

Empty batches return native success without shape or parameter validation,
including zero leading extents and invalid matrix-t df. The C/Python contracts
still validate requested extents, present array shapes, buffers, flags and state.
A nonempty batch with a zero leading extent returns native status 1 and has no
output elements. Neither empty calls nor rejected transport arguments advance
valid state.

| Status | Meaning | Output/state behavior |
| ---: | --- | --- |
| 0 | Success | Sample matrices and advanced state |
| 1 | Invalid extent/shape/flag; also matrix-normal nonfinite mean | Transport errors preserve buffers/state; native nonfinite mean gives NaNs |
| 2 | Native invalid parameter values or matrix-t df | Failed and later matrices NaN; no draws consumed |
| 3 | Native SPD/solve failure | Completed matrices retained; failed and later matrices NaN; consumed state exported |
| 4 | Invalid serialized RNG state | Output and state unchanged |

Both native batch procedures now define previously unwritten remainder slices
as NaN when a scalar draw fails. Completed matrices are retained. Matrix-normal
checks row/column SPD before drawing; its parameter/factorization failures do
not advance state. Matrix-t checks row SPD through the inverse-Wishart sampler,
then checks column SPD after generating that row covariance. A non-SPD column
matrix can therefore advance state by m*(m+1)/2 uniforms before returning 3.
Nonfinite matrix-t mean returns 2; native matrix-normal uses 1 for nonfinite mean.
Those existing status conventions are preserved. Tiny positive matrix-t df can
cause an inverse-Wishart generation failure after RNG consumption; status 3,
NaN remainder and advanced state are preserved. Singular row/column spreads are
not supported, and this interface does not replace the native numerical algorithms.

```c
void scifort_matrix_normal_rvs_array_f64(
    int m, int n, int n_samples, int64_t *state, double *samples,
    int has_mean, const double *mean, int has_rowcov, const double *rowcov,
    int has_colcov, const double *colcov, int *status);
void scifort_matrix_t_rvs_array_f64(
    int m, int n, int n_samples, int64_t *state, double *samples,
    int has_mean, const double *mean, int has_row_spread, const double *row_spread,
    int has_col_spread, const double *col_spread, double df, int *status);
```

Each flag is 0 absent or 1 present. All C pointers remain required, nonnull,
aligned, disjoint and fully sized even when a flag is zero or a batch is empty.
Absent inputs' buffer values are ignored. Matrices/tensors use column-major
storage and state uses six int64 components. No descriptors or derived-type
layout cross C. See [matrix-normal shape tables](../../docs/matrix_normal_shapes.md#matrix_normal_rvs_array)
and [matrix-t shape tables](../../docs/matrix_t_shapes.md#matrix_t_rvs_array).

## Matrix sampler validation (2026-11-26)

Validation uses GNU Fortran/GCC/G++ 17.0.0 20260510 (experimental), Python 3.13.3,
NumPy 2.2.6 and fparser 0.2.5 on Windows. Commands:

```bat
pushd tools\cabi_generator\tests
..\..\..\build\cabi-generator-venv\Scripts\python.exe -m unittest test_state_codecs test_rank_three_outputs -v
popd
python -m unittest discover -s bindings/python/tests -p test_matrix_sampling.py -v
python -m unittest discover -s bindings/python/tests -v
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json --check
```

The generator implementation is unchanged; seven focused state-codec/optional-
scalar and rank-three fixture tests check the capabilities reused here. Sampler
integration checks all eight optional-array combinations per family, genuine
absence and explicit zero/identity agreement, native/C/Python samples and final
state, split batches, native omitted df versus explicit df=1, and all three C
presence flags. Native clients use rectangular 2-by-3 samples. Moment checks
include the full vectorized row/column covariance, with covariance errors scaled
by marginal standard deviations. Error tests cover empty/zero extents, invalid
df/mean/SPD matrices, buffer ownership, nonfinite ignored C buffers, tiny-df
mixture failures and column-factorization failure after RNG consumption.

A deterministic later matrix-t failure (seed 2**40+123, ten initial uniforms,
one row/column, df=0.0001) checks retained completed output, NaN remainder and
advanced state in native/C/Python clients. Checked builds use
`-O0 -g -fcheck=all -fbacktrace`; optimized builds use `-O3`. Both use
`-std=f2018 -fimplicit-none -Werror=implicit-interface -shared -fPIC`.
C11 and C++17 checks use `-Wall -Wextra -Werror`.

The dependent native matrix-distribution regression compiles both edited modules
and reuses existing local FPM dependency objects:

```bat
gfortran -O0 -g -fcheck=all -fbacktrace ^
  -I build\gfortran_2806A042DAA8BBF7 -J build\multivariate-parameters-regression ^
  src\stats\matrix_normal.f90 src\stats\matrix_t.f90 test\test_checks.f90 ^
  test\test_matrix_distributions_reference.f90 test\test_matrix_distributions.f90 ^
  build\gfortran_3C159EB802C22FB9\scifort\libscifort.a ^
  -o build\multivariate-parameters-regression\matrix.exe
build\multivariate-parameters-regression\matrix.exe
```

These cache paths are local artifacts. The full FPM suite was not run: only the
two batch procedures' previously unwritten failure remainders change numerical
behavior, while focused fixtures, all configured binding integrations and the
native matrix-distribution regression cover the affected paths. The generator
and shared build/test process are unchanged. The full suite remains required
before merge or release. Other platforms were not newly validated. Successful
sampling algorithms and existing C prototypes are unchanged; the matrix-t SciPy
attribution and license remain intact. Generated coverage now includes 161
Python methods and 73 array adapters.

All seven focused generator tests, all six focused sampler tests and all 104
binding integration tests passed. The full binding run includes the added late-
failure native/C/Python case, checked and optimized libraries and strict C/C++
clients. The native dependent matrix-distribution regression passed. All 124
generated artifacts passed the read-only consistency check. Existing prototypes,
retained notices, bounded native changes, UTF-8/BOM, Fortran ASCII/line lengths,
added local links and whitespace checks passed. No numerical tests were repeated
for final documentation-only edits.

## Multivariate box CDF and normal log-CDF

```python
state = api.rng(123)
cov = [[1., .4], [.4, 2.]]
p = api.multivariate_normal.cdf(state, [.5, .3], cov=cov, maxpts=4000)
state.set_state(api.rng(123).state)
logp = api.multivariate_normal.logcdf(state, [.5, .3], cov=cov, maxpts=4000)
pt = api.multivariate_t.cdf(api.rng(123), [.5, .3], shape=cov, df=5, maxpts=4000)
```

The generated methods return float:

```python
api.multivariate_normal.cdf(state, x, mean=None, cov=None, allow_singular=False,
                           maxpts=0, abseps=1e-5, releps=1e-5, lower_limit=None)
api.multivariate_normal.logcdf(state, x, mean=None, cov=None, allow_singular=False,
                              maxpts=0, abseps=1e-5, releps=1e-5, lower_limit=None)
api.multivariate_t.cdf(state, x, loc=None, shape=None, df=1.0,
                      allow_singular=False, maxpts=0, lower_limit=None)
```

n=size(x) must be positive. x, location and lower_limit have shape `(n)`;
cov/shape has shape `(n,n)`. Optional arrays preserve genuine native absence:
zero location, identity matrix and negative-infinity lower bounds. The lower
matrix triangle is authoritative, but all matrix entries must be finite.
Location must be finite. Bounds may be infinite, but NaNs are invalid. PSD
singularity requires allow_singular=True; indefinite matrices are invalid.
Student t df must be positive or positive infinity, which delegates to normal
integration. Shape is not covariance. Native scalar parameters are required at
C and get explicit Python defaults. Python validates ranks, shapes, dtypes,
C-integer ranges and disjointness before C; no broadcasting is performed.

**SciFort interface convention:** state is explicitly required and comes first
in Python. Native Fortran still accepts absent state (fixed local seeds), and
SciPy uses its own RNG convention. The existing [stream API](#explicit-rng-streams-and-dirichlet-sampling)
provides snapshot/restore and library ownership. C always validates the six-int64
state, even when an analytic path consumes no random draws. One-dimensional
CDFs are analytic and leave state unchanged. Normal log-CDF with diagonal
covariance also leaves state unchanged and sums stable interval logarithms,
remaining finite when independent tail probabilities underflow. Multivariate
normal CDF generally still uses QMC for diagonal matrices. Correlated log-CDF
logs the ordinary probability estimate, so underflow can still produce -infinity.

Integration uses the existing randomized, tent-transformed Halton implementation,
not SciPy's lattice sequence. Native/C/Python values and final state agree for
the same method, controls and initial state. Snapshot/restore reproduces a call;
repeated calls on one advancing stream can differ. State advancement depends on
dimension and adaptive QMC batches, not a fixed number of draws per CDF call.
Native changes to stopping decisions or numerical platforms may change streams.
Matching SciPy seeds does not imply matching values or RNG state.

Nonpositive maxpts selects the native dimension-dependent default: 1000000*n for
normal, 1000*n for Student t, capped at the native integer maximum. Positive
maxpts requests a point budget. The backend floors that budget at ten points;
it is approximate and analytic paths need none. Normal abseps must be finite
and positive, but is an error target, not a certified bound. Student t CDF uses
the backend's default absolute target 1e-5. Normal releps must be positive and
not NaN (positive infinity is accepted), but currently does not affect stopping.
CDFs do not return the backend's error estimate, point count or a convergence
status. Exhausting a budget can still return a finite approximate result and
call status 0; select suitable budgets and validate numerical accuracy for the
intended calculation. Setting maxpts=0 follows native defaults rather than
adding a smaller wrapper budget.

Reversed bounds swap each reversed coordinate and multiply by one negative sign
per reversal. Thus CDF may return a negative signed box value. Normal log-CDF
returns NaN for negative mass and -infinity for zero mass. The API retains native
NaN/infinity results for invalid distributions or integration failures. Inputs
can be converted; arrays must not overlap one another or state. Scalar function
methods expose no out mapping.

```c
void scifort_multivariate_normal_cdf_f64(
    int n, const double *x, int has_mean, const double *mean,
    int has_cov, const double *cov, int allow_singular, int maxpts,
    double abseps, double releps, int has_lower_limit, const double *lower_limit,
    int64_t *state, double *value, int *status);
void scifort_multivariate_normal_logcdf_f64(
    int n, const double *x, int has_mean, const double *mean,
    int has_cov, const double *cov, int allow_singular, int maxpts,
    double abseps, double releps, int has_lower_limit, const double *lower_limit,
    int64_t *state, double *value, int *status);
void scifort_multivariate_t_cdf_f64(
    int n, const double *x, int has_loc, const double *loc,
    int has_shape, const double *shape, double df, int allow_singular, int maxpts,
    int has_lower_limit, const double *lower_limit, int64_t *state,
    double *value, int *status);
```

Every pointer must be valid, aligned, nonnull, fully sized and disjoint, including
absent optional-array placeholders. Flags are 0/1. State and matrices use the
same six-int64 and column-major conventions as sampling. Status 1 rejects invalid
extent/flags and 4 rejects state, preserving result/state. Status 0 means native
execution occurred; it does not guarantee a finite probability or convergence.
Native-consumed state is exported even when the result is NaN. Python raises
RuntimeError for nonzero C call status; its own validation raises TypeError or
ValueError before C. The public state setter validates snapshots, so corrupt
state normally arises only through raw C or private-buffer mutation. See the
[normal shape contracts](../../docs/multivariate_normal_shapes.md) and
[Student t shape contracts](../../docs/multivariate_t_shapes.md).


The finite-df Student t backend handles zero conditional variances as
point masses. Regression cases cover perfectly correlated and anticorrelated
variables, finite interval intersections, mixed deterministic coordinates and
zero-rank matrices, with inclusive deterministic endpoints.

Validation commands for this extension (CMD):

```bat
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -p test_stateful_functions.py -v
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
python -m unittest discover -s bindings/python/tests -p test_multivariate_cdf.py -v
python -m unittest discover -s bindings/python/tests -v
fpm test
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json --check
```

The binding tests build checked and optimized libraries and compile native,
C and C++ clients. The full Fortran suite is required here because the singular
Student t correction changes shared integration code.

## Real random matrix sampling

`api.ortho_group` and `api.special_ortho_group` expose `rvs(state, *, n=1)`
and `rvs_array(state, *, n=1, n_samples=1)`. Results are `ArrayResult` objects:
`outputs["sample"]` has shape `(n,n)` and `outputs["samples"]` has shape
`(n,n,n_samples)`. The batch axis is last, matching native Fortran, rather than
SciPy's leading batch axis. Dimensions are explicit keyword arguments and are
never inferred from output buffers. All nonnegative dimensions are permitted;
zero-dimensional matrices and empty batches succeed without consuming randomness.
Special-orthogonal samples have determinant +1; at n=1 they are always [[1]].
Orthogonal samples can have either determinant sign.

`api.random_correlation.rvs(state, eigs, tol=1e-13, diag_tol=1e-7)` infers
`n` from the one-dimensional eigenvalue array and returns `outputs["sample"]`
with shape `(n,n)`. It requires n >= 2, finite eigenvalues >= -tol, and
`abs(sum(eigs)-n) <= tol`. Both tolerances must be finite and nonnegative.
`tol` controls eigenvalue validation; `diag_tol` bounds the unit-diagonal error
before the native implementation sets the diagonal exactly to one. Eigenvalues
within tolerance can be slightly negative; they are not clipped. Numerical
accuracy still depends on conditioning, dimension, and floating-point rounding.
There is currently no correlation batch adapter.

```python
stream = api.rng(123)
q = api.ortho_group.rvs(stream, n=3)
assert q.status == 0
rotations = api.special_ortho_group.rvs_array(stream, n=3, n_samples=8)
corr = api.random_correlation.rvs(stream, [0.2, 0.8, 2.0])
assert corr.status == 0
```

Pass `out={"sample": buffer}` or `out={"samples": buffer}` to reuse a writable
float64 Fortran-contiguous array of the exact required shape. Inputs can be
copied into suitable storage. Outputs must not overlap any input, state, or other
output. Python rejects malformed dimensions, shapes, types and buffers before
calling C. Use an explicit `api.rng(seed)` stream belonging to the same API.
Repeated single draws and a batch with the same initial stream produce identical
matrices and final state. SciFort streams and algorithms do not reproduce SciPy's
bit sequences.

The generated C headers are included by `scifort.h`. Their symbols are
`scifort_ortho_group_rvs_f64`, `scifort_ortho_group_rvs_array_f64`,
`scifort_special_ortho_group_rvs_f64`,
`scifort_special_ortho_group_rvs_array_f64`, and
`scifort_random_correlation_rvs_f64`. Dimensions precede native arguments;
correlation retains the native argument order `(state,eigs,sample,status,tol,diag_tol)`.
C requires both tolerance scalars. Matrices use column-major layout, batches are
consecutive matrix slices, and state is a writable `int64_t[6]` buffer. All
buffers must be valid and nonoverlapping, including for zero extents; these thin
explicit-shape adapters do not check pointer validity or allocation capacity.

Status contracts:

| Status | Orthogonal groups | Random correlation |
| ---: | --- | --- |
| 0 | Success | Success |
| 1 | Invalid adapter dimension | Invalid adapter dimension (n < 2) |
| 2 | Native numerical failure | Invalid eigenvalues or tolerances |
| 3 | Unused | Native numerical failure |
| 4 | Invalid serialized RNG state | Invalid serialized RNG state |

Adapter validation errors preserve output and state, including invalid state on
an empty call. Correlation parameter errors zero the output without consuming
state. Numerical failures export all consumed RNG state. An orthogonal single
draw failure leaves a zero matrix; batch failure retains completed slices and
leaves later slices unspecified. Correlation numerical failure can leave a
partially transformed matrix. Always check status before using samples.

See the reviewed [orthogonal](../../docs/ortho_group_shapes.md),
[special-orthogonal](../../docs/special_ortho_group_shapes.md), and
[correlation](../../docs/random_correlation_shapes.md) native shape contracts.
Focused checks from the repository root (CMD):

```bat
python -m unittest discover -s bindings/python/tests -p test_random_matrix_sampling.py -v
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json --check
```

The focused test compiles the native/C agreement client
`test/test_random_matrix_c_api.f90`, runs a C11 client, checks C++17 headers,
and exercises Python with checked (`-O0 -g -fcheck=all -fbacktrace`) and `-O3`
shared libraries.

Local validation on 2026-11-27 used Python 3.13.3 and GNU Fortran/GCC/G++
17.0.0 20260510 (experimental) on Windows. All eight focused tests passed in
282.035 seconds, including both build configurations and the native/C/C++ clients.
The three native numerical implementations are unchanged apart from comments;
the full `fpm test` suite was not rerun for this bounded interface addition.

## Fixed-margin random tables

`api.random_table` provides `pmf(x, row, col)`, `logpmf(x, row, col)`,
`mean(row, col)`, `rvs(state, row, col)`, and
`rvs_array(state, row, col, *, n_samples=1)`. Probability methods return scalars
and `mean` returns a float64 matrix. Samplers return `ArrayResult` with an integer
`outputs["sample"]` matrix or `outputs["samples"]` cube.

Required shapes are row `(n_rows)`, col `(n_cols)`, x/sample/mean
`(n_rows,n_cols)`, and samples `(n_rows,n_cols,n_samples)`. Batch slices are last,
matching Fortran layout. Dimensions come from row and col; only n_samples is an
explicit output-size keyword. Probability, mean and single sampling require
nonempty margins. Empty batches also allow zero-length margins; use explicitly
typed integer arrays such as `np.array([], dtype=np.intc)`, since Python lists
with no elements do not identify an integer dtype.

```python
row, col = [3, 5], [2, 4, 2]
probability = api.random_table.pmf([[1, 1, 1], [1, 3, 1]], row, col)
mean = api.random_table.mean(row, col)
stream = api.rng(123)
result = api.random_table.rvs_array(stream, row, col, n_samples=10)
assert result.status == 0
samples = result.outputs["samples"]  # shape (2,3,10), dtype np.intc
```

Inputs must contain integers, with no floating, complex, or boolean values.
Python checks C-int range before conversion. Output buffers must be writable,
aligned, Fortran-contiguous `np.intc` arrays of the exact shape, passed through
`out={"sample": buffer}` or `out={"samples": buffer}`. The mean accepts a float64
output buffer through `out=buffer`. Inputs, output and state must not overlap.
The [native shape contracts](../../docs/random_table_shapes.md) also cover the
private marginal validator; its bounded-sum helper accepts any vector length.

Nonempty sampling requires nonnegative margins with equal totals. Native code
now checks sums before addition, requiring totals <= `huge(0)-1` so factorial
arguments `total+1` and `count+1` cannot overflow. This is 2147483646 for the
reviewed C-int build. The same domain applies to probability and mean calls.
Invalid margins produce NaN from pmf/logpmf/mean; negative table entries produce
NaN. Nonnegative tables whose row or column sums do not match the valid margins
have pmf 0 and logpmf -infinity, including tables whose sums exceed integer range.
Zero margins define an all-zero table with probability 1; its mean retains the
native NaN result for zero total. Avoid assuming uniformly high relative accuracy
for probabilities obtained by subtracting large log-factorials at large totals.

Status 1 denotes adapter dimension errors and preserves output/state. Invalid
serialized state returns 4, also preserving buffers, including on empty calls.
Native marginal errors return 2 with zeroed sample output and unchanged state;
status 3 denotes conditional-sampler failure. A numerical failure retains consumed
state and completed tables; the current table can be partial, and later tables
are zeroed. Empty batches skip native marginal validation but still validate
shapes, C-int representability and RNG state. Repeated single draws and split
batches match whole batches exactly for a shared initial stream. SciFort does
not reproduce SciPy's random bit sequences or its Boyett/Patefield method options.

The generated `scifort_generated_random_table.h`, included by `scifort.h`, exports
`scifort_random_table_{pmf,logpmf,mean,rvs,rvs_array}_f64`. Dimensions n_rows and
n_cols precede the native arguments; batch sampling additionally takes n_samples
before state. C uses contiguous column-major `int` table/count buffers, `double`
probability/mean outputs and caller-owned `int64_t[6]` state. All pointers must
be valid and non-NULL, even for empty arrays; sizes and nonoverlap are caller
responsibilities. Scalar native results retain NaN/infinity with status 0.
Integer-domain handling agrees with the native implementation.

Focused validation commands from the repository root (CMD):

```bat
python -m unittest discover -s bindings/python/tests -p test_random_table.py -v
build\cabi-generator-venv\Scripts\python.exe -m unittest discover -s tools/cabi_generator/tests -v
build\cabi-generator-venv\Scripts\python.exe tools\cabi_generator\pipeline.py bindings\cabi_config\pipeline.json --check
```

The focused test builds checked and optimized libraries, compares native/C/Python
values and RNG state, exercises a C11 client and C++17 headers, and compiles the
native random_table module with `-ftrapv` in its checked overflow client.

Local validation on 2026-11-27 used Python 3.13.3, SciPy 1.15.3 and GNU
Fortran/GCC/G++ 17.0.0 20260510 (experimental) on Windows. All 93 generator
tests passed in 599.194 seconds; all seven focused table tests passed in
275.567 seconds across checked and optimized builds. The pipeline checked all
140 artifacts. The numerical change is confined to random_table validation;
focused native/C/Python and overflow-trapping checks cover that behavior.
The full native `fpm test` suite was not rerun.
