# SciFort Python interface

This development interface exposes PDF, log-PDF, CDF, SF, log-CDF, log-SF, PPF,
and ISF for `normal`, `exponential`, `logistic`, `laplace`, `cauchy`, `weibull`,
`lognormal`, `gamma`, `beta`, and Student `t`. It uses NumPy arrays and the stable
binary64 C calling convention.
Its numerical calculations run in Fortran. Python 3.10+ and NumPy 1.24+ are
required; SciPy is needed only for comparisons/tests. The generated package does
not depend on fparser at runtime.

## Try it from the repository

In CMD, with gfortran on PATH and NumPy installed:

```text
python bindings\python\example.py
```

The example builds a cached, optimized, ten-family development shared library,
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

## Python contract

- Scalar input returns a Python float. Zero-dimensional NumPy input is scalar.
- Array-like real numeric input returns a native-endian float64 ndarray with the
  same shape. Lists, integer/float32 arrays, noncontiguous views, unaligned arrays
  and nonnative-endian arrays are converted when needed. Complex, string and
  object arrays are rejected. Empty arrays are supported and still validate
  uniform parameters.
- Location and scale default explicitly to 0 and 1 in the Python configuration.
  Weibull `c`, lognormal `s`, gamma `a`, beta `a, b`, and Student t `df` are required. All parameters must be
  scalar; broadcasting and per-element parameters are not implemented.
- Optional `out=` accepts a writable, aligned, C-contiguous, native float64 NumPy
  array of exactly the input shape. It is returned by identity. Wrong dtype,
  shape, layout, read-only storage and possible overlap with the input are
  rejected before the C call. Scalar input does not accept `out`.
- Scalar parameter errors retain native NaN results. Bulk call-level errors
  (invalid uniform parameters or rejected lengths) raise `ValueError`; output
  remains unchanged. Invalid individual points/probabilities retain native
  NaN/infinity/support behavior. An unknown ABI status raises `RuntimeError`.
- No numerical method is reimplemented in Python. Native behavior and intentional
  SciPy differences remain applicable. In particular SciPy 1.15.3 returns NaN
  for gamma/Weibull PDF/log-PDF at positive infinity, whereas SciFort returns the
  mathematical limits 0 and negative infinity.

## Generation and shared-library boundary

[Python policy](../cabi_config/python.json) references the existing C manifests.
It provides explicit defaults and mappings for the three normal methods already
maintained outside the C generator. The reusable
[generator/runtime](../../tools/cabi_generator/README.md) emits method signatures,
ctypes metadata, a runtime copy and small supplementary Fortran adapters.

```text
python tools\cabi_generator\generate_python.py bindings\cabi_config\python.json
python tools\cabi_generator\generate_python.py bindings\cabi_config\python.json --check
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
