# SciFort versus SciPy workload timings

See the [expanded 80-case Windows snapshot](scipy_speed_expanded_windows_20261122.md)
for current sample output. It was collected under concurrent CPU load; rerun
on an idle machine for performance conclusions.

The [historical nine-case Windows results from 2026-11-22](scipy_speed_windows_20261122.md)
preserve the original measurements, accuracy checks, environment details,
and individual samples.

That snapshot predates the direct inverse-normal kernel. See the
[normal quantile update](../validation/normal_quantile_2026-11-22.md) for its
validation and focused timing comparison; the historical snapshot is preserved.

Run from the repository root in CMD or a Unix shell:

```text
python -m pip install "numpy>=1.23.5" "scipy>=1.15"
python benchmarks/compare_scipy.py
```

Requires GNU Fortran on PATH. The script compiles only the required existing
SciFort modules and the benchmark driver, with `-O3` and no fast-math or runtime
checks. It caches the executable under `build/scipy-speed/` using a hash of the
compiler, flags, and source contents. It does not run FPM or the test suite and
does not modify the library. Compilation is excluded from timings.

## Workloads

All real calculations use binary64; FFTs use complex binary64. Two prepared
input sets alternate during each batch. Every output element in both sets is
checked against SciPy before timing that workload.

The suite has **80 workloads** (previously nine):

| Group | Cases | Default size | Compared operations |
| --- | ---: | --- | --- |
| Distribution arrays | 56 | 16,384 elements; normal CDF 65,536 | `pdf`, `logpdf`, `cdf`, `sf`, `logcdf`, `logsf`, `ppf`, `isf` for normal, exponential, logistic, Laplace, Cauchy, Weibull minimum, and lognormal |
| Scalar normal CDF | 1 | 256 calls | Fortran loop versus Python loop |
| Complex FFTs | 12 | 4,096, 4,050, or 4,093 points | Forward and normalized inverse transforms, each with reused plans or plan initialization included |
| Dense solves | 8 | 128 by 128 | LU factor-and-solve, reused LU, transposed LU, Cholesky factor-and-solve, reused Cholesky, lower/upper/transposed triangular solves |
| Barycentric interpolation | 3 | 4,096 queries, 33 nodes | Reused evaluation, construction plus evaluation, first derivative |

List exact case names and sizes without compiling:

```text
python benchmarks/compare_scipy.py --list
python benchmarks/compare_scipy.py --category linear_solve --list
```

Distribution array names have the form `FAMILY_METHOD_array`, for example
`weibull_logsf_array`. Standard location zero and scale one are used; Weibull
shape is 1.7 and lognormal shape is 0.8. CDF/density/tail arguments cover
[-8,8) for real-line families and [0.01,8.01) for positive-support families.
Quantile probabilities cover [1e-6,1-1e-6). These are throughput workloads on
ordinary finite inputs, not exhaustive numerical certification or extreme-tail
benchmarks. Existing distribution tests cover the broader numerical contract.

Normal CDF and PPF retain the direct `special.ndtr` and `special.ndtri` kernels
with reused output buffers. The other 54 distribution array cases call the
public `scipy.stats` methods, including their parameter validation, masking,
and output allocation. Native elemental calls reuse the output buffer. Thus
these rows compare application-facing APIs; they do not isolate the cost of
special-function arithmetic. Do not compare their ratios directly with the two
normal kernel ratios as if API overhead and allocation were identical.

Interpolation uses 33 Chebyshev nodes and samples of `1/(1+25*x*x)`, with
queries inside [-0.99,0.99]. Checks compare the same interpolating polynomial
and its derivative, not error relative to the sampled function. The setup row
includes constructing the interpolator each call on both sides; SciPy uses a
fixed RNG seed for weight construction. The derivative row includes ordinary
evaluation costs, with SciPy's derivative cache warmed before timing.

Linear systems use a symmetric positive-definite matrix
`A(i,j) = 1/(i+j) + n*delta(i,j)` with one-based indices. Triangular workloads
retain only the selected triangle. Both prepared right-hand sides have known
solutions `sin(i) + 0.1*k`, with `k=0,1`. Transposed solves use `A.T` for both
right-hand-side construction and residual checks. All solution elements,
known solutions, and scaled residuals are checked, rather than just checksums.

FFT lengths exercise power-of-two, mixed-radix, and prime-length paths.
Forward transforms are unnormalized; inverse transforms include division by
length, matching SciPy's `norm="backward"` convention.

## Timing and interpretation

Imports, process startup, input generation, validation, and file I/O are outside
the timers. Each native process stays alive while Python sends batch commands;
communication is outside the native timer. Timers measure elapsed wall time:
Fortran `system_clock` with int64 counts and Python `perf_counter`.

Both backends warm up. Each independently doubles its batch size until the batch
takes at least 0.1 seconds, then collects seven samples. Backend order alternates
between samples. Each timed call contributes to a checked checksum (including imaginary FFT
components); that small
bookkeeping cost and Python dispatch/loop overhead are included. Reported times
are milliseconds per operation (one full array, transform, solve, or query batch),
not per element. Ratios above one in `SciPy/SciFort` mean SciFort was faster.

The scalar case deliberately measures language-loop/call overhead and must not
be used to describe array throughput. Normal array CDF/PPF comparisons use SciPy's
direct special-function kernels rather than the higher-level `stats.norm` API.
Their output buffers are reused on both sides.

The FFT native plan operates in place, so copying the preserved input to its
work buffer is included. SciPy keeps its input and allocates an output. The
reuse cases exclude native plan construction; the setup case includes it.
SciPy's public `fft` API does not expose an equivalent reusable plan, and its
internal caching remains enabled. Compare these rows as application workloads,
not as identical allocation or planning strategies.

Linear-solve input matrices and right-hand sides are preserved. Both backends perform
their normal finite-input checks. SciPy uses its installed BLAS/LAPACK backend;
SciFort uses its native kernels. This compares the libraries available to the
application, not a language-only speed difference. Interpolator construction
is excluded except in the explicitly named setup case; evaluation includes
each library's ordinary allocation costs.

The script sets common BLAS/OpenMP thread environment variables to one before
importing NumPy/SciPy, and explicitly uses one FFT worker. Avoid running other
CPU-intensive jobs at the same time. Results depend on machine, compiler,
backends, input sizes, allocation, and numerical algorithms. Do not combine
these ratios into one overall library score or use timing thresholds in CI.
These are warm repeated-use measurements, not cold-start latency measurements.

Every successful run writes JSON and a sibling Markdown report, defaulting to
`build/scipy-speed/results.json` and `build/scipy-speed/results.md`. Markdown
contains a readable table and individual timing samples, or agreement checks
without invented timings in check-only mode. JSON contains raw
per-operation samples, medians, batch sizes, accuracy results, clock information,
compiler/version/flags, exact build commands, source hash, Python/NumPy/SciPy
versions, platform, thread settings, and NumPy/SciPy build configurations.
Generated binaries and reports stay in the ignored build directory.

## Focused checks and options

```text
python benchmarks/compare_scipy.py --check-only
python benchmarks/compare_scipy.py --quick --repeats 3 --target-seconds 0.02
python benchmarks/compare_scipy.py --case normal_cdf_array --case lu_reuse
python benchmarks/compare_scipy.py --category distribution
python benchmarks/compare_scipy.py --category fft --category interpolation
python benchmarks/compare_scipy.py --output build/scipy-speed/my-machine.json
```

Running all 80 workloads can take several minutes. Select a category or
specific cases for focused investigations.

`--case` and `--category` are repeatable; when both are supplied, their
intersection is selected. An empty intersection is an error. `--list` displays
the selection without invoking a compiler.

`--check-only` compiles with `-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow -fbacktrace`, checks all ordinary workload sizes,
and skips timing. `--quick` uses smaller workloads for smoke checks; those
timings are not representative of the default sizes. Change the compiler with
`--compiler PATH_TO_GFORTRAN`. Each new compiler/flag/source combination builds
its own executable. No compiler or Python packages are installed automatically.

The harness safeguards can be checked separately without a compiler:

```text
python benchmarks/test_compare_scipy.py
```

## Reference API documentation

- [SciPy normal CDF kernel](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ndtr.html)
- [SciPy normal quantile kernel](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ndtri.html)
- [SciPy FFT](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.fft.html)
- [SciPy LU solves](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lu_solve.html)
- [SciPy barycentric interpolation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.BarycentricInterpolator.html)

- [SciPy distribution methods and lognormal parameterization](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.lognorm.html)
- [SciPy Cholesky solves](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cho_solve.html)
