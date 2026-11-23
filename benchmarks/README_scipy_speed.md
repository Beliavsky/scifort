# SciFort versus SciPy workload timings

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

| Case | Size | Native SciFort call | SciPy call and timed work |
| --- | ---: | --- | --- |
| `normal_cdf_array` | 65,536 | elemental `normal_cdf` | `special.ndtr` with reused output |
| `normal_cdf_scalar` | 256 | Fortran loop calling `normal_cdf` | Python loop calling `special.ndtr` |
| `normal_ppf_array` | 16,384 | elemental `normal_ppf` | `special.ndtri` with reused output |
| `fft_power2_reuse` | 4,096 | `fft_plan%execute` with an existing plan | `fft.fft`, one worker |
| `fft_prime_reuse` | 4,093 | `fft_plan%execute` with an existing plan | `fft.fft`, one worker |
| `fft_prime_setup` | 4,093 | initialize a plan and execute it each call | `fft.fft`, one worker |
| `lu_factor_solve` | 128 by 128 | `lu_factor` then `lu_solve` | `linalg.lu_factor` then `linalg.lu_solve` |
| `lu_reuse` | 128 by 128 | `lu_solve` with existing factors | `linalg.lu_solve` with existing factors |
| `barycentric` | 4,096 queries | existing interpolator's `values` | existing `BarycentricInterpolator` evaluated at an array |

Normal CDF inputs cover [-8,8); quantile probabilities cover [1e-6,1-1e-6).
Interpolation uses 33 Chebyshev nodes and samples of `1/(1+25*x*x)`, with
queries inside [-0.99,0.99]. The comparison checks the same interpolating
polynomial, not its error relative to the sampled function. Linear systems
are well-conditioned dense matrices with two known solutions; the script also
checks their scaled residuals. Source tests supplying related correctness
coverage include `test_linalg_dense.f90`, `test_signal_fft_plans.f90`, and
`test_interpolate_barycentric.f90`.

## Timing and interpretation

Imports, process startup, input generation, validation, and file I/O are outside
the timers. Each native process stays alive while Python sends batch commands;
communication is outside the native timer. Timers measure elapsed wall time:
Fortran `system_clock` with int64 counts and Python `perf_counter`.

Both backends warm up. Each independently doubles its batch size until the batch
takes at least 0.1 seconds, then collects seven samples. Backend order alternates
between samples. Each timed call contributes to a checked checksum; that small
bookkeeping cost and Python dispatch/loop overhead are included. Reported times
are milliseconds per operation (one full array, transform, solve, or query batch),
not per element. Ratios above one in `SciPy/SciFort` mean SciFort was faster.

The scalar case deliberately measures language-loop/call overhead and must not
be used to describe array throughput. Array CDF/PPF comparisons use SciPy's
direct special-function kernels rather than the higher-level `stats.norm` API.
Their output buffers are reused on both sides.

The FFT native plan operates in place, so copying the preserved input to its
work buffer is included. SciPy keeps its input and allocates an output. The
reuse cases exclude native plan construction; the setup case includes it.
SciPy's public `fft` API does not expose an equivalent reusable plan, and its
internal caching remains enabled. Compare these rows as application workloads,
not as identical allocation or planning strategies.

LU input matrices and right-hand sides are preserved. Both backends perform
their normal finite-input checks. SciPy uses its installed BLAS/LAPACK backend;
SciFort uses its native kernels. This compares the libraries available to the
application, not a language-only speed difference. Interpolator construction
is excluded; evaluation includes each library's ordinary allocation costs.

The script sets common BLAS/OpenMP thread environment variables to one before
importing NumPy/SciPy, and explicitly uses one FFT worker. Avoid running other
CPU-intensive jobs at the same time. Results depend on machine, compiler,
backends, input sizes, allocation, and numerical algorithms. Do not combine
these ratios into one overall library score or use timing thresholds in CI.
These are warm repeated-use measurements, not cold-start latency measurements.

The JSON report defaults to `build/scipy-speed/results.json`. It contains raw
per-operation samples, medians, batch sizes, accuracy results, clock information,
compiler/version/flags, exact build commands, source hash, Python/NumPy/SciPy
versions, platform, thread settings, and NumPy/SciPy build configurations.
Generated binaries and reports stay in the ignored build directory.

## Focused checks and options

```text
python benchmarks/compare_scipy.py --check-only
python benchmarks/compare_scipy.py --quick --repeats 3 --target-seconds 0.02
python benchmarks/compare_scipy.py --case normal_cdf_array --case lu_reuse
python benchmarks/compare_scipy.py --output build/scipy-speed/my-machine.json
```

`--check-only` compiles with `-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow -fbacktrace`, checks all ordinary workload sizes,
and skips timing. `--quick` uses smaller workloads for smoke checks; those
timings are not representative of the default sizes. Change the compiler with
`--compiler PATH_TO_GFORTRAN`. Each new compiler/flag/source combination builds
its own executable. No compiler or Python packages are installed automatically.

## Reference API documentation

- [SciPy normal CDF kernel](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ndtr.html)
- [SciPy normal quantile kernel](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ndtri.html)
- [SciPy FFT](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.fft.html)
- [SciPy LU solves](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lu_solve.html)
- [SciPy barycentric interpolation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.BarycentricInterpolator.html)
