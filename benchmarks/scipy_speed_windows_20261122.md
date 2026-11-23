# SciFort versus SciPy: local Windows results, 2026-11-22

This is a readable snapshot of the saved benchmark run, not a new measurement.
All nine workloads passed accuracy and timed-checksum checks. Timings apply to
one full operation at the listed size, rather than one element. See the
[benchmark guide](README_scipy_speed.md) for inputs, APIs, timing boundaries,
and instructions for running the comparison.

The normal quantile row below uses the original bisection kernel. See the
[direct inverse-normal update](../validation/normal_quantile_2026-11-22.md)
for its subsequent accuracy checks and much faster local timings.

## Median elapsed times

Each median uses seven separately timed batches after warmup and calibration.
`SciPy/SciFort` is the ratio of median times: above one favors SciFort, below
one favors SciPy. Ratios are computed from unrounded measurements.

| Workload | Size | SciFort (ms) | SciPy (ms) | SciPy/SciFort |
| --- | ---: | ---: | ---: | ---: |
| Normal CDF array | 65,536 | 1.643592 | 2.627611 | 1.599 |
| Normal CDF scalar loop | 256 | 0.006380 | 0.249804 | 39.156 |
| Normal quantile array | 16,384 | 22.975663 | 0.403958 | 0.018 |
| FFT, power of two, plan reused | 4,096 | 0.079891 | 0.037425 | 0.468 |
| FFT, prime length, plan reused | 4,093 | 0.408775 | 0.182258 | 0.446 |
| FFT, prime length, plan setup included | 4,093 | 0.828479 | 0.164973 | 0.199 |
| LU factorization and solve | 128 x 128 | 0.303887 | 0.119786 | 0.394 |
| LU solve, factors reused | 128 x 128 | 0.025809 | 0.010034 | 0.389 |
| Barycentric interpolation, existing object | 4,096 | 1.154867 | 1.245400 | 1.078 |

The scalar CDF row includes Python loop/call overhead; use the array row for
array-throughput comparisons. FFT reuse excludes native plan construction,
while the setup row includes it; SciPy uses its normal public FFT API and
internal caching. LU results include each library's ordinary checks and use
SciPy's installed BLAS/LAPACK versus SciFort's native kernels. Interpolator
construction is excluded. These are workload comparisons, not a general
Fortran-versus-Python ranking. Small timing differences may be noise.

## Environment

| Item | Recorded value |
| --- | --- |
| Operating system | `Windows-11-10.0.26300-SP0` |
| Architecture | `AMD64` |
| Processor identifier | `Intel64 Family 6 Model 140 Stepping 1, GenuineIntel` |
| Logical CPU count | 8 |
| SciFort compiler | GNU Fortran (GCC) 17.0.0 20260510 (experimental) |
| Compiler flags | `-O3 -std=f2018 -fimplicit-none -Werror=implicit-interface` |
| Python | 3.13.3 |
| NumPy | 2.2.6 |
| SciPy | 1.15.3 |
| NumPy BLAS/LAPACK build | OpenBLAS 0.3.29, 64-bit integer interface |
| SciPy BLAS/LAPACK build | OpenBLAS 0.3.28, 32-bit integer interface |
| Thread settings | BLAS/OpenMP environment variables set to 1; FFT workers=1 |
| Target batch duration | At least 0.1 seconds during calibration |
| Run recorded (UTC) | `2026-11-23T01:25:06.416800+00:00` |

The UTC timestamp corresponds to the evening of November 22 in America/New_York.
Only the required SciFort modules were compiled, without FPM, fast-math, or
runtime checks. Compilation, process startup, imports, setup, validation, and
file I/O were excluded from timings. Native timing used `system_clock`; Python
used `perf_counter`. Common thread environment variables were set before
loading NumPy/SciPy; the report does not independently measure active threads.

## Accuracy checks

Every output element in both alternating input sets was checked before timing.
Absolute differences below describe this run; acceptance also uses the
workload-specific relative, absolute, or transform-scaled criteria in the
harness. The LU scaled residual was below 8.93e-16.

| Workload | Values compared | Maximum absolute difference |
| --- | ---: | ---: |
| Normal CDF array | 131,072 | 2.220446e-16 |
| Normal CDF scalar loop | 512 | 2.220446e-16 |
| Normal quantile array | 32,768 | 1.332268e-15 |
| FFT, power of two, plan reused | 8,192 | 7.635863e-11 |
| FFT, prime length, plan reused | 8,186 | 9.026592e-11 |
| FFT, prime length, plan setup included | 8,186 | 9.026592e-11 |
| LU factorization and solve | 256 | 2.220446e-15 |
| LU solve, factors reused | 256 | 2.220446e-15 |
| Barycentric interpolation, existing object | 8,192 | 1.998401e-15 |

Total: **197,620 values**.

## Timing samples

Each list contains the seven elapsed times per operation in milliseconds, in
measurement order. Batch counts are independent for the two backends. This
preserves the observed variation rather than just the median.

| Case | Backend | Calls per batch | Seven samples (ms per operation) |
| --- | --- | ---: | --- |
| `normal_cdf_array` | SciFort | 64 | 2.11279844, 1.68211563, 1.64359219, 1.60885, 1.56953125, 1.655375, 1.599425 |
| `normal_cdf_array` | SciPy | 64 | 3.04239531, 2.85905156, 2.61523438, 2.62761094, 2.56349219, 2.74220781, 2.53346406 |
| `normal_cdf_scalar` | SciFort | 16,384 | 0.00775369263, 0.00637967529, 0.00729484863, 0.00595914917, 0.00646707764, 0.00613326416, 0.00626757813 |
| `normal_cdf_scalar` | SciPy | 512 | 0.258397266, 0.260909766, 0.236207812, 0.239106641, 0.249804492, 0.241524805, 0.268476758 |
| `normal_ppf_array` | SciFort | 8 | 22.9756625, 22.537675, 24.9653, 22.131, 22.857, 24.3732375, 23.3252125 |
| `normal_ppf_array` | SciPy | 256 | 0.398973047, 0.413477344, 0.398909375, 0.407276953, 0.400016797, 0.412995703, 0.403957813 |
| `fft_power2_reuse` | SciFort | 2,048 | 0.0911542969, 0.0791108398, 0.079890625, 0.0771296387, 0.0754931152, 0.083627002, 0.0869914551 |
| `fft_power2_reuse` | SciPy | 4,096 | 0.0403893311, 0.0382054443, 0.0369816162, 0.0347127686, 0.0348431152, 0.0388219971, 0.0374247803 |
| `fft_prime_reuse` | SciFort | 256 | 0.415044922, 0.404567187, 0.397529297, 0.553105469, 0.583707422, 0.387232422, 0.408775391 |
| `fft_prime_reuse` | SciPy | 1,024 | 0.163931445, 0.163005762, 0.189086035, 0.283860742, 0.221933301, 0.182258105, 0.168311719 |
| `fft_prime_setup` | SciFort | 128 | 0.828478906, 0.792402344, 0.793365625, 0.80099375, 0.864180469, 0.861535938, 0.862236719 |
| `fft_prime_setup` | SciPy | 1,024 | 0.164973438, 0.159812891, 0.15556748, 0.157937598, 0.185149805, 0.165732031, 0.16561748 |
| `lu_factor_solve` | SciFort | 512 | 0.30815957, 0.302426562, 0.330520898, 0.29558125, 0.303886523, 0.292477539, 0.305625 |
| `lu_factor_solve` | SciPy | 1,024 | 0.123763379, 0.119786328, 0.12919043, 0.11864248, 0.11965957, 0.12401748, 0.111545898 |
| `lu_reuse` | SciFort | 8,192 | 0.0252354248, 0.024370105, 0.0244601807, 0.0258088745, 0.0268555908, 0.0277835938, 0.0265393188 |
| `lu_reuse` | SciPy | 16,384 | 0.00932518311, 0.00916712647, 0.0109070374, 0.0101392944, 0.0100338928, 0.0102394348, 0.0096585022 |
| `barycentric` | SciFort | 128 | 1.07132656, 1.17359141, 1.12684219, 1.15486719, 1.11076172, 1.34570078, 1.41069922 |
| `barycentric` | SciPy | 128 | 1.25539297, 1.244025, 1.191075, 1.2454, 1.24468359, 1.38444453, 1.39628203 |

## Reproduction and source record

The run used this command from the repository root:

```text
python benchmarks/compare_scipy.py --output build/scipy-speed/windows-20261122.json
```

The original JSON is stored locally at
`build/scipy-speed/windows-20261122.json`; the build directory is ignored by Git.
This Markdown snapshot preserves the timing samples and key environment and
accuracy information in the repository. It does not refresh automatically.

- Recorded Git HEAD: `28504b653cad9d3cfd52499f87d6d09bcdef432f`.
  This preceded committing the new benchmark files; the hashes below identify
  the measured harness and compiled sources.
- Python harness SHA-256: `2497928ea7c1993ed5cb5ecfee5f03e444901532c7152a5dacfba2dc06cc4da2`.
- Native build signature: `41748651e829dd2138b1bb4037f87d7a1eaecac9cae8853f9fc59812cce80807`.
  The build signature includes the compiler path/version, flags, and native source contents.
- Saved JSON SHA-256: `9b48176405a1f33da5afdf75032c93d16c0c1457a12b40f753a0418d5e35d3b9`.

Rerun on the target machine before making performance decisions. Compiler
versions, library backends, workload sizes, and background CPU activity can
change both timings and rankings.
