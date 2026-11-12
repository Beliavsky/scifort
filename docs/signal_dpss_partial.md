# Long-window DPSS eigenpairs

`dpss(m, nw, kmax, sym, norm, solver, ratio_method)` supports three
solver choices:

- `solver='auto'` (default): use the partial eigensolver for the internal
  tridiagonal size `n >= 256` and `kmax <= 32`; otherwise retain the original
  full-eigenvector QL method.
- `solver='partial'`: compute only the requested largest `kmax` eigenpairs
  using Sturm-count bisection and pivoted tridiagonal inverse iteration.
  Reorthogonalize the computed eigenvectors against preceding tapers.
- `solver='ql'`: force the existing implicit-shift tridiagonal QL algorithm.

The option is appended to the original argument list, so older positional and
keyword calls are compatible. It is a **SciFort extension**; SciPy does not
accept a `solver` keyword for its DPSS window function. The standard output
is unchanged: `out%tapers` stores samples by taper and `out%ratios` stores
band-concentration ratios. Normalizations `l2`, `approximate`, and
`subsample`, as well as symmetric and periodic windows, all work with both
solvers.

```fortran
use scifort_signal, only: dp, dpss, dpss_result

type(dpss_result) :: out
out = dpss(8192, 3.5_dp, kmax=5, solver='partial', ratio_method='fft')
if (out%status /= 0) error stop 'DPSS calculation failed'
```

## Numerical cost and limitations

The old QL path holds an `n`-by-`n` eigenvector matrix, while the partial
method uses an `n`-by-`kmax` array plus O(n) work vectors. Accordingly,
eigenvector workspace changes from **O(n^2)** to **O(n*kmax)**. Each selected
eigenvalue takes approximately 60 bisection iterations, and the vector is
refined using a pivoted tridiagonal factorization with residual-based
convergence checks. It returns an explicit `dpss_status_no_convergence` if
an eigenvector cannot be validated, rather than silently returning a
potentially unreliable result.

**Concentration-ratio methods:** the appended `ratio_method` keyword is a
SciFort extension, accepting `'auto'` (default), `'direct'`, or `'fft'`.
`direct` uses the original O(kmax*n^2) autocorrelation summation. `fft`
computes the *same* concentration Rayleigh quotient with zero-padded FFT
linear autocorrelations (O(kmax*n*log(n)) arithmetic and O(n) complex
scratch). The automatic method selects FFT at internal lengths `n >= 256`;
shorter lengths use the direct summation. Periodic DPSS concentration ratios
are always evaluated on the full length `m+1` eigenvectors before cropping.

To avoid a circular module dependency (`signal_spectral` itself imports
`signal_dpss`), DPSS uses a small private radix-2 transform, with no
external FFT library. Padded lengths are at least `2*n-1`; the autocorrelation
is linear, not circular. The direct method is retained as a verification
reference and can be forced when debugging extreme numerical cases.

The original `-O3` GNU Fortran 14.2 eigensolver benchmark of
`n=2048, kmax=6, nw=3.5` gave 0.0162 CPU seconds (partial) versus 6.54 CPU
seconds (QL), with 1.5 MB versus 34 MB peak process RSS. After adding FFT
ratios, a separate end-to-end six-taper test improved from 0.1969 s (direct)
to 0.0399 s (FFT) at n=8192, about 4.94x. These results are specific to the
compiler/hardware. Reproduce with
`benchmarks/bench_signal_dpss_partial.f90` and
`benchmarks/bench_signal_dpss_fft_ratios.f90`.

## Validation and provenance

`test/compare_signal_dpss_partial_scipy.py` verifies symmetric and periodic
windows, the three normalization modes, variable bandwidths, and lengths
through 16,384 with SciPy 1.17. The independently written eigensolver uses
classical symmetric tridiagonal Sturm counts and inverse iteration. No SciPy
or other third-party source is included or linked at runtime.

`test/test_signal_dpss_fft_ratios.f90` exercises forced FFT, direct and
automatic ratio selection for 150 cases, including small odd/even and
periodic windows. `test/compare_signal_dpss_fft_ratios_scipy.py` checks
256 independent SciPy DPSS concentration-ratio cases (both checked and
optimized GNU Fortran builds) with a maximum absolute difference of about
1.21e-14 in the optimized build. These tools use SciPy only for testing.
