# Gustafsson zero-phase forward/backward filtering

The `scifort_signal` module now implements the second initialization method
for `scipy.signal.filtfilt`: the default remains `method='pad'`, while
`method='gust'` determines forward and backward initial filter states so that
the order of the two passes agrees in a least-squares sense.

```fortran
use scifort_kinds, only : dp
use scifort_signal, only : filtfilt, filtfilt_gust
real(dp) :: b(3), a(3), x(100)
real(dp), allocatable :: y(:)
integer :: code
! Initialize coefficients and signal first.
y = filtfilt(b, a, x, method='gust', status=code)
y = filtfilt(b, a, x, method='gust', irlen=24, status=code)
y = filtfilt_gust(b, a, x, irlen=24, status=code)
```

`filtfilt` supports rank-1, rank-2, and rank-3 real or complex signals,
with real or complex coefficient arrays of the same type as the signal.
Use `axis=1`, `axis=2`, or `axis=3` as appropriate for rank-2/3 arrays;
axis indices are one-based, and the default is the last axis. For a rank-1
signal the axis argument is omitted. The standalone rank-1
`filtfilt_gust` also supports mixed real coefficients/complex samples and
complex coefficients/real samples.

- `method='pad'` (default): existing odd/even/constant/none extension and
  steady-state initialization are unchanged. `padtype` and `padlen` apply.
- `method='gust'`: no signal extension; `padtype` and `padlen` are ignored,
  including for signals shorter than the default padding length. The optional
  positive `irlen` limits the length of the impulse response used to estimate
  the start/end states. With `irlen` absent, the complete signal is used;
  truncation is used only when the signal length exceeds `2*irlen`.
- A length-1 signal and a zeroth-order transfer function are supported.
  Direct-form-II initial state is computed using a minimum-norm linear
  least-squares solution for underdetermined cases.
- Invalid coefficients, nonfinite input, unknown `method`, and nonpositive
  `irlen` return `filtfilt_status_invalid_input`, with an allocated quiet-NaN
  result. Failure of the Gustafsson numerical solve returns the new
  `filtfilt_status_numerical_failure` code (3).

**Assumed-rank API:** `filtfilt_nd` now also accepts the trailing optional
arguments `method='gust'` and `irlen=...`, for real/complex ranks 1-15 with
one-based `axis`. The output is caller-allocated and must match the rank
and shape of `x`:

```fortran
real(dp) :: x(4,2,100,3), y(4,2,100,3)
integer :: code
! Initialize input x and filter coefficients b,a first.
call filtfilt_nd(b, a, x, y, axis=3, method='gust', irlen=32, status=code)
```

The filtering axis may be shorter than the default padding length. Every
combination of real/complex numerator, denominator and input promotes to
complex output when necessary. `padtype` and `padlen` are ignored for
`method='gust'`, while `method='pad'` retains the original behavior.
Invalid methods and `irlen <= 0` return `filtfilt_status_invalid_input`;
invalid filters, nonfinite data or numerical failure produce a NaN-filled
output, including when a later channel fails. Empty outer dimensions are
allowed when filter parameters and the selected axis are valid. The rank-independent `sosfiltfilt_nd` also supports Gustafsson states,
implemented directly in the SOS cascade.

Gustafsson conditioning can be
poor for short high-order transfer-function filters; for high-order filtering,
second-order sections are usually preferable when the padding method suffices.

## Algorithm and numerical comparisons

The implementation constructs the observability/state-response matrix,
filters its time reversal, and solves the forward-backward consistency
condition in the least-squares sense (Gustafsson, *IEEE Transactions on
Signal Processing*, 46(4), 988-992, 1996). It uses an internal complex
Jacobi SVD for small systems and a pivoted, reorthogonalized QR method with
minimum-norm correction as a fallback. It does not use LAPACK or copy SciPy
source code. Persistent state is not required. Working memory scales as
`O(length*filter_order)` when `irlen` is absent; specifying a bounded
`irlen` limits the state-estimation matrices while retaining the output-sized
arrays.

Run `test/test_signal_filtfilt_gust.f90` with the standalone signal runner.
Build `tools/signal_reference/driver_signal_filtfilt_gust.f90` against the
public filter module, then run:

```text
python tools/signal_reference/compare_signal_filtfilt_gust_scipy.py <driver-executable>
```

The Python script includes direct comparisons with SciPy 1.17 for real
filters and signals (including complex-valued observations), and an
independent complex-valued least-squares reference for complex filter
coefficients. SciPy 1.17 uses a real observability matrix even for complex
coefficients and emits a `ComplexWarning`, so complex-coefficient cases are
not compared directly with that implementation.

## Validation record (2026-10-26)

Using GNU Fortran 14.2 and SciPy 1.17, the final comparison driver passed
**858 cases** and **41,897 sample values**: 758 real-coefficient cases
compared directly with SciPy, plus 100 complex-coefficient cases compared
with a separate complex least-squares implementation. The maximum absolute
difference was `2.26e-8` for a short, ill-conditioned, fourth-order filter;
the maximum absolute error divided by `1 + max(abs(reference))` was
`7.89e-11`. The largest complex-coefficient reference discrepancy was
`2.29e-15`.

The clean GNU Fortran checked build passed 5/5 relevant filtering regression
programs. Those five programs and the 858 reference comparisons also passed
with the modified filtering modules compiled at `-O3`. The full optimized
dependency rebuild was interrupted by the execution limit, and `fpm` was not
available.

## Rank-independent comparison and test record (2026-10-26)

`test/test_signal_filtfilt_nd_gust.f90` checks ranks 1, 5, and 15,
selected axes, all eight real/complex argument combinations, signals shorter
than the pad length, noncontiguous sections, empty channels and invalid
inputs. The previous `test/test_signal_filtfilt_nd.f90` padding test still
passes unchanged.

The driver `tools/signal_reference/driver_signal_filtfilt_nd_gust.f90` and
Python reference checker
`tools/signal_reference/compare_signal_filtfilt_nd_gust_scipy.py` compare
960 rank-1/2/5/15 cases (48,384 samples): 240 against `scipy.signal.filtfilt`
from SciPy 1.17 and 720 complex-coefficient cases against a separate complex
least-squares implementation. Both checked and `-O3` GNU Fortran 14.2
builds pass with maximum absolute deviations <= `2.24e-16`.

One way to reproduce (run from project root after building the minimal
`kinds`, `signal_lfilter`, `signal_sos`, `signal_filtfilt_gust`,
`signal_filtfilt`, and `signal_filtfilt_nd` modules) is:

```text
python tools/signal_reference/compare_signal_filtfilt_nd_gust_scipy.py DRIVER_EXECUTABLE
```

## SOS-native Gustafsson extension (2026-10-26)

SciFort now provides `sosfiltfilt(sos,x,method='gust',irlen=...)` for
rank-1/2/3 arrays and `sosfiltfilt_nd(sos,x,y,method='gust',irlen=...)`
for caller-supplied real/complex arrays through rank 15. A separate
`sosfiltfilt_gust(sos,x,irlen,status)` generic accepts all four rank-1
real/complex coefficient/sample combinations. Every `axis=` is one-based.
Existing calls retain their default `method='pad'` behavior.

```fortran
use scifort_signal, only : sosfiltfilt_nd, filtfilt_status_success
real(dp) :: sos(6,6), records(2,3,200), filtered(2,3,200)
integer :: code
! Initialize sos and records; each SOS row is [b0,b1,b2,1,a1,a2].
call sosfiltfilt_nd(sos,records,filtered,axis=3,method='gust', &
    irlen=48,status=code)
if (code /= filtfilt_status_success) error stop 'SOS filtering failed'
```

This SOS-native solver computes all state-response columns using the
second-order-section cascade, then uses the existing SciFort complex
least-squares solvers to choose initial states. It **never expands** the
cascade to the often ill-conditioned high-order numerator and denominator
polynomials. `irlen` bounds the portion of the signal used to estimate
initial states when the signal is longer than `2*irlen`; it does not remove
the need to calculate the signal-sized forward and backward outputs.
`padtype` and `padlen` are ignored for `method='gust'`, as with BA filtering.

**This is a SciFort extension, not a SciPy API option:** SciPy 1.17
`scipy.signal.sosfiltfilt` supports padding but has no `method='gust'`
parameter. The independent Python reference therefore builds the SOS
state-observability matrix with `scipy.signal.sosfilt` and solves for the
Gustafsson states using `scipy.linalg.lstsq`, without calling a nonexistent
SciPy API. High-order, short, untruncated systems may be numerically
underdetermined, even with SOS state coordinates. If the least-squares
consistency residual is too large, SciFort returns
`filtfilt_status_numerical_failure` (3), rather than inaccurate values;
outputs are quiet NaNs. Consider longer observations, bounded `irlen`, or
the default SOS padding method in such cases.

The new `test/test_signal_sos_gust.f90` checks ranks 1/2/3/5/15, mixed
real/complex coefficient/sample types, noncontiguous array sections,
method and input rejection, short signals, backward compatibility with
padding, and numerical-failure reporting. Numerical comparisons use
`tools/signal_reference/driver_signal_sos_gust.f90` and
`compare_signal_sos_gust_scipy.py`. The state-space construction and
Fortran test are original SciFort code; no new runtime dependency was
introduced.

## Under-determined high-order conditioning improvement (2026-10-26)

When the state matrix has at least as many unknowns as observations,
SciFort first applies column-pivoted *Householder QR* to its conjugate
transpose. For a full-row-rank design, this computes a minimum-norm
initial-state solution without squaring the condition number. Jacobi SVD
and the earlier reorthogonalized pivoted QR remain fallbacks when the
Householder system is rank-deficient or cannot be solved reliably.

The new solver resolves four of the eight previously rejected short,
high-order test cases from the 297-case SOS state-space reference suite.
Four 12th-order, 20-sample untruncated cases still return numerical
failure: their very small singular values are close to floating-point
resolution, and accepting the resulting samples would be misleading.
The SOS numerical-failure consistency check remains enabled.
