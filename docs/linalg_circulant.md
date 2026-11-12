# Circulant matrices and FFT-based solves (`scipy.linalg.circulant` / `solve_circulant`)

Module: `scifort_linalg_circulant` (MIT; independently authored for SciFort).

## Fortran API

`circulant(c, a, status)` constructs an `n x n` circulant matrix from its
first column `c(n)`. With one-based indices, `a(i,j) = c(modulo(i-j,n)+1)`.
A real `c` produces a real matrix; complex `c` produces a complex matrix.
The caller provides the output `a`. Empty input constructs an empty matrix.
Nonfinite input returns `linalg_status_no_convergence` and a zero matrix.

`solve_circulant(c, b, x, status, singular, tol, check_finite, baxis, outaxis)`
solves `circulant(c) x = b` using FFT diagonalization. `c` is always a
**one-dimensional vector**. `b` and `x` can be rank-one vectors or rank-two
arrays containing multiple right-hand sides. Output is caller-allocated and
set to zero if the solve fails. Real/complex and mixed-type inputs are
supported; a complex output is required if either input is complex.

- `singular` (optional): `'raise'` (default) reports a near-singular
  circulant matrix with `linalg_status_singular`, or `'lstsq'` returns the
  minimum-Euclidean-norm least-squares solution (spectral pseudoinverse).
  Other strings are rejected with `linalg_status_invalid_shape`.
- `tol` (optional): absolute eigenvalue cutoff, nonnegative and finite.
  The default is `n * epsilon(1.0_dp) * maxval(abs(fft(c)))`, matching SciPy.
  An eigenvalue is considered small if its absolute value is `<= tol`.
  `tol=0` treats exact zero eigenvalues as singular. A negative or nonfinite
  tolerance is rejected.
- `check_finite` (optional): `.true.` by default. NaN or infinity in `c`
  or `b` is rejected with `linalg_status_no_convergence`. With `.false.`,
  behavior for nonfinite input is not guaranteed.
- `baxis`, `outaxis` (optional, rank-two overloads only): axes holding the
  sample dimension in the input/output. These are **Fortran one-based**:
  `1` (default) for the first dimension, or `2` for the second. As in
  NumPy, `-1` selects the last axis and `-2` selects the first. For example,
  with `baxis=2` and `outaxis=1`, `b(nrhs,n)` becomes `x(n,nrhs)`.
  Output dimensions must agree with the selected axes.

The `status` argument uses the existing constants in `scifort_linalg`:
`linalg_status_success=0`, `linalg_status_invalid_shape=1`,
`linalg_status_no_convergence=4`, and `linalg_status_singular=5`.
Empty vectors and zero-column right-hand-side matrices are accepted when
shapes are consistent.

### Scope and SciPy compatibility

The vector/multiple-right-hand-side mathematical behavior, circulant
orientation, singular cutoff, and minimum-norm behavior agree with
SciPy 1.17.0, including real and complex inputs. This Fortran API uses
caller-owned output arrays and status codes in place of NumPy return arrays
and Python exceptions. Axes are one-based rather than zero-based.

SciPy's `caxis` and broadcasting among **batches of distinct circulant
coefficient vectors** are not implemented. SciFort accepts one `c(:)` per
call, and can solve it against several right-hand sides. Python's
arbitrary-rank broadcasting and array dtype promotion are not exposed in
the core Fortran interface. These restrictions are intentional and documented;
batch-level composition can be added in wrappers without changing the API.

## Example

```fortran
use scifort_kinds, only : dp
use scifort_linalg, only : linalg_status_success
use scifort_linalg_circulant, only : circulant, solve_circulant
real(dp) :: c(3), b(3), x(3), a(3,3)
integer :: status

c = [2.0_dp, 2.0_dp, 4.0_dp]
b = [1.0_dp, 2.0_dp, 3.0_dp]
call solve_circulant(c,b,x,status)
if (status /= linalg_status_success) error stop 'solve failed'
! x = [0.75, -0.25, 0.25]
```

See the [executable example](../example/linalg_circulant_example.f90).

## Algorithm and complexity

A circulant matrix has discrete Fourier eigenvectors. If `F` denotes the
unnormalized forward discrete Fourier transform, its eigenvalues are
`lambda = F(c)`. The system is solved via the elementwise spectrum
`F(b) / lambda`, followed by an inverse FFT. Small eigenvalues are either
rejected (`'raise'`) or given a zero pseudoinverse coefficient (`'lstsq'`).
This corresponds to the Moore-Penrose pseudoinverse because the Fourier
basis is unitary up to normalization and circulant matrices are normal.

For `n` rows and `nrhs` columns, the FFT path uses
`O(n*(1+nrhs)*log(n))` arithmetic and `O(n*nrhs)` total storage. The FFT
plan is constructed once per solve and reused for all right-hand sides.
SciFort's existing FFT plan supports radix-2, small direct transforms,
mixed radix, and Bluestein transforms for other lengths. No new FFT
implementation or external numerical dependency was added.

As with SciPy, numerically ill-conditioned (but not classified singular)
circulant systems can amplify floating-point errors. Explicit `tol` allows
callers to choose more conservative rank truncation. A finite status code
does not provide a condition-number certificate. Extremely large coefficients
may overflow the intermediate FFT; such cases return a numerical failure
status rather than a reliable solution.

## Provenance and validation

The solver was written independently from the standard Fourier
circulant diagonalization identity, not copied from SciPy or FFTPACK.
The existing `scifort_signal_spectral.fft_plan` is reused unchanged.

Public API specifications were verified against:

- https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.circulant.html
- https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_circulant.html

- [Native Fortran tests](../test/test_linalg_circulant.f90)
- [Independent SciPy comparison driver](../tools/linalg_reference/driver_circulant.f90)
- [Reference harness](../tools/linalg_reference/run_circulant_checks.py)
- [Validation report](linalg_circulant_validation_2026-11-09.md)
