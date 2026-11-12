# Matrix exponential (`scifort_linalg_expm`)

SciFort exposes a standalone `scipy.linalg.expm`-style routine:

```fortran
use scifort_kinds, only : dp
use scifort_linalg_expm, only : expm
real(dp) :: a(2,2), x(2,2)
integer :: status, scalings

a = 0.0_dp
a(1,2) = -1.0_dp
a(2,1) = 1.0_dp
call expm(a, x, status, scalings=scalings)
if (status /= 0) error stop 'matrix exponential failed'
```

### Interfaces

The generic `expm` has four overloads:

| Type | Input | Result |
|---|---|---|
| `real(dp)` | `(n,n)` | `real(dp) (n,n)` |
| `complex(dp)` | `(n,n)` | `complex(dp) (n,n)` |
| `real(dp)` | `(batch,n,n)` | `real(dp) (batch,n,n)` |
| `complex(dp)` | `(batch,n,n)` | `complex(dp) (batch,n,n)` |

All output arrays are caller-allocated. Two-dimensional routines accept the
optional integer `scalings` output, giving the number of squarings. Batch
routines use one **leading** batch dimension rather than arbitrary batch ranks.

Status values come from `scifort_linalg`:

- `linalg_status_success` (0): completed.
- `linalg_status_invalid_shape` (1): shape incompatibility.
- `linalg_status_no_convergence` (4): nonfinite input, an unsafe Pade/LU
  computation, nonfinite result, or scaling count beyond the supported range.

Output is initialized to zero on failure. Empty square matrices and empty
batches are supported. Zero matrices give identity matrices; diagonal matrices
use elementwise exponentiation, avoiding a full Pade solve.

### Algorithm and limitations

The implementation forms the fixed [13/13] Pade rational approximant to
`exp(A/2**s)` and applies `s` squarings. It selects the scaling count using the
matrix 1-norm and Higham's degree-13 backward-error threshold, and solves the
Pade denominator system using SciFort's pivoted LU implementation.

This is a **fixed-order**, dense algorithm. SciPy 1.17.0 uses adaptive Pade
orders and additional optimizations. The results are intended to agree
numerically, not to reproduce SciPy's floating-point operation order or speed.
It does not provide arbitrary batch-rank broadcasting, `float32`, or
`complex64` overloads. Very large exponentials may overflow `real(dp)`; the
routine returns a failure status rather than silently reporting success.

### Testing

From the repository root, with `gfortran`, NumPy and SciPy installed:

```shell
python tools/linalg_reference/run_expm_checks.py
```

The script builds a checked and an optimized executable in a temporary
directory, runs native tests and existing dense-LU regressions, and performs
an independent SciPy 1.17.0 comparison. No FPM installation is needed.
See `validation/linalg_expm_2026-11-10.md`.

### References

- Nicholas J. Higham, *Functions of Matrices: Theory and Computation*, SIAM,
  2008 (the standard Pade scaling-and-squaring algorithm).
- SciPy 1.17.0 API, `scipy.linalg.expm`:
  https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.linalg.expm.html
