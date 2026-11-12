# Spectral-norm estimates (`scipy.linalg.interpolative`)

New module: `src/linalg_spectral_norm.f90` (original MIT-licensed SciFort code).

## Algorithms and interfaces

The two SciPy analogues are
`scipy.linalg.interpolative.estimate_spectral_norm(A, its=20, rng=None)` and
`scipy.linalg.interpolative.estimate_spectral_norm_diff(A, B, its=20, rng=None)`.

SciFort estimates the matrix operator 2-norm using the power method:

1. Initialize a nonzero `x` and normalize it in the Euclidean norm.
2. For `its` steps, compute `y = A x`, normalize `y`, then compute
   `x = A^H y` and normalize `x`.
3. Return `||A x||_2` after the final step.

The difference routine applies the same algorithm to `A-B`. The adjoint is
**conjugate transpose** for complex matrices. Dense input is scaled before
iteration to avoid unnecessary over/underflow; near-cancelling differences
are subtracted before scaling when the subtraction is representable.

These are **approximate randomized estimators**, not SVDs or certified
upper/error bounds. The estimate is normally a lower bound up to rounding.
More iterations help when the two largest singular values are close.

### Dense matrix forms

```fortran
use, intrinsic :: iso_fortran_env, only : int64
use scifort_linalg_spectral_norm, only : estimate_spectral_norm, &
    estimate_spectral_norm_diff

real(dp) :: a(100,50), b(100,50), norm_a, norm_delta
integer :: status
! ... initialize a and b ...
call estimate_spectral_norm(a,norm_a,status,its=40,seed=123_int64)
call estimate_spectral_norm_diff(a,b,norm_delta,status,its=40,seed=123_int64)
```

Both generic interfaces accept either `real(dp)` or `complex(dp)` matrices.
For complex operands, the result `estimate` is real(dp). Matrix shapes must
agree for the difference interface. Rectangular and zero-sized matrices work.

### Matrix-free / sparse operator forms

For a real matrix-free operator, implement two procedures with the interfaces

```fortran
subroutine matvec(x,y,status)
    real(dp), intent(in) :: x(:)     ! n elements
    real(dp), intent(out) :: y(:)    ! m elements: A*x
    integer, intent(out) :: status  ! 0 on success
end subroutine
subroutine rmatvec(x,y,status)
    real(dp), intent(in) :: x(:)     ! m elements
    real(dp), intent(out) :: y(:)    ! n elements: transpose(A)*x
    integer, intent(out) :: status
end subroutine
```

and call

```fortran
call estimate_spectral_norm_operator(m,n,matvec,rmatvec, &
    estimate,status,its=40,seed=123_int64)
call estimate_spectral_norm_diff_operator(m,n,matvec_a,rmatvec_a, &
    matvec_b,rmatvec_b,estimate,status,its=40)
```

For `complex(dp)` callbacks, `rmatvec` **must return `A^H*x`**, not
`transpose(A)*x`. Because Fortran cannot distinguish generic procedures
solely by the type of their *procedure arguments*, the complex entry points
have explicit names:

- `estimate_spectral_norm_operator_complex`
- `estimate_spectral_norm_diff_operator_complex`

Both real and complex operator forms also accept an optional initial vector.
This design supports matrices in CSR format, structured FFT operators and
other implicit matrices **without converting them to dense arrays**. The
application is responsible for the callback implementation.

### Arguments and error handling

- `its` — optional number of complete power iterations, default `20`;
  can be zero for a single initial matrix-vector evaluation; must be >= 0.
- `seed` — optional `integer(int64)` local Park–Miller RNG seed; identical
  seeds reproduce SciFort results without changing the global Fortran RNG.
  Seed values intentionally do not match NumPy's random sequences.
- `start` — optional nonzero vector of length `n`; takes precedence over `seed`.
  Useful for deterministic testing and when a good initial guess is known.
- `estimate` — real(dp) scalar output, initialized to zero on errors.
- `status` — `linalg_status_success` (0), `linalg_status_invalid_shape` (1)
  for bad dimensions/iteration counts/starting vectors, or
  `linalg_status_no_convergence` (4) for invalid matrix values, callback
  failures or nonfinite iterates.

A zero operator or zero-dimensional input returns a zero norm. A nonzero
operator multiplied by a starting vector in its null space can also return
zero; selecting a different starting vector is recommended. For implicit
operators, overflow inside a user-supplied product cannot be prevented by
SciFort. Large matrices should use the operator API, which requires only
`O(m+n)` working storage; dense calls retain an `O(m*n)` scaled copy.

### Reproducible testing

```sh
# Run from project root, with gfortran and SciPy/NumPy installed
mkdir -p /tmp/scifort-spectral
cd /tmp/scifort-spectral
gfortran -O0 -g -fcheck=all -c /path/to/scifort/src/{kinds,linalg,linalg_spectral_norm}.f90
gfortran -O0 -g -fcheck=all /path/to/scifort/test/test_linalg_spectral_norm.f90 *.o -o native
./native
gfortran -O0 -g -fcheck=all /path/to/scifort/validation/reference_drivers/driver_linalg_spectral_norm.f90 *.o -o driver
python /path/to/scifort/test/compare_linalg_spectral_norm_scipy.py ./driver
```

The Python script explicitly supplies initial vectors to reproduce the same
power iterates and *also* compares independent SciPy randomized estimates and
exact NumPy SVD norms. Since SciPy's RNG and starting vectors differ, two
randomized estimates are **not expected to agree entry-for-entry**.

See [validation report](linalg_spectral_norm_validation_2026-11-09.md) and
[executable example](../example/linalg_spectral_norm_example.f90).
