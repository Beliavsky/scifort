# Hadamard, Leslie, and DFT matrices

SciFort implements the matrix constructors [`scipy.linalg.hadamard`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.hadamard.html), [`scipy.linalg.leslie`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.leslie.html), and [`scipy.linalg.dft`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.dft.html) in three dependency-free Fortran modules.

All routines use caller-allocated output arrays and return `status` (`0` for success and `1` for an invalid order, shape, or option). On failure, the output array is filled with zeros. They are `pure` and require only the existing `scifort_kinds` and `scifort_linalg` modules.

## `scifort_linalg_hadamard`: `hadamard`

```fortran
use, intrinsic :: iso_fortran_env, only : int32, int64
use scifort_kinds, only : dp
use scifort_linalg_hadamard, only : hadamard
integer(int64) :: h(8,8)
integer :: status
call hadamard(8,h,status)
```

- `n` must be a **positive power of two** (`n=1` is allowed).
- `h` must have shape `(n,n)`; overloads accept `integer(int32)`, `integer(int64)`, `real(dp)`, or `complex(dp)` outputs.
- Uses the Sylvester ordering: `h(i,j) = (-1)**popcnt(iand(i-1,j-1))`. Thus `matmul(h,transpose(h))=n*I` exactly for integer matrices of practical sizes.
- No temporary matrix is required. Construction time is `O(n^2)` with `O(n^2)` output storage and `O(1)` auxiliary storage.

## `scifort_linalg_leslie`: `leslie`

```fortran
use scifort_kinds, only : dp
use scifort_linalg_leslie, only : leslie
real(dp) :: a(4,4)
integer :: status
call leslie([0.1_dp,2.0_dp,1.0_dp,0.1_dp], &
    [0.2_dp,0.8_dp,0.7_dp],a,status)
```

The vector form takes fecundities `f(n)` and survival coefficients `s(n-1)`, with `n>=2`; it fills `a(1,:)=f`, `a(i+1,i)=s(i)`, and zero everywhere else.

- Real64, complex64, and homogeneous `integer(int32)` / `integer(int64)` vectors are supported. Mixed real64/complex64 inputs produce complex64 output.
- Rank-2 inputs `f(bf,n)` and `s(bs,n-1)` create a rank-3 output `a(b,n,n)` where `b=max(bf,bs)`. One of `bf` or `bs` may be `1` and will be broadcast to `b` matrices; otherwise they must match. This overload supports real, complex, and mixed real/complex inputs.
- The batch dimension comes **first** in the output array. This matches NumPy's logical shape, though the Fortran memory layout is column-major.
- SciFort currently supports one leading batch dimension; SciPy permits multiple batch dimensions. Integer inputs do not yet have a batched overload, and unlike NumPy, mixed integer/real inputs need explicit conversion by the caller.
- No range restriction is imposed on fecundity or survival coefficients, matching SciPy's constructor.

## `scifort_linalg_dft`: `dft`

```fortran
use scifort_kinds, only : dp
use scifort_linalg_dft, only : dft
complex(dp) :: a(8,8)
integer :: status
call dft(8,a,status,scale='sqrtn')
```

- `n >= 0`, with `a(n,n)`. The empty `n=0` case is accepted.
- Optional `scale` is `'none'` (default, representing SciPy's `None`), `'sqrtn'` (divide by `sqrt(n)`), or `'n'` (divide by `n`). Strings are case-sensitive.
- The matrix entries are `exp(-2*pi*1j*(i-1)*(j-1)/n)` times the optional scale, computed with modular integer phase reduction.
- `scale='sqrtn'` is unitary: `matmul(a,conjg(transpose(a)))=I` to roundoff.
- This is an `O(n^2)` matrix constructor; `matmul(a,x)` costs `O(n^2)` as well. Use SciFort's existing FFT APIs for an `O(n log n)` transform of a vector.
- Only double-precision complex output is presently provided.

## Running the example and independent tests

See [`example/linalg_hadamard_leslie_dft_example.f90`](../example/linalg_hadamard_leslie_dft_example.f90) and the [validation report](linalg_hadamard_leslie_dft_validation_2026-11-09.md).

From the project root:

```sh
mkdir -p build_hadamard_leslie_dft
cd build_hadamard_leslie_dft
gfortran -std=f2018 -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow \
  -c ../src/kinds.f90 ../src/linalg.f90 ../src/linalg_hadamard.f90 \
     ../src/linalg_leslie.f90 ../src/linalg_dft.f90
gfortran -std=f2018 -O0 -g -fcheck=all -o test_special \
  ../test/test_linalg_hadamard_leslie_dft.f90 *.o
./test_special
gfortran -std=f2018 -O0 -g -fcheck=all -o special_driver \
  ../reference/linalg_hadamard_leslie_dft_driver.f90 *.o
python ../tools/compare_linalg_hadamard_leslie_dft.py ./special_driver
```

The reference script requires NumPy and SciPy but no Python-to-Fortran bindings. Test results in the [validation report](linalg_hadamard_leslie_dft_validation_2026-11-09.md).
