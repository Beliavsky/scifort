# Clarkson–Woodruff (CountSketch) matrix sketching

`scifort_linalg_clarkson_woodruff` implements the CountSketch random projection
underlying `scipy.linalg.clarkson_woodruff_transform`. For an input `A(n,d)`,
the output `B(k,d)` is

\[
B_{h(i),:} \mathrel{+}= s(i) A_{i,:}, \qquad 1\leq i\leq n,
\]

where each source row independently chooses a destination `h(i)` from `1:k`
and a sign `s(i)` in `{+1,-1}`. The resulting sketch often approximately
preserves Euclidean norms. Accuracy is probabilistic, not guaranteed for an
individual sketch. See [SciPy CountSketch documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.clarkson_woodruff_transform.html).

## Interfaces

```fortran
use, intrinsic :: iso_fortran_env, only : int64
use scifort_kinds, only : dp
use scifort_linalg_clarkson_woodruff, only : countsketch_plan, &
    make_countsketch_plan, clarkson_woodruff_transform, clarkson_woodruff_csr

type(countsketch_plan) :: plan
real(dp) :: a(500,20), sa(80,20)
integer :: status

! Initialize a ...
call make_countsketch_plan(500, 80, plan, status, seed=12345_int64)
if (status /= 0) error stop 'Invalid sketch dimensions'
call clarkson_woodruff_transform(a, 80, sa, status, plan=plan)
if (status /= 0) error stop 'Sketch could not be computed'
```

The `clarkson_woodruff_transform` generic offers rank-2 `real(dp)` and
`complex(dp)` variants and rank-3 batch variants with shapes
`(batch,n,d)` to `(batch,k,d)`. Batched matrices all use the **same** random
mapping, matching the linear transformation across slices. Output arrays
are preallocated and the routine sets them to zero before processing.

`make_countsketch_plan(n,k,plan,status,seed)` creates a reusable mapping.
An optional `seed` uses an independent, deterministic 31-bit Park–Miller
pseudorandom generator. No seed uses the Fortran intrinsic RNG.
The same integer seed does **not** produce the same NumPy-generated mapping.
For exact SciPy/NumPy reproducibility, construct a `countsketch_plan` using
SciPy's one-based row buckets and `+1`/`-1` signs, then pass `plan=...`.
`plan` and `seed` cannot both be supplied in one transform call.

CSR sparse arrays use a separate generic:

```fortran
call clarkson_woodruff_csr(row_ptr,col_ind,values,ncols,k,sa,status,plan=plan)
```

Here `row_ptr` has length `n+1` and uses **one-based** offsets,
`col_ind` contains one-based column indices, and `values` has the same length
as `col_ind`. Duplicated columns and unsorted column indices are supported;
their contributions are summed. The sparse interface supports real and
complex data but does not depend on a particular SciFort sparse matrix type.

The implementation accumulates dense data in `O(nd)` time and CSR data in
`O(nnz(A))` time, after `O(kd)` zero initialization. Sketch plans require
`O(n)` integer storage. Large output magnitudes or nonfinite inputs return a
nonzero status and a zeroed output. For `n=0`, a `k=0` sketch is permitted;
for nonempty matrices `k` must be positive. No new dependencies are needed.

## Validation

- Native test: `test/test_linalg_clarkson_woodruff.f90`
- Comparison driver: `validation/reference_drivers/driver_linalg_clarkson_woodruff.f90`
- Independent test harness: `test/compare_linalg_clarkson_woodruff_scipy.py`
- Runnable example: `example/linalg_clarkson_woodruff.f90`

The Python harness gets signed mappings from SciPy's `cwt_matrix`, then
compares each Fortran output to the **public**
`scipy.linalg.clarkson_woodruff_transform` with the same NumPy Generator
seed. Consequently these tests check the actual operation, not agreement of
unrelated random-number streams. See
[`linalg_clarkson_woodruff_validation_2026-11-09.md`](linalg_clarkson_woodruff_validation_2026-11-09.md)
for the reproducible commands and results.

## Scope differences from SciPy

- No Python-style `rng` object is accepted; callers use an explicit plan,
  local int64 seed, or the intrinsic Fortran RNG.
- CSR data use arrays of offsets, indices, and values rather than SciPy
  sparse objects.
- One batch dimension is supported rather than arbitrary batch rank.
- Only real(dp) and complex(dp) are supported; integer input promotion is
  not included in this version.
