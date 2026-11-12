# Logarithmic-time maximum-length-sequence skip and rewind

`scifort_signal::max_len_seq` now accepts an optional final
`integer(int64) :: skip` argument. This is a **SciFort extension** to
SciPy's `signal.max_len_seq` convention.

- `skip > 0`: advance the specified Fibonacci register before generating bits.
- `skip < 0`: rewind it by `abs(skip)` register clocks before generating bits.
- `skip = 0` or omitted: preserve the pre-existing behavior.
- `length=0`: compute only the jumped register state, without sample allocation.
- `state=...`: specify any nonzero state, not just an all-ones seed.
- `taps=...`: choose the same custom tap positions used by SciPy.

```fortran
use, intrinsic :: iso_fortran_env, only : int64
use scifort_signal, only : max_len_seq, max_len_seq_result

type(max_len_seq_result) :: s, back
s = max_len_seq(8, length=20, skip=1000000_int64)
back = max_len_seq(8, state=s%state, length=0, skip=-1000020_int64)
! back%state is all ones; s%seq contains bits at offsets 1000000..1000019.
```

## Numerical method

Let the register be the binary column vector `s` and let the one-step
Fibonacci transition be `T`, defined by shifting every bit down and placing
the XOR feedback in the last bit. Then the register after `k` steps is
`T**k s` in GF(2). For a negative `k`, use the inverse transition: its
first bit is the final feedback bit XOR the tapped shifted bits, and its
remaining bits are the shifted register. The matrix is always invertible
because bit 1 is included in the feedback XOR.

Represent the 32-or-fewer columns of `T` as 64-bit integers, then compute
`T**k` by squaring the column representation, applying selected powers to
the packed seed. Time complexity is `O(nbits**2 log(abs(k)))`; storage is
`O(nbits)` packed integers. The minimum signed 64-bit `skip` is supported
by taking one reverse step before exponentiation, avoiding signed overflow.
No skipped sequence samples or large dense matrices are allocated.

The existing positional `max_len_seq(nbits,state,length,taps)` calls are
unchanged. Invalid widths, seeds, taps, and lengths retain their earlier
status codes. No new runtime dependencies or licenses are introduced.

## Validation

```bash
python tools/run_signal_sequence_skip.py
```

This builds the modified module with GNU Fortran in a checked `-fcheck=all`
mode and again with `-O3 -march=native`. It runs both native test programs
and checks 698 deterministic SciPy 1.17 comparisons per build (sample bits
and final register states). The reference constructs its transition matrix
from SciPy's one-sample output for every basis seed, inverts the matrix via
independent GF(2) Gauss-Jordan elimination, and raises it to large positive
or negative powers. The native tests also check exact direct-generation
prefixes, primitive full periods, custom nonprimitive/duplicate taps,
reversibility, continuation, and the most-negative 64-bit integer.
