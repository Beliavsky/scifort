# Standalone Fourier, cosine, and sine transforms

`scifort_fft` provides native binary64, one-dimensional and multidimensional transforms using the
existing SciFort mixed-radix/Bluestein FFT engine. Signal-analysis interfaces
are unchanged. Import `dp` from `scifort_kinds` or `scifort_fft`.

## Available API

The one-dimensional functions below return allocatable rank-one arrays.
Inputs are read-only and may be noncontiguous sections.

| Procedure | Input | Result |
| --- | --- | --- |
| `fft(x [, n, norm, status])` | Real or complex | Full complex forward DFT |
| `ifft(x [, n, norm, status])` | Real or complex | Full complex inverse DFT |
| `rfft(x [, n, norm, status])` | Real | Complex nonnegative-frequency half |
| `irfft(x [, n, norm, status])` | Real or complex half spectrum | Real full inverse |
| `hfft(x [, n, norm, status])` | Real or complex Hermitian half sequence | Real full forward transform |
| `ihfft(x [, n, norm, status])` | Real | Complex Hermitian half sequence |
| `dct`, `idct`, `dst`, `idst` | Real or complex | Same type as input |
| `fftfreq(n [, d, status])` | Positive length, sample spacing | Full signed frequency bins |
| `rfftfreq(n [, d, status])` | Positive length, sample spacing | Nonnegative frequency bins |
| `fftshift(x)`, `ifftshift(x)` | Real or complex | Shifted sequence of the same type |

Cosine/sine signatures are:

```text
dct(x [, type, n, norm, orthogonalize, status])
idct(x [, type, n, norm, orthogonalize, status])
dst(x [, type, n, norm, orthogonalize, status])
idst(x [, type, n, norm, orthogonalize, status])
```

The reusable `fft_plan` type is also re-exported, with its existing
`initialize`, `execute`, `inverse`, `length`, and `clear` methods. `execute`
is an unnormalized forward complex FFT; `inverse` divides by the plan length.
Plans contain mutable scratch space and must not be shared by concurrent
calls. Convenience functions use local plans with no global cache.

## Example

```fortran
program fft_example
    use scifort_fft, only : dp, rfft, irfft, dct, idct
    implicit none
    real(dp) :: x(5)
    real(dp), allocatable :: recovered(:), coefficients(:)
    complex(dp), allocatable :: spectrum(:)
    integer :: status

    x = [1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp, 3.0_dp]
    spectrum = rfft(x, status=status)
    if (status /= 0) error stop 'rfft failed'
    recovered = irfft(spectrum, n=size(x))
    coefficients = dct(x, norm='ortho')
    recovered = idct(coefficients, norm='ortho')
end program fft_example
```

## Length and normalization contracts

`n` is a positive transform length. Input is cropped or zero-padded at the
end before transformation. An empty input is allowed when an explicit
positive `n` supplies the length; otherwise it is invalid.

For `fft`, `ifft`, `rfft`, `ihfft`, and cosine/sine transforms, the default
length is the input size. For `irfft` and `hfft`, it is `2*(size(x)-1)`;
therefore a one-element half input needs an explicit positive `n`.
Supply `n` explicitly to recover an odd full length from a half spectrum.
A half result has `n/2+1` entries using integer division. Imaginary parts of
DC and the even-length Nyquist entry of a supplied Hermitian half are ignored.
An odd-length final half entry is not constrained to be real.

`norm` is `backward` by default:

| Normalization | Forward transform | Inverse transform |
| --- | --- | --- |
| `backward` | Unscaled | Divide by logical length |
| `forward` | Divide by logical length | Unscaled |
| `ortho` | Divide by square root of logical length | Same |

For Fourier transforms, logical length is `n`. For DCT-I it is `2*(n-1)`;
for DST-I, `2*(n+1)`; for types II-IV, `2*n`. An inverse cosine/sine call
uses the inverse of the requested type: types II and III interchange,
while types I and IV are self-inverse with the stated scaling.

Types I-IV are supported; type II is the default. DCT-I requires `n > 1`.
DST-I supports `n=1`, matching the tested SciPy 1.15.3 implementation.
`orthogonalize` defaults to true for `norm='ortho'`, false otherwise, and
may be overridden independently of `norm`. It applies SciPy's documented
endpoint weights for DCT-I/II/III and DST-II/III. DCT-IV and DST-I/IV need no
endpoint correction. Complex inputs follow the same definitions by linearity.
Cosine/sine types use FFT-based even/odd extensions with at most `8*n`
complex embedding entries, rather than quadratic summation. Requests with
`n > huge(default_integer)/8` are rejected before extension arithmetic.

## Bins and shifts

Sample spacing `d` defaults to 1 and must be finite and nonzero. Negative
spacing reverses frequency signs. Frequencies are in cycles per unit of `d`.
The zero-frequency entry is at Fortran index 1. For even `n`, `fftfreq` puts
Nyquist in the negative half; `rfftfreq` includes positive Nyquist.
The computation avoids the potentially overflowing product `n*d`.
Unrepresentable frequency bins reject the call with status 1.

`fftshift` moves zero frequency to the center; `ifftshift` reverses that
permutation. Their shifts differ for odd lengths. Empty sequences are valid,
and NaN/infinity entries are preserved by these permutations.

## Multidimensional transforms and axis selection

The following subroutines accept arrays of any rank from 1 through 15:

```text
call fftn(x, y [, s, axes, norm, status])
call ifftn(x, y [, s, axes, norm, status])
call rfftn(x, y [, s, axes, norm, status])
call irfftn(x, y [, s, axes, norm, status])
call hfftn(x, y [, s, axes, norm, status])
call ihfftn(x, y [, s, axes, norm, status])
call dctn(x, y [, type, s, axes, norm, orthogonalize, status])
call idctn(x, y [, type, s, axes, norm, orthogonalize, status])
call dstn(x, y [, type, s, axes, norm, orthogonalize, status])
call idstn(x, y [, type, s, axes, norm, orthogonalize, status])
```

Input and output must have the same rank. Allocate `y` to the expected shape
before calling; the routines do not allocate or resize it. Input and output
must not overlap. Strided sections are supported through compiler-managed
contiguous temporaries; large strided arrays can require additional memory.
The real/complex combinations match the one-dimensional table. DCT/DST input
and output have the same type.

`axes` contains distinct **one-based**, positive dimension numbers in the
requested order. With no `axes`, all dimensions are transformed, or the last
`size(s)` dimensions when `s` is supplied. With no `s`, lengths come from the
selected input dimensions, except for the inferred full length of
`irfftn`/`hfftn` described below. Every supplied length must be positive.
An empty `axes` list gives an identity for complex FFT and cosine/sine calls;
real/Hermitian pairs require at least one axis.

`s(i)` specifies the full transform length along `axes(i)`, independently
cropping or zero-padding that dimension. Unselected dimensions retain their
input extent. For `rfftn` and `ihfftn`, the **last selected axis** has output
length `s(last)/2+1`; other transformed extents have their full lengths.
For `irfftn` and `hfftn`, `s` gives full output lengths. Without `s`, the
last selected output extent is `2*(m-1)`, where `m` is that input extent.
Supply `s` for an odd full length or a length-one inverse.

The real-to-half axis is transformed first, followed by complex transforms
along the other selected axes. Half-to-real calls reverse these stages.
Endpoint reality is imposed after the other axes have been transformed,
so Hermitian symmetry couples all transformed dimensions correctly.
Each axis receives its own normalization; the total scale is their product.
Cosine/sine transforms apply the requested type and orthogonalization
independently along each selected axis.

```fortran
program fft_nd_example
    use scifort_fft, only : dp, rfftn, irfftn
    implicit none
    real(dp) :: image(3,5), recovered(3,5)
    complex(dp) :: spectrum(3,3)
    integer :: status

    image = 1.0_dp
    call rfftn(image, spectrum, axes=[1,2], norm='ortho', status=status)
    if (status /= 0) error stop 'rfftn failed'
    call irfftn(spectrum, recovered, s=[3,5], axes=[1,2], norm='ortho', status=status)
    if (status /= 0) error stop 'irfftn failed'
end program fft_nd_example
```

These routines use separability and the existing one-dimensional kernels.
Fourier plans are reused across all lines of a given axis within a call.
Cosine/sine calls currently set up their one-dimensional FFT extensions per
line. There is no public multidimensional plan or thread control.
Input, output, and conservative intermediate bounding shapes must fit the
default-integer element count. Nonfinite values are rejected even in cropped
regions. Zero input extents can be padded using explicit positive `s`;
unselected zero extents give an empty batch.

## Errors and compatibility limits

Optional `status` reports:

- `fft_status_success = 0`: successful computation.
- `fft_status_invalid_input = 1`: invalid length, normalization, transform
  type, nonfinite transform data, invalid spacing, or unrepresentable bins.
- `fft_status_setup_failure = 2`: caught allocation or FFT plan setup failure.

Rejected one-dimensional function calls return an allocated zero-length
result. Multidimensional calls reject invalid axes, mismatched lengths,
wrong output rank/shape, or invalid options with status 1. Representable
rank 1-15 outputs are filled with NaNs on failure. Transform inputs must
be finite, including entries beyond a requested crop; this intentionally
differs from SciPy's propagation of NaN/infinity. Shifts only permute data and
do not apply this restriction. Binary64 kernel range and rounding limits
still apply; extended precision is not used.

SciFort honors explicit `orthogonalize` for complex DCT/DST inputs. The local
SciPy 1.15.3 complex path ignores that argument; comparisons therefore apply
SciPy's real transforms separately to real and imaginary components.
For very large spacing, SciFort retains finite subnormal bins where the
tested SciPy/NumPy frequency helper's `n*d` intermediate overflows. Those
cases use independent high-precision references. Both differences are
deliberate and covered by validation.

Fast Hankel transforms, fast-length helpers, multidimensional shifts,
two-dimensional convenience names, worker controls, vendor backends, and
overwrite options remain open. Use `fftn` with two selected axes for a 2D FFT.
Negative or repeated axes and `s=-1` sentinel lengths are not accepted.
No C ABI is added. See the [gap inventory](scipy_gaps.md) and
[roadmap](../ROADMAP.md) for follow-up scope.

## Verification and provenance

`test/test_fft.f90` checks independent DFT and cosine/sine matrices, inverse
round trips, all normalizations and endpoint weights, orthonormal energy,
real/complex agreement, lengths 1-33, padding/cropping, strided inputs,
input preservation, Hermitian endpoints, frequency signs, shifts, nonfinite
input, invalid options, and integer-range rejection.

`validation/check_fft.py` compares with SciPy at smooth, prime, odd/even,
padded, and cropped lengths up to 1009, with complex and real inputs, all
four cosine/sine types, all normalizations, explicit/default orthogonalization,
and amplitudes `1e-100` and `1e100`. Frequency edge cases use mpmath.

`test/test_fft_nd.f90` adds an independent two-dimensional DFT, real and
Hermitian round trips, axis subsets/order, ranks 1/2/3/15, padding,
noncontiguous input/output, and invalid/nonfinite input checks.
`validation/check_fft_nd.py` compares rank-four batches against SciPy with
all transform families, selected and reordered axes, cropping/padding,
all normalization/type/orthogonalization settings, singleton/empty inputs,
and large/small amplitudes. See its
[validation record](../validation/fft_nd_2026-11-12.md).

See the [dated validation record](../validation/fft_2026-11-12.md) for results
and exact commands, and [provenance](../CODE_PROVENANCE.md) for sources.
Definitions follow the public SciPy specifications for
[FFT](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.fft.html),
[DCT](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.dct.html),
[DST](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.dst.html),
[Hermitian FFT](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.hfft.html),
and [frequency bins](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.fftfreq.html).
