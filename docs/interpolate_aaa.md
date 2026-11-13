# AAA rational approximation

`scifort_interpolate` exports `aaa`, `aaa_type`, and `complex_aaa_type` for
adaptive Antoulas-Anderson rational approximation of scalar sampled data.
The implementation module is `scifort_interpolate_aaa`.

```fortran
r = aaa(x, y, rtol=1e-10_dp, max_terms=100, clean_up=.true., &
    clean_up_tol=1e-13_dp, status=info)
```

`x(:)` and `y(:)` must be nonempty vectors of equal length. Coordinates must
be finite. Samples with nonfinite values are discarded; at least one finite
value must remain. Coordinates are sorted lexicographically, and duplicates
retain the first finite supplied value, as in SciPy. Both coordinates and
values use binary64. Real coordinates and values return `aaa_type`; either
complex input returns `complex_aaa_type`. No multichannel or broadcasting
interface is provided.

`rtol` is finite and nonnegative, default `epsilon(1.0_dp)**0.75`.
`max_terms` is a positive integer, default 100; the effective budget is capped
at the number of distinct retained samples. Zero tolerance is allowed, but
floating-point data may require many terms to attain exact interpolation.
`clean_up` defaults to true, and `clean_up_tol` is finite and nonnegative,
default `1e-13`. Invalid controls return an invalid object and status
`aaa_invalid_input`.

## Example

This complete program fits a function whose poles are `+i` and `-i`, residues
are both `1/2`, and only finite zero is zero.

```fortran
program aaa_example
    use scifort_interpolate, only : aaa, aaa_type, aaa_success
    use scifort_kinds, only : dp
    implicit none
    type(aaa_type) :: r
    real(dp) :: x(101), y(101)
    complex(dp), allocatable :: poles(:)
    integer :: i, status

    x = [(-2.0_dp+0.04_dp*real(i-1,dp), i=1,101)]
    y = x/(1.0_dp+x*x)
    r = aaa(x, y, status=status)
    if (status /= aaa_success) error stop 1
    print '(a,es24.16)', 'r(0.5): ', r%value(0.5_dp)
    poles = r%poles(status)
    if (status /= aaa_success) error stop 2
    do i=1,size(poles)
        print '(a,2es24.16)', 'Pole: ', real(poles(i),dp), aimag(poles(i))
    end do
end program aaa_example
```

## Object methods

| Method | Result |
| --- | --- |
| `r%value(q[,status])` | Scalar evaluation at a real or complex query. A real query on a real fit returns real; other combinations return complex. |
| `r%values(q[,status])` | Allocatable vector of evaluations, retaining query order. The optional status reports the first failed query. |
| `r%support_points()` | Copy of retained support coordinates, with the object's real/complex type. |
| `r%support_values()` | Copy of the original supplied values at the support points. |
| `r%weights()` | Copy of barycentric weights; a common complex phase or nonzero scale does not change the function. |
| `r%errors()` | Absolute maximum sample errors after successive greedy iterations, before cleanup. |
| `r%poles([status])` | Allocatable complex vector of finite denominator roots, in unspecified order. |
| `r%roots([status])` | Allocatable complex vector of finite numerator roots, in unspecified order. |
| `r%residues([status])` | Allocatable complex vector of simple-pole residues, paired with `poles()`. |
| `r%clean_up([cleanup_tol,status])` | Integer count of distinct supports removed for weak poles; updates the fit only after successful refitting. |

The following diagnostic methods are **SciFort extensions**:

| Method | Meaning |
| --- | --- |
| `r%valid()` | A usable representation was retained, even if tolerance was not met or cleanup failed. |
| `r%converged()` | The current maximum sample error meets the requested relative tolerance. |
| `r%fit_status()` | Fitting/convergence status, independent of cleanup failure. |
| `r%cleanup_status()` | Last cleanup result; zero also covers cleanup not attempted. |
| `r%max_error()` | Actual absolute maximum error on the retained samples, recomputed after successful cleanup. |

Accessors return independent copies. An invalid object has empty support and
history arrays and a NaN maximum error. Exact queries at retained support
points return the original supplied value, avoiding a divided-and-rescaled
copy that could differ by rounding.

For finite queries the implementation uses scaled barycentric sums to avoid
forming large reciprocals near support points. Query NaNs return NaN with
invalid-input status. Infinite queries use the leading-coefficient ratio
`sum(w*f)/sum(w)` when its denominator is nonzero. This ratio can be unreliable
when leading terms cancel, including polynomial limits. An unresolved zero
denominator, an unrepresentable result, or unresolved coordinate separation
returns NaN and numerical-failure status. In particular, exact poles do not
produce unchecked division by zero. These status-based failure conventions
are intentional differences from SciPy's warnings, exceptions, and possible
infinite results.

## Convergence and cleanup

AAA greedily selects the next sample with largest error, then obtains
barycentric weights from the smallest right singular vectors of the Loewner
matrix. The stopping criterion is the maximum sample error divided by the
maximum sample magnitude, bounded by `rtol`. This controls error at the
supplied samples, not between samples or under extrapolation. More sampling
or problem-specific analysis may be needed near singularities or sharp
features. Overly tight tolerances and noisy data can introduce weak poles.

Cleanup tests each finite pole's residue divided by its distance to the
nearest sample against `cleanup_tol` times the geometric mean of the retained
sample magnitudes. A sample with value zero makes this geometric mean zero,
so this criterion removes no poles, matching the upstream convention.
For each detected weak pole, the closest support point is removed and weights
are refitted. The native return count is the number of distinct support points
removed for detected weak poles; refitting can additionally discard exactly
zero-weight supports, which are not included in that count. A zero cleanup threshold is a successful no-op that skips spectral work.
No printing or warning emission occurs in library code.

Cleanup decisions are made entirely in normalized coordinates and values,
avoiding erroneous decisions when physical residues underflow. A failed
cleanup leaves the representation unchanged. Successful cleanup can loosen
the approximation: `max_error()`, `converged()`, and `fit_status()` describe
the resulting fit, while `errors()` retains the pre-cleanup history.

Always distinguish representation validity, convergence, and spectral
availability. Check constructor `status` together with the diagnostic methods
when automatic cleanup is enabled.

## Status codes and spectral limits

| Constant | Value | Meaning |
| --- | --- | --- |
| `aaa_success` | 0 | Requested operation completed. |
| `aaa_invalid_input` | 1 | Invalid data/control/query or uninitialized object. |
| `aaa_max_terms` | 2 | Tolerance not met: support budget exhausted, or successful cleanup loosened accuracy beyond tolerance. |
| `aaa_numerical_failure` | 3 | SVD/eigenvalue failure, unresolved geometry/arithmetic, or a failed spectral/refitting check. |
| `aaa_spectral_size_limit` | 4 | The retained support set exceeds the spectral backend's current capacity. |
| `aaa_setup_failure` | 5 | A caught workspace allocation failure. |
| `aaa_undefined` | 6 | Zeros requested for an identically zero function, whose zero set is not a finite list. |

Constructor status describes the final fit, or reports an automatic-cleanup
error if cleanup failed. A manual `clean_up` status describes the cleanup
operation; inspect `fit_status()` for the resulting approximation's tolerance.
Spectral failures return empty arrays and their status, retaining the fit.
A nonzero constant has no finite poles or zeros; the zero function has no
poles and reports undefined zeros.

**Poles, zeros, residues, and cleanup support up to 512 retained support
points**, using an internal arrowhead backend adapted from Reference-LAPACK's
`ZHGEQZ` and support kernels. The former 39-point ceiling is removed; the
default 100-point fitting budget is covered. Fitting and evaluation support larger
budgets. If automatic cleanup encounters a larger representation, constructor
status is 4 while `valid()` and `converged()` still describe the retained fit.
Set `clean_up=.false.` to omit automatic spectral work.

Pencils use normalized coordinates and coefficients. Projective roots with
`abs(beta) <= 64*epsilon*abs(alpha)` are treated as infinite and omitted.
Finite candidates are polished and checked against the normalized secular
equation. Exact numerator factors at zero-valued support points are snapped
to the corresponding support within a normalized distance of `1e-8`.
Very distant roots, close clusters, repeated poles, and severe degree
cancellation remain conditioning limits; the residue formula assumes simple
poles. The QZ backend does not remove these conditioning limits, and it is
restricted to normalized AAA pencils rather than arbitrary matrix pairs.
Binary64 underflow can also erase very small normalized samples or physical
residues after rescaling.

The backend reduces the arrowhead matrix to Hessenberg form by similarities
that leave `diag(0,I)` unchanged, then runs the licensed eigenvalue-only QZ
kernel. Its dense storage is quadratic and arithmetic cost is cubic in the
support count. Failed QZ iteration or an unresolvable secular root returns
status 3; the retained approximation remains valid. Status 4 now applies above
512 retained supports. This change adds no linked LAPACK/BLAS dependency and
does not change the order limits of other public matrix routines.

SciFort's native SVD can choose different supports, weights, degrees, and
cleanup paths from SciPy. The spectral backend is separately validated
against SciPy and analytic roots through 512 supports; larger public real and
complex fits also cover residues and actual cleanup removal. See the
[larger-spectrum validation record](../validation/aaa_large_2026-11-13.md).
Native status diagnostics replace Python warnings/exceptions. C ABI and
wrapper expansion remain deferred.

## References

- Y. Nakatsukasa, O. Sete, and L. N. Trefethen, *The AAA algorithm for rational
  approximation*, SIAM J. Sci. Comput. 40(3), A1494-A1522 (2018),
  [DOI](https://doi.org/10.1137/16M1106122), [preprint](https://arxiv.org/abs/1612.00337).
- [SciPy AAA API](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.AAA.html).
- [Provenance](../CODE_PROVENANCE.md), [retained licenses](../LICENSES/scipy-aaa-BSD-3-Clause.txt),
  and [validation record](../validation/aaa_2026-11-13.md).
- [Reference-LAPACK v3.12.1](https://github.com/Reference-LAPACK/lapack/tree/v3.12.1)
  and its [retained license](../LICENSES/lapack-BSD-3-Clause.txt).
