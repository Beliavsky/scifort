# Explicit-state spline continuation

Native `splrep` and `splprep` support `task=1` with a caller-owned
`fitpack_continuation_type`, exported by `scifort_interpolate`. Scalar and
parametric, ordinary and periodic fits support degrees 1-5, positive residual
targets and interpolation (`s=0`). No global cache or new dependency is used.

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: splrep,fitpack_tck_type,fitpack_continuation_type
implicit none
type(fitpack_tck_type) :: fit
type(fitpack_continuation_type) :: state
real(dp) :: u(61),y(61)
integer :: i,status

u=[(real(i-1,dp)/60.0_dp,i=1,61)]
y=cos(2*acos(-1.0_dp)*u)+0.08_dp*sin(14*acos(-1.0_dp)*u)
fit=splrep(u,y,s=0.05_dp,state=state,status=status)
if (status==0) fit=splrep(u,y,s=0.02_dp,task=1,state=state,status=status)
! Check status and fit%valid() before evaluation.
call state%reset()
```

## State contract

- `task=0` (default) with `state=` starts a new history and replaces that
  state's previous history only on success. A failed initialization preserves
  any earlier successful state.
- `task=1` requires a valid `state=` from the same tuple fitting family,
  degree, periodic flag, physical sites, active data, effective weights and
  bounds. All these inputs must match exactly; only `s` changes. Default
  weights match an explicitly supplied array of ones. Default degree and
  bounds retain their usual meanings; they are not inferred from state.
- `splprep` may generate chord parameters or accept `u`. Repeated automatic
  parameterization is deterministic for identical active data; if an earlier
  call supplied `u`, supply it again unless the generated parameters match
  exactly. The closing leg is included for periodic chord parameters.
- Periodic final data and weights remain ignored and may change, even to
  nonfinite values. All physical sites, including the final period endpoint,
  must still match. Active samples and weights retain their validation rules.
- `task=-1` with `state=` and any supplied `t` in `task=0/1` are invalid.
  Fixed-knot least squares does not initialize continuation.
- `state%valid()` reports whether a successful history is available. A default
  state is invalid; `state%reset()` releases its storage and is idempotent.
  Intrinsic assignment makes an independent deep copy. Distinct states can
  be interleaved without interference.
- A failed call returns an invalid tuple and preserves the supplied state.
  Status 1 denotes invalid input or incompatible/missing state; status 3
  denotes numerical/range/allocation failure. A successful call commits the
  metadata and numerical history together.

Existing calls without `state=` retain their contracts. A stateless
`task=1` remains invalid. Returned scalar coefficients are padded, parametric
coefficients unpadded, and parameter, endpoint, extrapolation, and residual
conventions remain those in the [tuple guide](interpolate_fitpack_tck_backend.md).
Positive residuals retain the relative `0.001` tolerance. Interpolation with
`state=` uses the attributed FITPACK kernels and their rounding checks;
ordinary zero-target bounds must equal the sample endpoints.

## Reused work and compatibility

The state saves **normalized knots, per-interval residual/count information,
and knot-refinement history** from the existing attributed `fpcurf`, `fppara`,
`fpperi`, or `fpclos` kernel. These drive the genuine `iopt=1` branches.
Observation QR factors and smoothing iterations are rebuilt; this interface
does not promise factorization reuse or a fixed speedup. Stored input copies
permit exact compatibility checks. Persistent storage scales with samples
times channels plus active knots, and assignment/call transactions copy it.

Continuation can retain more knots than a fresh fit at the same `s`. It can
restart from the polynomial/constant limit when appropriate. Thus comparison
with a fresh `task=0` fit alone is not a continuation correctness test.
The independent validation keeps the exact upstream work arrays across calls.

SciPy's public tuple interfaces describe continuation with internal stored
state. SciFort deliberately requires an explicit object for reentrancy and
independent histories. It does not reproduce Python's global-cache behavior.
The native closed-curve kernel also corrects an inspected upstream restart
capacity check: restarting after interpolation needs room for the minimal
seed, not room beyond the previous interpolation knot count.

Existing finite-range, rank and final physical-residual checks still apply.
Continuation does not cure ill-conditioning or guarantee every smoothing
request succeeds. The known moderate-target scalar failure with short cyclic
bases is [corrected](../validation/periodic_short_basis_2026-11-14.md) by ordered
coefficient extension. Broader periodic conditioning and numerical reliability
remain separate audits. C ABI and wrappers are deferred.

See the [validation record](../validation/continuation_2026-11-14.md),
[provenance](../CODE_PROVENANCE.md), and the public SciPy
[`splrep`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.splrep.html)
and [`splprep`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.splprep.html)
continuation contracts.
