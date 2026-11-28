# Censored likelihood array contracts

The initial C/Python family set is normal, exponential and minimum Weibull.
The adapters call the existing native log-likelihood and negative log-likelihood
(`nnlf`) procedures in [censored.f90](../src/censored.f90); their algorithms,
public names and endpoint behavior are unchanged.

## Array shapes

| Procedure | Array constraints | Empty behavior |
| --- | --- | --- |
| Native `make_censored_data` | optional `uncensored(nu)`, `left(nl)`, `right(nr)`, `interval(ni,2)`; independent group lengths | absent or zero-length groups contribute no observations |
| Six public procedures in `scifort_censored_bridge` | `uncensored(nu)`, `left(nl)`, `right(nr)`, `lower(ni)`, `upper(ni)`; endpoint vectors have equal lengths | all lengths may be zero; valid parameters give log-likelihood and NNLF zero |
| Private `pack_censored_samples` | the same five vector shapes; total `nu+nl+nr+ni` must fit the native default integer | copied groups, including empty groups, are local to the call |
| Native normal/exponential PDF, CDF, SF and interval callbacks | `params(2)` contains location and positive scale | length must be exactly two |
| Native Weibull PDF, CDF, SF, interval and log-hazard helpers | `params(3)` contains positive shape, location and positive scale | length must be exactly three |
| Generic native `censored_loglikelihood` and `tail_interval` | `params(:)` has the length required by the supplied callbacks | no universal parameter-vector length |

See the [generated bridge contracts](censored_bridge_shapes.md) and
[bridge source](../src/c_api/censored_bridge.f90). The native constructor copies
its inputs. No pointers, observations or derived-type objects are retained by the
new stateless adapters.

## Data and numerical conventions

Exact observations and left/right censoring limits must be finite. Interval
endpoints must be ordered and cannot contain NaNs. Equal finite endpoints
represent an exact observation, using a density rather than zero interval mass.
One endpoint may be infinite: `(-inf,b)` means left censoring and `(a,+inf)` means
right censoring. Both endpoints infinite are invalid. Observations below a
family's support are valid data and retain native zero-density/probability rules.

Location is finite; scale is finite and strictly positive. Weibull shape `c` is
finite and strictly positive. Normal scale is its standard deviation; exponential
scale is its inverse rate. Weibull uses `((x-loc)/scale)**c` as cumulative hazard
on its support. These are the existing native parameterizations.

Impossible observations can give log-likelihood `-inf` and NNLF `+inf`. Infinite
endpoint densities and indeterminate sums retain the native `+inf` or NaN result
with status 0. Status reports call validity, not finiteness of the answer.
Stable native interval algorithms are reused, including intervals whose rounded
CDF values coincide and far-tail intervals.

## C ABI

Each symbol is named `scifort_<family>_censored_<method>_f64`, with families
`normal`, `exponential`, `weibull` and methods `loglikelihood`, `nnlf`.
Argument order is:

```text
nu, nl, nr, ni, uncensored, left, right, lower, upper,
[c for Weibull], loc, scale, value, status
```

Counts are nonnegative C `int` values; scalar parameters are `double` values.
Inputs are explicit-length contiguous binary64 vectors. `value` and `status` are
caller-owned output pointers. All pointers must be valid, including zero-length
inputs; buffers are disjoint and output pointers non-null. No null-pointer checks
are supplied. The C ABI uses separate endpoint vectors, avoiding a matrix layout
convention for C callers. The generated Fortran wrappers use explicit shapes;
native numerical procedures retain assumed shapes and derived types.

| Status | Meaning and outputs |
| --- | --- |
| 0 | valid call evaluated; result may be finite, infinite or NaN |
| 1 | negative C extent, or native count/shape error; negative C extents preserve `value`; native errors set it to NaN |
| 2 | invalid distribution parameters or censored data; `value` is NaN |

The C caller supplies the same `ni` for both endpoint buffers; insufficient
storage cannot be detected. Native adapter calls additionally check equal endpoint
lengths. Python checks vector ranks and matching shapes before C entry.

## Python

`api.censored.normal_loglikelihood`, `normal_nnlf`, `exponential_loglikelihood`,
`exponential_nnlf`, `weibull_loglikelihood` and `weibull_nnlf` return Python floats.
The first argument is optional `uncensored`; keywords are `left`, `right`,
`interval`, `loc=0.0`, `scale=1.0`, plus required keyword `c` for Weibull.
Omitted groups are empty. `interval` has shape `(ni,2)` with lower/upper columns;
`[]` is also accepted as an empty interval group. Strided, read-only and overlapping
inputs are accepted through conversion/copying; inputs are not modified. No NumPy
broadcasting of parameters is provided. Invalid calls raise `ValueError` (invalid
numeric types may raise `TypeError`).

The lower-level `api.censored_bridge.<family>_<method>` methods take five input
vectors and return an `ArrayResult` containing `status` and scalar `value`.
