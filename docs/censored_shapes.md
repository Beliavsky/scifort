# Censored likelihood array contracts

The C/Python family set is normal, exponential, minimum Weibull, lognormal, gamma,
beta, logistic, Laplace, Cauchy and Rayleigh.
The adapters call the existing native log-likelihood and negative log-likelihood
(`nnlf`) procedures in [censored.f90](../src/censored.f90); their algorithms,
public names and endpoint behavior are unchanged.

## Array shapes

| Procedure | Array constraints | Empty behavior |
| --- | --- | --- |
| Native `make_censored_data` | optional `uncensored(nu)`, `left(nl)`, `right(nr)`, `interval(ni,2)`; independent group lengths | absent or zero-length groups contribute no observations |
| Twenty public procedures in `scifort_censored_bridge` | `uncensored(nu)`, `left(nl)`, `right(nr)`, `lower(ni)`, `upper(ni)`; endpoint vectors have equal lengths | all lengths may be zero; valid parameters give log-likelihood and NNLF zero |
| Private `pack_censored_samples` | the same five vector shapes; total `nu+nl+nr+ni` must fit the native default integer | copied groups, including empty groups, are local to the call |
| Native normal/exponential/logistic/Laplace/Cauchy/Rayleigh density, tail and interval callbacks, plus shared `location_scale_z` and Cauchy/Rayleigh distance/hazard helpers | `params(2)` contains location and positive scale | length must be exactly two |
| Native Weibull PDF, CDF, SF, interval and log-hazard helpers | `params(3)` contains positive shape, location and positive scale | length must be exactly three |
| Native lognormal PDF/CDF/SF, interval, standardization and log-coordinate helpers | `params(3)` contains positive log standard deviation, location and positive scale | length must be exactly three |
| Native gamma PDF/CDF/SF, interval and scaled-coordinate helpers | `params(3)` contains positive shape, location and positive inverse rate | length must be exactly three |
| Native beta parameter, distance/coordinate, density, tail and interval helpers | `params(4)` contains positive first/second shapes, location and positive support width | length must be exactly four |
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
on its support. Lognormal `s` is the standard deviation of `log(x-loc)` and
`scale` is `exp(mean(log(x-loc)))`. Gamma `a` is its shape and `scale` its inverse
rate. Beta `a,b` are its shapes and `scale` is its support width. All shape
parameters are finite and positive. These are the existing native parameterizations.
Logistic uses its standard location/scale convention; Laplace scale is the mean
absolute deviation about `loc`; Cauchy scale is its half width at half maximum.
Rayleigh support starts at `loc` and its cumulative hazard is
`0.5*((x-loc)/scale)**2`. None of these four families takes a shape parameter.
Beta additionally requires a finite rounded upper support edge: when `loc > 0`,
`scale <= huge(1.0_dp)-loc`. Violating this native limit returns status 2 even
for empty observations. No additional upper-edge representability requirement
is introduced by the wrappers.

Impossible observations can give log-likelihood `-inf` and NNLF `+inf`. Infinite
endpoint densities and indeterminate sums retain the native `+inf` or NaN result
with status 0. Status reports call validity, not finiteness of the answer.
Native numerical limitations also remain: for example, unsupported beta
scaled-underflow expansions can yield NaN with status 0. These are not reclassified
as invalid calls. Stable native interval algorithms are reused, including intervals whose rounded
CDF values coincide and far-tail intervals.

## C ABI

Each symbol is named `scifort_<family>_censored_<method>_f64`, with families
`normal`, `exponential`, `weibull`, `lognormal`, `gamma`, `beta`, `logistic`,
`laplace`, `cauchy`, `rayleigh` and methods `loglikelihood`, `nnlf`.
Argument order is:

```text
nu, nl, nr, ni, uncensored, left, right, lower, upper,
[family shape parameters], loc, scale, value, status
```

| Family | Shape parameters in C order | Required Python keywords |
| --- | --- | --- |
| Normal, exponential, logistic, Laplace, Cauchy, Rayleigh | none | none |
| Minimum Weibull | `c` | `c` |
| Lognormal | `s` | `s` |
| Gamma | `a` | `a` |
| Beta | `a, b` | `a`, `b` |

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

`api.censored.<family>_loglikelihood` and `api.censored.<family>_nnlf` return
Python floats for each family listed above.
The first argument is optional `uncensored`; keywords are `left`, `right`,
`interval`, `loc=0.0`, `scale=1.0`, plus the family shape keywords listed above.
Omitted groups are empty. `interval` has shape `(ni,2)` with lower/upper columns;
`[]` is also accepted as an empty interval group. Strided, read-only and overlapping
inputs are accepted through conversion/copying; inputs are not modified. No NumPy
broadcasting of parameters is provided. Invalid calls raise `ValueError` (invalid
numeric types may raise `TypeError`).

The lower-level `api.censored_bridge.<family>_<method>` methods take five input
vectors and return an `ArrayResult` containing `status` and scalar `value`.
