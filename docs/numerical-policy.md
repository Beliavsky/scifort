# Numerical policy

## Precision

The initial public API uses binary64 through `dp = kind(1.0d0)`. The C ABI uses
`double`. Support for additional real kinds requires a separate design review.

## Invalid inputs

Scalar functions return a quiet NaN for invalid distribution parameters or
probabilities outside `[0, 1]`. Valid support endpoints and probability
endpoints return their mathematically appropriate finite or infinite values.

## Tail calculations

Survival functions and logarithmic probabilities must use stable formulas.
Implementations should avoid cancellation such as:

```text
sf = 1 - cdf
logsf = log(1 - cdf)
```

unless a proof or regime check shows that the calculation is safe.

## Quantiles

Quantile implementations must be monotone and must define behavior at zero and
one. Fast approximations require tests across central, moderate-tail, and
extreme-tail regimes. Safeguarded iteration is preferred over unbracketed
iteration when convergence is not otherwise guaranteed.

## Floating-point exceptions

Expected underflow in extreme tails is acceptable when the mathematically
rounded binary64 result is zero. Library code must not deliberately trigger a
fatal exception. Tests should identify unexpected overflow, invalid, and
divide-by-zero behavior.

## Reference validation

Use several forms of evidence where practical:

- high-precision reference values;
- independent implementations;
- published tables;
- mathematical identities;
- SciPy or another mature library at a pinned version;
- property tests.

Reference data must record how and with which version it was generated.
