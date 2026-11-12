# Parametric FITPACK spline validation — 2026-11-01

## Scope

New independent SciFort `scifort_interpolate` features:

- `fitpack_parametric_tck_type` and `splprep` for nonperiodic parametric
  exact interpolation (`s=0`) and full-knot weighted least squares (`task=-1`)
- chord-length parameterization or supplied, strictly monotone parameter values
- degrees 1–5, real coordinate channels, first through third derivatives
- parametric scalar/vector `splev`, including FITPACK extrapolation policies

The implementation reuses native SciFort B-spline machinery and introduces no
third-party source code or new FPM dependencies. Positive-s smoothing,
periodic fitting, and FITPACK task=1 continuation are deliberately unsupported.

## Numerical checks

`tools/interpolate_reference/compare_fitpack_parametric_scipy.py` compares
randomized cases against `scipy.interpolate.splprep` and `splev` (SciPy 1.17.0):

| Test | Result |
|---|---|
| Main seed, checked GNU Fortran | 90/90 PASS |
| Main seed, GNU Fortran `-O3` | 90/90 PASS |
| Independent seed 8675309, checked | 90/90 PASS |
| Independent seed 8675309, optimized | 90/90 PASS |
| Maximum normalized numerical difference in the independent-seed run | `5.42e-12`, using `abs(a-b)/(1+abs(b))` |
| Maximum absolute difference in that run | `1.58e-7`, in high-order derivatives of large magnitude |
| Previous scalar FITPACK TCK tests | 45/45 PASS per checked/optimized build |
| Public interpolation native test suite, checked | 26/26 PASS |
| New parametric native test, optimized focused build | PASS |

The full optimized 26-program native suite **was not completed**: compiling
its whole dependency closure at `-O3` exceeded the execution limit. Likewise,
FPM, Windows/MinGW-w64, and the real optional `perazz/fitpack` dependency
have not been validated in this environment. The optional adapter's earlier
mock-interface test is not an upstream numerical test.

## Reproduction

```bash
python tools/interpolate_reference/compare_fitpack_parametric_scipy.py
python tools/interpolate_reference/compare_fitpack_parametric_scipy.py --seed 8675309
python tools/interpolate_reference/compare_fitpack_tck_scipy.py
python tools/interpolate_reference/run_all_native_interpolation.py --only checked
```

## Code and licensing

Original MIT-licensed SciFort code. No FITPACK source is included or copied.
The optional upstream project is `https://github.com/perazz/fitpack` under
BSD-3-Clause; its license notice is kept under `LICENSES/` for reference.
