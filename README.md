# SciFort

SciFort is a modern Fortran library for statistics and scientific computing,
with native Fortran APIs and a Fortran Package Manager (FPM) build. It provides
probability distributions, special functions, statistical methods, and numerical
tools for integration, interpolation, optimization, linear algebra, and signals.

The project follows SciPy conventions for shared functionality while keeping
Fortran interfaces. Coverage and numerical guarantees vary by routine; consult
the linked guides for supported options, validation, and limitations. SciFort
is under active development and is not a complete replacement for SciPy. Code
was written primarily by ChatGPT with Extra High thinking effort.

## What is included

| Area | Capabilities and documentation |
| --- | --- |
| Statistics and probability | Continuous and discrete distributions, descriptive statistics, hypothesis tests, random sampling, and fitting. Start with the [API reference](docs/api.md). |
| Special functions | Functions supporting distributions and numerical methods, including gamma, beta, Bessel, and related families. See the [API reference](docs/api.md). |
| Integration | Quadrature including [Lebedev sphere rules](docs/lebedev.md), initial-value solvers, and adaptive real/complex [boundary-value ODEs](docs/bvp.md); see the [examples](example/) and [validation records](validation/). |
| Interpolation | [Splines](docs/interpolate_bspline.md), [regular grids](docs/interpolate_regular_grid.md), [scattered data](docs/interpolate_scattered.md), [radial basis functions](docs/interpolate_rbf_interp1d.md), and [AAA rational approximation](docs/interpolate_aaa.md). |
| Optimization | Root finding, local minimization, least squares, and global optimization. Guides include [linear least squares](docs/optimize_linear_least_squares.md), [differential evolution](docs/optimize_differential_evolution.md), and [SHGO](docs/optimize_shgo.md). |
| Linear algebra | Dense and banded systems, matrix decompositions, matrix functions, and sparse solvers. See the [API reference](docs/api.md) and [banded-system guide](docs/linalg_banded.md). |
| Fourier transforms | 1D, 2D, and multidimensional FFTs and cosine/sine transforms, FFTLog fast Hankel transforms, axis selection, frequency bins, shifts, and fast-length helpers. See the [FFT guide](docs/fft.md). |
| Signals and images | Filtering, spectral analysis, streaming transforms, and multidimensional image operations. See [signal processing](docs/signal.md) and [image filtering and morphology](docs/ndimage.md). |
| Clustering and spatial methods | [Hierarchical clustering](docs/cluster_hierarchy.md), [vector quantization](docs/cluster_vq.md), distances, and geometric methods; see the [source modules](src/) and [examples](example/). |

## Build and run

Install Git, FPM, and a Fortran compiler with Fortran 2018 support. GNU Fortran
is used in the project's CI configuration for Linux, macOS, and Windows.
The first build fetches the QUADPACK, MINPACK, L-BFGS-B, and SLSQP dependencies
declared in [fpm.toml](fpm.toml).

From the repository directory, run these commands in CMD or a Unix shell:

```text
fpm build
fpm test
fpm run --example distributions_demo
```

`fpm test` rebuilds changed sources before running the tests, so a separate
`fpm build` is unnecessary when testing changes. Tests run without keyboard
input. The full suite is extensive and can take several minutes.

For an optimized build, use `fpm build --profile release`. To run the tests
with debug checks explicitly enabled, use `fpm test --profile debug`.
With GNU Fortran, array-temporary runtime warnings can be suppressed while
retaining the other debug checks:

```text
fpm test --profile debug --flag "-fcheck=no-array-temps"
```

Changing compilation flags triggers a separate build. Reference comparisons
driven by Python are documented under [validation](validation/) and are
separate from the native FPM test suite.

## Use from Fortran

Add SciFort to another project's `fpm.toml`:

```toml
[dependencies]
scifort = { git = "https://github.com/Beliavsky/scifort.git" }
```

This follows the default branch. Pin a known commit with `rev` when a
reproducible dependency version is needed.

For example, evaluate a normal CDF for scalars or arrays and compute a quantile:

```fortran
program normal_example
    use scifort_stats, only : dp, normal_cdf, normal_ppf
    implicit none

    real(dp) :: x(3)

    x = [-1.0_dp, 0.0_dp, 1.0_dp]
    print *, normal_cdf(1.96_dp)
    print *, normal_cdf(x)
    print *, normal_ppf(0.975_dp)
end program normal_example
```

Normal, exponential, Weibull, lognormal, gamma, beta, logistic, and both Gumbel distributions also support
[censored likelihoods and fitting](docs/censored_fitting.md).

See [runnable examples](example/) for other APIs and [numerical policy](docs/numerical-policy.md)
for parameter validation, endpoint behavior, precision, and error handling.

## Compatibility and language bindings

For shared distribution functionality, SciFort aims to follow SciPy's names,
parameterizations, parameter order, and mathematical behavior. Intentional
differences and restricted subsets are documented in the API guides.
Additional capabilities are labeled **SciFort extension**.

An experimental C ABI is declared in [include/scifort.h](include/scifort.h).
It provides distribution entry points using interoperable types, explicit
lengths, caller-owned output memory, and integer status codes. Its coverage
is smaller than the native Fortran API. Python, R, MATLAB, and Octave wrappers
are planned; see [language bindings](bindings/README.md).

SciFort is independent of SciPy, NumFOCUS, and the Python Software Foundation.
Development has used OpenAI and Anthropic language models; implementation
sources and adaptations are recorded in [code provenance](CODE_PROVENANCE.md).

## External packages and acknowledgments

SciFort builds on the following external packages. The tables distinguish
dependencies fetched by FPM from implementations adapted into SciFort's source.
Exact versions, source files, modifications, authorship, and retained notices
are recorded in [code provenance](CODE_PROVENANCE.md) and
[third-party licenses](THIRD_PARTY_LICENSES.md).

**FPM dependencies** (declared in [fpm.toml](fpm.toml)):

| Package | Use in SciFort |
| --- | --- |
| [Modern QUADPACK](https://github.com/jacobwilliams/quadpack) | One-dimensional numerical quadrature. |
| [Modern MINPACK](https://github.com/fortran-lang/minpack) | Nonlinear equation and least-squares kernels. |
| [L-BFGS-B](https://github.com/jacobwilliams/lbfgsb), tag `1.0.0` | Bound-constrained minimization. |
| [SLSQP](https://github.com/jacobwilliams/slsqp), tag `1.6.0` | Constrained optimization and least-squares adapters. |

**Adapted implementations, recurrences, and coefficients** included in SciFort:

| Package or source | Use in SciFort |
| --- | --- |
| [SciPy](https://github.com/scipy/scipy) | Selected optimization, differentiation, interpolation, sparse-solver, signal, spatial-transform, special-function, and FFTLog fast Hankel kernels. |
| [Chebfun](https://chebfun.org/) | AAA rational approximation, adapted through SciPy with Oxford/Chebfun notices retained. |
| [perazz/FITPACK](https://github.com/perazz/fitpack) (Paul Dierckx / Federico Perini) | Periodic/pole constraints for fixed-knot spherical splines and the banded solver and adaptive smoothing for [rectangular spherical splines](docs/interpolate_rect_sphere.md). |
| [Riemann_zeta-and-Gamma](https://github.com/Alexey-Kuznetsov-math/Riemann_zeta-and-Gamma) (Alexey Kuznetsov) | Complex zeta Euler-Maclaurin kernel and Bernoulli coefficients, adapted from its BSD-licensed Fortran implementation. |
| [PRIMA](https://github.com/libprima/prima) | COBYLA implementation adapted primarily through SciPy's `pyprima` translation. |
| [ODEPACK](https://github.com/jacobwilliams/odepack) | Selected LSODA operations, coefficients, and control logic. |
| [MATH77](https://netlib.org/math/) | Selected Bessel and sine/cosine-integral kernels. |
| [carlson-elliptic-integrals](https://github.com/jacobwilliams/carlson-elliptic-integrals) | Real Carlson symmetric elliptic-integral kernels. |
| [rklib](https://github.com/jacobwilliams/rklib) | Bogacki-Shampine and Dormand-Prince Runge-Kutta tableaux. |
| [DOP853](https://github.com/jacobwilliams/dop853) | Eighth-order Runge-Kutta coefficients, error estimator, and dense output. |
| [ForDiff](https://github.com/gha3mi/fordiff) | Scalar complex-step differentiation identity. |
| [roots-fortran](https://github.com/jacobwilliams/roots-fortran) | Bisection, Brent, and Ridders scalar root-finding kernels. |
| [fmin](https://github.com/jacobwilliams/fmin) | Bounded Brent scalar minimization. |
| [Fong-Saunders LSMR](THIRD_PARTY_LICENSES.md#lsmr-fong-saunders) | Iterative least squares, adapted through SciPy's LSMR implementation. |

**Optional backends**, built separately from the default project:

| Package | Integration |
| --- | --- |
| [perazz/FITPACK](https://github.com/perazz/fitpack) | Spline smoothing and fitting through [the optional FITPACK adapter](optional/fitpack_backend/). |
| [polyroots-fortran](https://github.com/jacobwilliams/polyroots-fortran) | Polynomial root finding through [the optional root backend](optional/polyroots_backend/). |

The provenance records also credit design, algorithmic, and cross-check
references, including
[nlesolver-fortran](https://github.com/jacobwilliams/nlesolver-fortran),
[kiss](https://github.com/ivan-pi/kiss),
[NumDiff](https://github.com/jacobwilliams/NumDiff),
[TRLIB](https://github.com/felixlen/trlib),
[Spectrum](https://github.com/jchristopherson/spectrum),
[PEAKS](https://github.com/jchristopherson/peaks), and
[Ian Morrall's differential-evolution example](https://github.com/ian-mmm/differential-evolution_f95).
These acknowledgments do not
imply that those packages are linked dependencies or that their source was
copied. SciPy and [mpmath](https://mpmath.org/) also supply independent numerical reference values
for many validation workflows.

## Documentation and development

- [API reference](docs/api.md) and [topic guides](docs/)
- [Architecture](docs/architecture.md) and [roadmap](ROADMAP.md)
- [Validation record](docs/validation.md) and [comparison reports](validation/)
- [Contributing](CONTRIBUTING.md) and [repository rules](AGENTS.md)
- [Changelog](CHANGELOG.md) and [historical README notes](README_HISTORY.md)

The historical notes preserve earlier feature announcements, implementation
details, and validation reports. They describe the project at the time each
entry was written; use the current API guides for supported behavior.

## License and provenance

Original SciFort code uses the [MIT License](LICENSE). Adapted third-party
components retain their own applicable notices and license terms. See
[code provenance](CODE_PROVENANCE.md), [third-party licenses](THIRD_PARTY_LICENSES.md),
and [license texts](LICENSES/) for their sources and obligations.
