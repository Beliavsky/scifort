# SciFort starter

SciFort is a proposed modern Fortran statistics library with an idiomatic
Fortran API, an FPM build, and a stable C ABI for future Python, R, MATLAB, and
Octave bindings.

This repository is an independent project. It is not affiliated with or
endorsed by SciPy, NumFOCUS, or the Python Software Foundation.

## Current scope

Version 0.1.0 contains dependency-free implementations of six continuous
probability distributions:

- normal
- uniform
- exponential
- Laplace
- logistic
- Cauchy

The unreleased development version adds the gamma and chi-square
distributions and, in `scifort_special`, the regularized incomplete gamma
functions `gammainc` and `gammaincc` with their inverses.

Each distribution currently provides:

- `pdf`
- `logpdf`
- `cdf`
- `sf`
- `logcdf`
- `logsf`
- `ppf`
- `isf`

All scalar functions are `pure elemental`, so the same routines operate on
scalars and conformable arrays. Invalid location or scale parameters return a
quiet NaN. Probability arguments outside `[0, 1]` also return a quiet NaN.

The normal CDF uses the standard `erfc` identity. The initial normal quantile
implementation uses safeguarded bisection against the CDF or log-CDF. It is a
correctness-first implementation intended to be replaced or supplemented by
a faster, carefully licensed approximation after accuracy and provenance
reviews.

An experimental C ABI exposes scalar and vector normal PDF, CDF, and PPF
routines. This provides the first vertical slice for future language bindings.
The Fortran API is the primary API in this starter release.

## Requirements

- A Fortran 2018 compiler
- FPM 0.12.0 or newer

The default FPM build produces a static library. FPM 0.13.0 or newer is needed
if the manifest is changed to build static and shared library targets together.

## Build and test

```text
fpm build
fpm test
fpm run
fpm run --example distributions_demo
```

For a debug build with current FPM:

```text
fpm test --profile debug
```

To clone the committed repository into a temporary directory and verify a
clean build and test run on Windows:

```text
test_clone.bat
test_clone.bat https://github.com/Beliavsky/scifort.git
```

With no argument, the script clones the local repository containing it. Local
working-tree changes and untracked files are not included in that clone.

## Fortran use

```fortran
program example
    use scifort_stats, only : dp, normal_cdf, normal_ppf
    implicit none

    real(dp) :: p
    real(dp) :: x(3)

    x = [-1.0_dp, 0.0_dp, 1.0_dp]
    p = normal_cdf(1.96_dp)

    print *, p
    print *, normal_cdf(x)
    print *, normal_ppf(0.975_dp)
end program example
```

To use SciFort from another FPM project:

```toml
[dependencies]
scifort = { git = "https://github.com/Beliavsky/scifort.git" }
```

This tracks the default branch. After a release is tagged, add
`tag = "vX.Y.Z"` to pin that version.

## Repository rules

Coding agents and contributors must read `AGENTS.md`. Claude Code must also
read `CLAUDE.md`. The most important rules are:

- preserve exact code provenance and original license notices
- do not import GPL or LGPL implementation code into the permissive core
- never claim tests were run when they were not
- prioritize numerical correctness, especially in tails and at boundaries
- keep the native Fortran API separate from the stable C ABI
- do not break public APIs or C symbols without an explicit versioning decision

See `CONTRIBUTING.md`, `CODE_PROVENANCE.md`, and the files under `docs/`.

## Validation

The starter has been compiled and tested directly with GNU Fortran 14.2.0,
including a shared-library C smoke test and selected numerical comparison with
SciPy. See `docs/validation.md` for exact checks and the explicit limitation
that FPM itself was not available in the archive-creation environment.

## C ABI

The experimental header is `include/scifort.h`. The initial ABI uses only
standard C scalar types, explicit array lengths, caller-owned output storage,
and integer status codes. It does not expose compiler-specific Fortran module
symbols or descriptors.

The ABI is intentionally small in version 0.1.0. Python, R, MATLAB, and Octave
wrappers should all call this common ABI rather than independently binding to
compiler-specific Fortran interfaces.

## License

Original SciFort code is licensed under the MIT License. Third-party code may
be included only when its license is compatible with distribution in this
project, and every imported or translated component must retain its applicable
notices. See `THIRD_PARTY_LICENSES.md` and `CODE_PROVENANCE.md`.
