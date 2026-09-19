# Architecture

## Layers

SciFort separates four concerns:

1. **Numerical core**: pure modern Fortran algorithms.
2. **Native API**: idiomatic module procedures for Fortran users.
3. **C ABI**: a small interoperable and versioned procedural interface.
4. **Language wrappers**: Python, R, MATLAB, and Octave conventions built over
   the common C ABI.

The numerical core must not depend on Python, R, MATLAB, Octave, Cython, MEX,
or compiler-specific Fortran descriptors.

## Source layout

```text
src/
    kinds.f90
    constants.f90
    math.f90
    special.f90
    special/
        <function family>.f90
    stats.f90
    stats/
        <distribution>.f90
    c_api/
        api.f90
```

## Native API

Distribution procedure names use a family prefix, for example:

```text
normal_pdf
normal_logpdf
normal_cdf
normal_sf
normal_logcdf
normal_logsf
normal_ppf
normal_isf
```

Scalar arguments are elemental. Array reductions, fitting routines, random
sampling, and algorithms requiring workspaces should use separate procedural
interfaces rather than overloading scalar distribution evaluation.

## C ABI

The C ABI uses:

- `double` for initial public real values;
- `size_t` for array lengths;
- caller-owned arrays;
- `int` status outputs;
- stable lower-case symbol names beginning with `scifort_`.

Future opaque state, workspace, fitted-model, or RNG objects must be represented
by `c_ptr` handles with explicit create and destroy functions.

## Wrapper policy

Language wrappers own language-specific behavior:

- Python: NumPy conversion, broadcasting, exceptions, wheel packaging.
- R: vector recycling decisions, `NA` policy, registration, package metadata.
- MATLAB: MEX arrays, column-major dimension handling, error translation.
- Octave: oct-file or compatible MEX interface and package metadata.

Wrappers must not independently reimplement numerical formulas.
