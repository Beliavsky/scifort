# Genuine `polyroots-fortran` integration validation (2026-10-31)

**Inputs:** `polyroots-fortran-master.zip` uploaded by the SciFort user, SHA-256
`d78040fb46302b5ab72fd629bcb68448b0708343cec68740fae518dcbf5db830`.
The archive contains Jacob Williams's `src/polyroots_module.F90` and upstream
`LICENSE.md`. Its label is **master**, not a confirmed v1.2.0 release snapshot.
No upstream source has been copied into the SciFort project.

**Toolchain:** GNU Fortran 14.2.0, Linux; no FPM installed.

**Compiler:** `gfortran -cpp -O2 -ffree-line-length-none -c
polyroots_module.F90` succeeded. The object directly references LAPACK's
`dgeev` and `zgeev`. Ordinary linking with `-llapack -lblas` succeeded;
linking the adapter **without** BLAS/LAPACK also succeeded with
`-ffunction-sections -fdata-sections` at compile time and
`-Wl,--gc-sections` at link time. These flags eliminate the unused
LAPACK-based root-solving routines for an `rpoly`-only executable.

**Checked build:** `-O0 -g -fcheck=all -fbacktrace`, real upstream
`test_optional_polyroots` passed; 120 independently constructed odd-degree
(5, 9, 13, 17, 21, 25) polynomials with three known real zeros and nonreal
conjugate pairs passed to absolute accuracy better than `1e-6`.

**Optimized build:** `-O3`, both tests passed as well. Both the normal
BLAS/LAPACK link and linker-garbage-collection modes passed.

**Actual upstream benchmark:** the opt-in program ran 100 repetitions at
degrees 5, 9, 13, and 17. Native and upstream callback root counts agreed;
maximum absolute root-location discrepancies respectively about
`3.33e-16`, `3.44e-15`, `8.58e-14`, and `3.33e-16`. Individual elapsed
CPU times were fractions of milliseconds and varied between runs, so they
are not reliable evidence of a large algorithmic speed advantage.

**Strict standard caveat:** the unmodified upstream module fails
`gfortran -std=f2018` because one LAPACK interface uses `COMPLEX*16`.
A temporary source copy with a standards-conforming `complex(real64)`
declaration and a local import of `iso_fortran_env:real64` compiled under
`-std=f2018`. The main upstream source in the user archive was untouched.

**Open items:** Full FPM builds, Windows/MinGW-w64, direct v1.2.0 tag
validation, and stress tests on highly clustered/multiple polynomial roots.
The optional backend remains opt-in; SciFort's native solver is still
the default. The original upstream `fpm.toml` declares BLAS/LAPACK links.

Reproduce from the complete SciFort archive:

```bash
python tools/interpolate_reference/check_optional_polyroots_local.py \
    /path/to/polyroots-fortran-master.zip
python tools/interpolate_reference/check_optional_polyroots_local.py \
    /path/to/polyroots-fortran-master.zip --link-mode lapack
```
