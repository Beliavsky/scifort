# Fortran-code-on-GitHub library review for dense/sparse matrix functions

Reviewed 2026-11-06. The catalog is a list of independent repositories: its
Unlicense license does not override the licenses of the listed source projects.

| Project | Upstream URL | Reported license | Possible use / qualification |
|---|---|---|---|
| NTPoly | https://github.com/william-dawson/NTPoly | MIT | Sparse Hermitian matrix functions using polynomial expansions. Good reuse candidate if MPI, BLAS, LAPACK, and CMake dependencies are acceptable; not drop-in dense `sqrtm`/`signm`. |
| PSFUN | https://github.com/Cirdans-Home/psfun | BSD-3-Clause | Matrix-function vector products, Krylov/quadrature methods. Suitable provenance/license for reuse with notices, but designed for PSBLAS/distributed arrays; upstream cautions that it is under development. |
| Reference LAPACK | https://github.com/Reference-LAPACK/lapack | Modified BSD | Schur, QR and Sylvester building blocks; easiest future route is an optional LAPACK backend rather than vendoring. Preserve notices and attribution. |
| Burkardt Fortran 90 collection | https://github.com/Beliavsky/Burkardt-Fortran-90-codes | Mixed by file: often MIT, some GPL/LGPL | `matrix_exponential` and related routines can be evaluated individually; **check the specific source-file license before copying**. Not approved for blanket copying. |
| MEPACK | https://github.com/mpimd-csc/mepack | GPL-3.0-or-later for main distribution | Useful design/reference for Sylvester/Lyapunov solvers, but unsuitable for copying into the permissively licensed SciFort source distribution without an alternative license. |
| FortranMatrix | https://github.com/KNAnagnostop/FortranMatrix | GPL-3.0 | Matrix wrappers; not a permissive-license source for SciFort incorporation. |

For the new `sqrtm` and `signm` routines, SciFort incorporates **no code from
these external projects**. The algorithms are original Fortran code using an
existing MIT-licensed SciFort complex QR/Schur routine; SciPy and NumPy are
used only to validate numerical outputs. If code from NTPoly, PSFUN, or LAPACK
is added later, retain exact copyright/license notices and record per-file
provenance in `CODE_PROVENANCE.md` and `THIRD_PARTY_LICENSES.md`.
