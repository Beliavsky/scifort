# Third-party licenses

Most SciFort numerical source is original MIT-licensed implementation.
`CODE_PROVENANCE.md` identifies the small amount of adapted third-party code.
Standard mathematical formulas, language intrinsics, and numerical identities
do not by themselves create a third-party source-code dependency.

## AAA arrowhead spectra (Reference-LAPACK, BSD-3-Clause)

`src/interpolate_aaa_spectrum.f90` adapts Reference-LAPACK v3.12.1's
`ZHGEQZ`, `ZLARTG`, `ZROT`, `ZSCAL`, `ZLADIV`, `DLADIV`/`DLADIV1`/`DLADIV2`,
and safe-scaling constants for an internal normalized AAA pencil backend.
Applicable University of Tennessee / Tennessee Research Foundation,
University of California Berkeley, and University of Colorado Denver
copyright notices, along with upstream NAG and ZLARTG author attribution,
are retained. Full terms are in
[`LICENSES/lapack-BSD-3-Clause.txt`](LICENSES/lapack-BSD-3-Clause.txt).
The exact source versions/hashes and modifications are recorded in provenance.
This is a source adaptation with no mandatory linked LAPACK/BLAS dependency.

## Curve smoothing and spherical spline solvers (perazz/FITPACK, BSD-3-Clause)

`src/interpolate_curfit_core.f90` adapts the scalar `fpcurf` automatic-knot
smoothing solver and its numerical helpers from the pinned perazz/FITPACK
source recorded in provenance. It also adapts the parametric `fppara` solver
and `fp_rotate_row_vec` / `fp_rotate_shifted_vec`, reusing those attributed
helpers. The existing
[`LICENSES/perazz-fitpack-BSD-3-Clause.txt`](LICENSES/perazz-fitpack-BSD-3-Clause.txt)
retains the full BSD terms and Federico Perini copyright notice. Original
SciFort normalization and interface code uses MIT; the combined file is
marked `MIT AND BSD-3-Clause`. No mandatory linked FITPACK package is added.
The same file also adapts periodic `fpperi`, `fpperi_reset_interp`, `fpbacp`,
and `fp_rotate_row_2mat`, covered by the same retained BSD license and notices.

`src/interpolate_sphere.f90` adapts the fixed-knot periodic coefficient and pole
constraint construction from `fpsphe` and `fprpsp` in Federico Perini's
[perazz/FITPACK](https://github.com/perazz/fitpack), based on Paul Dierckx's FITPACK.
The inspected commit is `1907c7bad9f5c4557970e906bdd4909eb52b2d4c`,
`src/fitpack_core.F90`. Copyright (c) 2022, Federico Perini. All rights reserved.
The adapted source retains BSD-3-Clause and the upstream attribution; full terms
are already retained in
[`LICENSES/perazz-fitpack-BSD-3-Clause.txt`](LICENSES/perazz-fitpack-BSD-3-Clause.txt).
`src/interpolate_sphere_grid_core.f90` additionally adapts the interpolation
and adaptive smoothing branches of `spgrid`/`fpspgr`, their knot-selection
and root-iteration helpers, and their banded numerical helpers from the same
pinned source. It retains the same copyright and BSD-3-Clause terms.
`src/interpolate_sphere_scattered_core.f90` also adapts the scattered-data
adaptive solver, rank reduction, coefficient expansion, observation ordering,
and rotation helpers from that source. It retains the same notices and terms.
SciFort uses its existing B-spline evaluation and, for scattered fitting,
dense SVD solver. No full FITPACK source distribution or mandatory package
dependency is added. Exact imported routines and modifications are recorded
in `CODE_PROVENANCE.md`.

## Complex Riemann zeta (Alexey Kuznetsov, BSD-3-Clause)

`src/special/complex_zeta.f90` adapts the Euler-Maclaurin kernel and Bernoulli
coefficients from Alexey Kuznetsov's
[Riemann_zeta-and-Gamma](https://github.com/Alexey-Kuznetsov-math/Riemann_zeta-and-Gamma),
commit `1f1ad8934044cfb20d6d70137868ddd0c74652c4`,
`Fortran/gamma_zeta_module.f90`. Copyright (c) 2025, Alexey-Kuznetsov-math.
The adapted file is BSD-3-Clause licensed; complete terms are retained in
[`LICENSES/kuznetsov-zeta-BSD-3-Clause.txt`](LICENSES/kuznetsov-zeta-BSD-3-Clause.txt).
Exact sources, hashes, and modifications are recorded in `CODE_PROVENANCE.md`.

## AAA rational approximation (Oxford/Chebfun and SciPy BSD-3-Clause)

`src/interpolate_aaa.f90` adapts SciPy v1.18.1's
[`scipy/interpolate/_bary_rational.py`](https://github.com/scipy/scipy/blob/v1.18.1/scipy/interpolate/_bary_rational.py),
which retains the 2017 University of Oxford/Chebfun copyright and license.
Both the full Oxford/Chebfun notice and the SciPy project BSD-3-Clause notice
are retained in
[`LICENSES/scipy-aaa-BSD-3-Clause.txt`](LICENSES/scipy-aaa-BSD-3-Clause.txt).
Original SciFort modifications are MIT licensed. Exact sources, modifications,
mathematical references, and the Fortran-source search are recorded in
`CODE_PROVENANCE.md`.

## Lebedev quadrature (getLebedevSphere / SciPy BSD-2-Clause)

`src/integrate/lebedev.f90` adapts SciPy v1.18.1's
[`scipy/integrate/_lebedev.py`](https://github.com/scipy/scipy/blob/v1.18.1/scipy/integrate/_lebedev.py),
Copyright (c) 2010, Robert Parrish. Complete terms are retained in
[`LICENSES/scipy-lebedev-BSD-2-Clause.txt`](LICENSES/scipy-lebedev-BSD-2-Clause.txt).
The source retains the Parrish / van Wuellen / Laikov attribution chain and
requested Lebedev-Laikov citation. Original SciFort modifications are MIT
licensed. Exact source, significant changes, and the unresolved direct CCL
Fortran import are documented in `CODE_PROVENANCE.md`.

## Boundary-value ODE solver (SciPy BSD-3-Clause)

`src/integrate/bvp.f90` adapts SciPy v1.15.3's
[`scipy/integrate/_bvp.py`](https://github.com/scipy/scipy/blob/v1.15.3/scipy/integrate/_bvp.py).
Full required notices and terms are retained in
[`LICENSES/scipy-bvp-BSD-3-Clause.txt`](LICENSES/scipy-bvp-BSD-3-Clause.txt).
Original SciFort modifications are MIT licensed. Exact source version,
modifications, mathematical references, and the reviewed Fortran alternative
are recorded in `CODE_PROVENANCE.md`.

## Fast Hankel transforms / FFTLog (SciPy BSD-3-Clause)

`src/fftlog.f90` and the coefficient-magnitude check in
`validation/check_fftlog.py` adapt SciPy v1.15.3's
[`scipy/fft/_fftlog_backend.py`](https://github.com/scipy/scipy/blob/v1.15.3/scipy/fft/_fftlog_backend.py).
The complete upstream BSD-3-Clause copyright, redistribution terms, and
disclaimer are retained in
[`LICENSES/scipy-fftlog-BSD-3-Clause.txt`](LICENSES/scipy-fftlog-BSD-3-Clause.txt).
Original SciFort modifications are MIT licensed. Attribution to Hamilton's
FFTLog method and the decision not to import unclear-license original
Fortran sources are recorded in `CODE_PROVENANCE.md`.

## MINRES symmetric Krylov recurrence (SciPy BSD-3-Clause)

The Paige--Saunders Lanczos/Givens recursion in
`src/sparse_linalg.f90` was adapted from the BSD-3-Clause SciPy 1.17
implementation at:
https://github.com/scipy/scipy/blob/v1.17.0/scipy/sparse/linalg/_isolve/minres.py

The SciPy copyright statement, three redistribution conditions, and full
warranty disclaimer are reproduced in
`LICENSES/scipy-minres-BSD-3-Clause.txt`. SciFort rewrites the numerical
recurrence in modern Fortran with one-based CSR storage, actual shifted
residual tests and its own validation/status semantics. The BSD-3-Clause
license is compatible with SciFort's MIT-licensed project.

## QMR sparse Krylov iteration (SciPy BSD-3-Clause)

The two-sided Lanczos, quasi-minimization, and solution-update recurrence
in `scifort_sparse_linalg.qmr` was adapted from SciPy v1.17.0's
`scipy/sparse/linalg/_isolve/iterative.py` (`qmr`):
https://github.com/scipy/scipy/blob/v1.17.0/scipy/sparse/linalg/_isolve/iterative.py

Retained SciPy BSD-3-Clause attribution and full license terms are recorded in
`LICENSES/scipy-qmr-BSD-3-Clause.txt`. SciFort provides its own canonical
CSR matrix/transpose, right-Jacobi preconditioner and verified true-residual
stopping logic. The newly added complex `qmr_complex` implementation in
`src/sparse_linalg_complex.f90` extends SciFort's real QMR recurrence
with conjugate-transpose shadow vectors and left preconditioning; it
retains the same BSD-3-Clause attribution. The package has no linked
SciPy dependency.

## TFQMR sparse Krylov iteration (SciPy BSD-3-Clause)

The right-preconditioned TFQMR iteration in `src/sparse_linalg.f90`
adapts SciFort's existing BSD-3-Clause matrix-free TFQMR implementation
in `src/optimize/krylov_iterative.f90`. That implementation follows the
SciPy v1.17.0 `scipy/sparse/linalg/_isolve/tfqmr.py` recurrence:
https://github.com/scipy/scipy/blob/v1.17.0/scipy/sparse/linalg/_isolve/tfqmr.py

The full upstream SciPy BSD-3-Clause license, copyright notice and disclaimer
are retained in `LICENSES/scipy-tfqmr-BSD-3-Clause.txt`. SciFort adds CSR
storage, consistent `A*M` right preconditioning, and true residual checks.

## ODEPACK

Selected LSODA persistent Nordsieck rescaling/prediction/correction operations, Adams
correction and error/order-selection coefficients, functional-corrector convergence logic,
switching constants, and control behavior in
`src/integrate/lsoda.f90` are adapted from the public-domain ODEPACK `DSTODA` and
`DCFODE` routines as modernized by Jacob Williams. Sources and notice inspected
2026-09-29:

- https://raw.githubusercontent.com/jacobwilliams/odepack/master/src/M_da1/dstoda.inc
- https://raw.githubusercontent.com/jacobwilliams/odepack/master/src/M_da1/dcfode.inc
- https://raw.githubusercontent.com/jacobwilliams/odepack/master/LICENSE

NOTICE

This work was produced at the University of California, Lawrence Livermore
National Laboratory (UC LLNL) under contract no. W-7405-ENG-48 (Contract
48) between the U.S. Department of Energy (DOE) and The Regents of the
University of California (University) for the operation of UC LLNL. The
rights of the Federal Government are reserved under Contract 48 subject
to the restrictions agreed upon by the DOE and University as allowed
under DOE Acquisition Letter 97-1.
The ODEPACK package has been declared to be in the Public Domain.

DISCLAIMER
This work was prepared as an account of work sponsored by an agency of
the United States Government. Neither the United States Government nor the
University of California nor any of their employees, makes any warranty,
express or implied, or assumes any liability or responsibility for the
accuracy, completeness, or usefulness of any information, apparatus,
product, or process disclosed, or represents that its use would not
infringe privately-owned rights. Reference herein to any specific
commercial products, process, or service by trade name, trademark,
manufacturer or otherwise does not necessarily constitute or imply
its endorsement, recommendation, or favoring by the United States
Government or the University of California. The views and opinions of
authors expressed herein do not necessarily state or reflect those of
the United States Government or the University of California, and shall
not be used for advertising or product endorsement purposes.

## Independently implemented Gaussian-spline and Kaiser formulas

`gauss_spline`, `kaiser_beta`, `kaiser_atten`, `kaiserord`, and the optional
parameterized `firwin` path were written independently from the published
Gaussian B-spline and Kaiser FIR design equations. No upstream SciPy source
is included or adapted for these helpers; SciPy 1.17 (BSD-3-Clause) is used
solely as an independent numerical reference. The new Fortran implementation
remains MIT-licensed. See CODE_PROVENANCE.md for bibliographic credit.

## SciPy cubic and quadratic spline-filter recurrences

`src/signal_spline.f90` adapts the recursive cubic/quadratic B-spline
coefficient construction, cubic smoothing recurrence, and mirror-boundary
evaluation conventions in SciPy 1.17.0's `scipy.signal._spline_filters`.
SciPy is distributed under BSD-3-Clause; its copyright attribution and
complete license notice appear in the SciPy section below. The Fortran
implementation is built into SciFort and does not link to Python or SciPy.

## Additional SciPy signal filter formulas

`src/signal_iir_extra.f90` now includes `gammatone` and `iircomb`
coefficient calculations adapted from SciPy 1.17.0
`scipy.signal._filter_design`. These follow the published formulas of
Malcolm Slaney (1993) and Sophocles Orfanidis (1996). The applicable
SciPy BSD-3-Clause copyright and redistribution terms appear under the
SciPy notice later in this file. The Fortran implementation does not link
against or redistribute SciPy.

## PRIMA

The COBYLA numerical kernel in `src/optimize/cobyla_prima.f90` is a
modern-Fortran adaptation of the PRIMA COBYLA design based primarily on SciPy
1.17.0's faithful `pyprima` translation of the modern PRIMA reference
implementation. PRIMA is available at <https://github.com/libprima/prima> and
<https://www.libprima.net>. The upstream license location is
<https://github.com/libprima/prima/blob/main/LICENCE.txt>. Notice inspected
2026-10-06:

Copyright (c) 2020--2025, Zaikun ZHANG (https://www.zhangzk.net)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## RotationSpline formulation

`src/spatial_rotation_spline.f90` implements the smooth attitude interpolation
formulation also used by `scipy.spatial.transform.RotationSpline`
(`scipy/spatial/transform/_rotation_spline.py`, SciPy v1.17.0). The SO(3)
Jacobian formulas, nonlinear knot-rate continuity equations and iteration
were informed by this SciPy BSD-3-Clause implementation and its published
mathematical reference. SciPy's existing BSD-3-Clause copyright, conditions
and disclaimer in the SciPy section below are retained for that attribution.
The Fortran code is self-contained and introduces no linked SciPy dependency.

## SciPy

The cyclic knot placement and endpoint derivative collocation scheme for
`make_interp_spline(...,bc_type='periodic')` in
`src/interpolate_bspline.f90` were adapted from SciPy 1.17.0
`scipy/interpolate/_bsplines.py` under the BSD-3-Clause license.
The SciPy copyright and conditions reproduced below apply to this
adaptation. The rest of the existing SciFort spline engine uses its own
independent implementation; SciPy code is not linked or vendored.

The new `chebwin`, `taylor`, and `kaiser_bessel_derived` code in
`src/signal_windows_extra.f90` adapts the corresponding SciPy 1.17.0
`scipy/signal/windows/_windows.py` formulas and synthesis steps, under
SciPy's BSD-3-Clause terms reproduced below. SciFort uses its own direct
DFT and log-scaled numeric evaluation; no upstream Python files are bundled.


`src/signal_cwt_peaks.f90` adapts the default Ricker-wavelet convolution,
ridge linking, percentile-based noise estimation, and peak-selection conventions
from SciPy 1.17.0 `scipy.signal._peak_finding` and
`scipy.signal._wavelets`. This implementation was written in modern Fortran;
no SciPy source files or SciPy runtime dependencies are bundled.
SciPy copyright and BSD-3-Clause terms appear in this section below.


The adaptive finite-difference stencil/weight/error/termination controller in
`src/differentiate.f90` is adapted from SciPy 1.17.0
`scipy/differentiate/_differentiate.py`.
The complex-valued CDF antiderivative in `src/stats/rel_breitwigner.f90` is
adapted from `scipy.stats.rel_breitwigner` in SciPy 1.17.0. The regime
selection, Durbin/Marsaglia-Tsang-Wang matrix calculation, and Pelz-Good
expansion in `src/stats/kstwo.f90` are adapted from `scipy.stats._ksstats` in
SciPy 1.17.0. The S1-to-S0 conversion, difficult-input handling, special cases,
and Nolan/Zolotarev piecewise formulas in `src/stats/levy_stable.f90` are
adapted from SciPy 1.17.0 `scipy.stats._levy_stable` and its `levyst.c` helper.
The conditional multivariate-normal/Student-t box transform and pivoted/permuted
Cholesky logic in `src/stats/multivariate_integration.f90` are adapted from
SciPy 1.17.0 `scipy/stats/_qmvnt.py` and `scipy/stats/_qmvnt_cy.pyx`; SciFort's
randomized tent-transformed Halton point generator is independently implemented.
The Shapiro-Wilk numerical kernel, finite-sample Cramer-von Mises correction,
and Anderson-Darling critical-value/public fitting conventions in
`src/goodness_of_fit.f90` are adapted from SciPy 1.17.0
`scipy/stats/_ansari_swilk_statistics.pyx`, `_hypotests.py`, and `_morestats.py`.
The support clipping, normalization, tail, and inverse-tail relationships in
`src/distribution_infrastructure.f90` are adapted from SciPy 1.17.0
`scipy/stats/_distribution_infrastructure.py` for the public
`make_distribution`/`truncate` compatibility layer.
The exact untied Kendall inversion-count recurrence in
`src/association_extended.f90` is adapted from SciPy 1.17.0
`scipy/stats/_mstats_basic.py`; the surrounding association/regression/trend
implementations are independently written.
The sequential hypergeometric sampling control flow in
`src/stats/multivariate_hypergeom.f90` and the inverse-Wishart/matrix-normal
sampling decomposition in `src/stats/matrix_t.f90` are also translated/adapted
from SciPy 1.17.0 `scipy/stats/_multivariate.py`; their probability formulas and
surrounding Fortran implementations are independently written. The Joe-Kuo Sobol
direction-number initialization data in `src/stats/qmc_sobol_data.f90` are generated
from SciPy 1.17.0's distributed `_sobol_direction_numbers.npz`. The `random-cd`
proposal loop and centered-discrepancy exchange update in `src/stats/qmc.f90`, and
the Lloyd relaxation/stopping semantics in `src/stats/qmc_lloyd.f90`, are adapted
from SciPy 1.17.0 `scipy/stats/_qmc.py`. The native Bowyer-Watson/Voronoi geometry
in `qmc_lloyd.f90` is original SciFort code and does not copy Qhull. The complex `log_softmax` max-shift branch convention in
`src/special/reductions_low.f90` is adapted from SciPy 1.17.0
`scipy/special/_logsumexp.py`. The irregular
Simpson formulas and built-in Newton-Cotes coefficient table in
`src/integrate/quadrature.f90` are adapted from SciPy 1.17.0
`scipy/integrate/_quadrature.py`. The RK45 quartic dense-output coefficient
matrix and polynomial convention in `src/integrate/ivp.f90` are adapted from
SciPy 1.17.0 `scipy/integrate/_ivp/rk.py`. The fifth-order Radau IIA stepper,
Newton convergence logic, cubic dense output, and dense finite-difference
Jacobian strategy in `src/integrate/radau.f90` are adapted from SciPy 1.17.0
`scipy/integrate/_ivp/radau.py` and `_ivp/common.py`. The variable-order BDF/NDF
stepper, Newton convergence logic, dense output, and dense finite-difference
Jacobian strategy in `src/integrate/bdf.f90` are adapted from SciPy 1.17.0
`scipy/integrate/_ivp/bdf.py` and `_ivp/common.py`. The Genz-Malik
embedded-rule formulas,
Gauss-Kronrod 15/21 tables, and global adaptive subdivision strategy in
`src/integrate/cubature.f90` are adapted from SciPy 1.17.0
`scipy/integrate/_rules/_genz_malik.py`, `_rules/_gauss_kronrod.py`, and
`_cubature.py`. The nonlinear solve controller, Broyden/Anderson/mixing updates,
matrix-free Krylov behavior, Broyden rank reduction, Armijo search, and
DCSRCH/DCSTEP strong-Wolfe search in `src/optimize/nonlin.f90` and
`src/optimize/wolfe.f90` are adapted from SciPy 1.17.0
`scipy/optimize/_nonlin.py`, `_linesearch.py`, `_dcsrch.py`, and `_root.py`.
The Nelder-Mead method including bound clipping/reflection, modified Powell
including bounded/one-sided line searches, BFGS, Polak-Ribiere+ conjugate-gradient,
Newton-CG inner iteration, finite-difference gradient/Hessian-product, dogleg,
trust-ncg, trust-exact, and trust-krylov trust-region controllers and
subproblems, and second strong-Wolfe/zoom logic in `src/optimize/minimize.f90` and `src/optimize/wolfe.f90`
are adapted from SciPy 1.17.0 `scipy/optimize/_optimize.py`, `_minimize.py`,
`_linesearch.py`, `_trustregion.py`, `_trustregion_dogleg.py`, `_trustregion_ncg.py`,
`_trustregion_exact.py`, `_trustregion_krylov.py`, and
`_trustregion_constr/equality_constrained_sqp.py`. The PRIMA-derived COBYLA
kernel in `src/optimize/cobyla_prima.f90` also uses SciPy 1.17.0
`scipy/_lib/pyprima/cobyla` as a translation/reference source; the underlying
PRIMA notice is reproduced separately above.
The dense TRF/dogbox nonlinear least-squares controllers, trust-region helpers,
robust-loss scaling, LinearOperator-style matrix-free TRF Jacobian path, and
`curve_fit` trust-region semantics in `src/optimize/least_squares.f90` are adapted from SciPy 1.17.0
`scipy/optimize/_lsq/trf.py`, `_lsq/dogbox.py`, `_lsq/common.py`,
`_lsq/least_squares.py`, and `_minpack_py.py`.
The applicable SciPy BSD-3-Clause license is reproduced below.

Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above
   copyright notice, this list of conditions and the following
   disclaimer in the documentation and/or other materials provided
   with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


## Modern MINPACK

SciFort depends through FPM on `https://github.com/fortran-lang/minpack.git`
for the multidimensional `root`/`fsolve` MINPACK kernels. The dependency is
not vendored in this repository. Upstream `LICENSE.txt` contains the permissive
modern MINPACK contributor license together with the historical University of
Chicago/Argonne MINPACK notice. See:

- https://github.com/fortran-lang/minpack/blob/main/LICENSE.txt
- https://github.com/fortran-lang/minpack

## L-BFGS-B 1.0.0

SciFort depends through FPM on Jacob Williams' `lbfgsb` release 1.0.0 for the
`scifort_optimize` L-BFGS-B numerical kernel. The dependency is not vendored.
The release manifest identifies the project as BSD-3 licensed and attributes the
copyright to the Regents of the University of California; the source is a modern
Fortran refactor of L-BFGS-B 3.0. The upstream release's `License.txt` is the
authoritative distributed license text (including its literal template fields):

- https://github.com/jacobwilliams/lbfgsb/tree/1.0.0
- https://github.com/jacobwilliams/lbfgsb/blob/1.0.0/License.txt

## SLSQP 1.6.0

SciFort depends through FPM on Jacob Williams' `slsqp` release 1.6.0 for the
`scifort_optimize` SLSQP numerical kernel. The dependency is not vendored. The
release is BSD-3-Clause, Copyright (c) 2016-2025 Jacob Williams, and its LICENSE
also retains the original SLSQP notice for Dieter Kraft (1988) and the
Association for Computing Machinery (1994). See the complete upstream text:

- https://github.com/jacobwilliams/slsqp/tree/1.6.0
- https://github.com/jacobwilliams/slsqp/blob/1.6.0/LICENSE

## nlesolver-fortran 1.1.0

SciFort cross-checks the dense Broyden rank-one update in the inexact-Newton
`scifort_optimize.root` methods against nlesolver-fortran 1.1.0 by Jacob Williams.
The project is not added as an FPM dependency. The reviewed upstream source is:

- https://github.com/jacobwilliams/nlesolver-fortran/tree/1.1.0

BSD 3-Clause License

Copyright (c) 2021-2022, Jacob Williams
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


## Quadpack2 / modern QUADPACK

SciFort depends through FPM on `https://github.com/jacobwilliams/quadpack.git`,
following its `master` branch as of 2026-09-30. The dependency is not vendored in
this repository. Its BSD-3-Clause license is reproduced here because
`scifort_integrate` directly uses its public double-precision QUADPACK interfaces.

BSD 3-Clause License

Copyright (c) 2021-2022, Jacob Williams
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The Quadpack2 LICENSE additionally states that the original QUADPACK from the
SLATEC Library is public-domain work of the United States government.


## SciPy ellipsoidal-harmonic kernel (BSD-2-Clause)

The first-kind Lame-polynomial class matrices, coefficient scaling, and
class-prefactor evaluation in `src/special/ellipsoidal_harmonics.f90` are
adapted from SciPy 1.17.0 `scipy/special/_ellip_harm.pxd`. That source carries
the following BSD-2-Clause notice, which applies to the adapted portions:

Copyright (c) 2012, Jaydeep P. Bardhan
Copyright (c) 2012, Matthew G. Knepley
Copyright (c) 2014, Janani Padmanabhan
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.


## Yoshimura Landau-distribution research reference

The left-tail asymptotic coefficient table used by `src/stats/landau.f90` is
attributed to the classical Koelbig-Schorr expansion as tabulated in:

Takuma Yoshimura, "Numerical Evaluation and High Precision Approximation
Formula for Landau Distribution" (2024), and the associated
`tk-yoshimura/LandauDistribution` research repository.

That research material is published under the Creative Commons Attribution
4.0 International (CC BY 4.0) license. SciFort does not copy or adapt the
repository's implementation source code; the Fortran implementation is original.
License: https://creativecommons.org/licenses/by/4.0/
Full unmodified license text: [LICENSES/CC-BY-4.0.txt](LICENSES/CC-BY-4.0.txt),
retrieved on 2026-09-23 from the research repository's license:
https://raw.githubusercontent.com/tk-yoshimura/LandauDistribution/main/LICENSE.
The canonical license is https://creativecommons.org/licenses/by/4.0/legalcode.txt.
The tabulated rational coefficients are expressed as Fortran constants; their
research attribution is retained in `src/stats/landau.f90` and
`CODE_PROVENANCE.md`.

## MATH77

`src/special/bessel_jy.f90` adapts the MATH77 `DBESJN`, `DBESYN`, and
`DBESPQ` routines, and `src/special/exponential_integrals.f90` adapts the MATH77
`DSI`, `DCI`, and `DCPVAL` kernels distributed by Netlib. The upstream license was
retrieved on
2026-09-23 from https://netlib.org/math/license.html and is reproduced below.

Copyright (c) 1996 California Institute of Technology, Pasadena, California.
ALL RIGHTS RESERVED. Based on Government Sponsored Research NAS7-03001.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer. Redistributions in binary
form must reproduce the above copyright notice, this list of conditions and
the following disclaimer in the documentation and/or other materials provided
with the distribution. Neither the name of the California Institute of
Technology (Caltech) nor the names of its contributors may be used to endorse
or promote products derived from this software without specific prior written
permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

## carlson-elliptic-integrals

`src/special/carlson_elliptic.f90` adapts the `DRC`, `DRD`, `DRF`, and `DRJ`
routines from Jacob Williams' `carlson-elliptic-integrals` release 1.0.0. The
upstream source and license were retrieved on 2026-09-24 from:

- https://raw.githubusercontent.com/jacobwilliams/carlson-elliptic-integrals/1.0.0/src/carlson_elliptic_module.f90
- https://raw.githubusercontent.com/jacobwilliams/carlson-elliptic-integrals/1.0.0/LICENSE

The upstream project notes that its original Fortran 77 SLATEC routines are
public domain. The modern refactoring is distributed under the following
BSD-3-Clause license.

carlson-elliptic-integrals: Carlson symmetric forms of elliptic integrals
https://github.com/jacobwilliams/carlson-elliptic-integrals

Copyright (c) 2016-2022, Jacob Williams
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* The names of its contributors may not be used to endorse or promote products
  derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


## rklib

`src/integrate/ivp.f90` adapts the Bogacki-Shampine 3(2) (`rkbs32`) and
Dormand-Prince 5(4) (`rkdp54`) Runge-Kutta tableaux from Jacob Williams'
`rklib`, branch `master`, retrieved on 2026-09-28 from:

- https://raw.githubusercontent.com/jacobwilliams/rklib/master/src/rklib_variable_steps.f90
- https://raw.githubusercontent.com/jacobwilliams/rklib/master/LICENSE

The upstream license is reproduced below.

BSD 3-Clause License

Copyright (c) 2022-2023, Jacob Williams
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## DOP853

`src/integrate/dop853.f90` adapts the DOP853 tableau, embedded error estimator,
and seventh-order dense-output stages/coefficients from Jacob Williams' modern
Fortran DOP853 implementation, branch `master`, retrieved
on 2026-09-28 from:

- https://raw.githubusercontent.com/jacobwilliams/dop853/master/src/dop853_module.F90
- https://raw.githubusercontent.com/jacobwilliams/dop853/master/LICENSE

The upstream notices are reproduced below.

Modern Fortran Edition of the DOP853 ODE Solver
https://github.com/jacobwilliams/dop853

Copyright (c) 2015-2022, Jacob Williams
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.
* The names of its contributors may not be used to endorse or promote products
  derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Original DOP853 License

Copyright (c) 2004, Ernst Hairer

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## ForDiff

`src/differentiate.f90` adapts the scalar complex-step derivative identity from
ForDiff 0.1.1 by Seyed Ali Ghasemi, retrieved from
<https://github.com/gha3mi/fordiff>. The vector complex-step routines are original
SciFort implementations. The upstream BSD-3-Clause license is reproduced below.

BSD 3-Clause License

Copyright (c) 2023, Seyed Ali Ghasemi

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## roots-fortran

`src/optimize/root_scalar.f90` adapts the bisection, `brentq`, `brenth`, and
Ridders numerical kernels from roots-fortran 1.0.0 by Jacob Williams, retrieved
from <https://github.com/jacobwilliams/roots-fortran/tree/1.0.0>. The upstream
BSD-3-Clause notice is reproduced below.

Roots-Fortran: Root solvers for modern Fortran
https://github.com/jacobwilliams/roots-fortran

Copyright (c) 2021, Jacob Williams
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
* The names of its contributors may not be used to endorse or promote products
  derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## fmin

`src/optimize/minimize_scalar.f90` adapts the bounded Brent minimization kernel
from fmin 1.1.2 by Jacob Williams, retrieved from
<https://github.com/jacobwilliams/fmin/tree/1.1.2>. The upstream BSD-3-Clause
notice is reproduced below.

BSD 3-Clause License

Copyright (c) 2021, Jacob Williams
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



## LSMR (Fong-Saunders)

`src/optimize/lsmr.f90` implements the LSMR iteration described and distributed by
David Fong and Michael Saunders. The implementation is adapted through SciPy 1.17.0's
LSMR code and was cross-checked against Jacob Williams's modern Fortran LSMR project.
The original LSMR software is distributed under the following BSD-2-Clause license.

Copyright (c) 2010-2014 David Fong, Michael Saunders, Stanford University
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.


## TRLIB

The native `trust_krylov` implementation follows the GLTR algorithm, Krylov
forcing rules, and hot-start behavior documented by TRLIB:
<https://github.com/felixlen/trlib>. SciFort does not vendor or link TRLIB; the
Fortran implementation stores its own explicit Krylov basis. TRLIB is licensed
under the MIT License.

Copyright (c) 2016--2017 Felix Lenders

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Spectrum (reference and algorithmic attribution)

The design of the 1D Welch/cross-spectrum estimation and segment accumulation
in `src/signal_spectral.f90` was informed by `spectrum`, by Jason Christopherson:
https://github.com/jchristopherson/spectrum

Source pages consulted 2026-10-15:
https://github.com/jchristopherson/spectrum/blob/main/LICENSE
https://raw.githubusercontent.com/jchristopherson/spectrum/main/src/spectrum_periodogram.f90
https://raw.githubusercontent.com/jchristopherson/spectrum/main/src/spectrum_windows.f90

No dependency on `spectrum` is required; SciFort's FFT and public APIs are
separately implemented. Its MIT license and copyright notice are included here
for attribution to the source that informed the design.

MIT License

Copyright (c) 2023-2026 Jason Christopherson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


## PEAKS (algorithmic inspiration, no copied source code)

The new `src/signal_peak_excursion.f90` independently implements a
minimum-excursion alternating peak/valley detector. Algorithmic attribution:
Eli Billauer, *peakdet*, as credited by Jason Christopherson's
MIT-licensed PEAKS project, https://github.com/jchristopherson/peaks.

PEAKS license (https://github.com/jchristopherson/peaks/blob/main/LICENSE):

MIT License

Copyright (c) 2023 Jason Christopherson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

No PEAKS source was copied, translated, vendored or linked: this license is
reproduced for transparent attribution to the algorithmic inspiration.

## SciPy Clough–Tocher interpolation formulas

`src/interpolate_clough_tocher.f90` adapts the gradient-estimation and
triangular cubic Bernstein–Bézier formulas from Pauli Virtanen's 2010
implementation in SciPy 1.17.0, `scipy/interpolate/_interpnd.pyx`.
The source is BSD-3-Clause; its notice, license conditions and disclaimer
are retained in `LICENSES/scipy-interpnd-BSD-3-Clause.txt`. No Qhull code
is incorporated. The surrounding SciFort Fortran API and mesh handling are
original work.

## FITPACK smoothing backend (optional, 2026-11-01)

The `optional/fitpack_backend/` manifest adds a separate version/tag-pinned
FPM dependency on Federico Perini's `perazz/fitpack` implementation of Paul
Dierckx's FITPACK. It is licensed BSD-3-Clause. SciFort's optional adapter
and dependency-free TCK compatibility layer are independently authored and
MIT-licensed; the optional adapter does not vendor the full upstream library.
The native spherical-constraint adaptation is described above. For clear
attribution and license compliance, the upstream BSD license is reproduced
at `LICENSES/perazz-fitpack-BSD-3-Clause.txt`. The pinned upstream package
itself was not available for a real integration build here; mock-interface
checks should not be confused with testing the upstream implementation.

## Davenport-angle conversion (SciPy BSD-3-Clause)

The quaternion half-angle decomposition used by
`rotation_type%as_davenport` follows the published generalized Euler
rotation algorithm of Shuster and Markley (2003), and Bernardes and
Viollet (2022), with the implementation cross-checked/adapted against
SciPy 1.17's `_rotation_xp.py`. The complete SciPy BSD-3-Clause license
and attribution for this small adaptation are preserved in
`LICENSES/scipy-spatial-transform-BSD-3-Clause.txt`.
The Fortran routines have no runtime dependency on SciPy or Python.

## Differential evolution algorithm attribution (Ian Morrall, MIT)

The independently written SciFort `scifort_optimize_differential_evolution`
module follows the DE/rand/1/bin evolutionary method and dithering as also
presented in Ian Morrall's MIT-licensed Fortran implementation:
https://github.com/ian-mmm/differential-evolution_f95 . The copyright (c) 2017
Ian Morrall and complete MIT permission and disclaimer are preserved in
`LICENSES/ian-morrall-differential-evolution-MIT.txt`. This does not add a
runtime dependency, and no original source file was vendored unchanged.
