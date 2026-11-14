## 2026-11-14 - Periodic native spline constructors

The `bc_type='periodic'` branches in `make_splrep` and `make_splprep` and
the private normalization/physical-validation wrapper in
`src/interpolate_fitpack_compat.f90` are original MIT SciFort code. They
reuse existing periodic native interpolation and the already attributed
`fpperi`/`fpclos` smoothing adaptations. No new third-party implementation
or linked dependency is introduced. The interpolation kernel retains its
existing SciPy 1.17.0 provenance recorded below; it is not rewritten here.

`src/interpolate_curfit_core.f90` gains trailing optional hard knot capacities
in its two periodic entry points and an integer-product guard for limited
parametric capacities. Default allocation, knot selection, QR factors,
penalties, root searches, and legacy tuple behavior are preserved. The exact
cached upstream `fpperi` and `fpclos` source and license were inspected again:
perazz/FITPACK commit `1907c7bad9f5c4557970e906bdd4909eb52b2d4c`,
[`src/fitpack_core.F90`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/src/fitpack_core.F90),
SHA-256 `02e413063f7df855c45eed6119130f7ab7a07cfa3e2f6ac0e42b6f496df414ec`.
The `MIT AND BSD-3-Clause` notice, Copyright (c) 2022 Federico Perini,
and [full retained license](LICENSES/perazz-fitpack-BSD-3-Clause.txt) remain.

The public mathematical reference is P. Dierckx, *Curve and Surface Fitting
with Splines*, Oxford University Press, 1993, periodic interpolation and
derivative-jump smoothing under a weighted residual constraint. Exact SciPy
v1.18.0 [`_fitpack_repro.py`](https://github.com/scipy/scipy/blob/v1.18.0/scipy/interpolate/_fitpack_repro.py),
its BSD-3-Clause license, and official constructor documentation were inspected
for interface contracts and differences; no new SciPy code was translated.
The new tests and validation driver/checker are original MIT code. See
[validation](validation/periodic_constructors_2026-11-14.md) for strict
periodic interpolation comparisons, smoothing-equation checks, and
automatic-path limitations.

## 2026-11-14 - Native parametric smoothing constructor

`make_splprep` in `src/interpolate_fitpack_compat.f90` is an original MIT
SciFort wrapper around native array interpolation and the existing attributed
parametric FITPACK solver in `src/interpolate_curfit_core.f90`. Its chord
parameter calculation reuses the existing SciFort stable coordinate-difference
helper. The new wrapper, tests, and validation programs are original MIT code.
No new third-party implementation or dependency is imported.

The existing `fppara` adaptation gains optional fixed-knot and polynomial-only
modes, full-knot polynomial representation, hard capacities, and factor checks
before back substitution. The new fixed mode uses reciprocal-mean triangular
diagonal initialization, as in the scalar solver; legacy automatic calls
retain their mean initialization and knot/root search. Existing normalization,
shared residuals, dynamic channel storage, and physical checks are reused.

The exact pinned perazz/FITPACK source and BSD license were inspected again:
commit `1907c7bad9f5c4557970e906bdd4909eb52b2d4c`,
[`src/fitpack_core.F90`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/src/fitpack_core.F90),
SHA-256 `02e413063f7df855c45eed6119130f7ab7a07cfa3e2f6ac0e42b6f496df414ec`.
Copyright (c) 2022 Federico Perini and the `MIT AND BSD-3-Clause` file notice
remain. Full terms remain in
[the retained license](LICENSES/perazz-fitpack-BSD-3-Clause.txt).
The mathematical reference is P. Dierckx, *Curve and Surface Fitting with
Splines*, Oxford University Press, 1993, derivative-jump smoothing with one
weighted residual constraint across coordinates.

The exact installed SciPy 1.15.3 `make_splprep` and `_make_splrep_impl`
sources and current official constructor documentation were inspected for
interface semantics and independent validation. They were not copied or
translated. Modern automatic knot selection is not adopted. Hard capacities,
limited knot forms, and unsupported periodic constructor boundaries are
documented. See [the validation record](validation/make_splprep_2026-11-14.md).

## 2026-11-13 - Native scalar smoothing constructor

`make_splrep` in `src/interpolate_fitpack_compat.f90` is an original MIT
SciFort wrapper around existing native interpolation and the attributed
scalar FITPACK kernel in `src/interpolate_curfit_core.f90`. No new
third-party implementation is imported. The existing `fpcurf` adaptation
is extended with optional fixed-knot and polynomial-only modes, retaining
its original defaults, derivative-jump penalty, banded QR, and root search.
Polynomial limits on supplied knots are represented by fitting values of
the already computed polynomial. A shifted-normalization option and
physical residual verification apply to the new constructor.

The adapted kernel remains attributed to perazz/FITPACK at commit
`1907c7bad9f5c4557970e906bdd4909eb52b2d4c`, source
[`src/fitpack_core.F90`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/src/fitpack_core.F90),
SHA-256 `02e413063f7df855c45eed6119130f7ab7a07cfa3e2f6ac0e42b6f496df414ec`.
The retained [BSD-3-Clause license](LICENSES/perazz-fitpack-BSD-3-Clause.txt),
Copyright (c) 2022 Federico Perini, and source `MIT AND BSD-3-Clause` notice
are unchanged. The mathematical reference is P. Dierckx, *Curve and Surface
Fitting with Splines*, Oxford University Press, 1993, derivative-jump
smoothing under a weighted residual constraint.

SciPy's official `make_splrep` documentation and the exact installed
SciPy 1.15.3 public constructor and helper source were inspected for
interface semantics and independent comparison, not copied or translated.
The modern automatic knot-selection implementation is not adopted; native
automatic fits retain the existing FITPACK path. Hard capacity, limited
knot forms, and unsupported periodic constructor boundaries are documented
compatibility limits. The new tests and validation programs are original
MIT code. See [the validation record](validation/make_splrep_2026-11-13.md).

## 2026-11-13 - Native periodic parametric spline fitting

`src/interpolate_curfit_core.f90` also adapts `fpclos`,
`fpclos_reset_interp`, and `fp_rotate_row_2mat_vec` from Federico Perini's
perazz/FITPACK, based on Paul Dierckx's FITPACK. The inspected source is
[`src/fitpack_core.F90`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/src/fitpack_core.F90),
commit `1907c7bad9f5c4557970e906bdd4909eb52b2d4c`, SHA-256
`02e413063f7df855c45eed6119130f7ab7a07cfa3e2f6ac0e42b6f496df414ec`.
The inspected [upstream license](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/LICENSE.md)
has SHA-256 `39fbeb302cce5182528712269832f82d356e684b9a3a11dc6949a79ef64781b8`.
Full BSD-3-Clause terms and Copyright (c) 2022 Federico Perini are retained
in [the existing license file](LICENSES/perazz-fitpack-BSD-3-Clause.txt).
The adapted file retains source attribution and `MIT AND BSD-3-Clause`.

Changes are binary64/default integer/logical kinds, lowercase free form,
explicit `implicit none`, documented dummy arguments, wrapped comments,
ordinary loops replacing `forall`, a channel-sized working vector replacing
the upstream fixed maximum, and existing periodic factor checks before
wrapped back substitution. Existing attributed basis, knot insertion, vector
Givens, wrapped solve, and root-iteration helpers are reused unchanged.
Original SciFort MIT wrappers validate active entries, compute scaled closed
chord lengths, normalize coordinates/weights/parameters, restore unpadded
physical tuples, and verify physical residuals. The public `splprep` signature
and existing nonperiodic behavior are preserved. No dependency is added.

Validation independently compiles the unchanged pinned original core and
compares 187 fresh fits, including seven fresh-fit cases using coordinates
read from upstream `test_closed_fit` in `test/fitpack_curve_tests.f90`
(inspected SHA-256
`0bb92548cb070e70e245b1f837fb931014e69b6f1c335f8cfdca49c514908847`).
Those upstream literals are read from an optional reference file at run time,
not copied into native tests or vendored. Their default-real rounding is
preserved before promotion to binary64. The full upstream object-oriented
and continuation sequence is not run. SciPy is a comparison reference only;
no new SciPy implementation is imported. See the
[validation record](validation/splprep_closed_2026-11-13.md).

Mathematical references: P. Dierckx, "Algorithms for smoothing data with
periodic and parametric splines", *Computer Graphics and Image Processing*
20 (1982), 171-184; and *Curve and Surface Fitting with Splines*, Oxford
University Press, 1993, periodic spline fitting.

## 2026-11-13 - Native periodic scalar spline fitting

`src/interpolate_curfit_core.f90` also adapts `fpperi`,
`fpperi_reset_interp`, `fpbacp`, and `fp_rotate_row_2mat` from Federico
Perini's perazz/FITPACK, based on Paul Dierckx's FITPACK. The exact source is
[`src/fitpack_core.F90`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/src/fitpack_core.F90),
commit `1907c7bad9f5c4557970e906bdd4909eb52b2d4c`. The inspected source
SHA-256 is `02e413063f7df855c45eed6119130f7ab7a07cfa3e2f6ac0e42b6f496df414ec`;
the inspected upstream `LICENSE.md` SHA-256 is
`39fbeb302cce5182528712269832f82d356e684b9a3a11dc6949a79ef64781b8`.
Full BSD-3-Clause terms and Copyright (c) 2022 Federico Perini are retained
in [the existing license file](LICENSES/perazz-fitpack-BSD-3-Clause.txt).
The combined file retains `MIT AND BSD-3-Clause` and source attribution.

Changes are binary64/default integer/logical kinds, lowercase free form,
explicit `implicit none`, documented dummy arguments, reuse of the existing
attributed basis/knot/Givens/root helpers, and factor-diagonal checks before
wrapped back substitution. Original MIT SciFort wrappers normalize inputs,
validate ignored/active entries separately, guard ranges and allocation,
restore physical tuples, and verify actual physical residuals. Fresh
interpolation and smoothing are supported; upstream continuation state and
fixed-knot validation are not imported. No build dependency is added.

The original core was compiled separately and compared on 181 fresh fits.
The sine/cosine geometry and analytic tolerance of upstream
`test_periodic_fit` were inspected in pinned `test/fitpack_curve_tests.f90`
(SHA-256 `0bb92548cb070e70e245b1f837fb931014e69b6f1c335f8cfdca49c514908847`)
and reproduced mathematically in the original validation checker. This does
not claim that the complete upstream object-oriented or continuation test
suite was run. SciPy 1.15.3 `fpperi.f` was inspected for comparison only;
no new SciPy source is imported. Adaptive paths can differ; see the
[validation record](validation/splrep_periodic_2026-11-13.md).
The mathematical reference is P. Dierckx, *Curve and Surface Fitting with
Splines*, Oxford University Press, 1993, periodic spline fitting.

## 2026-11-13 - Native parametric automatic-knot smoothing

`src/interpolate_curfit_core.f90` now also adapts `fppara`,
`fp_rotate_row_vec`, and `fp_rotate_shifted_vec` from Federico Perini's
perazz/FITPACK, based on Paul Dierckx's FITPACK:
[`src/fitpack_core.F90`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/src/fitpack_core.F90),
commit `1907c7bad9f5c4557970e906bdd4909eb52b2d4c`. The source and
[license](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/LICENSE.md)
were inspected; source SHA-256 is
`02e413063f7df855c45eed6119130f7ab7a07cfa3e2f6ac0e42b6f496df414ec`,
and license SHA-256 is
`39fbeb302cce5182528712269832f82d356e684b9a3a11dc6949a79ef64781b8`.
The existing [retained BSD license](LICENSES/perazz-fitpack-BSD-3-Clause.txt)
contains the full terms and Copyright (c) 2022 Federico Perini notice;
only blank-line formatting differs. The combined source retains
`MIT AND BSD-3-Clause` for original wrappers and adapted numerical kernels.

Adaptation changes include SciFort's binary64 kind and default integers/
logicals, lowercase free-form formatting, explicit `implicit none`, documented
dummy arguments, and one statement per line. The existing attributed private
B-spline, triangular-band, knot-insertion, Givens, and root-iteration helpers
are reused without changing their numerical implementations. The parametric
kernel retains shared-knot selection, a single penalty parameter across all
coordinate channels, its mean-diagonal initial guess, relative residual
tolerance 0.001, and a maximum of 20 smoothing iterations. The surrounding
upstream `parcur` routine was inspected; SciFort calls `fppara` directly
with task zero, explicit workspace, and supplied normalized parameters.

Original SciFort integration adds a normalized wrapper, coordinate-fastest
input packing/channel-major output unpacking, active-factor rank checks,
achieved-residual verification, finite/range checks, status mapping, and
scaled chord lengths that avoid overflowing the physical coordinate difference.
Chord differences are formed before scaling when representable, preserving
geometry at large common offsets, and before-subtraction scaling is used only
where needed to avoid overflow.
Origin-first parameter normalization preserves large-offset spacing, and a
second residual check validates physical knots and rounded/rescaled coefficients.
Unrepresentable physical spans or unresolved physical knot intervals fail explicitly.
The existing public signature, default s=0, unpadded coefficient shape,
interpolation/fixed-knot paths, periodic failure status, and arbitrary positive
channel-count commitment are preserved. More than ten channels are documented
as a SciFort extension to legacy SciPy splprep. There is no new package,
external linked dependency, or C ABI.

The inspected mathematical references include Dierckx, "Algorithms for
smoothing data with periodic and parametric splines", Computer Graphics and
Image Processing 20 (1982), 171-184, and *Curve and Surface Fitting with
Splines* (Oxford, 1993). Validation independently forms B-spline design and
kth-derivative-jump matrices to check the penalized normal equations. It also
compares the separately compiled original pinned core and SciPy's residual
contracts. SciPy 1.15.3's
[`fppara.f`](https://github.com/scipy/scipy/blob/v1.15.3/scipy/interpolate/fitpack/fppara.f)
was inspected to confirm its reciprocal-mean initial guess differs from the
adopted mean-diagonal guess; no SciPy implementation code was imported for
this change. Iteration-path and high-derivative differences are documented.

Optional upstream fixture input is read from `test_parametric_fit` in
[`test/fitpack_curve_tests.f90`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/test/fitpack_curve_tests.f90)
(SHA-256 `0bb92548cb070e70e245b1f837fb931014e69b6f1c335f8cfdca49c514908847`).
The checker preserves default-real literal rounding and exercises six fresh
fits drawn from that test's settings, rather than its stateful continuation
sequence. Upstream data are not vendored and the full upstream suite was not
run. New native tests, validation drivers, checker, and wrapper code are
original MIT-licensed SciFort contributions. See the
[validation record](validation/splprep_smoothing_2026-11-13.md).

## 2026-11-13 - Native scalar automatic-knot smoothing

`src/interpolate_curfit_core.f90` adapts scalar `fpcurf` and supporting kernels
from Federico Perini's perazz/FITPACK, based on Paul Dierckx's FITPACK:
[`src/fitpack_core.F90`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/src/fitpack_core.F90),
commit `1907c7bad9f5c4557970e906bdd4909eb52b2d4c`. The exact source and
[`LICENSE.md`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/LICENSE.md)
were inspected before adaptation. Source SHA-256:
`02e413063f7df855c45eed6119130f7ab7a07cfa3e2f6ac0e42b6f496df414ec`;
inspected license SHA-256:
`39fbeb302cce5182528712269832f82d356e684b9a3a11dc6949a79ef64781b8`.

Copyright (c) 2022, Federico Perini. All rights reserved. The BSD-3-Clause
terms and copyright are retained in
[`LICENSES/perazz-fitpack-BSD-3-Clause.txt`](LICENSES/perazz-fitpack-BSD-3-Clause.txt);
its text differs from the inspected license only in blank-line formatting.
The combined source is marked `MIT AND BSD-3-Clause`: the SciFort wrapper
is original MIT code; the numerical kernels remain attributed adaptations.

Adapted routines: `fpcurf`, `fpknot`, `fp_knots_to_add`, `fp_rank_deficient`,
`fp_rotate_row`, `fp_rotate_shifted`, `fpback`, `fpbspl`, `fpdisc`, `fpgivs`,
`fprota`, `fp_knot_interval`, `root_finding_iterate`, `fprati`, `equal`, and
`not_equal`. The last ten reuse already adapted, attributed copies from
SciFort's spherical grid core; they remain private copies here so scalar
smoothing does not depend on spherical fitting. No existing shared numerical
kernel is changed. The `curfit` workspace partition and source were inspected;
the wrapper calls `fpcurf` directly with task zero and fully allocated arrays.

Modifications include SciFort's binary64 kind, default integers/logicals,
explicit `implicit none`, documented dummy arguments, lowercase free-form
formatting, and a private module interface. FITPACK's relative tolerance
`0.001`, maximum 20 smoothing iterations, derivative-jump penalty, adaptive
knot selection, rank checks, and rational root iteration are preserved.
Original SciFort code normalizes coordinates/values/weights, rescales knots
and coefficients, rejects unrepresentable/range-limited inputs, verifies the
achieved residual, pads coefficients, and maps failure to the existing status 3.
Task-zero interpolation and fixed-knot fitting retain their former solvers.
Standalone script source lists include the new dependency; the package
manifest and external dependencies are unchanged.

References inspected in the upstream source include de Boor (1972),
"On calculating with b-splines", J. Approx. Theory 6, 50-62; Cox (1972),
"The numerical evaluation of b-splines", J. Inst. Maths Applics 10, 134-149;
and Dierckx (1993), *Curve and Surface Fitting with Splines*, Oxford.

Validation compares with SciPy 1.15.3 and the original, separately compiled
pinned Fortran core. The upstream `test_fpknot_crash` data are read from
[`test/fitpack_curve_tests.f90`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/test/fitpack_curve_tests.f90)
(SHA-256 `0bb92548cb070e70e245b1f837fb931014e69b6f1c335f8cfdca49c514908847`)
when that optional fixture is requested; no upstream test data are vendored.
This exercises its data at three explicit smoothing targets rather than
running the upstream spline-object test unchanged. The full upstream suite
was not run. New native tests and validation drivers/checker are original MIT
code. See the [validation record](validation/splrep_smoothing_2026-11-13.md).

## 2026-11-13 - Larger AAA arrowhead spectra

`src/interpolate_aaa_spectrum.f90` adapts the following files from
[Reference-LAPACK v3.12.1](https://github.com/Reference-LAPACK/lapack/tree/v3.12.1).
The exact source and project license were inspected before adaptation.

| Upstream file | Adapted routines/content | SHA-256 of inspected file |
| --- | --- | --- |
| `SRC/zhgeqz.f` | `ZHGEQZ` complex single-shift QZ | `8acad110608790e4348f69455f7cbed2f67a46b4a08bf8fb7c4a3e93fe14c6a4` |
| `SRC/zlartg.f90` | Safe complex plane rotations | `d25bf2ce6839de8c73fd6825a5cad0e1c1e154fb6a20c05b91a539868732e5b4` |
| `SRC/zrot.f` | Complex vector rotation | `1f0d5c47fc7c59dce2d8f57bf90cc1822c6489762cf51486052042350bcf4c59` |
| `BLAS/SRC/zscal.f` | Complex vector scaling | `e510a73d828efce4525349524e3ff82ac911cbff7039e68f41a7fb3837cc14fa` |
| `SRC/zladiv.f` | Scaled complex division wrapper | `4f7d0281958f53345cdd8d39a72be1a0f16fa71d0f152431d6c8fbe28f3dc28a` |
| `SRC/dladiv.f` | `DLADIV`, `DLADIV1`, `DLADIV2` | `e3ee0a30d5b573425aba201c9d9d494c23849d6046a061c39b92094728ce7f5b` |
| `SRC/la_constants.f90` | Safe-scaling constant definitions only | `1d9bdd0fad20e76f6c418b0d5405aa7a256d89cfe57eba2a36a9c301c5ce7d1b` |
| `LICENSE` | Full retained license | `978539245c405775bf347b1256c36df4daf8b69a64b099e94b1d53742201621e` |

Copyright (c) 1992-2023, University of Tennessee and the University of Tennessee
Research Foundation; (c) 2000-2023, University of California Berkeley;
(c) 2006-2023, University of Colorado Denver. All rights reserved. Applicable
notices and NAG attribution remain in source, along with ZLARTG's credit to
Weslley Pereira (2021) and Anderson's safe-scaling reference. Complete license
conditions, patent disclaimer, and warranty disclaimer are retained verbatim
in [`LICENSES/lapack-BSD-3-Clause.txt`](LICENSES/lapack-BSD-3-Clause.txt).
The adapted file is marked `MIT AND BSD-3-Clause` for original SciFort wrappers
and third-party kernels respectively. No linked LAPACK/BLAS package is required.

The new backend is restricted to normalized AAA arrowhead pencils
`A=[0,u;1,diag(s)]`, `B=diag(0,I)`, through 512 support points (pencil order
513). Native Hessenberg reduction fixes the first coordinate, so all its
Householder similarities commute with `B`; no separate generalized reduction
is needed. `ZHGEQZ` is invoked only with eigenvalues requested and no Schur
vectors. Existing public general-purpose QZ/eigenvalue/matrix-function routines
and their size limits are unchanged.

Significant modifications:

- Convert fixed-form source to lowercase free-form, use SciFort's binary64
  kind, explicit interfaces, `implicit none`, documented dummy intents, named
  function results, and structured `do` loops while preserving jump labels.
- Replace legacy double-specific intrinsics with kind-qualified standard
  intrinsics and statement functions with private module functions.
- Specialize private machine/option adapters to the internal binary64,
  lowercase-option calls. Scaling constants follow inspected `la_constants`;
  `INSTALL/dlamch.f` was also inspected for the E/S/O/B definitions.
- Replace `ZLANHS` calls with direct Frobenius norms. The wrapper bounds support
  coordinates/weights by two and count by 512, so these normalized matrix norms
  cannot overflow. The full general-purpose norm routine is not imported.
- Remove `XERBLA` printing; invalid calls return a code. Wrapper allocation,
  input, convergence, and finite-output checks return explicit statuses.
- Reuse SciFort's existing Hessenberg reducer without changing it. Retain
  AAA's existing projective infinity filter, secular polishing/checks, physical
  rescaling, and cleanup rules. Allocation failure maps to status 5, spectral
  numerical failure to status 3, and support counts above 512 to status 4.
- Track finite projective roots with an explicit logical mask rather than
  NaN sentinels, including in optimized builds with exception traps enabled.

The LAPACK `ZHGEQZ` reference is Moler and Stewart, "An Algorithm for
Generalized Matrix Eigenvalue Problems", SIAM J. Numer. Anal. 10 (1973),
241-256. `ZLARTG` cites Anderson, "Algorithm 978: Safe Scaling in the Level 1
BLAS", ACM Trans. Math. Softw. 44 (2017),
[DOI](https://doi.org/10.1145/3061665). The AAA algorithm and its original
SciPy/Oxford adaptation remain attributed in the existing AAA entry below.

Validation covers homogeneous finite/infinite spectra, clustered coordinates,
degree cancellation, analytic roots through 512 supports, public real/complex
spectra/residues, and actual large-representation cleanup. The full upstream
LAPACK suite was not run; this wrapper supports a restricted pencil family,
not the complete LAPACK interfaces. SciPy's independently compiled LAPACK
provides numerical comparisons alongside analytic oracles. Tests and validation
drivers/checkers are original MIT-licensed SciFort code. See the
[validation record](validation/aaa_large_2026-11-13.md). No full upstream
distribution, downloaded archive, binary artifact, dependency declaration,
or C ABI is added.

## 2026-11-13 - Adaptive scattered spherical smoothing

`src/interpolate_sphere_scattered_core.f90` adapts `fpsphe`, `fprank`,
`fprpsp`, `fporde`, `fp_rotate_row`, and `fp_rotate_shifted`, plus the
workspace partition from `sphere`, in Federico Perini's perazz/FITPACK:
[`src/fitpack_core.F90`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/src/fitpack_core.F90),
commit `1907c7bad9f5c4557970e906bdd4909eb52b2d4c`.
The exact source and license were inspected. Source SHA-256:
`02e413063f7df855c45eed6119130f7ab7a07cfa3e2f6ac0e42b6f496df414ec`.
Upstream `LICENSE.md` SHA-256:
`39fbeb302cce5182528712269832f82d356e684b9a3a11dc6949a79ef64781b8`.
Copyright (c) 2022, Federico Perini. All rights reserved. Upstream credits
Paul Dierckx's FITPACK. BSD-3-Clause notices remain in adapted source and full
terms remain in
[`LICENSES/perazz-fitpack-BSD-3-Clause.txt`](LICENSES/perazz-fitpack-BSD-3-Clause.txt).

The new core reuses existing adapted FITPACK basis, rotation, banded solve,
penalty, interval-search, and rational root-iteration helpers from
`src/interpolate_sphere_grid_core.f90`. Only their internal module visibility
changes; their numerical implementation is unchanged. The new constructor in
`src/interpolate_sphere.f90` validates inputs, guards workspace sizes using
bounded integer arithmetic, scales observations/weights and the residual
target, rejects unrepresentable coefficients, and returns the shared spline
object with physical residuals and explicit statuses. Original tests,
reference tools, and native wrapper additions are SciFort work; adapted source
retains its applicable third-party license.

Significant adaptation changes:

- Use SciFort's binary64 kind/constants, default integers, lowercase ASCII,
  explicit interfaces, `implicit none`, and documented dummy arguments.
- Allocate and initialize all solver state locally, including sufficient
  secondary workspace for rank reduction. No persistent cache or warm start
  is exposed. Optional knot capacities are SciFort extensions; defaults follow
  SciPy 1.15.3's `spherfit_smth` wrapper.
- Preserve `f(1:ncof)` in the first rank-deficient solve; the pinned modern
  source instead assigns zero to the transformed right-hand side. The legacy
  [`fpsphe.f`](https://github.com/scipy/scipy/blob/v1.15.3/scipy/interpolate/fitpack/fpsphe.f)
  was inspected to cross-check this correction, rather than imported.
- Correct `fprank` dummy intents for the factor and transformed right-hand
  side from `out` to `inout`, because they are read before being overwritten.
- After rank reduction, recompute the actual weighted residual from the
  physical tensor coefficients before knot selection or accepting a fit.
  The accumulated correction underestimated it in the equatorial fixture:
  approximately `0.00156929` versus the actual `0.003103635964350815`.
- Expose the final factor rank, with -1 for the initial polynomial where it
  is not estimated. Keep valid deficient fits with status 2, unlike SciPy
  1.15.3's high-level rejection of its rank-warning codes. Positive solver
  failures return an invalid fit; no partial solution is silently accepted.

An independent constrained tensor design and NumPy dense SVD verify the
rank-four minimum-norm fixture, including predictions away from observations;
the existing native dense SVD supplies another check. Full-rank fits are checked
against SciPy and polynomial limits against independent weighted least squares.
The upstream `mnsphe` fixture is exercised with fresh smoothing fits; no
upstream test source/data, full source distribution, archive, or binary is
redistributed. The full upstream package suite was not run. See the
[validation record](validation/smooth_sphere_2026-11-13.md).

Mathematical references, as cited by upstream: P. Dierckx, "Algorithms for
smoothing data on the sphere with tensor product splines", *Computing* 32
(1984), 319-342; and *Curve and Surface Fitting with Splines* (Oxford University
Press, 1993). The SciPy API/capacity reference inspected was
[`dfitpack.pyf`](https://github.com/scipy/scipy/blob/v1.15.3/scipy/interpolate/src/dfitpack.pyf).
No dependency or C ABI is added.

## 2026-11-13 - Adaptive rectangular spherical smoothing

`src/interpolate_sphere_grid_core.f90` now additionally adapts `fpspgr`,
`fpknot`, `fp_knots_to_add`, `new_knot_dimension`, `root_finding_iterate`, and
`fprati` from the same inspected perazz/FITPACK source as the interpolation
entry below: commit `1907c7bad9f5c4557970e906bdd4909eb52b2d4c`,
[`src/fitpack_core.F90`](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/src/fitpack_core.F90).
Source SHA-256 remains
`02e413063f7df855c45eed6119130f7ab7a07cfa3e2f6ac0e42b6f496df414ec`;
inspected `LICENSE.md` SHA-256 remains
`39fbeb302cce5182528712269832f82d356e684b9a3a11dc6949a79ef64781b8`.
Copyright (c) 2022, Federico Perini. All rights reserved. BSD-3-Clause source
notices and the full retained license in
[`LICENSES/perazz-fitpack-BSD-3-Clause.txt`](LICENSES/perazz-fitpack-BSD-3-Clause.txt)
continue to apply. Mathematical reference: Dierckx, *Curve and Surface Fitting
with Splines* (1993), Chapter 11, Section 11.2, as cited by upstream.

The interpolation path remains separate. Positive targets invoke adaptive
axis/knot selection followed by FITPACK's rational smoothing-parameter
iteration, with relative tolerance 0.001 and a 20-iteration limit. All mutable
cache state is local and initialized; warm starts and persistent knot reuse
are not exposed. Existing banded penalty calculations are now reachable.
New helpers use SciFort's binary64 kind, default integers, lowercase ASCII,
explicit interfaces, `implicit none`, and documented dummy arguments.

The axis selector alternates equal-cost choices according to the previous
axis, preserving the spherical FITPACK rule seen in SciPy v1.15.3's
[`fpspgr.f`](https://github.com/scipy/scipy/blob/v1.15.3/scipy/interpolate/fitpack/fpspgr.f).
The inspected modern generic selector always favors longitude on ties, which
can change the adaptive path. The previously documented south-only C1
cosine/sine initialization correction remains. The wrapper scales the
residual target in logarithms along with observations/poles, rejects scaled
targets below normal binary64 range, reports physical squared grid residuals,
and rejects unconverged or nonfinite fits with status 3. The former unsupported
smoothing constant remains exported for source compatibility.

Validation compares against SciPy C0 fits, independent constrained cubic
least squares, and a separately compiled full pinned modern FITPACK core with
the same two documented axis/pole corrections. This tests the extraction and
wrapper against the established implementation; it is not an independent
algorithmic reference for every C1 fit. SciPy 1.15.3 produced nonfinite or
call-dependent results for some tested C1 smoothing cases. The upstream
`mnspgr` spherical-grid fixture was also checked. See the
[validation record](validation/rect_sphere_smoothing_2026-11-13.md).
The full upstream package test suite was not run. Original reference drivers,
checkers, and tests are MIT licensed. No full upstream source, test data,
binary artifact, dependency, or C ABI is added to the repository.

## 2026-11-13 - Rectangular spherical grid interpolation

`src/interpolate_sphere_grid_core.f90` adapts the `s=0` branches of `spgrid`
and `fpspgr`, together with the following numerical helpers, from Federico
Perini's modern Fortran FITPACK implementation:
`equal`, `not_equal`, `fp_knot_interval`, `fp_rotate_2mat_stride`,
`fp_rotate_row_block`, `fp_rotate_row_stride`, `fpback`, `fpbacp`, `fpbspl`,
`fpcyt1`, `fpcyt2`, `fpdisc`, `fpgivs`, `fpgrsp`, `fpopsp`, `fprota`, `fpsysy`.
The source is
[src/fitpack_core.F90](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/src/fitpack_core.F90),
commit `1907c7bad9f5c4557970e906bdd4909eb52b2d4c`.
Source SHA-256: `02e413063f7df855c45eed6119130f7ab7a07cfa3e2f6ac0e42b6f496df414ec`.
Inspected upstream `LICENSE.md` SHA-256:
`39fbeb302cce5182528712269832f82d356e684b9a3a11dc6949a79ef64781b8`.
Copyright (c) 2022, Federico Perini. All rights reserved. The upstream source
credits Paul Dierckx's FITPACK. BSD-3-Clause notices remain in adapted source;
full terms are retained in
[`LICENSES/perazz-fitpack-BSD-3-Clause.txt`](LICENSES/perazz-fitpack-BSD-3-Clause.txt).
The retained and inspected license texts agree after whitespace normalization.

The new driver selects periodic longitude knots, estimates missing pole values
and derivatives using a coarse latitude fit, then constructs the exact
interpolating latitude knots and invokes the banded Givens solver. It does
not expose adaptive smoothing; shared helpers retain penalty branches that
are unreachable through the `s=0` driver. No full upstream package, downloaded
archive, binary artifact, or mandatory dependency is included.

Significant changes: SciFort binary64 kind and constants, default integers,
private module packaging, explicit interfaces and `implicit none`, documented
dummy arguments, and equivalent `do` loops replacing `forall`. Work arrays
are initialized. The public wrapper in `src/interpolate_sphere.f90` validates
inputs, bounds workspace and coordinate separation, scales grid and pole
data, checks physical coefficient range, and reuses the existing result type
and B-spline evaluator. Query longitude wrapping uses the grid origin.
The banded solver does not estimate rank; the public diagnostic returns -1.

The cosine/sine collocation initialization in `fpgrsp` now runs when either
pole requests C1 continuity (`iop0/=0 .or. iop1/=0`), fixing uninitialized
data in the south-only case. The pole-optimization step range includes both
supplied/estimated poles and all grid values. This differs from SciPy 1.15.3's
legacy initialization on a constant grid with an inconsistent opposite pole;
an independent constrained cubic least-squares oracle verifies the result.
SciPy v1.15.3 `spgrid.f`, `fpspgr.f`, `fpopsp.f`, and `fpgrsp.f` were inspected
to diagnose the difference, not copied or adapted into this addition.

Mathematical reference: P. Dierckx, *Curve and Surface Fitting with Splines*,
Oxford University Press, 1993, Chapter 11, Section 11.2, as cited by upstream.
The native tests and validation tools are original MIT-licensed SciFort code.
The [validation record](validation/rect_sphere_2026-11-13.md) reports analytic,
edge-case, independent least-squares, and SciPy numerical checks. The full
upstream `mnspgr` interpolation fixture was also checked (11-by-14 grid);
the full upstream package test suite was not run. No upstream test data or
test source is redistributed. C ABI expansion remains deferred.

## 2026-11-13 - Fixed-knot spherical spline constraints

`src/interpolate_sphere.f90` adapts the periodic coefficient mapping, pole
constraints, and cosine/sine spline collocation from `fpsphe` and `fprpsp` in
Federico Perini's modern Fortran FITPACK implementation:
[src/fitpack_core.F90](https://github.com/perazz/fitpack/blob/1907c7bad9f5c4557970e906bdd4909eb52b2d4c/src/fitpack_core.F90),
commit `1907c7bad9f5c4557970e906bdd4909eb52b2d4c`.
The exact source and `LICENSE.md` were inspected before adaptation.
Source SHA-256: `02e413063f7df855c45eed6119130f7ab7a07cfa3e2f6ac0e42b6f496df414ec`.
Upstream license SHA-256: `39fbeb302cce5182528712269832f82d356e684b9a3a11dc6949a79ef64781b8`.
Copyright (c) 2022, Federico Perini. All rights reserved. The source credits
Paul Dierckx's original FITPACK and retains BSD-3-Clause. The full license
terms are already retained, with whitespace differences, in
[`LICENSES/perazz-fitpack-BSD-3-Clause.txt`](LICENSES/perazz-fitpack-BSD-3-Clause.txt).

The established Fortran-source search identified perazz/FITPACK as suitable
for reusing this fixed-knot formulation. SciFort adapts the constraint
construction rather than importing the whole package or the large banded
fitter and its adaptive-knot/smoothing machinery. Longitude coefficients
repeat periodically; constant pole rows and adjacent rows expressed using
periodic cubic representatives of cosine/sine reproduce the upstream reduced
space. The result is converted to the standard longitude-fastest FITPACK tensor
coefficient layout.

Significant changes: use existing SciFort B-spline collocation/evaluation and
the native dense Jacobi SVD instead of FITPACK's banded Givens/`fprank` solver;
scale observations and weights; guard resource sizes and physical coefficient
rescaling; return a valid reduced-coordinate minimum-norm solution and explicit
rank warning when deficient. SVD cutoffs, rank-deficient extrapolations, dense
size limits, and wrapping query longitudes intentionally differ from FITPACK.
Public objects expose evaluation/derivatives and numerical diagnostics without
a C ABI. No adaptive knot placement, smoothing algorithm, or other upstream
routine is imported. Original tests and validation tools are MIT licensed.

Mathematical reference: P. Dierckx, *Curve and Surface Fitting with Splines*,
Oxford University Press, 1993, Chapter 11, Section 11.2, as cited by upstream
`sphere`. SciPy `LSQSphereBivariateSpline` is an independent numerical/API
reference. Existing SciFort numerical kernels and dependencies are unchanged.

## 2026-11-13 - Complex Bessel Y zeros

The complex-output overloads and private helpers in `src/special/bessel_zeros.f90`
are original MIT-licensed SciFort code, independently implemented from the
following public mathematical specifications:

- [DLMF 10.6](https://dlmf.nist.gov/10.6): `Y_0'=-Y_1` and
  `Y_1'=Y_0-Y_1/z`.
- [DLMF 10.2.1](https://dlmf.nist.gov/10.2.E1): the Bessel differential
  equation, used to evaluate `Y_1''` for Newton refinement of derivative roots.
- [DLMF 10.11.6](https://dlmf.nist.gov/10.11.E6): integer-order continuation
  across the negative-real cut, `Y_n(-w)=(-1)^n[Y_n(w)+2i J_n(w)]` on its upper side.
- [DLMF 10.17](https://dlmf.nist.gov/10.17): leading sine/cosine phases,
  combined with continuation to obtain real starting coordinates
  `-(k-1/4)*pi`, `-(k-3/4)*pi`, and `-(k-5/4)*pi` for the three root strings.
- [DLMF 10.21(ix)](https://dlmf.nist.gov/10.21#ix), especially 10.21.46:
  complex root geometry and the asymptote `Im(z)=log(3)/2`. The exceptional
  order-one roots start from coarse geometric guesses `-0.5+0.5i` and `0.5+1i`;
  these are not hard-coded reference roots or fitted numerical coefficients.

The solver reuses SciFort's existing complex `Y_0` and `Y_1` implementations,
with bounded Newton steps, backtracking, root-specific boxes, finite-value
checks, and both correction and residual convergence tests. No existing
numerical kernel was changed. The implementation search did not identify a
directly reusable established Fortran routine with verified suitable licensing
and behavior; no legacy specfun, SciPy, Zhang/Jin, or other zero-finder source
was copied or adapted. No dependency or third-party license was added.

SciPy's public API documentation and numerical values were inspected to check
the root sequence. Actual derivative values follow the DLMF identity, preserving
SciFort's established convention even where SciPy's numerical companion has
the opposite sign. Independent 80-digit mpmath solves start from mathematical
phase seeds, not SciPy/SciFort roots. The native tests and reference tools are
original SciFort code. C ABI expansion remains deferred.

## 2026-11-13 - Complex Riemann zeta

`src/special/complex_zeta.f90` adapts `zeta_Euler_Maclaurin` and its 25
Bernoulli coefficients from Alexey Kuznetsov's BSD-3-Clause Fortran implementation:
[Fortran/gamma_zeta_module.f90](https://github.com/Alexey-Kuznetsov-math/Riemann_zeta-and-Gamma/blob/1f1ad8934044cfb20d6d70137868ddd0c74652c4/Fortran/gamma_zeta_module.f90),
commit `1f1ad8934044cfb20d6d70137868ddd0c74652c4`.
The exact source and project license were inspected before adaptation.
Source SHA-256: `c4a8816cfb0f67bf709734000fa591b0a4bb5ec4008a1d3e4ec77749c7371315`.
License SHA-256: `9d515bd232c7caea94e38971f278bd3bdf4325fe0fd2c42218fe892336598820`.
Copyright (c) 2025, Alexey-Kuznetsov-math. Full terms are retained in
[`LICENSES/kuznetsov-zeta-BSD-3-Clause.txt`](LICENSES/kuznetsov-zeta-BSD-3-Clause.txt);
the adapted source retains BSD-3-Clause and upstream attribution.

Changes: replace quadruple precision with SciFort binary64; use compensated
direct summation instead of upstream `main_sums`; increase the shift padding
from 20 to 32; arrange the pole term separately to avoid an unnecessary large
intermediate; bound nonreal arguments and guard pole/result overflow. The
upstream high-frequency approximation, exponential-sum coefficients, Taylor
series, derivatives, gamma implementation, and STOP statements were not imported.
Continuation uses SciFort's existing complex log-gamma and a scaled logarithmic
sine, following [DLMF 25.4.2](https://dlmf.nist.gov/25.4.E2).
The adapted Euler-Maclaurin formula follows
[DLMF 25.2.10](https://dlmf.nist.gov/25.2.E10).

The established Fortran-source search identified this directly reusable
permissive kernel. SciPy's current XSF `zeta.h` was inspected for comparison
but was not adapted. SciPy and 90-digit mpmath values are independent test
references; all retained Bernoulli coefficients are verified against
`bernoulli(2*k)/factorial(2*k)` at high precision. Existing real zeta/Hurwitz
implementations remain unchanged. Tests and the generic forwarding overload
are original MIT-licensed SciFort code. No dependency or C ABI was added.

## 2026-11-13 - AAA rational approximation

`src/interpolate_aaa.f90` adapts SciPy **v1.18.1**
[`scipy/interpolate/_bary_rational.py`](https://github.com/scipy/scipy/blob/v1.18.1/scipy/interpolate/_bary_rational.py).
The exact tagged source, its Oxford/Chebfun BSD-3-Clause header, and the
[SciPy project license](https://github.com/scipy/scipy/blob/v1.18.1/LICENSE.txt)
were inspected before implementation. Downloaded source SHA-256:
`7d4cffe85836beec002e2889431370739620b78e47641b751dfaa3239c34c8ac`.
Copyright (c) 2017, The Chancellor, Masters and Scholars of the University
of Oxford, and the Chebfun Developers. All rights reserved.
Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
Both complete upstream BSD-3-Clause notices and terms are retained in
[`LICENSES/scipy-aaa-BSD-3-Clause.txt`](LICENSES/scipy-aaa-BSD-3-Clause.txt).
The adapted source carries MIT AND BSD-3-Clause and both copyright notices;
original SciFort modifications, tests, and reference drivers are MIT licensed.

Adapted operations: greedy support-point selection, Loewner-matrix SVD fits,
ill-conditioned-column equilibration, barycentric evaluation, arrowhead
generalized-eigenvalue pencils for finite poles/zeros, quotient-derivative
residues, and Froissart cleanup by removing nearby support points and refitting.
Mathematical reference: Y. Nakatsukasa, O. Sete, and L. N. Trefethen (2018),
"The AAA algorithm for rational approximation", SIAM J. Sci. Comput. 40(3),
A1494-A1522, [DOI:10.1137/16M1106122](https://doi.org/10.1137/16M1106122),
[author preprint](https://arxiv.org/abs/1612.00337).

Significant modifications: native real/complex Fortran objects and constructor
overloads, stable lexicographic input sorting, explicit status diagnostics,
finite-value filtering, affine coordinate/data scaling, guarded complex
division and rescaling, exact retention of supplied support values, scaled
barycentric evaluation, and normalized cleanup decisions. Reuse SciFort's
existing SVD and generalized-eigenvalue routines without new dependencies.
Drop projective roots within `64*epsilon` of infinity, polish finite roots,
and check normalized secular residuals. Respect the current order-40 eigenvalue
backend limit, reporting a size diagnostic beyond 39 retained support points.
Fit/evaluation remain available beyond that limit. Failed cleanup preserves
the fit; successful cleanup recomputes the actual sample error and convergence.
No raw upstream Python source, downloaded archive, or adaptation scratch files
are vendored.

A source search and the published
[rational-approximation software catalogue](https://rationalapproximation.github.io/)
did not identify a suitable established Fortran AAA implementation for reuse.
SciPy's exact permissively licensed implementation was selected; this is not
a claim that no other Fortran implementation exists. No third-party
nonpermissive implementation was used. The paper is cited as a mathematical
reference, not as an imported implementation source.

Native analytic tests and independent public SciPy comparisons are original
SciFort test code, not copied upstream test implementations. See
[`validation/aaa_2026-11-13.md`](validation/aaa_2026-11-13.md) and the
[API guide](docs/interpolate_aaa.md) for tested behavior and limitations.

## 2026-11-13 - Lebedev unit-sphere quadrature

`src/integrate/lebedev.f90` adapts all 32 coefficient tables and six symmetry
orbit expansions from SciPy **v1.18.1**
[`scipy/integrate/_lebedev.py`](https://github.com/scipy/scipy/blob/v1.18.1/scipy/integrate/_lebedev.py).
The exact tagged file and its complete BSD-2-Clause header were inspected.
The tagged [bundled-license inventory](https://github.com/scipy/scipy/blob/v1.18.1/LICENSES_bundled.txt)
explicitly identifies that file as BSD-2-Clause (`getLebedevSphere`).
Downloaded source SHA-256:
`9dfe2476c3cece2672694bbf9dabee87a3228eccaf123cbe15f8d184198dc4dc`.
Copyright (c) 2010, Robert Parrish. All rights reserved.
Complete required terms are retained in
[`LICENSES/scipy-lebedev-BSD-2-Clause.txt`](LICENSES/scipy-lebedev-BSD-2-Clause.txt).
The adapted source carries MIT AND BSD-2-Clause and source-file attribution;
original SciFort interfaces, tests, and validation drivers are MIT licensed.

Upstream attribution chain, preserved in the source and guide:
Robert Parrish's `getLebedevSphere` translated Christoph van Wuellen's Fortran
translation of Dmitri Laikov's original C routines. Requested reference:
V. I. Lebedev and D. N. Laikov, "A quadrature formula for the sphere of the
131st algebraic order of accuracy", Doklady Mathematics 59(3), 477-481 (1999).
The SciPy source also attributes its MATLAB origin to
[Parrish's published package](https://www.mathworks.com/matlabcentral/fileexchange/27097-getlebedevsphere).

Adaptation details: mechanically extracted the original decimal literals into
1,287 immutable orbit records; translated the explicit point expansion with
one-based indices while preserving its ordering and `4*pi` weight scaling.
Replaced the Python wrapper by a pure Fortran subroutine, allocatable outputs,
existing integration status codes, and allocation checks. No fitted, generated,
or rounded replacement coefficients were introduced. The Fortran translation
is lowercase, uses SciFort's `dp`, and has no I/O, global mutable state, or
runtime data-file dependency. Raw upstream source and build-time adaptation
scratch files are not committed. All coefficients are checked through rule
comparison and independent polynomial moments.

Established Fortran source reviewed first:
[CCL distribution](https://www.ccl.net/cca/software/SOURCES/FORTRAN/Lebedev-Laikov-Grids/index.shtml),
[exact source](https://www.ccl.net/cca/software/SOURCES/FORTRAN/Lebedev-Laikov-Grids/Lebedev-Laikov.F),
and [README](https://www.ccl.net/cca/software/SOURCES/FORTRAN/Lebedev-Laikov-Grids/README.TXT).
The README records permission to make the routines publicly available but does
not establish clear modification/redistribution terms for a new direct import.
That direct Fortran import is deferred pending explicit permissive permission;
no CCL source was directly imported. This implementation instead adopts the
explicitly BSD-2-Clause SciPy translation and retains its attribution chain.

Independent moment tests use the simplex beta integral in
[NIST DLMF 5.14.2](https://dlmf.nist.gov/5.14.E2), with squared Cartesian
coordinates giving Dirichlet parameters `(1/2,1/2,1/2)`. Native tests evaluate
finite product ratios; the Python oracle uses exact rational factorial ratios.
The installed SciPy 1.15.3 rules were also compared with the inspected 1.18.1
functions and were bitwise identical for all 32 rules. See
[`validation/lebedev_2026-11-13.md`](validation/lebedev_2026-11-13.md).

## 2026-11-12 - Boundary-value ODE collocation (SciPy)

`src/integrate/bvp.f90` adapts SciPy **v1.15.3**
[`scipy/integrate/_bvp.py`](https://github.com/scipy/scipy/blob/v1.15.3/scipy/integrate/_bvp.py).
The exact tagged source and [license](https://github.com/scipy/scipy/blob/v1.15.3/LICENSE.txt)
were inspected; installed SciPy 1.15.3 matched after line-ending normalization.
SHA-256 of the downloaded tagged source:
`d046b344a10268e12b61efcc0006e6811adb91d3a0497a709107f03cbc4fdf42`.
Copyright (c) 2001-2002 Enthought, Inc. 2003-2024, SciPy Developers.
Full BSD-3-Clause terms are retained in
[`LICENSES/scipy-bvp-BSD-3-Clause.txt`](LICENSES/scipy-bvp-BSD-3-Clause.txt).
The adapted source carries MIT AND BSD-3-Clause attribution; original SciFort
interfaces, modifications, tests, and validation drivers are MIT licensed.

Adapted operations are fourth-order Lobatto IIIA midpoint collocation,
analytic global Jacobian blocks, forward-difference Jacobians, damped Newton
line search, singular-term pseudoinverse projection, five-point Lobatto
relative-defect quadrature, and two-/three-way adaptive subdivision.
Mathematical reference: J. Kierzenka and L. F. Shampine (2001),
"A BVP Solver Based on Residual Control and the Matlab PSE", ACM Transactions
on Mathematical Software 27(3), 299-316,
[DOI:10.1145/502800.502801](https://doi.org/10.1145/502800.502801).
No article implementation code was imported.

Significant modifications: native real/complex callback and result types,
zero-length parameter conventions, direct CSC block assembly, row scaling,
existing SciFort sparse LU/SVD reuse, Jacobian refresh on every Newton step
(up to eight steps per mesh), real squared-magnitude complex line search,
explicit argument/nonfinite/size/stagnation checks and statuses, independent
cubic dense storage and derivative evaluation. No upstream Python source file
is vendored. Validation uses public mathematical problems and independent
SciPy calls, not copied upstream test implementation.
See `validation/bvp_2026-11-12.md`.

Established Fortran alternative reviewed before adoption:
[Paul Muir's BVP_SOLVER page](https://cs.smu.ca/~muir/BVP_SOLVER_Webpage.shtml),
[exact BVP_M-2.f90](https://cs.smu.ca/~muir/BVP_SOLVER-2_Files/BVP_M-2.f90),
and [linear-algebra license notes](https://cs.smu.ca/~muir/BVP_LA_LicenseInfo.txt).
The BVP module itself includes BSD-3-Clause terms. Its mutable module state,
separated boundary-condition interface, and separately sourced COLROW/LINPACK
linear algebra would require substantially more refactoring and additional
license resolution. It was not imported; this is not a claim that its own
BSD license is incompatible. SciFort's existing sparse LU avoids importing
that separate linear-algebra package.

## 2026-11-12 - Fast Hankel transforms (FFTLog)

`src/fftlog.f90` and the coefficient-magnitude calculation in
`validation/check_fftlog.py` adapt SciPy **v1.15.3**
[`scipy/fft/_fftlog_backend.py`](https://github.com/scipy/scipy/blob/v1.15.3/scipy/fft/_fftlog_backend.py).
The exact installed 1.15.3 source was inspected, together with its
[tagged license](https://github.com/scipy/scipy/blob/v1.15.3/LICENSE.txt).
The installed source matched the tagged upstream file after line-ending
normalization. SHA-256 of the downloaded tagged source:
`520a1e3f08683282f1bd0e6ca125199155af5964cfef2d71595003f4196b7496`.
SciPy copyright: Copyright (c) 2001-2002 Enthought, Inc. 2003-2024,
SciPy Developers. The full BSD-3-Clause terms are retained in
[`LICENSES/scipy-fftlog-BSD-3-Clause.txt`](LICENSES/scipy-fftlog-BSD-3-Clause.txt).
Adapted files carry MIT AND BSD-3-Clause SPDX identifiers and attribution;
original SciFort modifications and tests are MIT licensed.

Adapted operations: logarithmic input/output biasing, complex log-gamma
Mellin coefficients, real Nyquist projection, reciprocal-conjugate inverse
weights, reversal, and nearest low-ringing offset. Mathematical specification:
A. J. S. Hamilton (2000), MNRAS 312, 257-284,
[Appendix B](https://arxiv.org/abs/astro-ph/9905191), also published as
[FFTLog sections 6-9](https://jila.colorado.edu/~ajsh/FFTLog/).
The original algorithm is attributed there to J. D. Talman (1978),
J. Comput. Phys. 29, 35-48.
The removable simultaneous-pole limit uses the gamma residues specified in
[DLMF 5.2(i)](https://dlmf.nist.gov/5.2#i): their ratio gives the signed
factorial ratio used by the zero-frequency coefficient.

Significant changes: use SciFort's existing FFT plans and complex log-gamma
kernels; normalize spectra and restore amplitudes in logarithmic form;
catch setup, singularity and range errors without library output or STOP;
evaluate simultaneous gamma poles through their residue ratio; preserve
well-defined zero-DC directions; reject singular directions instead of
dropping their constant mode; reject a nearly zero even Nyquist projection;
add rank 1-15 batches with selectable axes and per-call plan/weight reuse.
No NCAR FFT, Ooura gamma, Python runtime, or new dependency is imported.

Validation uses public SciPy calls, 90-digit direct Mellin/DFT sums, and the
four historical Fortran-output fixtures in installed SciPy's
[`test_fftlog.py`](https://github.com/scipy/scipy/blob/v1.15.3/scipy/fft/tests/test_fftlog.py).
The fixtures are read from the installed test source rather than vendored.
Wide bias spans amplify errors in the sequence of scalings and FFTs;
validation reports both raw errors and errors relative to a factor bound.

### Documented issue: established Fortran FFTLog import deferred

Inspected `emsig/fftlog` commit
`8475c0fb40f6b0f55fa492b2db6a786464f7635e`:
[`src/fftlog.f`](https://github.com/emsig/fftlog/blob/8475c0fb40f6b0f55fa492b2db6a786464f7635e/src/fftlog.f),
[`README.rst`](https://github.com/emsig/fftlog/blob/8475c0fb40f6b0f55fa492b2db6a786464f7635e/README.rst),
and [LICENSE](https://github.com/emsig/fftlog/blob/8475c0fb40f6b0f55fa492b2db6a786464f7635e/LICENSE).
The README applies CC0 to wrapper additions and describes permission to
distribute the original Fortran with that package. It does not clearly grant
general adaptation rights for the original Hamilton implementation.
Import was therefore stopped; no Fortran source from that repository is
included. Future adoption requires a clear permissive license covering the
original code and any retained dependencies. The SciPy BSD adaptation above
is the implementation source for this change.

## 2026-11-12 - FFT helpers, multidimensional shifts, and 2D names

The additions to `src/fft.f90`, `test/test_fft_helpers.f90`,
`validation/fft_nd_reference_driver.f90`, and
`validation/check_fft_helpers.py` are original MIT-licensed SciFort code.
No upstream implementation, coefficients, dependency, binary, or archive
was imported. The FFT implementation and existing notices are unchanged.

Public specifications inspected 2026-11-12:

- [SciPy next_fast_len](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.next_fast_len.html)
  and [prev_fast_len](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.prev_fast_len.html):
  bounding smooth-number searches and backend-dependent radix policies.
- [SciPy fftshift](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.fftshift.html)
  and [ifftshift](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.ifftshift.html):
  selected-axis half-space permutations and odd/even inverse relationships.
- [SciPy fft2](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.fft2.html)
  and [hfft2](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.hfft2.html):
  two-axis convenience conventions, inherited N-D definitions and half-axis rules.

Fast lengths are independently enumerated as nonnegative powers of 2/3/5/7,
matching the exact radices inspected in SciFort's existing
`signal_spectral.f90:mixed_radix_compatible`. The real flag deliberately
selects the same policy because both transform paths use the same backend.
Native tests independently factor each integer for small-range searches;
the Python oracle builds the full representable smooth set with a heap.
Shift source indices follow modular coordinate addition in Fortran element
order. The 2D wrappers delegate to the existing multidimensional kernels.

## 2026-11-12 - Multidimensional FFT, DCT, and DST APIs

The multidimensional additions to `src/fft.f90`,
`test/test_fft_nd.f90`, `validation/fft_nd_reference_driver.f90`, and
`validation/check_fft_nd.py` are original MIT-licensed SciFort code.
They implement separable transforms using the existing SciFort 1D kernels;
no upstream implementation code, coefficients, binary, archive, or dependency
was imported. Existing kernel attribution and license notices remain intact.

Public specifications inspected 2026-11-12:

- [SciPy fftn](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.fftn.html):
  separability, axis selection, shape, and normalization conventions.
- [SciPy rfftn](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.rfftn.html)
  and [irfftn](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.irfftn.html):
  half-spectrum dimension, stage order, and inferred full length.
- [SciPy hfftn](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.hfftn.html)
  and [ihfftn](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.ihfftn.html):
  real-spectrum Hermitian pair and normalization.
- [SciPy dctn](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.dctn.html):
  separable application of the 1D cosine transforms.
  Sine transforms use the previously cited DST sums along each selected axis.

The native tests evaluate a double DFT sum independently. External checks
use public SciPy transforms, with real-component comparisons for complex
cosine/sine transforms as described below. The native interface deliberately
uses one-based positive axes and caller-owned output arrays.

## 2026-11-12 - Standalone FFT, DCT, and DST APIs

`src/fft.f90`, `test/test_fft.f90`, `validation/fft_reference_driver.f90`, and
`validation/check_fft.py` are original MIT-licensed SciFort code, implemented
from public mathematical specifications. No implementation code, numerical
coefficients, archive, binary, or dependency was imported. The existing
`scifort_signal_spectral.fft_plan` is reused unchanged; its mixed-radix,
Bluestein, and radix-2 provenance and notices remain applicable.

Specifications inspected 2026-11-12:

- [SciPy FFT](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.fft.html):
  defining DFT, padding/cropping, inverse sign, and normalization conventions.
- [SciPy DCT](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.dct.html)
  and [DST](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.dst.html):
  type I-IV sums, logical sizes, inverse relationships, and endpoint weights.
- [SciPy Hermitian FFT](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.hfft.html):
  half-sequence interpretation, inferred length, and odd/even endpoint rules.
- [SciPy frequency bins](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.fftfreq.html):
  signed bin ordering and cycles-per-spacing units.

Even/odd FFT extensions are derived algebraically from those sums. DCT-I
uses length `2*(n-1)` and DST-I `2*(n+1)`; types II/III use length `4*n`, and
type IV uses `8*n`. Cosine extensions are even and sine extensions odd.
For type II, input positions are `2*j+1`; for type III they are `j` (cosine)
or `j+1` (sine), with the endpoint half-weight required by the defining sum.
Type IV samples odd output frequencies from odd input positions. Symmetric
endpoints and orthogonalization weights follow the definitions, not copied
upstream code. Complex transforms extend the same linear real definitions.
Inverse Fourier transforms use conjugation and the documented normalization.

The installed SciPy 1.15.3 `_pocketfft.realtransforms._r2r` source was inspected
only to confirm that its complex path omits the `orthogonalize` argument.
No source was adopted. Reference checks therefore call the public real
transform independently for each complex component. mpmath supplies bin
references when the tested frequency helper's product `n*d` overflows;
other external transform references use SciPy. Native tests independently
evaluate the defining matrices. Scope expansion to a standalone FFT API is
recorded in the roadmap in response to the user's explicit FFT request.

## 2026-11-12 - Censored Gumbel extension

The Gumbel additions to `src/censored.f90` and `src/fit.f90`, their tests,
and reference scripts are independent MIT-licensed SciFort implementations
from the mathematical definitions below. No upstream implementation source,
coefficients, dependency, archive, or binary was imported. Existing SciFort
kernels retain their notices and provenance.

Specifications inspected 2026-11-12:

- [SciPy gumbel_r](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gumbel_r.html)
  and [gumbel_l](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gumbel_l.html):
  defining densities, support, and location/scale parameterizations.
- [SciPy CensoredData](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.CensoredData.html):
  exact, left, right, and interval observations and endpoint conventions.

Integrating the right-skewed density gives `F(z)=exp(-exp(-z))`.
For standardized endpoints `za < zb`, algebra gives
`F(b)-F(a)=F(b)*(1-exp(-delta))`, where
`delta=exp(-za)*(1-exp(-(b-a)/scale))`.
SciFort evaluates this identity with the existing log-width and exponential
mass helpers, keeping the increment in the log domain. Left-skewed kernels
use reflection, reversing endpoints and interchanging tails without copying
or modifying observations. The prior private `logistic_z` helper is renamed
`location_scale_z` and reused with its numerical algorithm unchanged.
The censored fitting overloads reuse the existing bounded local search.

Reference likelihoods use independently evaluated mpmath defining formulas
and high-precision tail subtraction; SciPy is used only for external fit
comparisons. No public Fortran API is removed or renamed, and no C ABI symbol
is added. Optimizer convergence does not establish global optimality,
identifiability, or existence of an unconstrained finite MLE.

## 2026-11-12 - Censored logistic extension

The logistic additions to `src/censored.f90` and `src/fit.f90`, their tests,
and reference scripts are independent MIT-licensed SciFort implementations
from the mathematical definitions below. No upstream implementation source,
coefficients, dependency, archive, or binary was imported. Existing kernels
retain their provenance and notices.

Specifications inspected 2026-11-12:

- [SciPy logistic distribution](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.logistic.html):
  standardized density and location/scale transformation.
- [SciPy CensoredData](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.CensoredData.html):
  exact, left, right, and interval observation semantics, including infinite
  endpoints and equal endpoints representing exact observations.

Integrating the defining density gives `F(z)=1/(1+exp(-z))`.
Algebra then gives `F(b)-F(a)=F(b)*S(a)*(1-exp(-(b-a)/scale))`.
The implementation evaluates this identity with log tails and the existing
log-width exponential-mass helper, without subtracting nearly equal CDFs.
Standardization rescales displacements near subtraction overflow and divides
bounded significands with a checked binary exponent. The bounded fitter reuses SciFort's existing local search.
Reference values use independent mpmath defining formulas; SciPy is used only
as an external fit comparison. No guarantee of a global optimum, finite MLE,
or parameter identifiability follows from optimizer convergence.

## 2026-11-12 - Censored beta extension

The beta additions to `src/censored.f90` and `src/fit.f90`, with tests and
reference checks, are MIT-licensed SciFort implementations from mathematical
specifications. They reuse existing incomplete-beta, log-beta, gamma-difference,
and scaled gamma-ratio kernels with their notices unchanged. No third-party
implementation source, dependency, binary, or archive was imported.
Specifications inspected 2026-11-12:

- Beta density and parameter convention:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.beta.html
- Defining integral, symmetry, and series: https://dlmf.nist.gov/8.17
- Gauss quadrature: https://dlmf.nist.gov/3.5.v
- Stirling gamma-ratio limit: https://dlmf.nist.gov/5.11
- FastTwoSum error-free sum, Algorithm 1.1 in T. Ogita, S. M. Rump, and
  S. Oishi, "Accurate Sum and Dot Product", SIAM Journal on Scientific
  Computing 26(6), 1955-1988 (2005):
  https://www.tuhh.de/ti3/paper/rump/OgRuOi05.pdf
- Censored likelihood specification:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rv_continuous.fit.html

The implementation of the upper-edge correction follows the cited arithmetic
specification with operands ordered by magnitude. Endpoint distances are
normalized independently. Narrow intervals integrate the lower-endpoint
factored density with log ratio `(a-1)*log1p(t*ql)+(b-1)*log1p(-t*qr)`.
Three-point Gauss-Legendre nodes and weights follow polynomial exactness
through degree five. Conservative relative-distance/variation bounds of .001
select this rule. Other intervals use log tails. At an underflowed coordinate,
`I_z(a,b) = z**a/(a*B(a,b)) * (1 + O(max(1,b)*z))` is evaluated in logs,
with a far-shape/coordinate bound before using the approximation. Existing
gamma-difference and scaled-ratio identities retain tiny-shape normalization.
For a tiny first shape and a far shape above `1/epsilon`, the Stirling
limit `log(Gamma(b+a)/Gamma(b)) = a*log(b)+O(a/b)` avoids powers of an
enormous argument. The neglected term is below roundoff under that branch's
`a < sqrt(epsilon)` and underflow-tail conditions.
The associated one-line fix in `src/special/log_gamma.f90` evaluates negative
powers as positive powers of the reciprocal, avoiding intermediate overflow
in compiler integer-power routines without changing the Stirling coefficients.
No C ABI symbols are added.

## 2026-11-12 - Censored gamma extension

Gamma additions in `src/censored.f90` and `src/fit.f90`, with tests and reference
checks, are original MIT-licensed SciFort code. No upstream implementation code
or dependencies were imported. Existing gamma, incomplete-gamma, logarithmic
Gamma, and arithmetic helpers are reused with their notices unchanged.
Public mathematical specifications inspected 2026-11-12:

- Gamma density, shape/scale convention, and exponential identity:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gamma.html
- Defining incomplete-gamma integrals: https://dlmf.nist.gov/8.2
- Small-argument series: https://dlmf.nist.gov/8.7
- Gauss quadrature: https://dlmf.nist.gov/3.5.v
- Stirling form used for large-shape small-argument kernels: https://dlmf.nist.gov/5.11
- Censored likelihood contributions:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rv_continuous.fit.html

The narrow-interval rule factors out the lower-endpoint gamma density and
integrates `exp((a-1)*log1p(t*q)-t*d)` for `t` in `[0,1]`. Three-point
Gauss-Legendre nodes/weights follow from polynomial exactness through degree
five. Conservative bounds `q*max(1,abs(a-1)) <= .001`, `d <= .001` limit
endpoint separation and density variation. At scaled arguments below the
smallest normal real, `log(P(a,z)) = a*log(z)-log(Gamma(a+1))+O(z)` is used;
the omitted relative correction is bounded by `z`. Other intervals reuse
existing log tails. These algebraic evaluations are independent derivations,
not adaptations of a third-party implementation. No C ABI symbols are added.

## 2026-11-12 - Censored lognormal extension

The lognormal additions to `src/censored.f90` and `src/fit.f90`, their native
tests, and high-precision reference checks are original MIT-licensed SciFort
code derived from the public mathematical specifications:

- https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.lognorm.html
- https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rv_continuous.fit.html

Both were inspected 2026-11-12. No upstream implementation was imported.
The change of variable `y=log(x-loc)` has normal mean `log(scale)` and
standard deviation `s`; exact densities include Jacobian `1/(x-loc)`.
Interval probabilities reuse SciFort's three-point normal quadrature, with
width `log1p((b-a)/(a-loc))` evaluated through log distances. Existing notices
and licenses are retained. No dependency or C ABI change is introduced.

## 2026-11-12 - Censored likelihoods and distribution fitting

`src/censored.f90`, the censored fitting additions in `src/fit.f90`, extensions
to the existing `censored_data` in `src/empirical_distribution.f90`, and their
tests and documentation are original SciFort MIT-licensed implementations.
No third-party source was copied, translated, or imported, and no dependency
was added. Existing notices on reused distribution kernels remain unchanged.

Mathematical and behavioral specifications:

- SciPy's public censored log-likelihood definition, inspected 2026-11-12:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rv_continuous.fit.html
- SciPy's public observation/endpoint conventions:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.CensoredData.html

The analytic fixed-location exponential estimate follows by differentiating
`-d*log(scale)-exposure/scale`, giving `scale=exposure/d`. Interval calculations
use algebraic log-difference and cumulative-hazard identities. Narrow normal
intervals use three-point Gauss-Legendre quadrature, specified by NIST DLMF
https://dlmf.nist.gov/3.5.v . The nodes `0,+/-sqrt(3/5)` and weights
`8/9,5/9,5/9` follow by exact integration of polynomials through degree five.
The standardized width bound of `0.01/max(1,abs(z))` makes the sixth-order
remainder negligible at binary64 precision. Large Weibull hazard-ratio differences omit a relative
term below `exp(-50) < 2e-22`; small-increment shortcuts use an epsilon-based
bound. These are original algebraic evaluations, not upstream adaptations.

`tools/generate_censored_reference.py` evaluates independent mathematical
formulas with mpmath and compares fits with SciPy. Initial native reference
constants were generated with mpmath 1.3.0 (100 to 700 decimal digits) and
SciPy 1.15.3. This numerical reference use does not import implementation code.

## 2026-11-11 - Repository refresh from the sibling working tree

The refreshed project files were copied from `../scifort-main`, based on Git
commit `fca1668e08c7278263ce63fa4365013e6c795e29`, including that working tree's
uncommitted streaming spectrogram default-window fix. Existing implementation
attribution, copyright notices, full license texts, and dependency declarations
were retained. This refresh introduces no separately sourced implementation.

Local adjustments replace non-ASCII Frechet/Matern text in a comment in
`src/linalg_expm_frechet.f90` and a comment/status message in
`test/test_interpolate_rbf.f90`. Numerical expressions and license notices in
these files are unchanged. The refresh is also recorded in `CHANGELOG.md`.

## 2026-11-11 - Streaming spectrogram default window

The omitted-window branch in `src/signal_stft_stream.f90` now supplies alpha
0.25 to SciFort's existing Tukey window implementation, matching the existing
batch spectrogram. The public default is documented by SciPy at
https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.spectrogram.html.
This is an original SciFort default-selection correction, with native and
independent SciPy regression checks. No upstream implementation code was
copied, no coefficients were introduced, and existing license notices remain.

## 2026-11-11 - Reference-driver relocation

Twenty comparison producers formerly under `test/` now live under
`validation/reference_drivers/` to keep `fpm test` unattended. This is a
source-preserving relocation: existing authorship, attribution, license notices,
and input/output protocols are unchanged. The directory README identifies the
driver groups and their comparison workflows. No implementation code or new
third-party dependency was imported.

## 2026-11-10 — Mixed-integer linear programming (`scipy.optimize.milp`)

The module `src/optimize/milp.f90`, its native test, standalone Python
validation scripts, example, and documentation are newly written original
MIT-licensed SciFort contributions. They implement the textbook mathematical
branch-and-bound method using lower bounds inherited from parent LP
relaxations, best-bound node selection, and disjunctions for integer and
semi-continuous variables. SciFort's existing independently implemented
`scifort_optimize_linprog` two-phase simplex code is reused, with a small
original minimum-ratio tie-tolerance bug fix for large RHS magnitudes.
No source was copied from SciPy, HiGHS, or an external MILP project.

SciPy's public API, variable type meanings, and result conventions were
consulted at <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.milp.html>.
SciPy 1.17.0/HiGHS was used only for independent solution comparisons.
Additional verification enumerates all assignments in bounded integer
reference models. This is not a sparse or performance-equivalent replacement
for SciPy's production HiGHS solver; full FPM integration is still pending.

## 2026-11-10 — Quadratic assignment (FAQ and 2-opt)

The new module `src/optimize/quadratic_assignment.f90`, Fortran/Python tests,
example program, and documentation are original SciFort MIT-licensed code.
They independently implement a Frank-Wolfe relaxation with quadratic line
search, final linear-assignment rounding, and best-improvement permutation
swaps with analytically calculated O(n) swap deltas. The existing SciFort
Hungarian `linear_sum_assignment` is reused; no upstream SciPy source code
was copied, adapted, or vendored. No new external dependency or third-party
license is required.

Mathematical reference: Vogelstein et al., *Fast Approximate Quadratic
Programming for Graph Matching* (2015), DOI 10.1371/journal.pone.0121002.
Interface and independent numerical reference: SciPy's published
`quadratic_assignment` FAQ/2-opt documentation:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.quadratic_assignment.html
(consulted 2026-11-10). The code differs intentionally in 2-opt
initialization and iteration accounting; Fortran indices are one-based.
SciPy 1.17.0 served only for regression comparisons. No provenance statement
here implies global optimality, full API parity, or successful FPM linkage.

## 2026-11-10 — Dual annealing (`scipy.optimize.dual_annealing`)

New independent MIT-licensed SciFort source modules:
`src/optimize/dual_annealing.f90` and
`src/optimize/dual_annealing_minimize.f90`. The generalized simulated-annealing
mathematical procedure derives from Tsallis (1988), Xiang, Sun, Fan and Gong
(1997), and Xiang and Gong (2000). SciPy's public API and mathematical
expressions were consulted at
https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.dual_annealing.html
and source https://github.com/scipy/scipy/blob/main/scipy/optimize/_dual_annealing.py .
SciPy is BSD-3-Clause licensed. All new Fortran implementation source is
independently written; **no SciPy or other external source code was copied**.
The existing SciFort L-BFGS-B and Powell wrappers are reused without editing
or adding dependencies. SciPy and SciFort PRNG streams are deliberately
independent. Reference tests are in
`test/optimize_dual_annealing_reference.f90` and
`validation/compare_optimize_dual_annealing_scipy.py`. The mock optimizer in
`validation/optimize_dual_annealing_test_minimize.f90` is for an isolated
compiler contract test, not part of the distributed library numerical source.

## 2026-11-10 — Basin hopping (`scipy.optimize.basinhopping`)

`src/optimize/basinhopping.f90` and `src/optimize/basinhopping_minimize.f90`
are independently implemented MIT-licensed SciFort code. The mathematical
basin-hopping/Metropolis procedure traces to Wales and Doye (1997), and
SciPy's documented public behavior is used as an API and numerical testing
reference: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html .
No SciPy or external Fortran source code was copied. The implementation
reuses SciFort's existing local optimization module through a separate
adapter. The reference driver, tests, executable example and documentation
are original MIT-licensed contributions. The test-only substitute module in
`validation/optimize_basinhopping_test_minimize.f90` is not an FPM dependency
and must never be linked into the production library.

## 2026-11-10 — Dense linear programming

`src/optimize/linprog.f90` and its tests, example and comparison harness are
original MIT-licensed SciFort contributions implementing the standard
textbook two-phase primal tableau simplex method with Bland pivot selection.
No source was copied from SciPy, HiGHS, or any external Fortran repository.
We inspected the MIT-licensed `thchang/DualSimplex` project as a potential
backend: https://github.com/thchang/DualSimplex . That implementation
requires LAPACK/BLAS and targets a specific dense dual formulation; it was
not incorporated. SciPy 1.17.0 with HiGHS is used solely as an independent
numerical and status reference, with HiGHS presolve disabled in the randomized
comparison script to avoid a demonstrated presolve false-infeasibility
classification on a feasible, unbounded floating-point problem.

## 2026-11-10 — Bounded and nonnegative least squares

- `src/optimize/linear_least_squares.f90`: independently implemented active
  set solver inspired by the published Lawson–Hanson/BVLS mathematical method.
  It reuses SciFort's original Jacobi-SVD least-squares code. MIT; no foreign
  source copied.
- `src/optimize/linear_least_squares_slsqp.f90`: original SciFort adapter to
  `bvls_module::bvls` provided by the existing FPM dependency
  `https://github.com/jacobwilliams/slsqp`, pinned to tag `1.6.0` in
  `fpm.toml`. The upstream package is BSD-licensed and contains work by
  C. L. Lawson, R. J. Hanson, and Jacob Williams. It is dynamically resolved
  as an FPM dependency rather than vendored.
- The adapter source itself was compiled with a source-compatible mock only;
  actual dependency verification awaits a full FPM build.

## 2026-11-10 — Differential evolution in `scipy.optimize`

`src/optimize/differential_evolution.f90`, the native tests, reference driver,
SciPy comparison harness, example, and documentation are new MIT-licensed
SciFort code. The standard Storn–Price DE/rand/1/bin donor formula and
per-generation dithering were informed by Ian Morrall's MIT-licensed Fortran
example at https://github.com/ian-mmm/differential-evolution_f95 ; see his
copyright notice and full original MIT text in
`LICENSES/ian-morrall-differential-evolution-MIT.txt`. No upstream source was
copied verbatim. SciFort provides local reproducible random draws, correct
distinct donor sampling (the upstream implementation is not used directly),
12 standard donor/crossover strategies, per-coordinate bounds repair,
initialization, evaluation accounting, and SciPy-style convergence reporting.
StochOptim, https://github.com/keurfonluu/StochOptim (MIT, Keurfon Luu), was
also reviewed as an alternative, but its MPI-oriented framework was not
incorporated. SciPy 1.17.0 was used as an independent optimization benchmark,
not as an implementation source. No new runtime libraries or FPM dependencies.

## 2026-11-10 — Standalone matrix exponential

`src/linalg_expm.f90`, its native tests, reference driver, independent
SciPy comparison script, executable example, and documentation are original
MIT-licensed SciFort contributions. The method is the standard
[13/13] Pade rational approximation with scaling and squaring (Higham,
*Functions of Matrices*, 2008) and reuses the existing SciFort pivoted LU
implementation. No SciPy source code, third-party package, or binary library
was copied. SciPy 1.17.0 was used as an independent numerical oracle:
https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.linalg.expm.html .

## 2026-11-10 — Matrix-free randomized interpolative SVD

`src/linalg_interpolative_svd_operator.f90`, its tests, Python comparison
harness, reference driver, example and documentation are original
MIT-licensed SciFort contributions. The new module composes the already
implemented SciFort matrix-free ID with the standard QR/SVD identity for
converting a skeleton/interpolation matrix into singular triplets. No
SciPy, LAPACK or third-party code was copied or vendored. SciPy's public
`scipy.linalg.interpolative.svd(LinearOperator, ...)` API and independently
computed NumPy/SciPy numerical results were used for validation only:
https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.linalg.interpolative.svd.html .

## 2026-11-10 — Matrix-free randomized interpolative decomposition

`src/linalg_interpolative_operator.f90`, its native test, reference driver,
Python comparison harness, example, and documentation are original
independent implementations using randomized sketching, pivoted Householder
QR, and least-squares projection identities. Existing MIT-licensed SciFort
QR, SVD and interpolative-decomposition primitives are reused. SciPy's
public API and mathematical conventions were consulted for compatibility;
no SciPy source code or external third-party source was copied or vendored.

## 2026-11-10 — Randomized interpolative decomposition

`src/linalg_interpolative_random.f90`, the new tests, driver, example, and
API/validation documents are original MIT-licensed SciFort contributions.
They implement standard Gaussian row sketching, pivoted QR, and least-squares
interpolation using existing SciFort numerical modules. No SciPy or external
Fortran implementation source was incorporated. SciPy's public API at
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.interpolative.interp_decomp.html
and SciPy 1.17.0 numerical results were used for validation only. This
addition does not add runtime dependencies.

## 2026-11-10 — Randomized interpolative SVD

`src/linalg_interpolative_svd_random.f90`, its test drivers, comparison
script, example, and documentation are original MIT-licensed SciFort
contributions. The implementation uses standard Gaussian random projection,
column-pivoted QR, an original-matrix least-squares correction, and the
existing SciFort `id_to_svd` function. No SciPy or third-party implementation
source is incorporated. Only SciPy's public randomized-ID API and numerical
outputs were used for comparison:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.interpolative.svd.html .
SciPy 1.17.0 and NumPy 2.3.5 are test-only dependencies. The new module
has no runtime dependencies beyond SciFort itself.

## 2026-11-09 — Interpolative SVD (`scipy.linalg.interpolative.svd`)

`src/linalg_interpolative_svd.f90`, its native and independent reference tests,
example, and documentation are original MIT-licensed SciFort contributions.
The algorithm is the mathematical composition of the existing SciFort
`interp_decomp` and `id_to_svd` routines; no SciPy or third-party
implementation source was incorporated. SciPy's public deterministic
SVD interface was consulted at
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.interpolative.svd.html .
SciPy 1.17.0 and NumPy 2.3.5 are independent test references, not runtime
dependencies. The new interface does not implement randomized sampling.

## 2026-11-09 — Spectral-norm estimation from matrix-vector products

`src/linalg_spectral_norm.f90`, associated tests, example and documentation
are original MIT-licensed SciFort contributions. The numerical method is the
standard power iteration on `A^H A`, with vector rescaling between products
and one final `||A*x||` evaluation. No SciPy source code or external package
was incorporated. Public SciPy mathematical and API conventions were checked
against its documentation at
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.interpolative.estimate_spectral_norm.html
and
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.interpolative.estimate_spectral_norm_diff.html.
SciPy 1.17.0 is used solely in an independent reference-test script.

## 2026-11-09 — Clarkson-Woodruff (CountSketch) matrix sketching

`src/linalg_clarkson_woodruff.f90` and its tests and example are original
MIT-licensed SciFort implementations of the standard signed-hashing matrix
sketch (`S(h(i),i)=sign(i)`), without importing or adapting SciPy or any
third-party implementation code. Public SciPy API behavior and mathematical
conventions were checked against:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.clarkson_woodruff_transform.html .
The test harness uses SciPy 1.17.0 as an independent numerical oracle and
reads SciPy's private `cwt_matrix` only to obtain identical random maps for
cross-language testing. This test-only use introduces no runtime dependency
in the Fortran library.

## 2026-11-09 — Interpolative decomposition (`scipy.linalg.interpolative`)

`src/linalg_interpolative.f90` and its tests, example, and documentation are
original MIT-licensed SciFort contributions. The implementation derives from
the standard column-pivoted QR identity `A Pi = Q [R11 R12; 0 R22]`,
solves `R11*T = R12` for the interpolation coefficients, and converts
`Aapprox = B*P` to SVD through QR(B) and SVD of the small core.
No SciPy implementation source or external Fortran package was copied.
Behavior and interfaces were checked against the public SciPy documentation:
https://docs.scipy.org/doc/scipy/reference/linalg.interpolative.html .
Independent numerical comparison uses SciPy 1.17.0 only as a reference.

## 2026-11-09 — Hermitian pseudoinverse (`pinvh`)

`src/linalg_pinvh.f90` implements the standard spectral pseudoinverse
identity `A+ = Q diag(1/lambda for retained lambda) Q^H` using the already
existing self-contained SciFort `eigh`. The implementation and all tests
are original MIT-licensed project contributions, without copying SciPy
source code or incorporating external Fortran code. SciPy API conventions
were checked against the public documentation at
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.pinvh.html .

## 2026-11-09 — Real block eigenvalue conversion (`cdf2rdf`)

`src/linalg_cdf2rdf.f90`, its two native tests, independent reference
protocol and comparison script, example and documentation are original
MIT-licensed SciFort code. The conversion is independently derived by
separating real and imaginary parts of the eigenvalue equation
`A(a+ib)=(alpha+i beta)(a+ib)` for a real matrix A, yielding the
real invariant basis `[-b, a]` and its real 2x2 eigenvalue block.

The behavior and shape conventions were checked against the public
SciPy documentation:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cdf2rdf.html .
SciPy 1.17.0 was used as an independent validation oracle. No SciPy,
LAPACK, Netlib, or other third-party implementation source was copied
or adapted, and no new license or runtime dependencies were introduced.

## 2026-11-09 — `linalg_lu`, `linalg_det_inv`

Independent implementation of standard partial-pivoting Gaussian elimination,
permutation reconstruction, determinant-from-LU, and inverse-by-LU-solves.
The public interfaces and edge-case conventions were cross-checked against
SciPy 1.17.0. No external source code was copied or vendored; MIT license.

## 2026-11-09 — RQ decomposition

`src/linalg_rq.f90`, its tests, independent reference harness,
example and documentation are original MIT-licensed SciFort additions.
They use the elementary relationship between RQ and QR obtained by
reversing the rows/columns of a conjugate-transposed matrix; no SciPy
implementation source or third-party code was copied. Public SciPy
API documentation was consulted for output shapes, modes and batch cases:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.rq.html.
SciPy 1.17.0 served only as an independent validation reference.

## 2026-11-09 — QR multiplication and reusable Householder application

`src/linalg_qr_multiply.f90`, its native tests, reference protocol and
comparison script, example and documentation are original MIT-licensed
SciFort additions using the elementary Householder identity
`H = I - tau*v*v^H` and the existing SciFort QR factors. They neither
copy SciPy's source code nor vendor third-party implementation code.
Public SciPy documentation was consulted to confirm the side, pivoting,
conjugation, and output shape contracts:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.qr_multiply.html.
SciPy 1.17.0 served only as an independent validation oracle.

## 2026-11-09 — Matrix structure and SVD diagonal helpers

`src/linalg_structure.f90`, `src/linalg_diagsvd.f90`, their native tests,
independent Fortran/Python reference harness, example and documentation
are original SciFort MIT-licensed contributions. The structure helpers
use direct loops and the elementary definitions of transpose, conjugate
transpose, and diagonal bandwidth. The SVD diagonal helper performs
caller-allocated diagonal construction. No SciPy source or third-party
kernel was copied or vendored. Public SciPy documentation was consulted
only to confirm conventions and edge cases:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.bandwidth.html,
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.issymmetric.html,
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.ishermitian.html,
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.diagsvd.html.
SciPy 1.17.0 was used solely as the independent test oracle.

## 2026-11-09 — Block-diagonal and Khatri–Rao matrix operations

`src/linalg_block_diag.f90`, `src/linalg_khatri_rao.f90`, their native test,
independent Fortran/Python reference harness, runnable example, and API docs
are original MIT-licensed SciFort contributions. `block_diag` uses elementary
block copying; `khatri_rao` uses the defining columnwise Kronecker formula.
No SciPy source, third-party software, or numerical kernels were copied or
vendored. The public API documentation was used only to confirm input shapes,
empty-input and broadcasting conventions:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.block_diag.html
and https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.khatri_rao.html.
SciPy 1.17.0 was used exclusively as an independent validation oracle.

## 2026-11-09 — Hadamard, Leslie and DFT constructors

`src/linalg_hadamard.f90`, `src/linalg_leslie.f90`, `src/linalg_dft.f90`,
`test/test_linalg_hadamard_leslie_dft.f90`, the independent reference
Fortran driver and Python comparison script, executable example, and
API/validation docs are original MIT-licensed SciFort contributions.
Hadamard uses the standard bit-parity identity for Sylvester matrices;
Leslie uses its defining first-row/subdiagonal structure; DFT uses the
negative-exponent roots-of-unity formula with modular phase reduction.
No SciPy implementation source or third-party code was copied or vendored.
Public SciPy documentation was consulted only for its APIs and conventions:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.hadamard.html,
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.leslie.html,
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.dft.html.
SciPy 1.17.0 was used exclusively as a numerical test oracle.

## 2026-11-09 — Pascal and inverse Pascal matrices

`src/linalg_pascal.f90`, `test/test_linalg_pascal.f90`,
`reference/linalg_pascal_driver.f90`, `tools/compare_linalg_pascal.py`,
`example/linalg_pascal_example.f90`, and the documentation are original
MIT-licensed SciFort contributions. The formulas use Pascal's triangle,
`P = L L^T`, and `P^-1 = (L^-1)^T L^-1`, with the signed binomial formula for
`L^-1`. The decimal arithmetic uses elementary positional arithmetic in
base 10^9. No third-party implementation, big-integer package, or SciPy
source was copied or translated. Public SciPy documentation was consulted
for argument and result conventions:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.pascal.html
and https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.invpascal.html.
SciPy 1.17.0 was used solely as an independent numerical test oracle.

## 2026-11-09 — Helmert, Hilbert and inverse Hilbert constructors

`src/linalg_helmert.f90` and `src/linalg_hilbert.f90`, native and independent
reference tests, example and documentation are newly authored original
MIT-licensed SciFort contributions. Helmert uses the elementary orthonormal
contrast formula; inverse Hilbert uses the classical closed-form binomial
coefficient identity, evaluated with real arithmetic, bounded int64, or
independent elementary prime-factorized large-integer multiplication.
No SciPy implementation source or third-party big-integer code was copied,
linked or vendored. The public SciPy API documentation for `helmert`,
`hilbert` and `invhilbert` was consulted for conventions and dtype choices:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.helmert.html,
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.hilbert.html,
and https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.invhilbert.html.
SciPy 1.17.0 was used solely as an independent test oracle.

## 2026-11-09 — Hankel matrices and FFT multiplication

`src/linalg_hankel.f90`, `src/linalg_matmul_hankel.f90`, and their
new tests, reference driver, documentation, and example are original
MIT-licensed contributions. The constructor uses the defining Hankel
anti-diagonal indexing identity. The fast product uses the elementary
reversal identity `H = T*J` and the project's preexisting FFT-based
Toeplitz multiplication. No external code was copied or vendored.
SciPy's public API documentation was consulted for its `hankel(c,r)`
conventions: optional last row, ignored `r[0]`, and complex output:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.hankel.html.
`matmul_hankel` is a SciFort extension, not a SciPy function.

## 2026-11-09 — FFT-based rectangular Toeplitz multiplication

`src/linalg_matmul_toeplitz.f90`, native and SciPy reference tests,
example, and documentation were independently authored for SciFort under
its MIT license. They implement the standard circulant-embedding / Fourier
convolution identity using SciFort's preexisting `fft_plan` implementation.
SciPy's public `matmul_toeplitz` API documentation was consulted for its
first-column / optional-first-row convention, ignored `r(1)`, and rank-one
or rank-two multiplication semantics:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.matmul_toeplitz.html.
SciPy 1.17.0 and NumPy 2.3.5 were used exclusively as independent numerical
test oracles; no SciPy source code or external FFT implementation was
copied, translated, bundled, or linked. No new license notice or runtime
dependency was added.

## 2026-11-09 — Circulant matrices and FFT-based linear solves

`src/linalg_circulant.f90`, the accompanying test programs, SciPy comparison
harness, example and documentation are original SciFort MIT contributions.
The mathematical algorithm is the standard Fourier diagonalization of
circulant matrices and spectral Moore-Penrose pseudoinverse, implemented
independently. It reuses the existing `scifort_signal_spectral.fft_plan`
without copying or modifying third-party code. The public SciPy documentation
for `scipy.linalg.circulant` and `scipy.linalg.solve_circulant`
(https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.circulant.html
and https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_circulant.html)
was consulted for matrix orientation, singular handling and default tolerance.
SciPy 1.17.0 and NumPy 2.3.5 were used exclusively as independent test
oracles; no SciPy, LAPACK or other implementation source was imported.
No new third-party license or runtime dependency is introduced.

## 2026-11-09 — Toeplitz construction and bordering solves

`src/linalg_toeplitz.f90`, native tests, SciPy/NumPy reference harness,
example and documentation are independently authored MIT-licensed SciFort
contributions. The O(n^2) general Toeplitz solver follows elementary
block-matrix Schur complement and Toeplitz reversal identities, described
in `docs/linalg_toeplitz.md`. Its pivoted dense fallback reuses SciFort's
existing MIT `scifort_linalg_dense` routines. SciPy's public API contracts
were consulted for first-column/first-row conventions, optional conjugate
row, ignored `r(1)`, and RHS shapes:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.toeplitz.html
and https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_toeplitz.html.
SciPy 1.17.0 and NumPy 2.3.5 are used only as independent numerical test
oracles; no SciPy, LAPACK, Netlib, or other external implementation code
was copied or translated. No new license or runtime dependencies.

## 2026-11-09 — Symmetric/Hermitian LDL-based linear-system solvers

`src/linalg_indefinite_solve.f90` and its native tests, SciPy/NumPy reference
harness, example, and documentation are original, MIT-licensed SciFort
contributions. The implementation reuses the original SciFort
`scifort_linalg_ldl` Bunch–Kaufman factorization and implements the standard
mathematical three-step solve through a permuted unit triangular factor,
1x1/2x2 block-diagonal factor, and the factor's transpose or adjoint.
SciPy's public `solve` API was consulted for the symmetric/Hermitian solver
semantics: https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve.html.
SciPy 1.17.0 and NumPy are used only for independent numerical comparisons.
No LAPACK, SciPy, Netlib, or other third-party implementation source was
copied or translated; no additional runtime dependencies or license
obligations were introduced.

## 2026-11-08 — Generalized nonsymmetric eig/eigvals

`src/linalg_generalized_eigen.f90`, the native tests, independent SciPy
reference harness, example and documentation are original MIT-licensed
SciFort contributions. They reuse SciFort's own `scifort_linalg_qz` and
implement standard generalized Schur-pencil triangular eigenvector
back/forward substitution, homogeneous eigenvalue representation, and
normalization/residual checks. Mathematical basis: Golub and Van Loan,
*Matrix Computations*, 4th ed., chapters 7–8. Public SciPy API documentation
was consulted (https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eig.html).
SciPy and NumPy are used only for independent reference tests; no third-party
implementation source was copied. No new mandatory dependencies were added.

## 2026-11-08 — Two-sided complex QZ iteration

The iterative QZ routines added to `src/linalg_qz.f90`, new native tests,
example, updated SciPy comparison harness, and documentation are independently
authored, MIT-licensed SciFort contributions. The code uses classical
Hessenberg–triangular reduction and implicit two-sided single-shift QZ
iteration via complex Givens rotations, with generalized 2-by-2 shifts.
Mathematical/algorithmic references include Moler and Stewart, "An Algorithm
for Generalized Matrix Eigenvalue Problems," SIAM J. Numer. Anal. 10 (1973),
241–256, and the documented LAPACK `ZHGEQZ` mathematical contract:
https://math.utah.edu/software/lapack/lapack-z/zhgeqz.html.
The public SciPy API and independent SciPy 1.17.0 numerical results were used
for validation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.qz.html.
No LAPACK, SciPy, Netlib, Eigen, or other external implementation source was
copied, translated, or bundled. The older inverse-based SciFort method remains
available as a compatibility fallback; no additional library dependencies
or license notices are introduced.

## 2026-11-08 — Real-to-complex Schur factor conversion

The `rsf2csf` procedure in `src/linalg_schur.f90`, native/reference tests,
example and documentation implement the standard local 2-by-2 complex
unitary similarity obtained from an eigenvector of each real Schur block.
The implementation is a separately authored MIT-licensed SciFort addition
composed from standard linear-algebra operations; it uses the existing
SciFort Schur module and does not copy external source code. The public
SciPy API and independent SciPy 1.17.0 numerical behavior were consulted
for interface and reference validation:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.rsf2csf.html.
No new runtime dependency, vendored third-party code, or license notice
was introduced.

## 2026-11-08 — Real quasi-triangular Schur decomposition

The new real-output `schur` overload in `src/linalg_schur.f90`, native tests,
reference driver, Python comparisons, example, and documentation are original
MIT-licensed SciFort contributions. It reuses SciFort's existing complex
Schur routine and complex adjacent-unitary-swap routine, then constructs a
nested real orthonormal invariant-subspace basis using real/imaginary parts
of the reordered Schur vectors and two-pass Gram-Schmidt. A real orthogonal
similarity yields the quasi-upper-triangular form. This is an independently
written composition of classical algorithms; no SciPy, LAPACK, Netlib, or
third-party implementation source was copied or translated. SciPy 1.17.0
is used only for independent numerical checks and publicly documented API
semantics: https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.schur.html.
No external library dependency, imported license, or attribution requirement
was added.

## 2026-11-08 — Complex Schur and sorted Schur

`src/linalg_schur.f90`, its native tests, comparison driver, Python reference,
example, and documentation are independently authored MIT-licensed SciFort
contributions. The implementation reuses SciFort's existing MIT-licensed
`dense_complex_schur` and original neighboring unitary Schur swaps from
`src/linalg_funm.f90`. Mathematical basis: Golub and Van Loan,
*Matrix Computations*, 4th ed. (2013), chapters 7–8. The public
`scipy.linalg.schur` specification was consulted for the sorting conventions
(https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.schur.html).
SciPy 1.17.0 was used solely as an independent numerical reference.
No SciPy, LAPACK, Netlib, or other third-party implementation source was
copied or translated; no new library or license obligations arise.

## 2026-11-08 — Ordinary nonsymmetric eig/eigvals

`src/linalg_general_eigen.f90`, its native tests, SciPy comparison harness,
example and documentation are independently authored MIT-licensed SciFort
contributions. The code reuses SciFort's existing MIT-licensed
`scifort_linalg_matrix_functions.dense_complex_schur` and applies standard
complex Schur decomposition and triangular eigenvector substitution; see
Golub and Van Loan, *Matrix Computations*, 4th ed. (2013), chapters 7–8.
SciPy 1.17.0 provides only an independent numerical reference and public
API conventions (https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eig.html).
No LAPACK, SciPy, Netlib or other third-party implementation source was
copied, translated, or bundled. No new library or license obligations arise.

## 2026-11-08 — Banded full-spectrum inverse iteration

The new real-arithmetic `selected_banded_real_vectors` implementation and
full-spectrum `method='inverse'` routing in `src/linalg_banded_eigen.f90`
are independently authored MIT-licensed SciFort contributions. The work
reuses SciFort's existing banded Sturm bisection and partial-pivoted LU
routines with inverse iteration and two-pass modified Gram-Schmidt, classical
algorithms covered by Golub and Van Loan, *Matrix Computations*, 4th edition
(2013), chapter 8. The public SciPy `eig_banded` API was consulted for the
Hermitian band storage and eigenvector conventions; SciPy 1.17.0 served solely
as an independent numerical reference. No SciPy/Netlib/LAPACK implementation
code was copied, adapted or bundled; no third-party license obligations or
external dependencies were added.

## 2026-11-08 — Selected tridiagonal inverse iteration

The changes to `src/linalg_banded_eigen.f90`, `test/test_linalg_tridiagonal_vectors.f90`,
`tools/linalg_reference/driver_tridiagonal_vectors.f90`,
`tools/linalg_reference/run_tridiagonal_vectors_checks.py`, the example and
associated documentation are independently authored MIT-licensed SciFort
work. The implementation combines classical Sturm bisection and inverse
iteration of shifted symmetric tridiagonal matrices, with local row pivoting
and modified Gram-Schmidt reorthogonalization. The algorithmic basis is
Golub and Van Loan, *Matrix Computations*, 4th ed. (2013), chapter 8;
see also the public LAPACK `DSTEIN` documentation for its mathematical role:
https://www.netlib.org/lapack/explore-html/d6/d71/group__stein_ga387a6a45bb4b2072faf764c7afc33a72.html
No LAPACK/Netlib/SciPy implementation code was copied, translated, or bundled.
SciPy 1.17.0 was used solely as an independent numerical reference. No
third-party numerical library or new license obligation was introduced.

## 2026-11-08 — Eigenvalue selection by half-open interval

The `subset_by_value`/`num_selected` changes in
`src/linalg_dense_eigen.f90`, native tests, reference driver, example, and
new API documentation are independently authored MIT-licensed SciFort code.
They select from the sorted eigenvalues produced by SciFort's existing
Householder/QL solver and generalized-Cholesky reduction. The public
`scipy.linalg.eigh`/`eigvalsh` documentation supplies the `(low,high]`
interval convention; local SciPy 1.17.0 serves only as an independent
numerical reference. No SciPy/LAPACK implementation code is copied or
bundled, no third-party library is added, and no extra license obligations
arise.

## 2026-11-08 — Generalized Hermitian definite eigenproblems

The optional `b` and `problem_type` extensions in
`src/linalg_dense_eigen.f90`, their native tests, reference harness,
example and documentation are independently authored MIT-licensed SciFort
code. The algorithms are the classical Cholesky reductions of types 1, 2,
and 3 generalized Hermitian definite eigenvalue problems, as described in
Golub and Van Loan, *Matrix Computations*, 4th edition, 2013, chapter 8.
They reuse the preexisting SciFort `scifort_linalg_dense.cholesky` and
`solve_triangular` routines and the standard SciFort Householder/QL
eigensolver. Public SciPy `scipy.linalg.eigh`/`eigvalsh` documentation
(https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html)
provided API conventions and SciPy 1.17.0 served as an independent numerical
reference only. No SciPy, LAPACK or third-party source code was copied,
translated or bundled. No new dependencies or license notices are required.

## 2026-11-07 — Dense symmetric/Hermitian eigensolvers

`src/linalg_dense_eigen.f90`, native tests, Python/SciPy comparison scripts,
example and documentation are independently authored MIT-licensed SciFort
code. They compose existing SciFort MIT-licensed
`scifort_linalg_hessenberg.hessenberg` Householder reduction with
`scifort_linalg_banded_eigen.eigh_tridiagonal` implicit-shift QL iteration.
The complex phase transformation is standard unitary diagonal similarity
mathematics (see Golub and Van Loan, *Matrix Computations*, 4th ed., 2013,
chapters 7-8). The public `scipy.linalg.eigh` and `eigvalsh` APIs were
consulted for compatibility and numerical comparison only; no SciPy,
LAPACK or other third-party implementation code was copied, translated,
or bundled. No new numerical dependency or third-party license obligation
was introduced.

## 2026-11-07 — Tridiagonal and Hermitian band eigenvalue solvers

`src/linalg_banded_eigen.f90`, native tests, reference comparison, example and
API documentation are independently authored MIT-licensed SciFort work. The
scaled implicit-shift QL recurrence is adapted conceptually from the
preexisting SciFort MIT-licensed `src/signal_dpss.f90` tridiagonal QL
algorithm. Hermitian band handling reuses SciFort's MIT-licensed
`scifort_linalg_hessenberg` Householder reduction with diagonal-unitary
phase corrections. Both use the classical mathematical algorithms described
in Golub and Van Loan, *Matrix Computations*, 4th ed. (2013), chapters 7–8.
`scipy.linalg.eigh_tridiagonal`, `eigvalsh_tridiagonal`, `eig_banded`, and
`eigvals_banded` are used as independent numerical references and API
conventions only. No SciPy/LAPACK/third-party source was copied or bundled;
no new numerical dependencies or third-party license obligations were added.

## 2026-11-07 — Dense matrix balancing

`src/linalg_balance.f90`, its tests, example and reference harness are
original MIT-licensed SciFort work. The algorithm is an independent
implementation of the mathematical approach described by B. N. Parlett and
G. W. Reinsch, "Balancing a Matrix for Calculation of Eigenvalues and
Eigenvectors", *Numerische Mathematik* **13** (1969), 293--304,
DOI: 10.1007/BF02165404. SciPy's documented `matrix_balance` API defines
the intended similarity contract and supplied independent numerical test
values. No SciPy, LAPACK `xGEBAL`, or other third-party implementation
source was copied, adapted, translated, or bundled. No new dependencies or
third-party license obligations were introduced.

## 2026-11-07 — Real and complex Householder Hessenberg reduction

The implementation in `src/linalg_hessenberg.f90`, its native tests, example,
comparison harness and documentation are independently authored,
MIT-licensed SciFort work. The mathematical algorithm uses standard
Householder two-sided similarity transformations as described in Gene H.
Golub and Charles F. Van Loan, *Matrix Computations*, 4th edition (2013),
chapter 7 (nonsymmetric eigenvalue problem). The public
`scipy.linalg.hessenberg` contract supplies compatibility requirements and
an independently evaluated numerical reference; no SciPy, LAPACK or other
third-party implementation source has been copied, translated, linked or
bundled. No new third-party licenses or runtime dependencies were added.

## 2026-11-07 — Schur-based Sylvester and Lyapunov matrix equations

The changes in `src/linalg_matrix_equations.f90`, native tests, reference
harness, example and documentation are independently authored MIT-licensed
SciFort work. The classical Bartels–Stewart reduction and triangular
recurrences follow standard matrix-equation mathematics (see R. H. Bartels
and G. W. Stewart, "Solution of the Matrix Equation AX + XB = C",
*Communications of the ACM* 15 (1972), 820–826). Schur decomposition uses
SciFort's own MIT-licensed `scifort_linalg_matrix_functions.dense_complex_schur`;
no external library, SciPy, or LAPACK implementation is copied, translated,
linked or bundled. SciPy 1.17.0 is used solely for independent reference
comparisons against the public `solve_sylvester`, `solve_continuous_lyapunov`,
and `solve_discrete_lyapunov` APIs. The original direct Kronecker solver is
retained unchanged, except for dispatch and an explicit method option.

## 2026-11-07 — Orthogonal/unitary Procrustes (SciPy-compatible)

`src/linalg_procrustes.f90`, its tests, reference driver, example and
associated documentation are independently authored, MIT-licensed SciFort
work based on the standard singular-value decomposition of the Procrustes
cross-product. The mathematical solution is described by Peter H.
Schoenemann, "A generalized solution of the orthogonal Procrustes problem",
*Psychometrika* **31** (1966), 1–10, DOI: 10.1007/BF02289451.
The implementation uses SciFort's own `scifort_linalg_svd.svd` and scaled
arithmetic; it does not reproduce, translate or bundle any SciPy or LAPACK
source. The installed SciPy 1.17.0 `scipy.linalg.orthogonal_procrustes`
was used only as an independent reference for numerical validation and
public API semantics. No external numerical library, third-party source,
new license text, or dependency is introduced.

## 2026-11-07 — Values-only Jacobi SVD and spectral diagnostics

The changes in `src/linalg_svd.f90`, accompanying tests, reference harness,
benchmark, example, and documentation are independently authored MIT-licensed
SciFort work. The Jacobi rotations follow the same classical one-sided Jacobi
SVD mathematics as the existing `tall_svd` implementation; see Golub and
Van Loan, *Matrix Computations*, 4th edition (2013), chapter 8. The 2-norm
condition-number formula `sigma_max/sigma_min` and singular-value rank
threshold follow standard SVD theory. `scipy.linalg.svdvals` and
`numpy.linalg.cond`/`matrix_rank` were used for independent numerical
comparison only. No SciPy/NumPy/LAPACK source was copied, translated, linked,
or bundled; no new third-party licenses or dependencies were introduced.

## 2026-11-07 - Full-SVD basis completion, svdvals and diagsvd

The additions to `src/linalg_svd.f90`, full-size SVD example, native tests,
and reference harness are independently authored MIT-licensed SciFort work.
They use classical twice-applied modified Gram-Schmidt and greedy coordinate
pivot selection to complete the Jacobi singular-vector bases, as in standard
numerical linear algebra treatments (e.g., Golub and Van Loan, *Matrix
Computations*, 4th ed., chapters 5 and 8). The mathematical definition of
full singular-value decomposition and rectangular diagonal matrices follows
standard SVD theory. SciPy 1.17.0 public `svd`, `svdvals` and `diagsvd`
functions are used only as independent numerical reference outputs; their
implementations were not copied, translated, or bundled. No new third-party
license, source, or dependency is introduced. See `docs/linalg_svd.md` and
`docs/linalg_svd_full_validation_2026-11-07.md`.

## 2026-11-07 — SVD range/null-space bases and principal angles

`src/linalg_subspaces.f90`, its tests, reference driver, example, and associated
numerical documentation are original MIT-licensed SciFort code. The range and
null-space formulas follow the singular value decomposition; the principal
angle singular-value identities and complementary projection-residual method
are described by A. V. Knyazev and M. E. Argentati, "Principal Angles between
Subspaces in an A-Based Scalar Product: Algorithms and Perturbation Estimates",
SIAM J. Sci. Comput. 23 (2002),
pp. 2008–2040 (article DOI 10.1137/S1064827500377332). Twice-applied modified
Gram-Schmidt and canonical completion use standard orthogonalization formulas.
The SciPy `scipy.linalg.orth`, `null_space`, and `subspace_angles` APIs provide
independent numerical reference results; no SciPy, LAPACK, or other third-party
implementation source was translated, copied, linked, or bundled. The module
reuses SciFort's existing independently authored economy SVD. No new license
or runtime dependency is introduced. See `docs/linalg_subspaces.md` and its
validation record.


## 2026-11-07 — Economy Jacobi SVD and rank-deficient minimum-norm solves

`src/linalg_svd.f90` and its native tests, example, SciPy comparison harness
and documentation are MIT-licensed SciFort code. The Jacobi sweeps are based
on the standard one-sided Jacobi SVD mathematical algorithm as described in
G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed. (2013),
chapter 8. The rotation formulas and basis completion reuse concepts from
SciFort's own `src/linalg_polar.f90` Jacobi implementation, with changed
outputs, economy-wide handling, scale normalization and rank filtering.
Minimum-norm solutions and pseudoinversion use the standard truncated SVD
formulas. SciPy 1.17.0 `svd`, `pinv`, and `lstsq` are used only for independent
numerical reference; no SciPy, LAPACK or other upstream implementation source
was copied, translated, linked or bundled. No new runtime dependency or
third-party license was added. See `docs/linalg_svd.md` and its validation log.

## 2026-11-07 — Real/complex Householder QR and least squares

`src/linalg_qr.f90` and its native tests, standalone example, Python/SciPy
reference comparison and documentation are independently authored MIT-licensed
SciFort code. The algorithms use the standard Householder QR transformation,
explicit trailing-column norm pivoting, and full-rank triangular least-squares
and minimum-norm solving formulas, described in Golub and Van Loan,
*Matrix Computations*, 4th edition (2013), chapters 5 and 6. SciPy 1.17.0
is used solely for independent numerical comparison via its public `qr`
and `lstsq` interfaces. No SciPy, LAPACK or other implementation source was
copied, translated, linked or bundled. No new license or dependency is added.
See `docs/linalg_qr.md` and its validation record.

## 2026-11-07 — Reusable dense Cholesky factorization and solve

The new real/complex `cholesky`, `cho_factor`, and `cho_solve` routines in
`src/linalg_dense.f90`, along with tests, reference driver, example, and
associated documentation are independently authored MIT-licensed SciFort
code. They use the standard unblocked Cholesky factorization and triangular
substitution, as described in Golub and Van Loan, *Matrix Computations*,
4th edition (2013). `scipy.linalg.cho_factor`, `cholesky`, and `cho_solve` are
used as independent numerical references only, without copying, translating,
linking, or bundling SciPy/LAPACK source. The unused triangle is zeroed (unlike
SciPy's `cho_factor`, where the unused half is not defined), and the module
checks all entries for finiteness. The old dense LU implementation also had
two pairs of empty no-op loops removed; the mathematical algorithm is unchanged.
No new third-party code, notices, licenses, or runtime dependencies are added.
See `docs/linalg_dense_cholesky_validation_2026-11-07.md`.

## 2026-11-07 — Dense triangular and reusable LU solvers

`src/linalg_dense.f90`, its native tests, executable example, Python/SciPy
comparison harness and documentation are original MIT-licensed SciFort code.
The mathematical algorithms are standard forward/back substitution and dense
Gaussian elimination with partial pivoting as described in G. H. Golub and
C. F. Van Loan, *Matrix Computations*, 4th ed. (2013). SciPy 1.17.0 is used
only for independent numerical reference outputs and API terminology; no
SciPy/LAPACK source was copied, translated, linked, or bundled. Input arrays
are preserved, LU pivot indices are one-based, and the module adds no external
dependencies. See `docs/linalg_dense.md` and the corresponding validation
record.

## 2026-11-07 — Banded LU factor reuse and transpose systems

The `lu_factor_banded` and `lu_solve_banded` additions to
`src/linalg_banded.f90` and their new native and reference tests, driver,
example, and docs are independently authored MIT-licensed SciFort code.
They implement standard Gaussian elimination with partial pivoting in
expanded-band storage, and the corresponding forward/backward and
transpose/conjugate-transpose substitutions. LAPACK's GBTRF/GBTRS concepts
are used for public-facing storage terminology only. The SciPy 1.17.0
public reference solves are used strictly as independent test outputs;
no SciPy/LAPACK source was copied, translated, included, or linked.
No new runtime dependencies or third-party licenses were introduced.
See `docs/linalg_banded.md` and
`docs/linalg_banded_lu_validation_2026-11-07.md`.

## 2026-11-07 — Dense banded Gaussian elimination and Cholesky

The `src/linalg_banded.f90` module, native tests, example, comparison driver,
Python/SciPy validation harness and documentation are original MIT-licensed
SciFort work. The packed diagonal indexing follows the public SciPy API;
the pivoted Gaussian elimination and Cholesky formulas are standard numerical
linear algebra. SciPy 1.17.0 is used solely to generate independent reference
outputs. No SciPy, LAPACK or other third-party source code has been copied,
translated, vendored or linked. See `docs/linalg_banded.md` and
`docs/linalg_banded_validation_2026-11-07.md`.

## 2026-11-07 — Extended adaptive Taylor evaluation of `funm` diagonal blocks

The optional higher-order derivative evaluation, convergence tests, scaled
nilpotence screening, native tests, Python reference generator and
documentation in this release are original MIT-licensed SciFort work.
The Taylor-series mathematical identity and Schur--Parlett framework are
standard; see N. J. Higham, *Functions of Matrices: Theory and Computation*,
SIAM (2008). The SciPy 1.17.0 public matrix function interfaces are used only
to compute independent numerical reference values. No SciPy source or other
third-party implementation was copied or imported, and no dependency or
license has changed. See
`docs/linalg_funm_adaptive_taylor_validation_2026-11-07.md`.

## 2026-11-07 — Block Schur--Parlett with derivative-free semisimple clusters

`src/linalg_funm.f90` now uses an original MIT-licensed triangular Sylvester
solver for inter-block matrix-function equations. Its formulas follow the
standard block Schur--Parlett relation derived from `T*F = F*T` and triangular
back substitution. Global unitary Schur reordering is applied even when no
scalar derivative is available; diagonalizable repeated-eigenvalue blocks
need only function values. The mathematical framework is described in
N. J. Higham, *Functions of Matrices: Theory and Computation*, SIAM (2008).
The new native tests and Python/SciPy validation generator are original;
SciPy is used for numerical reference outputs only. No external source has
been copied, translated, or vendored, and no runtime dependency or license
change is required. See `docs/linalg_funm_block_parlett_validation_2026-11-07.md`.

## 2026-11-07 — Schur reordering for repeated dense `funm` eigenvalues

The improvements to `src/linalg_funm.f90`, native test, standalone reference
comparison script, and supporting documentation are original MIT-licensed
SciFort code. They use the standard mathematical construction of a unitary
2x2 Schur swap from an eigenvector of an upper triangular 2x2 matrix, and
stable grouping of close eigenvalues as described in the Schur--Parlett
framework (N. J. Higham, *Functions of Matrices: Theory and Computation*,
SIAM, 2008). The existing SciFort complex Schur algorithm is reused.
Independent numerical comparisons use installed SciPy 1.17 specialized matrix
functions and do not copy, translate, or incorporate SciPy source. No new
license or runtime dependency is introduced. Details and limits are in
`docs/linalg_funm.md` and its reordered-validation record.

## 2026-11-07 — General analytic dense matrix functions

`src/linalg_funm.f90`, `test/test_linalg_funm.f90`, the example,
Python comparison harnesses, documentation and validation records are
original MIT-licensed SciFort code. The algorithm follows the published
Schur--Parlett commutator recurrence and Taylor expansion for repeated
Schur eigenvalues, as described in N. J. Higham, *Functions of Matrices:
Theory and Computation* (SIAM, 2008), and Golub and Van Loan,
*Matrix Computations*, 4th edition (2013). The existing MIT-licensed
SciFort complex Schur helper is reused without copying any new third-party
implementation. SciPy 1.17.0 public `funm` and specialized matrix functions
were consulted as numerical *references only*. No SciPy source, LAPACK or
other third-party source code has been translated or imported. No new runtime
or licensing dependency is introduced. See `docs/linalg_funm.md`.

## 2026-11-07 — Dense exponential condition numbers

`src/linalg_expm_condition.f90`, its native test, Fortran example,
comparison drivers, documentation, and validation checks are original
MIT-licensed SciFort code. The mathematical method is based on the
standard Frobenius norm of a Fréchet derivative, the fact that the adjoint
of `L_exp(A,.)` is `L_exp(A^H,.)`, and fully reorthogonalized Lanczos
iteration for the largest eigenvalue of `L* L`. Diagonal matrix divided
differences and scalar-shift invariance are standard exponential-function
identities. Independent outputs are compared against SciPy 1.17.0; no
SciPy, LAPACK, or other external source was translated, copied, or
bundled. Scientific background: A. H. Al-Mohy and N. J. Higham,
"Computing the Fréchet Derivative of the Matrix Exponential, with an
Application to Condition Number Estimation", *SIAM Journal on Matrix
Analysis and Applications* 30(4), 1639-1657 (2009),
DOI: 10.1137/080716426. No new runtime dependency or license is added.
See `docs/linalg_expm_condition.md`.

## 2026-11-07 — Dense trigonometric and hyperbolic matrix functions

`src/linalg_trigonometric.f90`, its standalone native tests, Fortran example,
Python reference driver/comparison harness and documentation are original
MIT-licensed SciFort work. The implementations use the standard Taylor
series of sin/cos/sinh/cosh, argument scaling and double-angle identities,
and SciFort's own MIT-licensed pivoted complex inverse from
`src/linalg_matrix_functions.f90`. Mathematical background: N. J. Higham,
*Functions of Matrices: Theory and Computation*, SIAM (2008), matrix
analytic-function identities. Independent test values come from SciPy 1.17.0;
no SciPy source is translated or incorporated. No new third-party code,
license or runtime dependency is introduced. See
`docs/linalg_trigonometric.md`.

## 2026-11-07 — Dense matrix logarithm and fractional powers

`src/linalg_logarithm.f90`, `test/test_linalg_logarithm.f90`,
`example/linalg_logarithm_example.f90` and the Python reference harness are
original MIT-licensed SciFort source. The algorithm composes the existing
MIT-licensed SciFort `sqrtm` and `expm_frechet` implementations with standard
inverse scaling-and-squaring and a matrix-atanh Taylor expansion. Its pivoted
dense solve and integral exponentiation are original implementations of
standard linear-algebra techniques. Public API semantics and independent
numerical outputs are compared with SciPy 1.17.0; no SciPy code is copied,
translated or bundled. Scientific background: A. H. Al-Mohy and N. J.
Higham, "Improved Inverse Scaling and Squaring Algorithms for the Matrix
Logarithm", SIAM Journal on Scientific Computing 34(4), C152-C169 (2012),
and N. J. Higham and L. Lin, "A Schur-Pade Algorithm for Fractional Powers
of a Matrix", SIAM Journal on Matrix Analysis and Applications 32(3),
1056-1078 (2011). The algorithm here is a simpler independent implementation,
not an adaptation of the Schur-Pade source. No additional license notices or
dependencies are required; limitations are recorded at
`docs/linalg_logarithm.md`.

## 2026-11-06 — Dense polar decomposition and exponential Fréchet derivative

`src/linalg_polar.f90`, `src/linalg_expm_frechet.f90`, the Fortran native test,
example, and Python comparison harness are original MIT-licensed SciFort
code. The one-sided complex Jacobi polar method and differentiated Taylor
scaling-and-squaring method are standard published numerical ideas; no SciPy,
LAPACK, or external repository source has been copied or bundled. SciPy 1.17
and NumPy were used only as independent numerical reference implementations.
There is no new runtime library or licensing dependency. See
`docs/linalg_polar_frechet.md` and the associated validation record.

## 2026-11-06 — Dense `sqrtm` and `signm` matrix functions

`src/linalg_matrix_functions.f90`, the native test, example and comparison
harness are SciFort MIT-licensed code. The QR/Schur helper adapts existing
SciFort MIT-licensed `src/sparse_linalg_complex_eigs.f90`; it does **not**
copy code from SciPy, NTPoly, PSFUN, or LAPACK. The algorithms are standard
Schur triangular square-root recurrence and scaled Newton matrix sign.
SciPy/NumPy are numerical testing references only. Candidate permissively
licensed Fortran libraries were reviewed, but **none were incorporated**;
see `docs/linalg_external_library_review_2026-11-06.md`.

## 2026-11-06 — Dense algebraic Riccati equation solvers

`src/linalg_riccati.f90` and its associated native test, example, Fortran
reference driver, and Python comparison harness are original MIT-licensed
SciFort work. The continuous matrix-sign/Newton and discrete structured
Doubling methods are established numerical algorithms; this code is newly
written and does not incorporate SciPy or third-party implementation source.
The public SciPy ARE equations and independent reference solutions were used
only for mathematical/API and numerical checks. The implementation reuses
SciFort's existing MIT-licensed dense Sylvester/Lyapunov kernels and introduces
no new external numerical dependencies. See `docs/linalg_riccati.md`.

## 2026-11-06 — Dense Sylvester and Lyapunov equations

`src/linalg_matrix_equations.f90`, the native tests, example, and reference
harness are original MIT-licensed SciFort work. The implementation uses the
standard vec/Kronecker identities for the Sylvester and continuous/discrete
Lyapunov equations and ordinary complex Gaussian elimination with partial
pivoting. No SciPy/LAPACK source was translated, copied, or bundled. SciPy's
public mathematical API and independently computed solutions are used as
references only. In contrast to SciPy's Schur/Bartels–Stewart methods, this
implementation constructs a dense Kronecker system with an explicit
400-unknown limit; see `docs/linalg_matrix_equations.md`. The new implementation
has no additional runtime or license dependencies.

## 2026-11-06 — Complex sparse exponential actions

`src/sparse_linalg_complex_expm.f90`, its tests, reference driver and Python
comparison harness are original MIT-licensed SciFort work. They adapt the
algorithmic design of SciFort's original MIT-licensed real `expm_multiply`
for complex CSR/CSC and matrix-free operations, including the standard
scaling-and-Taylor series of the exponential and a geometric remainder
bound. SciPy and NumPy are independent numerical references only. No SciPy
source code is embedded and no new third-party dependency is introduced.

## 2026-11-06 — Matrix-free complex Arnoldi and shift-cache extension

The matrix-free callback interface, reusable factor-cache type, and refactoring
of `src/sparse_linalg_complex_eigs.f90` are original MIT-licensed SciFort
work extending the existing original MIT-licensed complex Arnoldi and sparse
pivoted LU modules. The standard transformed eigenvalue mapping and true
residual checks are well-known algorithms. The new tests, streamed Fortran
driver, and randomized NumPy reference harness are original. Neither SciPy nor
NumPy source is embedded. No new third-party license or runtime dependency is
introduced. See `docs/sparse_complex_eigs_operator.md`.

## Complex shift-invert eigensolver extension (2026-11-06)

The optional `sigma` branch in `src/sparse_linalg_complex_eigs.f90` is
original MIT-licensed SciFort code. It composes SciFort's existing complex
pivoted LU solver with its original complex Arnoldi/Rayleigh--Ritz solver.
The transformed spectral mapping and residual check follow standard
shift-invert linear algebra. The accompanying Fortran regression test and
Python reference comparison script are original; SciPy and NumPy are
numerical test references only. No external code was copied or runtime
library dependency added. See `docs/sparse_complex_eigs_shift.md`.

## Complex nonsymmetric sparse eigenvalues (2026-11-06)

`src/sparse_linalg_complex_eigs.f90` uses original SciFort MIT-licensed
Arnoldi/Rayleigh--Ritz logic. Its dense projected complex QR/Schur helper was
adapted from SciFort's existing original MIT-licensed real-matrix
`eigs_arnoldi` implementation. The tests and reference script are original
SciFort code; SciPy and NumPy are independent numerical reference packages,
not sources of incorporated implementation code. No third-party source or
new runtime dependency was added. See `docs/sparse_complex_eigs.md`.

## 2026-11-06 — Complex sparse singular-value decomposition

`src/sparse_linalg_complex_svds.f90`, the Fortran test and the reference
comparison harness are original MIT-licensed SciFort source, adapting
SciFort's existing MIT-licensed real `svds` design. The mathematics follows
standard Golub--Kahan bidiagonalization and complex one-sided Jacobi SVD.
No SciPy, ARPACK, PROPACK, LAPACK or other external numerical source was
copied; SciPy and NumPy are independent numerical test references. No new
link or package dependency is introduced.

## 2026-11-06 — Complex sparse LSMR

`src/sparse_linalg_complex_lsmr.f90` adapts the established real LSMR
recurrence in `src/optimize/lsmr.f90` for complex arithmetic with the
Hermitian adjoint, SciFort's CSR/CSC types, and true residual diagnostics.
The original real LSMR carries BSD-3-Clause SciPy/Fong--Saunders provenance,
so this complex adaptation retains **BSD-3-Clause**, rather than describing
it as an original MIT implementation. The retained notice is at
`LICENSES/scipy-lsmr-BSD-3-Clause.txt`. The test harness and wrapper were
written for SciFort; SciPy is used as a numerical reference.

## 2026-11-06 — Complex sparse LSQR

`src/sparse_linalg_complex_lsqr.f90` and its native and SciPy/NumPy comparison
harnesses are original MIT-licensed SciFort code. The mathematics follows the
published Golub–Kahan bidiagonalization and Paige–Saunders LSQR recurrences
(C. C. Paige and M. A. Saunders, "LSQR: An Algorithm for Sparse Linear
Equations and Sparse Least Squares," ACM Transactions on Mathematical
Software 8(1), 43-71, 1982).
SciFort's existing real `lsqr` is used as an implementation reference; no
external LSQR, SciPy or NumPy source code was copied or incorporated. No new
runtime, build-system or licensing dependencies were introduced.

## 2026-11-06 -- Complex sparse pivoted LU and direct solves

`src/sparse_linalg_complex_direct.f90`, its native regression test, reference
Fortran driver and Python SciPy comparison harness are original MIT-licensed
SciFort work. The numerical method is classical Gaussian elimination with
partial row pivoting and sparse fill-in. The existing SciFort real sparse LU
provides style and API conventions; no SciPy or SuperLU source was copied,
translated or bundled. The documented SciPy public `splu`/`spsolve` API and
SciPy/NumPy computations are used solely as behavior and numerical references.
No external numerical dependency or third-party source notice is introduced.

## 2026-11-05 -- Complex ILU(0), BiCG and QMR

`src/sparse_complex_ilu.f90` is an original MIT-licensed implementation of
standard unpivoted zero-fill LU factorization and forward/adjoint triangular
substitution on the pre-existing canonical one-based CSR pattern. The complex
BiCG and ILU-preconditioned BiCGSTAB/GMRES extensions are original MIT code
based on textbook Krylov recurrences, not imported from SciPy/SuperLU.

The complex QMR solver adapts the *existing SciFort real QMR recurrence* in
`src/sparse_linalg.f90`, which itself adapts the BSD-3-Clause SciPy QMR
algorithm. The changes include complex-conjugate shadow directions,
complex-scaled quasi-minimal updates, and a left-preconditioned operator
whose Hermitian adjoint applies the conjugate-transposed ILU inverse. As a
conservative provenance measure, `src/sparse_linalg_complex.f90` now carries
`SPDX-License-Identifier: MIT AND BSD-3-Clause`; the upstream notice is
already retained at `LICENSES/scipy-qmr-BSD-3-Clause.txt` and listed in
`THIRD_PARTY_LICENSES.md`. The numerical algorithm was originally developed
by R. W. Freund and N. M. Nachtigal, "QMR: a quasi-minimal residual method
for non-Hermitian linear systems", *Numerische Mathematik* 60, 315-339
(1991), DOI: 10.1007/BF01385726.

SciPy 1.17 and NumPy are used solely as numerical testing oracles; the new
complex ILU factorization is independently checked against a small Python
zero-fill factor and both direct and adjoint triangular solves. The native
regression suite, Python comparison script, and input driver are original
SciFort test code. No external numerical dependency was added.

## 2026-11-05 -- Complex compressed sparse storage and Krylov methods

`src/sparse_complex.f90`, `src/sparse_linalg_complex.f90`, the native tests,
and the independent SciPy/NumPy comparison harness are original MIT-licensed
SciFort implementations. CSR/CSC formats, modified Gram-Schmidt Arnoldi,
complex Givens rotations, CG and BiCGSTAB are standard publicly documented
numerical techniques. Existing local MIT-licensed COO duplicate aggregation
is reused through `scifort_sparse_coo_nd_complex`. No SciPy source code was
copied or incorporated; SciPy and NumPy are external test references only.
No additional Fortran dependencies or third-party license notices.

## 2026-11-05 — Complex N-dimensional COO arrays and tensor operations

`src/sparse_coo_nd_complex.f90`, `src/sparse_coo_nd_complex_ops.f90`,
new tests and comparison scripts are MIT-licensed SciFort work derived from
SciFort's existing MIT-licensed real N-D COO modules. The new complex
arithmetic and Hermitian operations use standard Fortran intrinsic operations
and the same sparse sort/merge methods as the local real implementation.
The SciPy and NumPy public numerical APIs provide independent test references;
no third-party source code was incorporated. There are no new runtime or
Fortran-package dependencies.

## 2026-11-05 -- N-D sparse Boolean broadcasting and reductions

`src/sparse_coo_nd_bool_reduce.f90`, its native regression, streaming
reference driver, and NumPy comparison script are original MIT-licensed
SciFort work. Sparse coordinate projection/aggregation, broadcast shape
matching and Boolean truth tables are general techniques. NumPy dense
arrays are used only as independent reference results; no third-party
source is copied or incorporated. No new dependencies or notices.

## 2026-11-05 -- Complement-aware sparse Boolean masks

`src/sparse_coo_nd_bool.f90`, its native tests, streaming comparison driver,
and independent NumPy/SciPy comparison harness are original MIT-licensed
SciFort code. The implementation uses an implicit Boolean background and a
sorted coordinate exception set, standard sparse set-merge algorithms, and
existing SciFort COO duplicate aggregation; no third-party source is copied.
SciPy COO and NumPy Boolean arrays serve as independent numerical/structural
reference oracles. No new third-party dependencies or license obligations.

## 2026-11-05 -- Sparse COO N-D predicates and coordinate masks

`src/sparse_coo_nd_filter.f90`, its native Fortran regression tests, streaming
comparison driver and Python SciPy/NumPy reference harness are original
MIT-licensed SciFort work. The algorithms use conventional stable COO duplicate
aggregation, sorted sparse coordinate lookup and binary search. The public
SciPy 1.17.0 API was used for numerical/structural reference comparisons;
no SciPy source was copied, adapted or imported. No new dependencies or
third-party license obligations were introduced.

## 2026-11-05 -- N-dimensional COO boolean selection and indexed updates

`src/sparse_coo_nd_edit.f90`, its Fortran regression, reference driver and
Python comparison script are original MIT-licensed SciFort work. The code
uses standard stable merge sorting, COO aggregation, and sparse coordinate
selection techniques. SciPy 1.17.0 and NumPy serve as independent numerical
reference oracles; no SciPy source code was incorporated. The mutable COO
scatter API is Fortran-specific and introduces no additional dependency.

## 2026-11-05 -- N-dimensional COO indexing and slicing

`src/sparse_coo_nd_index.f90`, its native regression tests, streaming driver,
and Python comparison script are original MIT-licensed SciFort code. Public
SciPy 1.17 indexing behavior and NumPy's mathematical selection semantics were
used as independent numerical references; no SciPy or NumPy source code was
copied, ported or incorporated. The one-based inclusive-slice and orthogonal
multi-axis selection APIs are documented Fortran-specific choices. No new
third-party dependency or license obligation is introduced.

## 2026-11-05 — N-dimensional sparse reductions and tensor operations

`src/sparse_coo_nd_ops.f90`, `test/test_sparse_coo_nd_ops.f90`, and the
reference driver/comparison script are original MIT-licensed SciFort work.
They use standard COO coordinate aggregation, NumPy-style dimension alignment,
stable merge sorting and a sparse sort-merge contraction join. SciPy 1.17.0
and NumPy are independent reference oracles, not incorporated source code.
SciPy's published `coo_array` API was consulted for input/duplicate semantics;
its source implementation was not copied or translated. No additional
third-party dependencies or notices are required.

## 2026-11-05 — N-dimensional COO arrays

`src/sparse_coo_nd.f90` and its native/SciPy comparison tests are original,
MIT-licensed SciFort implementations. The storage model and expected public
behavior were informed by the documented `scipy.sparse.coo_array` API, but
no SciPy numerical source code was copied, translated or incorporated. Stable
bottom-up mergesort provides C-order canonicalization, and 64-bit checked
strides support large sparse shapes without dense allocation. No new
third-party dependencies or license notices are required.

## 2026-11-05 — Incremental DOK and LIL sparse formats

`src/sparse_dok_lil.f90` is original MIT-licensed SciFort Fortran code for
open-addressed hash-table DOK and dynamically grown, sorted row-list LIL
sparse matrices. It interoperates with the existing SciFort canonical
CSR/CSC and COO modules and uses textbook data-structure techniques.
SciPy 1.17.0 public documentation and its installed numerical behavior were
consulted for API semantics, especially explicit-zero assignment and
conversion behavior. SciPy is used as an independent regression oracle;
no SciPy source implementation was copied or translated.
No third-party code or additional dependency was introduced.

## 2026-11-05 — COO sparse storage and sparse.find

`src/sparse_coo.f90` is an original SciFort MIT-licensed implementation of
2D coordinate sparse matrix storage, duplicate coalescing, CSR/CSC
interoperation and stored-nonzero extraction. It uses only existing SciFort
canonical CSR/CSC kernels; implementation was not copied or translated from
SciPy source. SciPy 1.17.0 documentation is consulted for API semantics, and
SciPy is used as an independent oracle for numerical/structural regression
tests. No additional third-party code, attribution obligation, or dependency.

## 2026-11-05 — Generalized matrix-function Krylov actions

`scifort_sparse_linalg.funm_multiply_krylov` in `src/sparse_linalg.f90`
is an original SciFort MIT-licensed implementation of restarted
Arnoldi/Lanczos projection, block-Hessenberg restart coupling and
matrix-function correction updates described mathematically in:

- M. Eiermann and O. G. Ernst, "A Restarted Krylov Subspace Method for the
  Evaluation of Matrix Functions," *SIAM J. Numerical Analysis* 44 (2006),
  2481--2504, DOI: 10.1137/050633846.
- M. Afanasjew, M. Eiermann, O. G. Ernst and S. Guettel, "Implementation of
  a Restarted Krylov Subspace Method for the Evaluation of Matrix Functions,"
  *Linear Algebra and Its Applications* 429 (2008), 2293--2314,
  DOI: 10.1016/j.laa.2008.06.029.

SciPy 1.17.0's public `funm_multiply_krylov` interface is referenced for
API behavior. The new Fortran numerical code and tests are independently
written, not a line-by-line copy or translation of the SciPy source.
SciPy is also used as a second independent test oracle, alongside dense
SciPy matrix functions. No new third-party numerical source or license
is imported, and existing third-party notices remain unchanged.

## 2026-11-05 — N-dimensional grid Laplacians and sparse structural queries

`src/sparse_laplacian_nd.f90` implements a finite-difference grid Laplacian
using the standard Kronecker-sum stencil and separated analytic eigenpairs
for Dirichlet, Neumann, and periodic boundary conditions. The formulas are
classical discrete sine, cosine and periodic Fourier bases. For reference:
G. Strang, *The Discrete Cosine Transform*, SIAM Review 41 (1999), 135–147.
`src/sparse_linalg.f90` adds independent original CSR/CSC bandwidth and
structural triangularity scans. These are original SciFort MIT-licensed
implementations, not translations or adaptations of SciPy source code.
SciPy 1.17 is used only for behavioral comparisons in the regression tools.
No additional third-party license notices or FPM dependencies are required.

## 2026-11-05 — Augmented/recycled sparse GMRES (LGMRES, GCROT(m,k))

`lgmres`, `gcrotmk` and their shared `recycled_arnoldi_cycle` helper in
`src/sparse_linalg.f90` are original SciFort MIT-licensed implementations
of the mathematical augmented/recycled flexible-GMRES framework:

- A. H. Baker, E. R. Jessup, and T. Manteuffel, "A Technique for
  Accelerating the Convergence of Restarted GMRES", *SIAM J. Matrix
  Anal. Appl.* 26, 962 (2005), for LGMRES.
- E. de Sturler, "Truncation Strategies for Optimal Krylov Subspace
  Methods", *SIAM J. Numer. Anal.* 36, 864 (1999), and
  J. E. Hicken and D. W. Zingg, "A Simplified and Flexible Variant of
  GCROT for Solving Nonsymmetric Linear Systems", *SIAM J. Sci.
  Comput.* 32, 172 (2010), for GCROT(m,k).

The code follows these standard recurrences and uses SciFort's existing
CSR matvec, Jacobi/ILU(0) application, and status conventions. Neither
SciPy's LGMRES nor its GCROT Python numerical source was copied, translated,
or adapted. SciPy 1.17.0 is used only as an independent behavioral
reference in `tools/sparse_reference/compare_sparse_recycling_scipy.py`.
Tests and documentation are original SciFort MIT work. No additional
third-party source-code notices or FPM dependencies are required.

## 2026-11-05 — Sparse QMR and TFQMR iteration

`scifort_sparse_linalg.qmr` in `src/sparse_linalg.f90` adapts the
real two-sided Lanczos and quasi-minimal residual update recurrence from
SciPy 1.17.0 `scipy/sparse/linalg/_isolve/iterative.py`, function `qmr`,
BSD-3-Clause:
https://github.com/scipy/scipy/blob/v1.17.0/scipy/sparse/linalg/_isolve/iterative.py.
The full applicable license and SciPy copyright notice are retained in
`LICENSES/scipy-qmr-BSD-3-Clause.txt`; `THIRD_PARTY_LICENSES.md` lists the
adaptation. The SciFort code replaces SciPy's operator callbacks with a
one-based canonical CSR matrix and a materialized CSR transpose, uses a
right-sided Jacobi inverse, and independently verifies `||b-A*x||` after
every step. QMR is a method of R. W. Freund and N. M. Nachtigal (1991), also
described in Barrett et al., *Templates for the Solution of Linear Systems*
(SIAM, 1994): https://www.netlib.org/templates/templates.pdf.

`scifort_sparse_linalg.tfqmr` is a derivative of SciFort's existing
**BSD-3-Clause** `src/optimize/krylov_iterative.f90`, `solve_tfqmr`,
whose file header credits the SciPy 1.17.0 BSD-3-Clause
`scipy/sparse/linalg/_isolve/tfqmr.py` recurrence. It uses the opt-in
mathematically consistent right-preconditioned `A*M` variant, now specialized
to real CSR and existing Jacobi/ILU(0) inverse operations. The CSR adapter
performs actual residual checks every half-step. The full upstream SciPy
BSD-3-Clause notice for this adapted recurrence is retained in
`LICENSES/scipy-tfqmr-BSD-3-Clause.txt`. SciPy remains the independent
numerical reference for testing, but this code also has a SciPy-derived
algorithmic provenance, not MIT-only provenance.
The shared mathematics is documented in R. W. Freund, *A Transpose-Free
Quasi-Minimal Residual Algorithm for Non-Hermitian Linear Systems*, SIAM
J. Sci. Comput. 14(2), 470-482 (1993).

New `test/test_sparse_qmr_tfqmr.f90` and the independent comparison driver
in `tools/sparse_reference/` are original SciFort MIT-licensed tests. This
change adds no new FPM dependency. Both new sparse solvers retain the
applicable SciPy BSD-3-Clause notices.

## 2026-11-05 — BiCG/CGS iterations and Kronecker sum

`bicg` and `cgs` in `src/sparse_linalg.f90` are original SciFort
MIT-licensed implementations of well-established iterative methods, following
standard mathematical recurrences and preconditioning specifications described
in R. Barrett et al., *Templates for the Solution of Linear Systems: Building
Blocks for Iterative Methods*, SIAM (1994):
https://www.netlib.org/templates/templates.pdf (BiCG and CGS chapters).
These routines reuse SciFort's existing CSR, true residual, Jacobi and ILU(0)
operations. Neither SciPy's implementation nor the Templates source code was
copied, translated, or adapted. SciPy is used solely as a behavioral oracle.

`csr_kronsum` in `src/sparse.f90` is an independent CSR row-merge
implementation of the mathematical identity
`A (plus) B = I_n (tensor) A + B (tensor) I_m`, specified in
https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.kronsum.html.
It is not translated or copied from SciPy. All added tests and driver source
are SciFort MIT-licensed; no new runtime dependencies or imported notices.

## 2026-11-05 — Sparse constructors (`csr_random`, `csr_rand`, `csr_diags_array`)

`src/sparse_constructors.f90` is original SciFort MIT-licensed implementation
using standard Floyd sampling, open-addressed hashing, Box–Muller Gaussian
variates, SciFort's pre-existing MRG32k3a generator, and the existing CSR
construction routines. The code was not imported, copied, or translated from
SciPy or NumPy. SciPy 1.17 provides independent numerical and API-behavior
reference results only. All existing third-party notices are unchanged.

## 2026-11-05 — CSR deletion, reshape and sparse-to-sparse assignment

`csr_delete_entries`, `csr_delete_rows`, `csr_delete_columns`, `csr_reshape`,
and the revised sparse-source `csr_assign` in `src/sparse.f90` are original,
MIT-licensed SciFort implementations using standard sorted-index merging
and CSR row-pointer construction. No source code was copied from SciPy,
NumPy, or another repository. SciPy/NumPy are used only as independent
numerical references; all earlier third-party notices are preserved.

## 2026-11-05 — Two-dimensional sparse blocks and scatter assignment

The implementations of `csr_bmat`, `csr_block`, `csr_set`,
`csr_assign_coo`, and the generic `csr_assign` in `src/sparse.f90` are
original MIT-licensed SciFort code. They use direct block-row assembly and
standard sorted-coordinate merge algorithms. No code was copied, translated,
or imported from SciPy, NumPy, or another third-party source. SciPy 1.17 and
NumPy are used only as independent numerical references. All pre-existing
third-party attributions and license files are unchanged.

## 2026-11-05 — CSR assembly and indexing operations

The `csr_hstack`, `csr_vstack`, `csr_block_diag`, `csr_kron`, `csr_get`,
`csr_select`, and `csr_slice` additions in `src/sparse.f90` are original
MIT-licensed SciFort code. Implementations operate directly on one-based
canonical sparse rows and were not copied, adapted, or translated from SciPy
or any other external library. SciPy 1.17 is used solely for independent
numerical comparison tests. Existing upstream copyright notices and licenses
are preserved unchanged. See `docs/sparse.md` for API semantics.

## 2026-11-05 — DIA/BSR and sparse norms

The new `dia_matrix`, `bsr_matrix`, and related conversion and matrix-vector
routines in `src/sparse.f90`, together with `norm` and `norm_axis` in
`src/sparse_linalg.f90`, are original MIT-licensed SciFort implementations
of standard sparse matrix storage and numerical norm definitions. The code
was not copied or translated from SciPy, NumPy, or an external repository.
Existing SciFort `svds` is reused for the approximate spectral 2-norm.
SciPy 1.17 serves only as an independent test oracle, and pre-existing
third-party copyright notices and BSD license terms remain unchanged.

## 2026-11-05 — Sparse matrix powers and explicit inverse

`matrix_power` and `inv` in `src/sparse_linalg.f90` are independently
implemented MIT-licensed SciFort code using standard exponentiation by
squaring and the project's existing pivoted sparse LU/triangular solver.
These are mathematical operations also exposed by SciPy's
`scipy.sparse.linalg.matrix_power` and `scipy.sparse.linalg.inv`.
No SciPy, NumPy, or third-party source code has been copied, translated,
or imported for these routines. Randomized NumPy/SciPy computations are
used as test references only. Attribution and licenses of existing
BSD-origin SciFort components are preserved.

## 2026-11-05 — Locally optimal block eigensolver (`lobpcg`)

`scifort_sparse_linalg.lobpcg`, its validation driver and regression tests
are independently written MIT-licensed SciFort code implementing the
published LOBPCG algorithm (A. V. Knyazev, “Toward the Optimal Preconditioned
Eigensolver: Locally Optimal Block Preconditioned Conjugate Gradient Method,”
SIAM J. Sci. Comput. 23(2), 517–541, 2001). The implementation uses the
existing independently implemented projected symmetric Jacobi eigensolver.
SciPy/NumPy are used for numerical comparisons only; no SciPy or another
repository's LOBPCG source is copied or imported. Older BSD-licensed
components and their notices are preserved unchanged.

## Sparse triangular and direct factorization (2026-11-05)

The `spsolve_triangular`, `splu`, `splu_solve`, and `spsolve` additions in
`src/sparse_linalg.f90` are original MIT-licensed SciFort implementations.
They implement standard forward/back substitution and sparse Gaussian
elimination with partial row pivoting, retaining fill-in through a growable
sorted sparse-row representation. `src/integrate/sparse_lu.f90` already uses
related mathematical techniques for SciFort's ODE Jacobians, but the general
CSR/CSC APIs and factorization object are independent and do not import
solver code from SciPy, SuperLU, UMFPACK, or other external packages.
The Python reference comparisons invoke SciPy 1.17 solely as an independent
numerical oracle and reconstruct `P*A=L*U` with NumPy to verify factors.
New source, regression, and documentation files are SciFort contributions.

## 2026-11-04 — 1-norm estimation and CSR matrix exponential

The routines `onenormest` and `expm` added to `src/sparse_linalg.f90`
are independently written MIT-licensed SciFort code. `onenormest` uses the
published norm-estimation principle of repeated sign-vector application and
transposed matrix-vector products (Higham and Tisseur, 2000), not a copy of
SciPy's estimator implementation. It returns a lower bound with an explicit
normalized witness `(v,w=A*v)`. Exact-width probes use sparse column sums.

`expm` reuses SciFort's own scaling-and-Taylor `expm_multiply` to apply the
exponential to blocks of coordinate vectors, then packs the result into
canonical CSR form. It differs from SciPy's Padé-based sparse `expm`
algorithm. No SciPy source or other third-party code is imported. SciPy
1.17 and `scipy.linalg.expm` are used only as independent test references.

Reference: N. J. Higham and F. Tisseur, “A Block Algorithm for Matrix 1-Norm
Estimation, with an Application to 1-Norm Pseudospectra,” SIAM J. Matrix
Anal. Appl. 21(4), 1185–1201 (2000). The existing BSD-origin components
retain their attribution and license notices.

## 2026-11-04 — sparse exponential action (`expm_multiply`)

The new action routines in `src/sparse_linalg.f90`, their native regression,
reference driver, and comparison scripts are original MIT-licensed SciFort
code. They use the published scaling/Taylor matrix-exponential-action method
of Al-Mohy and Higham (2011), with independently implemented trace shifting,
a conservative induced-infinity-norm tail bound, and a bidirectional time-grid
strategy. No SciPy, EXPOKIT, Krylov, BLAS, or other library source has been
copied or translated. SciPy 1.17 `scipy.sparse.linalg.expm_multiply` and SciPy's
independent dense matrix exponential are test-only numerical references.

Reference: A. H. Al-Mohy and N. J. Higham, “Computing the Action of the Matrix
Exponential, with an Application to Exponential Integrators,” *SIAM Journal
on Scientific Computing* **33**(2), 488–511 (2011), DOI: 10.1137/100788860.
The existing BSD-licensed MINRES and LSMR components are unchanged.

## 2026-11-04 — partial sparse SVD (`svds`)

The implementation of `svds`, its independent-seed Golub--Kahan search, and
one-sided Jacobi projected singular-vector extraction in
`src/sparse_linalg.f90` are independently written MIT-licensed SciFort code.
They use standard mathematical Golub--Kahan bidiagonalization, full
reorthogonalization, and Jacobi singular-value rotations. No ARPACK,
PROPACK, LAPACK, SciPy or third-party library implementation was copied or
adapted. SciPy 1.17 `scipy.sparse.linalg.svds` and NumPy dense SVD are used
solely as *numerical reference implementations* for the randomized tests.
See `docs/sparse_svds_validation_2026-11-04.md`.

## 2026-11-04 — block-restarted nonsymmetric Arnoldi

The optional multiple-seed Arnoldi expansion and residual-block restart in
`src/sparse_linalg.f90`, and the new native regression and SciPy reference
scripts, are independently written MIT-licensed SciFort code. They apply
standard block Krylov, Rayleigh--Ritz and thick-restart mathematical ideas
within the existing SciFort Arnoldi/Schur framework. No external source from
ARPACK, SciPy, SLEPc, PRIMME, or other solver packages was imported or
translated; SciPy 1.17's dense eigenvalue solver is used only as an
independent test reference. The pre-existing BSD notices for MINRES and
LSMR are retained. See `docs/sparse_arnoldi_block_validation_2026-11-04.md`.

## 2026-11-04 — Schur invariant-subspace locking for nonsymmetric Arnoldi

The new opt-in Schur-locking path in `src/sparse_linalg.f90`, its tests and
reference comparison scripts are original SciFort MIT-licensed implementations
of established Arnoldi, Rayleigh--Ritz, Schur decomposition and invariant
subspace deflation concepts. They reuse SciFort's existing small complex QR
Schur decomposition and complex shifted GMRES. NumPy/SciPy are used **only**
as independent numerical references; no ARPACK, SciPy, LAPACK, or external
solver source has been copied. The earlier MINRES and LSMR BSD notices remain
unchanged. See `docs/sparse_arnoldi_schur_validation_2026-11-04.md`.

## 2026-11-04 — hard locking and block-seeded symmetric sparse eigenspaces

The `eigsh_lanczos_deflated` and `orthogonalize_deflated` implementations
in `src/sparse_linalg.f90` are new MIT-licensed SciFort code. They apply
two-pass reorthogonalization and orthogonal projection against an
accepted approximate symmetric invariant subspace. The optional multi-seed
first cycle uses deterministic block-Krylov expansion to represent
repeated eigenspaces. No SciPy, ARPACK, PRIMME, TRLan, or third-party
source is imported, translated, or reproduced. The general mathematical
background is Kesheng Wu and Horst Simon, "Thick-Restart Lanczos Method
for Large Symmetric Eigenvalue Problems," *SIAM Journal on Matrix
Analysis and Applications* **22**(2), 602--616 (2000), DOI:
[10.1137/S0895479898334605](https://doi.org/10.1137/S0895479898334605).
Source and tests remain independently authored; SciPy/NumPy are used only
as external numerical comparison oracles. Preexisting BSD licenses for
MINRES and LSMR remain unchanged.

## 2026-11-04 — shifted restarted sparse eigenvalue extensions

The added normal-mode shift-invert paths in `eigsh_lanczos` and
`eigs_arnoldi`, their dispatcher forwarding and standalone regression and
SciPy/NumPy comparisons are new MIT-licensed SciFort code. They compose
SciFort's existing shifted MINRES and complex restarted-GMRES inner solvers
with thick-restarted Krylov/Rayleigh--Ritz projection. No third-party
implementation was imported or adapted. The mathematical shifted ranking
`mu=1/(lambda-sigma)` follows the public SciPy `eigsh`/`eigs` normal-mode
specifications. For a real matrix with complex shift, SciFort uses the full
complex inverse rather than SciPy's default real-matrix `OPpart='r'`.
Prior BSD-3-Clause MINRES and LSMR notices remain in place.

## 2026-11-04 — symmetric thick-restarted Lanczos and solver method selection

`src/sparse_linalg.f90` includes the new independently written
`eigsh_lanczos` and `lanczos_orthogonalize` procedures. They implement
standard symmetric Krylov/Rayleigh--Ritz methods with two-pass full
reorthogonalization, explicit Ritz-vector retention, and soft convergence
locking. The restart design follows general thick-restart eigensolver
principles; this is not ARPACK's implicit QR restarted Lanczos or a
hard-deflated invariant-subspace method. No source was copied from SciPy,
ARPACK, PRIMME, or any external eigensolver.

The optional `method` arguments on the pre-existing `eigsh` and `eigs`
interfaces select the existing implementations or these two restarted
solvers; default behavior remains unchanged. Independent test-only oracles
use NumPy/SciPy 1.17; these are not runtime dependencies. Existing SciPy
MINRES and SciFort LSMR BSD notices remain preserved.

## 2026-11-04 — thick-restarted Arnoldi with retained Ritz vectors

`src/sparse_linalg.f90` now contains `eigs_arnoldi` and
`arnoldi_orthogonalize`. These are new MIT-licensed SciFort implementations
of two-pass Arnoldi expansion, Rayleigh--Ritz extraction, and explicit Ritz
vector retention with soft locking; existing SciFort small projected
complex-Schur routines are reused. The algorithm was developed from standard
mathematical descriptions rather than adapted source code. The relevant
restart motivation is Ronald B. Morgan, "On Restarting the Arnoldi Method for
Large Nonsymmetric Eigenvalue Problems," *Mathematics of Computation* 65
(1996), 1213--1230, DOI:10.1090/S0025-5718-96-00745-4.

No source was imported from ARPACK, SciPy, Morgan's publications, or
third-party eigensolvers. SciPy 1.17 `scipy.linalg.eigvals` and
`scipy.sparse.linalg.eigs` appear in **test-only** comparison scripts.
This implementation is not Morgan's or Sorensen's implicitly restarted
QR algorithm; converged Ritz vectors are retained but are not hard-deflated.

## 2026-11-04 — complex GMRES shift-invert and transformed eigenvalue selection

The complex restarted-GMRES inner solver and shift-invert extensions of
`eigs`/`eigsh` in `src/sparse_linalg.f90` are new MIT-licensed SciFort code
implementing standard Arnoldi/Givens and subspace iteration techniques.
They do not copy ARPACK, SciPy, LAPACK, kiss, or other upstream eigenvalue
implementations. Independent Python test oracles use NumPy dense eigenvalues
and SciPy 1.17 `eigs`/`eigsh`; these tests are not runtime dependencies.
The preexisting BSD-3-Clause MINRES attribution and LSMR provenance remain
unchanged. This update does **not** implement implicitly restarted Arnoldi.

## 2026-11-04 — nonsymmetric sparse eigenpairs and symmetric shift-invert

The `eigs` implementation in `src/sparse_linalg.f90`, its small complex
Householder/Hessenberg/shifted-QR Schur routines, tests in
`test/test_sparse_eigs.f90`, and `tools/sparse_linalg_reference/` compare
programs are independently written MIT-licensed SciFort source implementing
standard QR/Rayleigh--Ritz numerical methods. No ARPACK, LAPACK, SciPy, kiss,
or other external eigenvalue solver source was copied or linked. The expanded
`eigsh` shift-invert option calls the existing `minres` routine; the MINRES
BSD-3-Clause attribution and license remain preserved in
`LICENSES/scipy-minres-BSD-3-Clause.txt`. NumPy and SciPy 1.17 are used as
test oracles only, not build/runtime dependencies. This is not a translation
of ARPACK's implicit restart algorithm, and its limitations are documented in
`docs/sparse_linalg.md`.

## 2026-11-04 — symmetric MINRES and real sparse Ritz eigenpairs

The new Fortran `minres` entry point in `src/sparse_linalg.f90` implements
Paige and Saunders, "Solution of Sparse Indefinite Systems of Linear Equations",
*SIAM Journal on Numerical Analysis* 12 (1975). Its streamed Lanczos/Givens
recurrences were written for SciFort's one-based CSR storage with SciPy's
BSD-3-Clause `scipy/sparse/linalg/_isolve/minres.py` (v1.17.0) consulted as a
procedural reference. As a conservative attribution measure, the SciPy
copyright, conditions and disclaimer are retained in
`LICENSES/scipy-minres-BSD-3-Clause.txt` and `THIRD_PARTY_LICENSES.md`.
Changes include explicit CSR storage, shifted true-residual convergence,
positive Jacobi validation, and independent status/validation behavior.

The `eigsh` subset, its small projected Jacobi eigensolver, native regression
tests and Python comparison drivers are new original SciFort MIT-licensed
implementations. The symmetric sparse eigenvalue implementation uses
orthogonal block subspace iteration and Rayleigh--Ritz projection.
SciPy 1.17 APIs and eigenvalues are used as a numerical compatibility oracle;
**no ARPACK/kiss source code was imported**. The
public `eigsh` subset deliberately does not claim ARPACK-equivalent iterative
restarting or convergence for clustered eigenvalues. The existing LSMR BSD
attribution and linker caveat remain unchanged.

## 2026-11-04 — independently implemented LSQR and incomplete LU(0)

`src/sparse_linalg.f90` was extended with an independent SciFort MIT-licensed
implementation of LSQR's Golub--Kahan bidiagonalization recurrences and a
zero-fill incomplete LU factorization for canonical CSR matrices. The
optional ILU(0) BiCGSTAB/GMRES preconditioning paths, native regression
`test/test_sparse_lsqr_ilu0.f90`, and files in
`tools/sparse_linalg_reference/` are also original SciFort code. The primary
mathematical LSQR reference is C. C. Paige and M. A. Saunders, "LSQR: An
Algorithm for Sparse Linear Equations and Sparse Least Squares", *ACM
Transactions on Mathematical Software* **8** (1982). Standard incomplete
LU(0) recurrences are used with no fill or pivoting. The SciPy 1.17 public
`lsqr` interface and numerical outputs were consulted for compatibility,
but **no SciPy, jacobwilliams/LSQR, kiss, or other third-party source code was
copied or adapted**. ILU(0) results were independently verified against a
separate Python implementation of the textbook no-fill row updates, not
against SciPy's more sophisticated SuperLU `spilu`. No new copyright/license
notices, mandatory dependencies, or BLAS/LAPACK linking were required.

## 2026-11-04 — sparse CSR iterative solvers and LSMR bridge

`src/sparse_linalg.f90`, `test/test_sparse_linalg.f90`,
`tools/sparse_linalg_reference/`, and `docs/sparse_linalg.md` are new SciFort
components. CG, BiCGSTAB, and GMRES follow the published mathematical
recurrences, implemented independently in MIT-licensed source. Algorithmic
references: Hestenes & Stiefel, "Methods of Conjugate Gradients for Solving
Linear Systems", *J. Research NBS* **49** (1952); H. A. van der Vorst,
"Bi-CGSTAB: A Fast and Smoothly Converging Variant of Bi-CG", *SIAM J.
Scientific and Statistical Computing* **13** (1992); Y. Saad & M. H. Schultz,
"GMRES: A Generalized Minimal Residual Algorithm for Solving Nonsymmetric
Linear Systems", *SIAM J. Scientific and Statistical Computing* **7** (1986).

`lsmr` is a CSR callback adapter of SciFort's **existing**
`src/optimize/lsmr.f90` kernel (BSD-3-Clause; copyrights and earlier SciPy /
Fong--Saunders attribution preserved there and in `THIRD_PARTY_LICENSES.md`).
No external implementation has been newly copied into this update. Public
SciPy 1.17 `scipy.sparse.linalg` serves as an independent numerical oracle.
The BSD-3-Clause project https://github.com/ivan-pi/kiss was reviewed as a
possible later reuse source; its `fpm.toml` currently specifies BLAS/LAPACK
linkage, so it was **not** imported. No new build dependency or license text
is necessary for this update.

## 2026-11-04 — CSC storage, diagonal construction, triangular extraction

The extensions to `src/sparse.f90`, `test/test_sparse_csc_diagonals.f90`,
`tools/sparse_reference/compare_sparse_csc_scipy.py`, and associated driver
are independent MIT-licensed implementations of standard CSC compressed-column
storage, prefix-sum CSR/CSC conversion, packed diagonal construction, sparse
triangular extraction, and diagonal replacement. SciPy 1.17 serves exclusively
as a numerical test oracle. No code was copied, adapted or translated from
FSPARSE, kiss, LSQR, LSMR, LUSOL, SPOOF, SPRAL, or any other external source.
Consequently there is no new vendored license, build dependency or linked
runtime library. Relevant upstream projects and their licensing are surveyed
separately in `docs/sparse_reuse_repositories.md` for prospective reuse.

## 2026-11-03 — generic sparse matrix arithmetic

`src/sparse.f90`, `test/test_sparse_matrix.f90`,
`tools/sparse_reference/`, and `docs/sparse.md` are independent MIT-licensed
implementations of conventional compressed sparse row operations. COO
construction uses row grouping, heapsort of column keys and aggregation;
CSR matrix multiplication uses a two-pass row-wise sparse accumulator.
SciPy 1.17's public CSR numerical results serve as independent reference
values in randomized tests; no SciPy implementation or other third-party
source code was copied. No added runtime dependencies or licenses.

## 2026-11-03 — Yen K-shortest simple paths

The new public `yen` routine and private spur-search/candidate-heap utilities
in `src/sparse_csgraph.f90` implement the standard Yen loopless-path
algorithm independently in modern Fortran. Johnson reweighting is reused
from SciFort's existing Bellman–Ford potential helper and binary Dijkstra
heap; no SciPy source code was copied. Reference tests use the public
`scipy.sparse.csgraph.yen` API in SciPy 1.17, with independent exhaustive
DFS enumeration as a separate correctness oracle. New files include
`test/test_sparse_csgraph_yen.f90`,
`tools/csgraph_reference/csgraph_yen_driver.f90`, and
`tools/csgraph_reference/compare_csgraph_yen_scipy.py`.
The MIT license and Fortran runtime dependency set are unchanged.

## 2026-11-13 - Sparse weighted matching compatibility entry point

`min_weight_full_bipartite_matching` in `src/sparse_csgraph.f90` is an original
MIT-licensed forwarding procedure to the existing independently written sparse
assignment solver described below. No upstream implementation was imported or
adapted. The API contract was checked against the SciPy public specification:
<https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csgraph.min_weight_full_bipartite_matching.html>.
The optional total-cost output, stored zero edges, and parallel-edge selection
are documented SciFort extensions. The new exhaustive tests independently
enumerate matchings; the existing randomized SciPy comparison driver now uses
the compatibility entry point. No algorithm, dependency, or license changed.

## 2026-11-03 — sparse weighted matching and CSR Laplacians

The new implementation in `src/sparse_csgraph.f90`, native regression
`test/test_sparse_csgraph_weighted.f90`, and the driver/comparison scripts
in `tools/csgraph_reference/` was written independently for SciFort.
It implements a generic primal-dual successive-shortest-augmenting-path
method for weighted assignment and a row-wise sparse accumulator for graph
Laplacians. SciPy 1.17's public numerical behavior serves as an external
reference, not as copied or translated code. No new runtime dependencies,
vendored code or licenses were added.

## 2026-11-03 — graph ordering, flow and bipartite matching

`src/sparse_csgraph.f90`, `test/test_sparse_csgraph_advanced.f90`,
`tools/csgraph_reference/csgraph_advanced_driver.f90`, and
`tools/csgraph_reference/compare_csgraph_advanced_scipy.py` extend the
independently implemented MIT-licensed SciFort graph module. Reverse
Cuthill–McKee ordering uses structural BFS from pseudo-peripheral roots,
Dinic maximum flow uses level graphs and residual arcs, and iterative
Hopcroft–Karp matching uses alternating paths. SciPy 1.17 and its public
numerical API are used as behavioral references and independent test oracles;
no SciPy source or third-party implementation is copied. No new runtime
packages are introduced.

## 2026-11-03 — compressed sparse graph algorithms

`src/sparse_csgraph.f90`, `test/test_sparse_csgraph.f90`,
`tools/csgraph_reference/`, and `docs/sparse_csgraph.md` are independently
written MIT-licensed numerical algorithms and tests. The module uses standard
CSR storage, iterative graph traversal, Kosaraju's strongly connected
components algorithm, binary-heap Dijkstra, Bellman-Ford edge relaxation,
and Kruskal's union-find minimum-spanning-forest algorithm. The seven named
linkage wrappers in `src/cluster_hierarchy.f90` reuse SciFort's existing
Lance-Williams implementation. The public mathematical/API documentation of
SciPy's `scipy.sparse.csgraph` and `scipy.cluster.hierarchy` served as a
behavioral reference, and SciPy 1.17 as an independent numerical oracle. No
SciPy implementation or third-party code was copied or vendored. No new
licenses or runtime dependencies were introduced.

## 2026-11-03 — hierarchy visualization data and Gaussian K-means initialization

The MIT-licensed additions to `src/cluster_hierarchy.f90` and
`src/cluster_vq.f90`, associated tests and updated documentation were
independently implemented from the public SciPy mathematical and API
specifications. Numeric dendrogram coordinates are computed by iterative
postorder traversal and the standard ten-unit horizontal spacing.
`is_isomorphic` uses bidirectional equivalence-class checking with stable
merge sorting; MATLAB linkage conversion reconstructs cluster leaf counts.
Gaussian K-means initialization samples an exact empirical-covariance
Gaussian distribution by projecting standard normal variates through the
centered observation matrix. No SciPy or other third-party source has been
copied, translated, linked or vendored; SciPy 1.17 serves only as a numerical
reference. No licenses or external dependencies are added.

## 2026-11-03 — optimal hierarchical leaf ordering and tree helpers

The additions to `src/cluster_hierarchy.f90`,
`test/test_cluster_hierarchy_ordering.f90`, and
`tools/cluster_reference/{cluster_hierarchy_ordering_driver.f90,
cluster_hierarchy_leaders_driver.f90,compare_hierarchy_ordering_scipy.py,
compare_hierarchy_leaders_scipy.py}` are an independent MIT-licensed
implementation of endpoint-pair dynamic programming for optimal leaf
orientation, array-based tree materialization, and bottom-up subtree-label
propagation. SciPy's public `optimal_leaf_ordering`, `to_tree`, `leaders`,
and `correspond` API contracts are used for feature and numerical
comparisons. No third-party implementation source was copied, converted,
or included. The solver stores separators and discards intermediate costs
to bound working storage by O(n^2). Documentation and license requirements
remain unchanged.

## 2026-11-03 — hierarchical diagnostics, monotone cuts and cut trees

The extensions in `src/cluster_hierarchy.f90`,
`test/test_cluster_hierarchy.f90`,
`tools/cluster_reference/cluster_hierarchy_diagnostics_driver.f90`, and
`tools/cluster_reference/compare_hierarchy_diagnostics_scipy.py` were
independently implemented in MIT-licensed Fortran. The mathematical algorithms
use bottom-up subtree statistics, depth-limited sample standard deviation,
monotone merge criteria, and union-find partition formation, following the
public mathematical API specification of `scipy.cluster.hierarchy`. SciPy
1.17 is a numerical validation oracle, not a source-code dependency. No
third-party source code was copied or translated, and no licenses/dependencies
were added. See `docs/cluster_hierarchy.md` for API and numerical differences.

## 2026-11-03 — independent vector quantization and Lloyd K-means

`src/cluster_vq.f90`, `test/test_cluster_vq.f90`,
`tools/cluster_reference/{cluster_vq_driver.f90,run_native_vq.py,compare_vq_scipy.py}`,
and `docs/cluster_vq.md` are new MIT-licensed SciFort implementation,
tests and documentation. The algorithms were implemented independently from
Euclidean nearest-centroid assignment, population-standard-deviation
normalization, Lloyd's iterative centroid updates and the K-means++
squared-distance initialization rule. Design references:

- S. P. Lloyd (1982), "Least squares quantization in PCM," *IEEE
  Transactions on Information Theory* 28(2), 129-137.
- D. Arthur and S. Vassilvitskii (2007), "k-means++: The Advantages of
  Careful Seeding," *Proceedings of SODA 2007*, 1027-1035.
- SciPy `scipy.cluster.vq` mathematical API documentation:
  https://docs.scipy.org/doc/scipy/reference/cluster.vq.html .

SciPy 1.17 is used as an independent numerical test oracle only. No SciPy,
Burkardt, FSML or other external library source was copied, translated or
adapted; all randomness uses existing `scifort_random` APIs. No additional
FPM, binary or license dependency is introduced.

## 2026-11-03 — independent agglomerative hierarchical clustering

`src/cluster_hierarchy.f90`, `test/test_cluster_hierarchy.f90`,
`tools/cluster_reference/` and `docs/cluster_hierarchy.md` are original
MIT-licensed SciFort implementation, tests and documentation. The algorithms
implement the published Lance-Williams cluster-distance recurrences and
SciPy 1.17 `linkage`/`fcluster` mathematical API specifications:
https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.cluster.hierarchy.linkage.html .
For general context, see Murtagh and Contreras (2012), *Algorithms for
hierarchical clustering: an overview*, WIREs Data Mining and Knowledge
Discovery, 2, 86-97. SciPy is used only as an independent numerical test
oracle; no SciPy or other third-party source was imported, translated or
adapted. The implementation reuses existing SciFort pairwise distances and
adds no external numerical dependencies.

## 2026-11-03 — independent arbitrary-rank transform broadcasting

`src/spatial_transform_nd.f90`, `test/test_spatial_transform_nd.f90`,
`tools/spatial_reference/{rotation_nd_driver.f90,compare_rotation_nd_scipy.py}`,
and `docs/spatial_transform_nd.md` are original SciFort MIT-licensed code.
The implementation reuses SciFort's own scalar SO(3) and SE(3) kernels; its
shape-aware broadcasting and flattened component representation were authored
independently. SciPy 1.17 is used only as an independent numerical reference
for orientation/transform composition, point application, and batched shapes.
No SciPy or external Fortran source was incorporated and no dependencies
were added.

## 2026-11-03 — generalized Davenport-angle rotations

`src/spatial_transform.f90`, `src/spatial_rotation_stack.f90` and
`src/spatial_rotation_grid.f90` implement the established generalized
Euler/Davenport-angle algorithm. The quaternion half-angle extraction was
cross-checked/adapted against SciPy 1.17 `_rotation_xp.py`, BSD-3-Clause,
with the copyright/license notice in
`LICENSES/scipy-spatial-transform-BSD-3-Clause.txt`; the Fortran
implementation, input validation, and batched wrappers were written here.
See Shuster & Markley (2003), *Generalization of the Euler Angles*, and
Bernardes & Viollet (2022), *Quaternion to Euler angles conversion*.
Native tests and SciPy comparisons are in `test/test_spatial_rotation_davenport.f90`,
`tools/spatial_reference/rotation_davenport_driver.f90`, and
`tools/spatial_reference/compare_rotation_davenport_scipy.py`.

## 2026-11-03 — independent rank-two SE(3) batch operations

`src/spatial_rigid_transform_grid.f90`,
`test/test_spatial_rigid_transform_grid.f90`,
`tools/spatial_reference/{rigid_grid_driver.f90,compare_rigid_grid_scipy.py}`
and `docs/spatial_rigid_transform_grid.md` are new original MIT-licensed
SciFort implementations and tests. They reuse the existing scalar SE(3)
transforms and the SciFort quaternion chordal mean, implementing standard
component-last batch array operations and NumPy-compatible two-axis
broadcasting. Independent numerical comparisons use only the public
SciPy 1.17 `RigidTransform` API as a reference:
https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.spatial.transform.RigidTransform.html .
No SciPy code or third-party numerical implementation was incorporated.

## 2026-11-02 — independent SE(3) rigid transformations

`src/spatial_rigid_transform.f90`, `test/test_spatial_rigid_transform.f90`,
`tools/spatial_reference/{rigid_transform_driver.f90,rigid_twist_driver.f90,
compare_rigid_transform_scipy.py}`
and `docs/spatial_rigid_transform.md` are original MIT-licensed SciFort code.
The implementation uses the standard SE(3) group law, SO(3) left Jacobian
and its inverse, Hamilton dual-quaternion arithmetic, and a weighted
translation mean combined with SciFort's previous rotation mean. Numerical
comparisons use the public SciPy 1.17.0 `RigidTransform` API as a reference;
no SciPy code was copied and no third-party Fortran implementation used.
https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.RigidTransform.html

## 2026-11-02 — native finite rotation groups and symmetry-aware reduction

`src/spatial_rotation_group.f90`, `test/test_spatial_rotation_group.f90`,
`tools/spatial_reference/{rotation_group_driver.f90,compare_rotation_group_scipy.py}`
and `docs/spatial_rotation_group.md` are original MIT-licensed SciFort source.
Cyclic/dihedral groups are composed from standard axial rotations and
perpendicular half turns. Tetrahedral, octahedral and icosahedral groups use
the classical 24-cell/600-cell unit-quaternion construction with antipodal
identification, implemented independently. Reduction uses standard SO(3)
quaternion composition and minimum principal-angle selection. SciPy 1.17.0
is used exclusively as an external numerical oracle, based on documented
`Rotation.create_group` and `Rotation.reduce` behavior. No SciPy or third-party
implementation source was copied or adapted; no new dependency was added.

## 2026-11-02 — random rotations and rank-two spatial rotation arrays

`src/spatial_rotation_grid.f90`, extensions to `src/spatial_rotation_stack.f90`,
`test/test_spatial_rotation_random_grid.f90`,
`tools/spatial_reference/{rotation_grid_driver.f90,compare_rotation_grid_scipy.py}`,
and `docs/spatial_rotation_grid.md` are original MIT-licensed SciFort code.
The Shoemake uniform-quaternion construction and Park–Miller integer PRNG
are established mathematical algorithms implemented independently, without
copying any external Fortran or SciPy source. Comparisons use independently
computed mathematical references and the public SciPy 1.17.0 interface.

## 2026-11-02 — independently implemented batched rotation stacks

`src/spatial_rotation_stack.f90`, `test/test_spatial_rotation_stack.f90`,
`tools/spatial_reference/{rotation_stack_driver.f90,compare_rotation_stack_scipy.py}`
and `docs/spatial_rotation_stack.md` are original, MIT-licensed SciFort additions.
The stack wraps the existing native quaternion type and implements standard
SO(3) operations and 1D broadcasting semantics without copying SciPy's
implementation. SciPy 1.17.0 is used only as an independent numerical oracle.
Public behavior was checked against the SciPy Rotation documentation:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.html .
No external libraries or third-party implementation code were added.

## 2026-11-02 — smooth SO(3) rotation spline

The new `src/spatial_rotation_spline.f90`, native regression and independent
SciPy comparison driver implement the mathematical cubic SO(3) attitude spline
formulation. The angular-rate Jacobian, nonlinear angular-acceleration
correction, and iterative knot-rate equations follow the published *Smooth
Attitude Interpolation* method and SciPy 1.17.0's BSD-3-Clause
`scipy/spatial/transform/_rotation_spline.py`:
https://github.com/scipy/scipy/blob/v1.17.0/scipy/spatial/transform/_rotation_spline.py .
The Fortran implementation, interface, block-tridiagonal solver and tests
are written for SciFort; SciPy serves as the external numerical oracle.
SciPy's original licensing and copyright notices are preserved in
`THIRD_PARTY_LICENSES.md` under the existing SciPy BSD-3-Clause section.
No SciPy/Python runtime or new third-party binary is required.

## 2026-11-02 — independently implemented rotation alignment and means

The `rotation_align_vectors` and `rotation_mean` additions to
`src/spatial_transform.f90`, the native
`test/test_spatial_rotation_estimation.f90`, and the independent comparison
harness `tools/spatial_reference/{rotation_estimation_driver.f90,
compare_rotation_estimation_scipy.py}` are original MIT-licensed SciFort work.
The mathematical constructions use Davenport's quaternion formulation of
Wahba's least-squares rotation problem, the standard shortest-arc alignment
of a pair of vectors, a constrained axial optimization for a fixed primary
vector, and the weighted quaternion outer-product / Markley chordal mean.
A small symmetric Jacobi eigensolver was implemented locally without adapting
third-party numerical source code. SciPy 1.17.0 was used only as an independent
numerical reference, with the documented contracts available at
https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.align_vectors.html
and
https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.mean.html .
The addition introduces no third-party implementation, dependency, or license.

## 2026-11-02 — independently implemented Euler/MRP rotation conversions

Additions to `src/spatial_transform.f90`, the native regression
`test/test_spatial_rotation_euler_mrp.f90`, and the independent reference
`tools/spatial_reference/{rotation_euler_mrp_driver.f90,compare_rotation_euler_mrp_scipy.py}`
are original MIT-licensed SciFort code. Euler extraction uses standard active
rotation-matrix identities for proper Euler and Tait–Bryan angles, using the
reverse-axis relation between intrinsic and extrinsic conventions. Modified
Rodrigues parameters are derived from the standard quaternion stereographic
projection formula. No implementation was copied or adapted from SciPy or
third-party libraries. SciPy 1.17 serves exclusively as a numerical oracle.
Public semantics were checked against
https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.spatial.transform.Rotation.from_euler.html
and
https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.spatial.transform.Rotation.as_euler.html .

## 2026-11-02 — independent spatial 3D rotations and Slerp

`src/spatial_transform.f90`, `test/test_spatial_rotation.f90`,
`test/test_spatial_rotation_slerp.f90`, `docs/spatial_transform.md`, and the
`tools/spatial_reference/rotation*` and `compare_rotation_scipy.py` test files
are original MIT-licensed SciFort code. The implementation follows standard
Hamilton quaternion multiplication, axis/angle rotation geometry, and the
shortest-path unit-quaternion interpolation described by Ken Shoemake,
"Animating rotation with quaternion curves," SIGGRAPH 1985,
https://doi.org/10.1145/325334.325242 .

The public API semantics were checked against SciPy 1.17 reference docs:
https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.spatial.transform.Rotation.html
and https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.spatial.transform.Slerp.html .
SciPy 1.17.0 was used only for independently comparing numerical outputs;
its implementation code was not adapted, copied or vendored. This addition
introduces no third-party source code, new licenses, or external dependencies.

## 2026-11-02 — 3D convex-hull validation and near-coplanar recovery

The modifications to `src/spatial_convex_hull.f90`, the new native
`test/test_spatial_convex_hull_nearcoplanar.f90`, the data fixture, and
`tools/spatial_reference/compare_convex_hull_nearcoplanar_scipy.py` are
independently written, MIT-licensed SciFort code. The methods use standard
convex-polytope supporting-plane tests, Euler triangulation identities,
permutation-independent convex hull definitions, and planar monotone-chain
triangulation, not source imported or adapted from Qhull or SciPy.
SciPy's `scipy.spatial.ConvexHull` (Qhull backend) served only as an
independent numerical reference for facet geometry, area, volume, and
near-coplanar cases. The compatibility reference is
https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html .
No third-party implementation source or new external dependency was added.

## 2026-11-02 — native bounded halfspace intersection geometry

`src/spatial_halfspace.f90`, `test/test_spatial_halfspace.f90`,
`tools/spatial_reference/{halfspace_driver.f90,compare_halfspace_scipy.py}`,
and `docs/spatial_halfspace.md` are independently written MIT-licensed
SciFort code. The mathematical construction uses the standard polar duality
between a bounded intersection of halfspaces and a convex hull, with candidate
vertices obtained by solving pairs/triples of active affine constraints and
checked by direct inequality substitution. The existing SciFort planar convex
hull code triangulates dual faces. Public API behavior was checked against
https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.HalfspaceIntersection.html .
SciPy 1.17 and Qhull served only as external numerical references. No source
code from SciPy, Qhull, or third-party Fortran libraries was copied; no
additional runtime dependency or license was introduced.

## 2026-11-02 — original 2D/3D convex hull implementation

`src/spatial_convex_hull.f90`, `test/test_spatial_convex_hull.f90`, and
`tools/spatial_reference/{convex_hull_driver.f90,compare_convex_hull_scipy.py}`
are independently written, MIT-licensed SciFort source and tests. The 2D
algorithm uses the standard monotone-chain construction; the 3D implementation
extends the geometric incremental hull approach already used internally by
SciFort's spherical Voronoi routine. No third-party implementation has been
copied. The compatibility contract was checked against the public SciPy
`ConvexHull` reference documentation:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html .
SciPy 1.17.0 was used only as an independent numerical oracle. This addition
requires no new runtime libraries or license notices.

## 2026-11-02 — original spherical Voronoi geometry

`src/spatial_voronoi.f90`, `test/test_spatial_voronoi.f90`,
`tools/spatial_reference/{voronoi_driver.f90,compare_voronoi_scipy.py}` and
`docs/spatial_voronoi.md` are new, independent MIT-licensed SciFort work.
Construction follows standard convex-hull/spherical-Voronoi duality;
implementation uses an original incremental outward-oriented triangular
hull builder and edge-adjacency ordering, with spherical triangle solid
angles using the geometric `2*atan2` formula. Public SciPy API semantics
and 2D hull-edge area behavior were checked against
https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.SphericalVoronoi.html .
SciPy 1.17.0 serves solely as an independent validation oracle; no SciPy,
fishman, Qhull or other third-party implementation code was copied or
vendored, and no external build dependencies were added.

## 2026-11-02 -- Procrustes similarity fits and orthogonal transformations

`src/spatial_procrustes.f90`, `test/test_spatial_procrustes.f90`, and the
independent oracle files `tools/spatial_reference/{procrustes_driver.f90,
compare_procrustes_scipy.py}` are newly written MIT-licensed SciFort work.
The mathematical objective and orthogonal Procrustes polar-factor formula
follow Peter H. Schonemann, "A generalized solution of the orthogonal
Procrustes problem," *Psychometrika* 31 (1966),
https://doi.org/10.1007/BF02289451 . The one-sided Jacobi iteration,
scaling/centering safeguards, deterministic null-space completion, status
handling, and Fortran interfaces were implemented independently for SciFort.

The public compatibility contracts were checked against
https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.procrustes.html
and
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.orthogonal_procrustes.html .
SciPy 1.17.0 was used solely as a numerical test oracle, with no SciPy or
other third-party implementation source copied, adapted, or vendored. No
new license notices or build dependencies were necessary.

## 2026-11-02 — Hausdorff point-cloud and spherical interpolation algorithms

`src/spatial_hausdorff.f90`, `src/spatial_slerp.f90`,
`test/test_spatial_{hausdorff,slerp}.f90` and the spatial geometry comparison
harness are newly written, MIT-licensed SciFort code. The Hausdorff
implementation follows the mathematical max-min Euclidean definition with
safe scaled norms and deterministic early abandonment; see A. A. Taha and
A. Hanbury, "An efficient algorithm for calculating the exact Hausdorff
distance," IEEE TPAMI 37(11), 2015, DOI:10.1109/TPAMI.2015.2408351
(early break motivation). The spherical interpolation implementation uses
the standard great-circle / geometric SLERP formula, rather than copying
SciPy's implementation. See Ken Shoemake, "Animating rotation with
quaternion curves," SIGGRAPH 1985, DOI:10.1145/325334.325242.

SciPy 1.17.0 is used as a **numerical test oracle only**. Its documented
API behavior and references are found at
https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.directed_hausdorff.html
and https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.geometric_slerp.html .
No implementation source or other third-party code was incorporated; no
new license obligations or dependencies were introduced.

## 2026-11-02 — weighted Boolean `scipy.spatial.distance` metrics

Changes to `src/spatial_distance.f90`,
`test/test_spatial_distance_boolean.f90`, and
`tools/spatial_reference/compare_boolean_distance_scipy.py` are independent
MIT-licensed implementations and test harnesses for the standard Boolean
contingency-table definitions of Dice, Rogers--Tanimoto, Russell--Rao,
Sokal--Sneath, and Yule dissimilarities. Mathematical conventions and
interface behavior were checked against the public SciPy 1.17 documentation:
https://docs.scipy.org/doc/scipy/reference/spatial.distance.html
SciPy's installed API was used exclusively as a numerical oracle; no SciPy
implementation source or third-party code was incorporated. The existing
library dependencies and provenance requirements remain unchanged.

## 2026-11-02 — covariance and Jensen--Shannon spatial metrics

Changes to `src/spatial_distance.f90`, `test/test_spatial_distance_covariance.f90`,
`tools/spatial_reference/covariance_distance_driver.f90`, and
`tools/spatial_reference/compare_covariance_distance_scipy.py` are original
MIT-licensed SciFort implementations of the standardized Euclidean,
Mahalanobis, and Jensen--Shannon mathematical definitions. Sample variance,
correlation-scaled Cholesky covariance inversion, and probability normalization
are implemented independently; no SciPy or third-party implementation code
was copied. SciPy 1.17.0 serves only as the numerical validation reference.
Public formula/API references:
https://docs.scipy.org/doc/scipy/reference/spatial.distance.html
No new build dependency or third-party license is introduced.

## 2026-11-01 — weighted and multiple-radius KD-tree neighbor counts

`src/interpolate_kdtree.f90`, `test/test_spatial_kdtree_weighted_count.f90`,
and `tools/spatial_reference/{kdtree_weighted_count_driver.f90,compare_kdtree_weighted_count_scipy.py}`
contain new, independently implemented MIT-licensed SciFort code. The
vector-radius difference-array traversal, subtree-product short-circuits,
site-weight aggregation and Fortran overloads were developed independently.
SciPy 1.17.0's `scipy.spatial.cKDTree.count_neighbors` serves as a test
oracle only. No upstream implementation source, license, or dependency has
been copied or introduced.

## 2026-11-01 — Periodic boxes and Minkowski metrics for KD-trees

Changes to `src/interpolate_kdtree.f90` implement independent, MIT-licensed
periodic/toroidal coordinate geometry and scaled Minkowski distance bounds.
The scalar/vector box constructors, exact search/pruning logic, native tests
in `test/test_spatial_kdtree_periodic.f90`, and Python/Fortran comparison
drivers in `tools/spatial_reference/` are original SciFort code.  The SciPy
1.17.0 `cKDTree` API serves solely as a numerical reference; no SciPy, FLANN,
or other third-party KD-tree implementation source was copied. All existing
dependency and licensing requirements are unchanged.

## 2026-11-01 — Native KD-tree pair and cross-tree spatial search

`src/interpolate_kdtree.f90` extends SciFort's original MIT-licensed
balanced KD-tree with independent `query_pairs`, `query_ball_tree`,
`count_neighbors`, and `sparse_distance_matrix` implementations.
`src/spatial.f90` exposes the type through the public spatial facade.
Counting uses an independent dual-tree bounding-box algorithm with subtree
cardinalities; the enumeration APIs use native closed-ball searches and
CSR/COO storage. `test/test_spatial_kdtree_pairs.f90` and
`tools/spatial_reference/compare_kdtree_pairs_scipy.py` perform original
native and SciPy 1.17 reference checks. No third-party source code was
copied, and SciFort's MIT licensing and dependencies are unchanged.

## 2026-11-01 — `scipy.spatial.distance` scalar and pairwise metrics

`src/spatial_distance.f90`, `src/spatial.f90`,
`test/test_spatial_distance.f90` and `tools/spatial_reference/` are new
independent MIT-licensed SciFort implementations. Formulas for norms,
weighted cosine and correlation, Boolean Hamming/Jaccard distances,
Canberra and Bray--Curtis dissimilarities, and condensed row-pair indexing
are standard mathematical definitions described by the SciPy distance
API specification (https://docs.scipy.org/doc/scipy/reference/spatial.distance.html).
No SciPy, Qhull or third-party implementation source was copied, adapted,
or vendored; no new runtime or build dependency was introduced. The native
implementation is compared numerically with SciPy 1.17.0 using the included
Python reference runner, including documentation of intentional differences
in invalid-input handling.

## 2026-11-01 — periodic interpolating B-splines

The periodic knot-extension construction and endpoint-derivative collocation
convention in `src/interpolate_bspline.f90` were adapted from
SciPy 1.17.0's `scipy/interpolate/_bsplines.py`, specifically the
`_periodic_knots` and `_make_interp_per_full_matr` mathematical scheme:
https://github.com/scipy/scipy/blob/v1.17.0/scipy/interpolate/_bsplines.py
SciPy is BSD-3-Clause; copyright and redistribution notices are retained in
`THIRD_PARTY_LICENSES.md`. SciFort uses its pre-existing spline basis
and pivoted linear solver and implements its own Fortran array/channel
layout, validation, evaluation, and integration. No upstream source file is
vendored; no SciPy or LAPACK dependency is introduced. Independent numerical
comparison to SciPy 1.17.0 is performed by
`tools/interpolate_reference/compare_periodic_bspline_scipy.py`.

## 2026-11-01 — bivariate TCK differentiation and knot insertion

`src/interpolate_fitpack_bivariate.f90`,
`test/test_interpolate_bivariate_calculus.f90` and
`tools/interpolate_reference/{bivariate_calculus_driver.f90,compare_bivariate_calculus_scipy.py}`
add original MIT-licensed SciFort routines implementing the standard tensor
B-spline derivative identity and Boehm knot insertion through SciFort's
existing one-dimensional knot-insertion implementation. SciPy 1.17.0's
`BivariateSpline` and `BSpline` are used as numerical reference oracles only;
no SciPy or FITPACK implementation source is imported. No new dependency or
third-party licensing terms are introduced.

## 2026-11-01 — rectangular tensor spline interpolation and exact surface integration

`src/interpolate_fitpack_bivariate.f90`,
`test/test_interpolate_rect_bivariate.f90`, and
`tools/interpolate_reference/{rect_bivariate_driver.f90,compare_rect_bivariate_scipy.py}`
contain independently implemented, MIT-licensed SciFort additions. The new
`rect_bivariate_spline` constructor applies SciFort's previously implemented
1D B-spline collocation algorithm successively to a rectangular grid,
using the standard tensor-product interpolation identity. Exact rectangle
integration contracts successive one-dimensional B-spline antiderivatives.
These are standard numerical spline constructions; no new third-party
implementation code, FITPACK routines or SciPy source is incorporated.

SciPy 1.17.0's `RectBivariateSpline` serves only as an independent test
reference for knots, coefficients, derivatives and integrals. The optional
`perazz/fitpack` backend and its existing license/provenance are unchanged.

## 2026-11-01 — Complete-knot `make_lsq_spline` factory

`src/interpolate_bspline.f90` independently implements SciPy-style
`make_lsq_spline` using the existing SciFort B-spline evaluator and its
banded online Givens QR solver, shared with `lsq_univariate_spline`.
SciPy 1.17 serves strictly as a numerical reference for randomly generated
least-squares problems, not as the source of the numerical implementation.
No upstream Python, C, Cython, or Fortran code has been copied; there is no
new package dependency or separate third-party license.

## 2026-10-31 — FITPACK-style least-squares spline API

- `src/interpolate_bspline.f90` adds original banded Givens-rotation QR code
  for weighted least-squares spline coefficients with specified interior
  knots. It reuses the project's own B-spline basis evaluation.
- `interpolated_univariate_spline` reuses SciFort's prior not-a-knot spline
  implementation. No FITPACK, SciPy, or other third-party source was copied.
- SciPy 1.17 `LSQUnivariateSpline` and `InterpolatedUnivariateSpline` serve as
  numerical reference implementations only. The FPM dependencies are unchanged.

## Natural cubic smoothing splines (2026-10-31)

`src/interpolate_smoothing.f90` and its native and Python reference tests
are independently developed MIT-licensed SciFort code. The finite-difference
`Q` / roughness matrix `R` formulation, banded Cholesky factorization and
selected-inverse GCV trace are implementations of standard published
smoothing-spline mathematics, not copied from SciPy's B-spline-basis
implementation or the older Woltring Fortran package. SciPy 1.17 is
used only as an independent numerical test reference. No upstream code
was incorporated, and no dependency was added.

## Accelerated 3D tetrahedron lookup (2026-10-31)

The new AABB/BVH construction and candidate search in
`src/interpolate_scattered.f90`, `test/test_interpolate_tetra_bvh.f90`,
`benchmarks/bench_interpolate_tetra_bvh.f90`, and
`tools/interpolate_reference/benchmark_tetra_bvh.py` were independently
implemented for SciFort under the MIT license. The median-partition BVH and
bounding-box pruning are standard computational geometry techniques.
The existing tetrahedralization and barycentric formulas are unchanged.
SciPy/Qhull are used solely as independent reference implementations; no
third-party code has been copied and no additional external library is
required.

## 3D Delaunay tetrahedralization and scattered interpolation (2026-10-31)

`src/interpolate_scattered.f90`, `test/test_interpolate_scattered_3d.f90`, and
`tools/interpolate_reference/{compare_scattered_3d_scipy.py,delaunay_3d_driver.f90,check_delaunay_3d_scipy.py}`
add independently developed SciFort MIT code for the standard incremental
Bowyer–Watson Delaunay algorithm in three dimensions. No Qhull, SciPy,
third-party Fortran source, or external numerical library was copied.
SciPy's `LinearNDInterpolator`/`griddata` and Qhull's 3D Delaunay and convex
hull computations are used only as independent numerical and geometric test
references. Main-package FPM dependencies are unchanged.

## Cached Clough–Tocher patches and spatial indexing (2026-10-31)

The interpolation's existing numerical patch and gradient formulas remain
adapted from SciPy's BSD-3-Clause Cython implementation with the original
license notice retained. Precomputation of the 19 contributing cubic
Bernstein coefficients, CSR uniform-bin indexing, coordinate normalization
roundoff treatment, benchmark, and additional regression are original
SciFort work; they do not copy new upstream source or change required
external dependencies. The code still uses the independently implemented
native Delaunay triangulator rather than Qhull.

## 2D Clough–Tocher cubic scattered interpolation (2026-10-31)

`src/interpolate_clough_tocher.f90` implements a Fortran translation of the
minimum-curvature vertex-gradient algorithm and affine-invariant cubic Bézier
patch formulas from SciPy's `scipy/interpolate/_interpnd.pyx` (v1.17.0;
Copyright (C) Pauli Virtanen, 2010), distributed under BSD-3-Clause.
SciFort retains the SciPy copyright, redistribution conditions, and disclaimer
in `LICENSES/scipy-interpnd-BSD-3-Clause.txt`; this adaptation is **not**
represented as an independent MIT-only implementation. SciFort's triangulation,
mesh checks, data types and API wrappers remain original contributions.
The `griddata(method='cubic')` public dispatch and test/reference drivers are
new SciFort code. The only upstream source reused in this increment is the
identified BSD-3-Clause numerical algorithm; no Qhull, BLAS/LAPACK, or external
Fortran package is copied or linked.

## Balanced KD-tree interpolation index (2026-10-31)

`src/interpolate_kdtree.f90`, `test/test_interpolate_kdtree.f90`, and
`tools/interpolate_reference/{kdtree_driver.f90,compare_kdtree_scipy.py,benchmark_kdtree.f90,run_kdtree_benchmark.py}`
are independently developed under SciFort's MIT license. The median
partitioning, axis-aligned bounding boxes, max-heap k-nearest selection,
and radius-query pruning are standard algorithms. The module is used by
`src/interpolate_scattered.f90` and `src/interpolate_rbf.f90` and publicly
exposed through `src/interpolate.f90`. SciPy's `cKDTree` is a reference for
numerical outputs only; no SciPy, Qhull, or third-party Fortran implementation
has been copied. No new mandatory runtime dependency was introduced.

## Piecewise extension, degree-selected Hermite interpolation and root finding review (2026-10-31)

`src/interpolate_piecewise.f90`,
`test/test_interpolate_piecewise_extend_orders.f90`, and
`tools/interpolate_reference/{piecewise_extend_driver.f90,bernstein_orders_driver.f90,compare_piecewise_extend_scipy.py,compare_bernstein_orders_scipy.py}`
are new or independently developed SciFort MIT contributions. Power-basis
promotion uses zero-padding, Bernstein-basis promotion uses standard degree
elevation, and Hermite constraints are solved via endpoint derivative identities.
SciPy 1.17 is used for comparison only. For ambiguous upstream SciPy cases, an
independent Hermite linear solve serves as a numerical reference.

`docs/root_finding_upstream_reuse.md` is a license and architecture assessment
of potential upstream packages: it is **not** evidence of third-party source
incorporation. No source from `polyroots-fortran`, `roots-fortran`, or other
catalogued solvers is included in this increment. If such a package is used
later, retain its license and any algorithm-level third-party copyright and
attribution notices and pin the tested FPM revision.

## Piecewise polynomial real root isolation, B-spline conversion and complex support (2026-10-31)

`src/interpolate_piecewise.f90`, `src/interpolate.f90`,
`test/test_interpolate_piecewise_roots_complex.f90`, and the new
`tools/interpolate_reference/ppoly_roots_driver.f90`,
`ppoly_from_spline_driver.f90`, `complex_piecewise_driver.f90`, and matching
Python comparison scripts are original MIT-licensed SciFort contributions.
Real polynomial root isolation uses the mathematical fact that the roots of
an ordinary derivative partition the real axis into monotone intervals;
bisection finds sign-changing zeros and direct evaluation tests repeated
roots at stationary points. B-spline to power-basis conversion uses the
Taylor coefficients `f^(j)(x_i)/j!` of each polynomial segment. Complex
PPoly/BPoly evaluation and calculus use elementary Horner/de Casteljau and
basis-conversion identities already documented in this repository. These
are direct mathematical implementations, not imports or translations of
SciPy or another upstream interpolation library. SciPy 1.17 is an optional
**external numerical reference only**. No third-party source or coefficients
were copied and no new dependency or license notice is required.

## Tensor-product regular-grid B-splines and partial derivatives (2026-10-31)

`src/interpolate_regular_grid.f90`, `test/test_interpolate_tensor_grid.f90`,
`tools/interpolate_reference/tensor_grid_driver.f90`, and
`tools/interpolate_reference/compare_tensor_grid_scipy.py` were extended or
created independently under SciFort's MIT license. The implementation uses
standard not-a-knot knot placement, Cox–de Boor basis recurrences, the
Kronecker/separable collocation identity, per-axis pivoted LU, and analytic
basis differentiation. It calls SciFort's own already-established
`bspline_design_matrix` and requires no external code. The SciPy 1.17
`RegularGridInterpolator` implementation serves only as the independent
numerical reference, with `scipy.sparse.linalg.spsolve` to avoid default
iterative solver tolerance. No source or numerical tables have been copied
from SciPy, `bspline-fortran`, `finterp`, or other upstream packages.

Jacob Williams' BSD-style [`bspline-fortran`](https://github.com/jacobwilliams/bspline-fortran)
remains an excellent optional, pinned FPM backend candidate, especially for
large grid-performance improvements. A normal dependency retains its upstream
license; direct source reuse would require preserving its copyright, license
and disclaimers. No such dependency or source copy was made in this increment.

## Native BSpline and make_interp_spline (2026-10-31)

`src/interpolate_bspline.f90`, `test/test_interpolate_bspline.f90`,
`tools/interpolate_reference/bspline_driver.f90`,
`tools/interpolate_reference/compare_bspline_scipy.py`, and
`docs/interpolate_bspline.md` are independently written, MIT-licensed SciFort
implementations and verification materials. The numerical methods are the
standard Cox–de Boor basis recurrence, de Boor's local spline evaluation,
Boehm knot insertion, coefficient recurrences for derivatives and
antiderivatives, and a pivoted collocation solve for interpolating splines.
The SciPy 1.17 reference package is used only for randomized validation;
**no SciPy, bspline-fortran, or other third-party source was copied or linked**.
The upstream BSD-style [`bspline-fortran`](https://github.com/jacobwilliams/bspline-fortran)
project was reviewed as a possible revision-pinned FPM backend for later
1D–6D tensor-product splines. Its license and notices would need to accompany
any future source copying. Current SciFort B-spline functionality uses only
standard Fortran and `scifort_kinds`.

## CubicSpline, vector-valued ordinate support and interpolation licensing review (2026-10-31)

`src/interpolate.f90` adds original MIT-licensed cubic-spline continuity and
boundary-condition solvers (tridiagonal and periodic cyclic), based on the
standard mathematical spline interpolation identities, not on copied source.
`test/test_interpolate_cubic_spline.f90` and the new reference drivers and
Python comparison scripts in `tools/interpolate_reference/` are original test
implementations. SciPy 1.17.0 was used as an independent reference only.
No third-party Fortran interpolation source was imported, linked, or adapted.
`docs/interpolation_upstream_reuse.md` lists libraries with potentially
compatible licenses, with a caution for mixed/other licenses; these are
future integration candidates, **not current dependencies**. The required
notices will be retained if their code is incorporated in a later revision.

## Piecewise cubic PCHIP, Akima, Makima and Hermite (2026-10-31)

`src/interpolate.f90`, `test/test_interpolate_cubic.f90`,
`tools/interpolate_reference/pchip_akima_driver.f90`,
`tools/interpolate_reference/compare_pchip_akima_scipy.py`, and
`docs/interpolate_piecewise_cubic.md` are original MIT-licensed SciFort
implementations and verification materials. PCHIP slope selection follows the
published Fritsch–Butland (1984, DOI:10.1137/0905021) algorithm, and Akima
slope selection follows Akima (1970, DOI:10.1145/321607.321609). The
modified Akima weights follow Cleve Moler and Cosmin Ionita (2019),
https://blogs.mathworks.com/cleve/2019/04/29/makima-piecewise-cubic-interpolation/.
The cubic Hermite polynomial and its derivatives/antiderivatives are
implemented directly from their elementary interpolation identities. The
zero-weight Akima fallback and default extrapolation policies were checked
against publicly documented SciPy behavior. No SciPy implementation code,
tables or third-party source was imported. `scipy.interpolate` and NumPy are
used by an optional independent Python numerical regression only. The Fortran
production module uses only `scifort_kinds` and language intrinsics.

## Native fixed-taper exact-band spectrogram and DPSS VAR(p) study (2026-10-31)

`src/signal_fixed_taper_spectrogram.f90`, its native/public regression tests, the exact-band Fortran/Python verifier and study drivers in `tools/signal_reference/`, `tools/run_signal_varp_multitaper_coverage.py`, `docs/signal_fixed_taper_spectrogram.md` and `docs/signal_varp_multitaper_coverage.md` are original MIT-licensed SciFort code and documentation. The mathematical implementation follows the finite Fourier transform outer-product definition and the standard analytic integral of a complex exponential over an interval (sinc-band kernel). No third-party source code, lookup tables, numerical coefficients or Fortran implementations were copied. NumPy and SciPy provide only test references: SciPy periodic discrete prolate spheroidal windows and a discrete Lyapunov solve for stationary VAR(p) covariance, plus independent quadrature and SVD tests. The core library depends only on other SciFort modules. Exact-band fixed-weight estimation and the coverage study are SciFort extensions, not SciPy public APIs; simulated coverage is not a guarantee of valid inference.

## Fixed-weight VAR(p) multitaper spectral/canonical targets (2026-10-31)

`src/signal_varp_multitaper_target.f90`, its native and public tests,
`tools/run_signal_varp_multitaper.py`, the Fortran/Python independent reference
pair under `tools/signal_reference/`, and
`docs/signal_varp_multitaper_targets.md` are original MIT-licensed SciFort
contributions. The mathematical method follows linearity of expectation and
integration applied to fixed-weight averages of finite-window spectral outer
products, together with the existing documented discrete Lyapunov VAR(p)
stationary autocovariance identities and Schur-complement/Cholesky/SVD
canonical coherence. No code was imported or translated from a third-party
library, and no third-party numerical coefficients were copied. NumPy and
SciPy (`scipy.linalg.solve_discrete_lyapunov`,
`scipy.signal.windows.dpss` and SciPy/NumPy linear algebra) are used only to
construct independent test references, not by the Fortran implementation.
The routines are explicitly SciFort extensions, not SciPy APIs; fixed-taper
spectral population targets are not inferential confidence intervals.

## Conditional multimode VAR(p) coverage diagnostics (2026-10-30)

`tools/signal_reference/study_canonical_varp_conditional.py`,
`tools/signal_reference/canonical_varp_coverage_driver.f90`,
`tools/signal_reference/test_study_canonical_varp_conditional.py`,
`tools/run_signal_varp_conditional_coverage.py`, and
`docs/signal_varp_conditional_coverage.md` are original MIT-licensed SciFort
simulation and verification contributions. They implement standard Gaussian
VAR(p) state-space simulation and finite-window Fourier-kernel integration
from mathematical definitions, without copying any external source code.
SciPy's discrete Lyapunov solver, NumPy's linear algebra and SciPy/NumPy
canonical-SVD references are **test tools only**. No production source, C ABI,
Fortran numerical algorithm, or third-party runtime dependency was changed.
The empirical 95% coverage statistics are Monte Carlo diagnostics, not a
validated frequentist guarantee.

## Exact finite-window VAR(p) and conditional canonical targets (2026-10-30)

`src/signal_varp_taper_spectrum.f90`, `src/signal_varp_canonical_target.f90`,
new native/public regression tests, `tools/signal_reference/verify_varp_*`,
`tools/run_signal_varp_population.py`, and
`docs/signal_varp_population_targets.md` are original MIT-licensed SciFort
contributions. The implementation uses standard VAR companion-state
construction, the discrete Lyapunov identity, Smith doubling, exact integrals
of finite Fourier exponentials, and the existing SciFort canonical solver.
No third-party source code or numerical coefficients were copied. SciPy
`solve_discrete_lyapunov`, NumPy and SciPy singular value routines are used
only as independent test references; the new Fortran source does not link
against BLAS/LAPACK or require external libraries. These methods are
SciFort extensions, not SciPy signal functions, and are not claimed to
calibrate confidence intervals.

## Finite-window VAR(1) population spectra and bootstrap coverage (2026-10-30)

`src/signal_var1_taper_spectrum.f90`, the new native/public tests, NumPy/SciPy
independent verification scripts and the finite-window VAR coverage study are
original SciFort MIT-licensed contributions. The population spectral-moment
formula follows the standard stationary VAR autocovariance identities, and
the doubling solution of the discrete Lyapunov equation is independently
implemented from its algebraic fixed-point identity. This is commonly
called Smith doubling. No external source code or numerical coefficients
were copied. SciPy `solve_discrete_lyapunov` and
NumPy are used solely for independently calculating test references;
Fortran code adds no third-party dependency. The additional Monte Carlo
coverage study remains an exploratory diagnostic, not evidence that any
bootstrap procedure achieves calibrated intervals.

## Block-aware BCa canonical moving-block bootstrap (2026-10-30)

The BCa acceleration helpers, opt-in `interval_method='bca'`, conditional
mode/group acceleration output arrays, new native/public tests, NumPy/SciPy
reference driver and VAR coverage diagnostic, and
`docs/signal_canonical_bca.md` are original MIT-licensed SciFort changes.
The BCa endpoint formula and jackknife skewness definition follow Bradley
Efron (1987), *Better Bootstrap Confidence Intervals*, JASA 82(397),
171-185, https://doi.org/10.1080/01621459.1987.10478410. Dependence-aware
block resampling background: Künsch (1989), *The Jackknife and the Bootstrap
for General Stationary Observations*, Annals of Statistics 17(3), 1217-1241.
These are mathematical citations; no external source implementation has been
copied. Unequal deletion block weighted pseudovalues are a first-order
heuristic, not established higher-order BCa calibration. NumPy/SciPy are
used solely as numerical references; no new Fortran dependency was added.

## Studentized canonical moving-block bootstrap (2026-10-30)

The `studentized` option, block pseudovalue variance helper, bootstrap-t
endpoint helper and opt-in standard-error fields in
`src/signal_canonical_bootstrap.f90`, plus the native/public regression tests,
independent Python/Fortran NumPy reference and documentation, are original
MIT-licensed SciFort contributions. The studentized interval algebra
follows textbook bootstrap-t formulas (Efron and Tibshirani, 1993,
*An Introduction to the Bootstrap*); block-resampling background follows
Künsch (1989), *Annals of Statistics* 17(3), 1217–1241. These are cited as
mathematical references only: no external implementation was copied. The
weighted pseudovalue extension for incomplete blocks is heuristic. NumPy
and SciPy are used only for validation; no new Fortran dependencies or
license obligations were added. This is not an API of `scipy.signal` and
is not asserted to achieve nominal frequentist coverage.

## Bias-aware canonical bootstrap interval construction (2026-10-30)

The new trailing `interval_method` options and the `normal_cdf`,
`normal_quantile`, and `interval_bounds` helpers in
`src/signal_canonical_bootstrap.f90`, along with their Fortran test drivers,
independent Python references and documentation, are original MIT-licensed
SciFort contributions. The basic/reverse and bias-corrected (BC) percentile
interval identities are based on published mathematical formulas; see
Bradley Efron (1987), *Better Bootstrap Confidence Intervals*, *Journal of the
American Statistical Association* 82(397), 171-185, and Efron and Tibshirani
(1993), *An Introduction to the Bootstrap*. The serial-window resampling
background is H. R. Künsch (1989), *The Jackknife and the Bootstrap for
General Stationary Observations*, *Annals of Statistics* 17(3), 1217-1241.
No algorithm source code or numerical coefficients were copied from an
external implementation. NumPy and SciPy are used for independent verification
and exploratory coverage simulations only; no additional Fortran runtime
libraries or third-party licensing obligations were introduced. The new
methods are SciFort extensions, not SciPy signal APIs, and are not claimed
to provide calibrated confidence levels.

## Automatic canonical coherence bootstrap block-length pilot (2026-10-30)

`src/signal_canonical_block_selection.f90`, the optional block-length interface
in `src/signal_canonical_bootstrap.f90`, new native/public tests, the independent
NumPy pilot and AR(1) study and `docs/signal_canonical_auto_block_length.md`
are original MIT-licensed SciFort contributions. The heuristic combines
window-matrix Frobenius autocovariances, a simple fixed pilot-lag cutoff,
Bartlett tapering and an exploratory variance-oriented cube-root scale.
Künsch (1989), *The Jackknife and the Bootstrap for General Stationary
Observations*, *Annals of Statistics* 17(3), 1217–1241 is cited only for
background on dependent-data block bootstrap. No code was copied, and
this method is not a scipy.signal API or a claim of optimal block selection.
NumPy and SciPy are used in verification scripts only; no new Fortran runtime
dependencies were introduced.

## Moving-block bootstrap for canonical spectral coherence (2026-10-30)

`src/signal_canonical_bootstrap.f90`, its public/native regression tests,
independent NumPy/SciPy reference driver and validation scripts, plus
`docs/signal_canonical_moving_bootstrap.md`, are original MIT-licensed
SciFort contributions. The block resampling scheme follows standard
mathematical moving-block-bootstrap descriptions such as H. R. Künsch
(1989), *The Jackknife and the Bootstrap for General Stationary
Observations*, *Annals of Statistics* 17(3), 1217-1241, and Efron and
Tibshirani (1993), *An Introduction to the Bootstrap*. These are
mathematical background sources only; no implementation code was copied.
The independent numerical reference uses NumPy complex linear algebra
and SciPy interpolation/trapezoidal integration solely during validation.
The Fortran implementation reuses SciFort's existing canonical coherence,
multitaper and frequency-band integration APIs; no external Fortran
runtime dependencies, imported library code or license obligations were
introduced. This is a SciFort extension, not a scipy.signal function.

## Contiguous-block jackknife for canonical spectral coherence (2026-10-30)

The equal-block deletion, metadata, optional Fortran wrapper interfaces,
independent NumPy/SciPy reference driver, tests and documentation are original
MIT-licensed SciFort contributions in `src/signal_canonical_jackknife.f90`
and associated files. The implementation uses the standard equal-group
jackknife variance and bias identities, with contiguous groups of complete
spectral windows, and reuses the existing SciFort canonical spectral estimator.
Background: H. R. Künsch (1989), *The Jackknife and the Bootstrap for General
Stationary Observations*, Annals of Statistics 17(3), 1217–1241; Efron and
Tibshirani (1993), *An Introduction to the Bootstrap*. These are mathematical
references, not copied implementations. NumPy Cholesky/SVD and SciPy
interpolation/integration are used only in independent test calculations.
This is a SciFort extension, not a SciPy `signal` function; there are no new
third-party licensing obligations or runtime dependencies.

## Nonoverlapping-window canonical spectral jackknife (2026-10-30)

`src/signal_canonical_jackknife.f90`, the public and native tests,
Python/Fortran independent reference driver and documentation are original
MIT-licensed SciFort contributions. The method uses standard delete-one
jackknife and Fisher-transformed correlation identities, as described in
Efron, B. and Tibshirani, R. J. (1993), *An Introduction to the Bootstrap*,
Chapman & Hall/CRC, and reuses SciFort
canonical coherence, multitaper and band spectral integration code.
NumPy complex Cholesky/SVD and SciPy piecewise-linear interpolation and
trapezoid integration are independently written test references only.
No NumPy/SciPy implementation code or new external Fortran dependency was
incorporated. The function is an original SciFort extension, not a
`scipy.signal` API, and the intervals are exploratory rather than
calibrated.

## Covariance-metric spectral subspace matching (2026-10-30)

The conditioned covariance output in `src/signal_canonical_coherence.f90`,
midpoint-metric Cholesky whitening and comparison in
`src/signal_canonical_subspace_tracking.f90`, test drivers and mathematical
documentation are original MIT-licensed SciFort contributions based on
standard Hermitian positive-definite inner-product, Cholesky, Schur-complement,
principal-angle and unitary Procrustes identities. Existing SciFort principal-
angle, Hungarian assignment and small-matrix eigenvalue code is reused.
NumPy/SciPy routines serve as independent validation oracles only; no source
was imported from NumPy, SciPy or a third-party BLAS/LAPACK implementation.
This is a SciFort-specific extension, not a `scipy.signal` API. No additional
Fortran runtime dependency or third-party licensing obligations were added.

## Canonical spectral subspace tracking (2026-10-30)

`src/signal_canonical_subspace_tracking.f90`, its native/public test programs,
Python/Fortran reference driver, test scripts, and the corresponding
mathematical documentation are original MIT-licensed SciFort contributions.
The small Hermitian Jacobi eigensolver is adapted from SciFort's existing
`src/signal_canonical_coherence.f90` in the same source tree, with the
squared-singular-value clipping modified for general positive-semidefinite
Gram matrices. The principal-angle, complex polar/Procrustes, and global
matching procedures implement standard linear-algebra identities independently.
NumPy SVD and SciPy's `linear_sum_assignment` are test oracles only; no source
code from SciPy, NumPy, or external BLAS/LAPACK was imported. This is a
SciFort extension, not a public `scipy.signal` function. There are no new
third-party Fortran dependencies.

## Frequency-continuous canonical mode tracking (2026-10-30)

`src/signal_canonical_tracking.f90`, its native and public regressions, the
reference driver, documentation, and test harness are original MIT-licensed
SciFort contributions based on standard complex inner-product phase matching
and the classical Hungarian shortest-augmenting-path assignment principle.
The implementation was independently written for this project; no SciPy,
NumPy, or third-party source code was copied or incorporated. NumPy and SciPy's
`scipy.optimize.linear_sum_assignment` are independent verification oracles
only. This is a SciFort extension rather than a public `scipy.signal` routine.
No new Fortran runtime dependency or third-party code was added.

## Canonical variate coefficient vectors (2026-10-29)

The optional complex canonical coefficient vectors, Hermitian triangular
solve, accumulated Jacobi eigenvectors, zero-root orthogonal complement,
phase convention, direction-stability heuristic, native and Python/NumPy
reference tests, and documentation are original MIT-licensed SciFort
contributions based on standard Cholesky/SVD canonical-correlation identities.
NumPy SVD/Cholesky is used as an independent test oracle only: no NumPy,
SciPy, LAPACK or other third-party implementation was copied or imported.
This is a SciFort extension, not a top-level `scipy.signal` API; the Fortran
runtime has no new dependencies.

## Canonical spectral coherence modes (2026-10-29)

The conditional Schur-complement whitening, complex lower-triangular solves,
cyclic Hermitian Jacobi eigenvalue routine, sorted squared canonical coherence
modes, public interfaces, tests and documentation in
`src/signal_canonical_coherence.f90` are original MIT-licensed SciFort
contributions, based on standard canonical-correlation and Hermitian Jacobi
linear-algebra identities. Existing SciFort group-coherence regularization
and determinant summaries are reused without source adaptation. Independent
NumPy Cholesky/SVD calculations are used as verification references only;
no NumPy, SciPy, LAPACK or third-party implementation code was copied. The
feature is a **SciFort extension** (not a public `scipy.signal` function) and
requires no new Fortran runtime dependencies.

## Conditional spectral coherence between channel groups (2026-10-29)

The determinant group-coherence estimator, its Hermitian Cholesky logdet
subroutine, conditional Wilks-ratio outputs, spectral/band/time-frequency
interfaces, regression programs and NumPy reference driver are original
MIT-licensed SciFort contributions. The methods use standard algebraic
identities for Schur complements, determinants and canonical correlations.
NumPy singular-value decomposition and linear algebra are used solely for
independent verification; no source code from NumPy or SciPy was copied.
The addition is a **SciFort extension**, not a public `scipy.signal` function;
no new Fortran dependencies or third-party sources were introduced.

## Eigenvalue-weighted automatic spectral shrinkage (2026-10-29)

The fixed-weight equivalent taper count, metadata, public utility,
shrinkage integration and associated tests are original MIT-licensed
SciFort contributions based on standard weighted complex Gaussian
covariance moment identities. No external numerical source code was
copied, and there are no new runtime dependencies. The algorithm is
an explicitly labeled SciFort extension, not a `scipy.signal` function.

## Automatic Gaussian-Wishart spectral shrinkage (2026-10-29)

The automatic diagonal-target shrinkage and its reported coefficient arrays,
new native tests, public interface checks, NumPy reference comparisons and
documentation are original MIT-licensed SciFort contributions. The estimator
uses standard complex Wishart moment identities for diagonal-target covariance
regularization. It is not an exact oracle, classical Ledoit-Wolf estimator,
or public SciPy API. NumPy spectral-matrix inversions and SciPy 1.17 older
regression suites are verification references only; neither library's source
was copied into SciFort. No new Fortran runtime dependencies were introduced.

## Complex partial coherency and diagonal shrinkage (2026-10-29)

The optional complex coherency, individual-channel diagonal-target spectral
shrinkage, Fortran tests, and comparison scripts are original MIT-licensed
SciFort contributions implementing standard inverse spectral-matrix formulas.
NumPy linear algebra and SciPy 1.17 DPSS/FFT references are used solely for
independent verification and were not incorporated into the library.
No new Fortran dependencies or third-party source files were added.

## Time-frequency conditional multitaper coherence (2026-10-29)

`src/signal_partial_coherence.f90` now applies its existing original
MIT-licensed Hermitian Cholesky precision-matrix algorithm independently to
sliding cross-spectral matrices from `scifort_signal_multitaper`.
`partial_coherence_tf_matrix`, `multitaper_tf_partial_coherence`,
`partial_coherence_tf_result`, the public re-exports, test programs and
reference comparison harness are original SciFort contributions. These are
SciFort extensions and not public `scipy.signal` functions.
SciPy 1.17 DPSS windows and NumPy FFT/inversion algorithms are used solely
for independently computed test references; no external source code has been
incorporated. No new Fortran dependencies were added.

## Multichannel partial spectral coherence (2026-10-29)

The complex Hermitian Cholesky solver, inverse-spectral conditional coherence,
relative diagonal loading, public interfaces, tests, and documentation in
`src/signal_partial_coherence.f90` and associated files are original
MIT-licensed SciFort contributions implementing standard spectral precision
matrix formulas. NumPy `linalg.inv` and SciPy 1.17 `signal.windows.dpss` are
used solely in an independent reference driver; no source from either project
was copied and no new Fortran dependency was added. `multitaper_partial_coherence`
is a SciFort extension, not a standard `scipy.signal` API.

## Irregular-grid Savitzky-Golay local regression (2026-10-29)

The `savgol_irregular` procedures in `src/signal_savgol.f90`, public exports,
native tests and independent numerical references are original MIT-licensed
SciFort contributions. They implement standard weighted local-polynomial
least squares via centered, scaled, twice-reorthogonalized QR. NumPy
`lstsq` and SciPy 1.17 `signal.savgol_filter` are comparison references
only, not copied source code or runtime dependencies. This feature is a
SciFort extension, not a public SciPy irregular-grid API.

## Prepared weighted multiharmonic least-squares plans (2026-10-29)

The cached harmonic QR factors, channel-wise fit and uncertainty interfaces,
new native and public tests, and independent reference harness are original
MIT-licensed SciFort contributions. The algorithm is a direct use of
weighted modified Gram-Schmidt with twice-reorthogonalization and standard
least-squares triangular solves; NumPy weighted QR and SciPy 1.17
`signal.lombscargle` are used only for numerical verification. The new public
interfaces are SciFort extensions, not public SciPy APIs. No external
Fortran code was copied and no new Fortran dependencies were introduced.
See `docs/signal_harmonic_plan.md`.

## Harmonic coefficient and curve uncertainty (2026-10-29)

The weighted covariance calculation and `lombscargle_harmonic_predict_std`
interface in `src/signal_irregular.f90`, public export, tests, and
reference runner are original MIT-licensed SciFort contributions. They use
standard weighted regression formulas and triangular solves based on the
existing independently implemented QR factor. NumPy's independent full-design
weighted least-squares and matrix inverse are used solely in comparison
scripts; no SciPy or NumPy source code was copied into the library. No new
Fortran dependencies or license changes. See
`docs/signal_harmonic_uncertainty.md`.

## Multiharmonic coefficient extraction and prediction (2026-10-28)

The `lombscargle_harmonic_fit` and `lombscargle_harmonic_predict` interfaces,
weighted QR coefficient recovery, native regression tests, and independent
reference harness are original MIT-licensed SciFort code. The algorithms use
standard Fourier-series weighted least squares and twice-reorthogonalized
modified Gram–Schmidt with triangular back substitution. NumPy `lstsq` and
SciPy 1.17 `signal.lombscargle` are used only for numerical comparisons; no
upstream source was copied, and no Fortran dependency was introduced. These
interfaces are SciFort extensions, not functions in SciPy's public signal API.
See `docs/signal_harmonic_fit.md`.

## Weighted multiharmonic Lomb–Scargle (2026-10-28)

The weighted multi-sinusoid projection and twice-reorthogonalized modified
Gram–Schmidt algorithm in `src/signal_irregular.f90` are original
MIT-licensed SciFort code implementing standard weighted Fourier-series
least squares. The rank-one `nterms=` extension, regression tests, and
comparison driver were written for SciFort. Public SciPy 1.17
`signal.lombscargle` is used solely to verify one-harmonic outputs;
NumPy's generic weighted least-squares solver independently verifies
multiple harmonics. No third-party Fortran or SciPy source was copied,
no vendored routines were added, and Fortran runtime dependencies are
unchanged. See `docs/signal_multiharmonic.md`.

## Fast skip-ahead and rewind for Fibonacci maximum-length sequences (2026-10-28)

The `skip=` argument of `max_len_seq` and the packed GF(2) transition
exponentiation in `src/signal_sequences.f90` are original MIT-licensed
SciFort implementations of standard binary linear algebra and Fibonacci
LFSR feedback. This is a SciFort extension, not a public SciPy
`max_len_seq` parameter. The reference tests use public SciPy 1.17
`signal.max_len_seq` only to obtain independent one-clock transition
columns; Python bitwise Gaussian elimination and exponentiation verify
long positive and negative jumps. No third-party source code was imported,
no license notices were removed, and no new Fortran dependencies were added.
Details: `docs/signal_sequence_skip.md`.

## Parallel digital SOS filter banks (2026-10-28)

The new `dlsim_sos_bank` overloads and `discrete_sos_bank_result` type in
`src/signal_discrete_simulation.f90`, with the public exports, tests and
comparison scripts, are original MIT-licensed SciFort work. The algorithm
applies the existing independently written transposed direct-form-II
per-section recurrence separately to each input channel, using a shared
or per-channel SOS cascade. SciPy 1.17 `signal.sosfilt` is used only as
an independent numerical reference for filter outputs and state histories.
No upstream SciPy source was imported and no runtime dependencies were added.
See `docs/signal_dlsim_sos_bank.md` for conventions and numerical checks.

## SOS-native arbitrary-input discrete simulation (2026-10-28)

The `dlsim_sos` function and `dlsim` generic overloads in
`src/signal_discrete_simulation.f90`, top-level `scifort_signal` exports,
tests, and reference programs are original MIT-licensed SciFort code.
The routine uses the standard transposed direct-form-II second-order-section
recurrence already implemented in SciFort's `sosfilt`, with no high-order
polynomial expansion or imported third-party routines. See
`docs/signal_dlsim_sos.md` for equations, state conventions, validation and
references. SciPy 1.17's `signal.sosfilt` and NumPy serve only as independent
numerical verification; no SciPy source was copied, and no runtime license
or dependency was added.

## Direct digital SOS step/impulse responses (2026-10-28)

The `dstep_sos` and `dimpulse_sos` overloads in
`src/signal_discrete_simulation.f90` are original MIT-licensed SciFort
wrappers using the existing independently implemented SOS time-domain
recurrence. No SciPy source code or external Fortran routines were copied.
SciPy 1.17 `signal.sosfilt`, and for low orders `signal.dstep`/
`signal.dimpulse`, were used as independent numerical references only.
No new third-party dependencies or license changes. See
`docs/signal_sos_responses.md` and `tools/run_signal_sos_responses.py`.

## Stable digital ZPK group delay (2026-10-28)

The `group_delay_zpk` overloads and the root-wise delay implementation in
`src/signal_design.f90`, public export, tests, and scripts are original
MIT-licensed SciFort contributions. The mathematical method is the
logarithmic derivative of the standard digital ZPK transfer function.
No SciPy source or third-party Fortran implementation was copied.
SciPy 1.17 and NumPy are used only to verify the numerical results via
public `signal.group_delay` and `signal.zpk2tf` calls in a Python test.
`docs/signal_group_delay_zpk.md` records the ZPK convention, sample-delay
units, singularities and reference tests. There are no new runtime
libraries or license changes.

## Stable group delay from second-order sections (2026-10-28)

The new `group_delay_sos` functions in `src/signal_design.f90` and the public
export in `src/signal.f90` are original MIT-licensed SciFort contributions.
They derive group delay by differentiating the phase of the standard
second-order-section transfer function, using the additive property of
logarithmic derivatives. The derivation is documented in
`docs/signal_group_delay_sos.md`; see Julius O. Smith III,
*Introduction to Digital Filters with Audio Applications*, "Numerical
Computation of Group Delay", and Richard G. Lyons, *Understanding Digital
Signal Processing*, 3rd edition, p. 830 (references appearing in the SciPy
1.17 `signal.group_delay` documentation).
No upstream SciPy code or third-party Fortran code was copied. The native
regressions and reference input/output drivers are original SciFort test
code; NumPy and SciPy 1.17 are used only for independent comparisons and
are not runtime dependencies. This is a SciFort extension, not an existing
SciPy SOS API.

## Complex-wavelet magnitude ridge detection (2026-10-28)

The new `find_peaks_cwt_complex` and `find_peaks_cwt_magnitude` interfaces
in `src/signal_cwt_peaks.f90`, public exports, native tests, and standalone
reference runner are original MIT-licensed SciFort integration work.
They reuse SciFort's preexisting real CWT ridge tracker and independently
implemented complex Morlet/convolution kernels. The original tracker already
records its SciPy-compatible ridge/noise conventions and BSD-3-Clause
attribution; those notices are retained. The new explicit modulus criterion
is a SciFort extension, not a direct implementation of a SciPy public complex
peak finder. NumPy and SciPy 1.17 are used exclusively in reference tests:
SciPy's public `signal.convolve` constructs a CWT, and its private
`_identify_ridge_lines` and `_filter_ridge_lines` provide an independent
real-valued magnitude-map oracle. These test-only APIs may change.
No new runtime dependencies or third-party Fortran source were added.

## Complex CWT and Morlet analysis (2026-10-28)

The complex wavelet callback, normalized Morlet generator and conjugated
wavelet-convolution implementation in `src/signal_cwt_peaks.f90`, and the
public exports in `src/signal.f90`, are original MIT-licensed SciFort
contributions using standard wavelet formulas. The independent SciPy 1.17
`signal.convolve` numerical oracle, native Fortran tests and standalone
validation runner are original test code. SciPy's historical `_cwt`
conjugate-reversal and same-mode geometry are used as compatibility
conventions; no SciPy source code is copied into the Fortran implementation.
NumPy and SciPy remain test-only dependencies. The previously recorded SciPy
BSD-3-Clause attribution for real CWT peak finding remains unchanged.

## Mixed-type rank-one convolution (2026-10-28)

The mixed real/complex overloads in `src/signal_convolve.f90` and
`src/signal_fftconvolve.f90` are original MIT-licensed SciFort adapter code.
They reuse preexisting independent Fortran complex convolution, correlation,
FFT and overlap-add implementations; no SciPy source code or third-party
Fortran code was copied. The native tests, reference driver and validation
runner are original project contributions. SciPy 1.17 and NumPy are used
only for independent numerical verification. No runtime dependencies or
licenses have changed.

## ShortTimeFFT user-procedure detrending (2026-10-28)

New procedure overloads and per-slice dispatch in
`src/signal_short_time_fft.f90`, the public callback interface reexports,
native tests, reference output driver, and standalone numerical comparison
scripts are original MIT-licensed SciFort contributions. They use the
independent Fortran STFT implementation already present in the project. No
third-party Fortran code or SciPy source was copied. SciPy 1.17 is used
solely to test numerical results and coordinate conventions, including
callable transforms along time segments. FPM is not required for the
standalone regression runner.

## Mixed cross-spectral estimators (2026-10-28)

The rank-one mixed real/complex `csd`, `coherence`, and ShortTimeFFT
cross-spectrogram overloads in `src/signal_spectral_complex.f90` and
`src/signal_short_time_fft.f90`, and exports in `src/signal.f90`, are original
MIT-licensed SciFort glue for the existing independent Fortran Welch and
ShortTimeFFT kernels. No external Fortran code or SciPy source was copied.
The reference driver, native tests and comparison tool are original. SciPy
1.17 and NumPy are test-only dependencies. For the legacy CSD/coherence tests,
the real input is explicitly cast to `complex128` to bypass an order-dependent
mixed-dtype issue in SciPy's own dispatcher. ShortTimeFFT tests use direct
mixed-type SciPy calls. See `docs/signal_mixed_cross.md`.

## Complex local Wiener filtering (2026-10-28)

The rank-1/2/3 complex `wiener` extensions and shared flat-array kernel in
`src/signal_denoise.f90`, together with the test, numerical reference driver
and validation scripts, are original MIT-licensed SciFort contributions.
They reuse the existing independent SciFort zero-padded neighborhood walker,
without introducing any third-party Fortran code or runtime dependencies.
The complex pseudo-variance semantics were checked using the public SciPy
1.17 `signal.wiener` API (https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.wiener.html).
SciPy and NumPy are test-only dependencies; no SciPy source was incorporated.

## Callable spectral detrending (2026-10-28)

The in-place callback interfaces, segment forwarding, error validation and
multidimensional dispatch in `src/signal_spectral.f90`,
`src/signal_spectral_complex.f90`, `src/signal_spectral_nd.f90`,
`src/signal_timefreq_nd.f90`, and the public exports in `src/signal.f90`
are original MIT-licensed SciFort work. They compose the preexisting SciFort
FFT, window-generation and rank-independent lane-walking routines. No SciPy
implementation source or third-party Fortran library was incorporated.
Behavior was independently tested against SciPy 1.17's documented public
callable detrend option and Fourier/spectral routines, especially
https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.welch.html .
The new native Fortran tests, SciPy comparison driver, and reproducible runner
are original project contributions. NumPy and SciPy are test-only
requirements. No new runtime dependencies or licenses were introduced.

## Parameterized Fourier resampling (2026-10-28)

The new trailing Fortran window-shape arguments, propagation through
`src/signal_fourier_resample.f90` and `src/signal_multirate_nd_into.f90`,
and accompanying tests and documentation are original MIT-licensed SciFort
work. The code composes existing independently implemented SciFort FFT and
window-generation routines. No SciPy implementation source or third-party
Fortran library was copied. Numerical behavior was checked independently
against the public SciPy 1.17 `signal.resample` API and public
`signal.windows.general_cosine` for its unsupported tuple dispatch case:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.resample.html .
SciPy and NumPy are test-only requirements, not runtime dependencies.

## Complex second-order-section decimation (2026-10-28)

The complex SOS overloads in `src/signal_resample.f90` and
`src/signal_multirate_nd_into.f90`, their Fortran tests, Python reference
harness and documentation are original SciFort contributions under MIT.
They compose existing SciFort complex `sosfilt` and `sosfiltfilt` routines;
no SciPy implementation source or third-party Fortran code was copied.
The independent numerical reference uses SciPy 1.17's public
[`sosfilt`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.sosfilt.html)
and [`sosfiltfilt`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.sosfiltfilt.html)
interfaces, with decimation by indexing. SciPy and NumPy are test-only
requirements, not SciFort runtime dependencies.

## Transfer-function `decimate` filters (2026-10-28)

The implementation of optional `b`/`a` transfer-function coefficient arrays
in `src/signal_resample.f90` and `src/signal_multirate_nd_into.f90`, the
regeneration safeguard in `tools/generate_signal_multirate_nd.py`, and
associated reference and regression tests are original SciFort code under
the existing MIT license. The implementation composes preexisting independent
SciFort `lfilter`, `filtfilt`, `resample_poly` and `upfirdn` algorithms; no
SciPy source code or external Fortran dependency was incorporated. Parameter
semantics and numerical values were compared independently against SciPy
1.17's public `signal.dlti` and `signal.decimate` APIs. NumPy and SciPy are
used only by the test scripts.

## User-supplied FIR/SOS filters in decimate (2026-10-28)

The implementation and native regression tests for the trailing `taps` and
`sos` options in `src/signal_resample.f90` and
`src/signal_multirate_nd_into.f90` are original MIT-licensed SciFort code.
They reuse the existing independently implemented FIR polyphase and IIR SOS
kernels. The numerical output is compared against the documented public
`scipy.signal.decimate` API using SciPy 1.17 `dlti` objects as independent
references: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.decimate.html .
SciPy/NumPy are used only in verification scripts; no external source code,
library dependency, Fortran package or license was copied or added.

## Parameterized `resample_poly` windows (2026-10-28)

The changes to `src/signal_resample.f90` and
`src/signal_multirate_nd_into.f90` and their accompanying tests are original
MIT-licensed SciFort code. They reuse SciFort's existing symmetric named
window generator, FIR prototype, polyphase filtering and multidimensional
wrappers. The parameter conventions and numerical outputs were independently
compared to SciPy 1.17's public `scipy.signal.resample_poly` interface:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.resample_poly.html .
For the general-cosine case, the independent test constructs the symmetric
window using SciPy's public `windows.general_cosine` and passes a custom FIR
taper, due to a SciPy 1.17 argument-dispatch limitation. SciPy/NumPy are
comparison-only dependencies; no SciPy implementation source or third-party
Fortran library was copied or imported. No new runtime dependency or license
was introduced.
## Legacy spectrogram output modes (2026-10-28)

The new output-mode handling in `src/signal_spectral.f90`,
`src/signal_spectral_complex.f90`, and `src/signal_timefreq_nd.f90` is
original MIT-licensed SciFort code. It reuses the existing independently
implemented SciFort spectral window, detrending and FFT kernels; it does not
copy SciPy implementation source. The option names, amplitude normalization,
frequency-axis phase unwrapping and signed-zero branch conventions were
independently compared against the documented public SciPy 1.17
`scipy.signal.spectrogram` interface:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.spectrogram.html .
The native and Python-to-Fortran reference tests are new code for SciFort;
NumPy/SciPy are test-only dependencies. No licenses or runtime dependencies
were added.

## Scalar/logical and odd-edge `firwin` design (2026-10-27)

The generic overloads and generalized alternating multiband logic in
`src/signal_design.f90` are original MIT-licensed SciFort code, reusing the
existing FIR sinc/window/normalization implementation. They implement the
publicly documented cutoff and passband conventions from
https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.firwin.html .
Native Fortran regression tests and a seeded independent SciPy coefficient
comparison were authored for SciFort; they use SciPy only as a reference,
not as a runtime dependency. No SciPy source or third-party code was copied,
and no new licenses or external Fortran dependencies were introduced.

## FIR transition-width Kaiser and symmetric shaped `firwin2` windows (2026-10-27)

The modifications to `src/signal_design.f90`, `src/signal_fir_advanced.f90`,
and the optional symmetric-window argument of `spectral_named_window` in
`src/signal_spectral.f90` are original SciFort MIT-licensed integration code.
They reuse the existing original FIR design, Kaiser formulas, general
cosine/window families, and FFT kernels. Parameter naming, the normalized
transition-width convention, and symmetric-window conventions were checked
against public `scipy.signal.firwin`, `firwin2`, `kaiser_beta`, `kaiser_atten`
and `get_window` behavior (SciPy 1.17):
https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.firwin.html
and https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.firwin2.html .
The new native and independent Fortran/Python verification programs were
written for SciFort and do not copy SciPy's implementation source.
`window_values` in `firwin2` is a SciFort-only compatibility extension.
SciPy and NumPy are reference-only test dependencies; no third-party code,
new licenses, runtime libraries or exported C-ABI symbols were imported.

## Two-sided legacy real-input spectral estimators (2026-10-27)

The optional `return_onesided` extension of `src/signal_spectral.f90`,
`src/signal_spectral_complex.f90`, `src/signal_spectral_nd.f90`,
`src/signal_timefreq_nd.f90`, and `scripts/gen_signal_spectral_nd.py` is
original MIT-licensed SciFort integration code reusing the existing FFT and
spectral kernels. The frequency ordering and normalization were independently
validated against the public SciPy periodogram, Welch, CSD, STFT and
spectrogram interfaces (https://docs.scipy.org/doc/scipy/reference/signal.html)
and NumPy's FFT frequency convention. `test/test_signal_twosided.f90` and
`tools/{reference_signal_twosided.f90,compare_signal_twosided.py,
run_signal_twosided.py}` are new independently authored tests. SciPy and NumPy
are comparison-only dependencies, not required by the Fortran library. No
third-party implementation source was copied or linked.

## Parameterized legacy inverse STFT windows (2026-10-27)

The changes in `src/signal_istft.f90`, `src/signal_istft_nd.f90`, and the
common spectral-window validator in `src/signal_spectral.f90` are original
SciFort MIT-licensed integration code. They reuse the project's existing
FFT, overlap-add, DPSS and window-generation implementations. Parameter and
output conventions follow SciPy's documented public `stft`, `istft`,
`get_window`, `check_COLA` and `check_NOLA` interfaces:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.istft.html .
The new native regression and independent SciPy/NumPy-to-Fortran comparison
programs are original test code. SciPy/NumPy are validation-only dependencies;
no SciPy source was copied and no runtime dependencies or license obligations
were added.

## Parameterized and caller-provided legacy spectral windows (2026-10-27)

The additions to `src/signal_spectral.f90`,
`src/signal_spectral_complex.f90`, `src/signal_spectral_nd.f90`,
`src/signal_timefreq_nd.f90`, and `scripts/gen_signal_spectral_nd.py`
are original SciFort Fortran and Python integration code, MIT-licensed
with the surrounding project. Window conventions and public expected
outputs were independently checked against the SciPy signal API reference:
https://docs.scipy.org/doc/scipy/reference/signal.html .
The newly authored native test, Fortran reference-output driver,
Python comparison script, and standalone runner are independent test code.
SciPy and NumPy are required only for the optional numerical comparisons,
not for compilation or execution of the Fortran library. No SciPy source
was copied or translated, and no new third-party licenses were introduced.

## Legacy STFT boundary extensions (2026-10-27)

The boundary-extension helpers added to `src/signal_spectral.f90` are original,
MIT-licensed Fortran code implementing the standard even reflection, odd
reflection about the endpoint, and constant-endpoint signal extension
formulas. The mathematical definitions and naming conventions are documented
in the SciPy STFT reference:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.stft.html .
The existing FFT and spectral-estimation kernels are reused; no third-party
implementation code has been copied or translated. `src/signal_spectral_complex.f90`
reuses the same helper. `test/test_signal_stft_boundaries.f90` and the new
Fortran/Python comparison programs are original test harnesses; SciPy is
used only as an independent numerical reference and is not a runtime
library dependency. No new third-party license obligations were introduced.

## Arbitrary-rank Hilbert analytic signals (2026-10-27)

`src/signal_hilbert_nd.f90` is an original MIT-licensed Fortran 2018 adapter
that applies the established SciFort FFT and one-dimensional Hilbert masks
along independently selected, explicitly strided tensor axes. It uses local
reusable SciFort FFT plans and staged allocation rather than external code.
The mathematical definitions correspond to the public SciPy Hilbert
transform documentation:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.hilbert.html
and https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.hilbert2.html .
`test/test_signal_hilbert_nd*.f90` and the independent SciPy comparison driver
and Python harness are new original validation code; SciPy/NumPy are used for
independent numerical comparison only and are not linked or included as
implementation dependencies. No SciPy source or third-party source was
copied, and no license obligations were added.

## Rank-independent relative-extrema masks (2026-10-27)

The new assumed-rank adapters in `src/signal_relative_extrema.f90` are original
MIT-licensed Fortran implementations of neighbor comparisons along a selected
array axis. Their semantics follow the publicly documented `order`, `axis`,
and `mode` definitions of `scipy.signal.argrelextrema`:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.argrelextrema.html .
The new logical-mask return convention is a SciFort extension and is not the
SciPy coordinate-tuple API. The implementation was written independently and
no SciPy or third-party source was copied. Independent SciPy 1.17 reference
checks appear in `tools/signal_reference/verify_extrema_nd_scipy.py` and its
Fortran test driver; native regressions are in `test/test_signal_extrema_nd.f90`.
No external libraries or license obligations were added.

## Arbitrary-axis CZT and zoom FFT adapters (2026-10-27)

`src/signal_czt_nd.f90` is an original MIT-licensed Fortran adapter
built around the previously implemented `scifort_signal_czt` plan kernels.
It uses ordinary Fortran 2018 assumed-rank selection, axis/lane packing,
reusable plans, and output staging; it does not import or translate SciPy
implementation source. The mathematical definition and SciPy argument
conventions are described in the existing chirp-Z provenance entry and
at https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.czt.html
and https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.zoom_fft.html .
`test/test_signal_czt_nd.f90`, `reference/signal_czt_nd_reference_driver.f90`,
and `test/compare_signal_czt_nd_scipy.py` are original regression/reference
files. SciPy 1.17 and NumPy are used only as external numerical reference
calculations. No new third-party code or runtime dependency is introduced.

## Anisotropic and elliptical two-dimensional FIR design (2026-10-27)

`src/signal_fir2d.f90` extends SciFort's independently written separable
FIR implementation using an outer product of two conventional windowed-sinc
filters with independent frequency edges and sampling frequencies. The
elliptical option applies an elementary rotation and axis scaling to the
normalized coordinate grid before one-dimensional linear interpolation of
the oversampled radial filter. These options are SciFort extensions, not
SciPy 1.17 `firwin_2d` API options, and are based on standard geometric and
FIR windowing formulas rather than copied external implementation code.

`test/test_signal_fir2d_anisotropic*.f90` and
`tools/signal_reference/fir2d_anisotropic_driver.f90` are original tests.
`tools/signal_reference/compare_signal_fir2d_anisotropic_scipy.py` uses
SciPy's `firwin` and NumPy's `outer`, `meshgrid`, `hypot`, and `interp` only
as independent reference calculations. No third-party source has been
imported or modified, and no new dependencies were added to the library.

## DPSS half-sample normalization (2026-10-27)

`src/signal_dpss.f90` adds an original Fortran evaluation of the midpoint
value of the periodic bandlimited interpolant of the first even-length
Slepian taper. The finite Fourier-series correction is expressed using
`1 + 2 sum_{k=1}^{n/2} cos(k theta) = sin((n+1)theta/2)/sin(theta/2)`,
which eliminates the FFT from the peak-normalization calculation. This
implements the numerical convention documented for `norm='subsample'` at
https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.dpss.html .
No SciPy implementation source or third-party Fortran code was copied.
SciPy 1.17 is used exclusively as an independent comparison reference in
`test/compare_signal_dpss_subsample_scipy.py`; a native regression in
`test/test_signal_dpss_subsample.f90` checks parity, scaling, defaults,
errors, and the public module. No new license or build dependency is added.

## FFT-plan reuse in convolution and multitaper estimation (2026-10-26)

The changes in `src/signal_fftconvolve.f90` and `src/signal_multitaper.f90`
reuse SciFort's own `fft_plan` within its existing separable FFT and DPSS
spectral-estimation kernels. They introduce no external FFT implementation,
code copying, or Fortran dependency. The new native test and benchmark are
original SciFort programs; the original overlap-add SciPy comparison script
was repaired to compile modules in dependency order. SciPy is used only as
an independent numerical reference. The plan is mutable and owned per
operation or streaming instance; simultaneous threads require separate plans.

## Mixed-radix Cooley-Tukey FFTs (2026-10-26)

The new mixed-radix path in `src/signal_spectral.f90` is an original
MIT-licensed Fortran implementation of the standard decimation-in-time
Cooley-Tukey decomposition, selecting factors 2, 3, 5 and 7 and reusing
an original-length table of complex roots of unity. Its mathematical basis
is J. W. Cooley and J. W. Tukey, "An algorithm for the machine calculation
of complex Fourier series," *Mathematics of Computation* **19** (1965),
297-301. No third-party implementation code was imported or adapted.
`test/test_signal_fft_mixed_radix.f90`,
`tools/signal_reference/compare_signal_fft_mixed_radix.py`, and
`benchmarks/bench_signal_fft_mixed_radix.f90` are original SciFort test
and benchmark code. SciPy is used only as an independent numerical
reference. No new library dependencies or third-party licenses are added.

## Reusable Bluestein FFT plans (2026-10-26)

The new original MIT-licensed code in `src/signal_spectral.f90` factors the
previous SciFort Bluestein convolution into an explicit reusable `fft_plan`,
using the mathematical chirp identity described in the Bluestein FFT section
below. No external FFT library or source has been imported. Calls from
`src/signal_stft_stream.f90`, `src/signal_istft_stream.f90`, and internal batch
spectral estimators reuse their own plan. The new regression, SciPy comparison
harness, and benchmarks are original SciFort test/measurement code; SciPy is
only an independent numerical reference.

## Bluestein FFT acceleration (2026-10-26)

The new non-power-of-two FFT path in `src/signal_spectral.f90` is original
MIT-licensed Fortran code implementing the chirp-convolution identity of
L. Bluestein, "A linear filtering approach to the computation of discrete
Fourier transform," *IEEE Transactions on Audio and Electroacoustics* 18(4),
451-455 (1970), DOI: https://doi.org/10.1109/TAU.1970.1162132 .
No paper source text or third-party implementation code was copied.
The original SciFort radix-2 FFT is reused for the padded convolution.
`test/test_signal_fft_bluestein.f90`,
`tools/signal_reference/driver_signal_fft_bluestein.f90`,
`tools/signal_reference/compare_signal_fft_bluestein.py`, and
`benchmarks/bench_signal_fft_non_power_two.f90` are new project code.
SciPy 1.17 is used only as an independent numerical reference. No new
library license or dependency was introduced.

## Callable Fourier resampling windows (2026-10-26)

`src/signal_fourier_resample.f90` and
`src/signal_multirate_nd_into.f90` add original Fortran frequency-response
callbacks and one-shot per-axis FFT-bin evaluation. These APIs follow the
numerical contract described by SciPy's publicly documented `resample`
function; no SciPy source was copied into SciFort. The new native regression
`test/test_signal_resample_callable_window.f90` and independent reference
`test/compare_signal_resample_callable_window_scipy.py` with driver
`tools/signal_reference/resample_callable_window_driver.f90` were authored
for SciFort. SciPy 1.17 is used exclusively as a numerical test reference.
No additional external Fortran dependencies are required.

## Frequency-windowed Fourier resampling (2026-10-26)

`src/signal_fourier_resample.f90` and
`src/signal_multirate_nd_into.f90` extend the original MIT-licensed SciFort
FFT resampling and assumed-rank lane adapters to accept caller-supplied real
or complex spectral windows. FFT-bin pairing uses the mathematical real-FFT
symmetry convention; no SciPy source code was copied or incorporated.
`test/test_signal_resample_frequency_window.f90`,
`tools/signal_reference/resample_frequency_window_driver.f90` and
`test/compare_signal_resample_frequency_window_scipy.py` are original test
code. They use SciPy 1.17 only as an independent runtime numerical reference.
No extra library dependency is introduced.

### Windowed and streaming multitaper colored-noise backgrounds (2026-10-26)

The new SciFort interfaces in `src/signal_multitaper.f90` reapply the
original MIT-licensed local-median background estimator independently along
the frequency dimension of per-window, per-channel multitaper PSDs.
New tests `test/test_signal_multitaper_noise_floor_window_stream.f90`,
`test/compare_signal_multitaper_noise_floor_window.py` and
`reference/signal_multitaper_noise_floor_window_reference_driver.f90` are
original SciFort test code. The independent verification uses SciPy 1.17
DPSS windows, `signal.detrend`, and `ndimage.median_filter` and NumPy FFTs;
no third-party source was imported into the library. SciPy and NumPy are
comparison-only dependencies. The feature is explicitly a SciFort extension,
not a port of a SciPy public API.

## Whole-record and streaming multitaper detrending (2026-10-26)

`src/signal_multitaper.f90` extends the original MIT-licensed SciFort
multitaper implementation with per-channel least-squares detrending over a
full record or over streaming windows. It reuses the centered-coordinate
linear regression added to windowed time-frequency routines; no external
Fortran code or libraries were copied or introduced. New native tests and
`tools/signal_reference/driver_signal_multitaper_whole_stream_detrend.f90`
plus `compare_signal_multitaper_whole_stream_detrend.py` are original.
Validation uses SciPy 1.17 `signal.detrend` and `signal.windows.dpss`,
combined independently with NumPy FFTs for whole-record and streaming
complex spectral matrices. Adaptive taper weighting is cross-checked
against SciFort's independently computed batch results.

## Gustafsson Householder solve and multitaper linear detrending (2026-10-26)

The rank-revealing pivoted Householder QR routine in
`src/signal_filtfilt_gust.f90` and windowed per-segment linear least-squares
trend removal in `src/signal_multitaper.f90` are newly authored MIT-licensed
SciFort code. No third-party Fortran code is copied. Standard numerical
methods are used: conjugate-transpose QR for full-row-rank minimum-norm
least squares, and centered least-squares regression for segment trends.
Independent reference verification uses SciPy 1.17 `signal.sosfilt`,
`scipy.linalg.lstsq`, `signal.detrend` and `signal.windows.dpss` with NumPy
FFTs. Verification scripts are
`tools/signal_reference/compare_signal_sos_gust_scipy.py` and
`tools/signal_reference/compare_signal_multitaper_linear_detrend.py`.

## SOS-native Gustafsson filtering (2026-10-26)

Original MIT-licensed SciFort code in `src/signal_filtfilt_gust.f90`,
`src/signal_filtfilt.f90` and `src/signal_filtfilt_nd.f90` constructs
Gustafsson observability matrices from SOS direct-form-II states and uses
SciFort's existing Jacobi-SVD/pivoted-QR least-squares solvers, without
polynomial expansion. The implementation does not copy SciPy source and adds
no dependency. Native tests and the Python reference harness are original;
`scipy.signal.sosfilt` and `scipy.linalg.lstsq` are used *only* for independent
numerical validation, because `scipy.signal.sosfiltfilt` does not offer a
Gustafsson method. Reference: F. Gustafsson (1996), "Determining the initial
states in forward-backward filtering," IEEE Transactions on Signal Processing
46(4), 988-992.

## Assumed-rank Gustafsson filtering (2026-10-26)

`src/signal_filtfilt_nd.f90` extends its original MIT-licensed rank-independent
axis adapters to reuse the existing pure-Fortran Gustafsson 1D kernel and
its status semantics for all real/complex operand combinations. The axis
remapping and error validation are original SciFort code; no SciPy source
was copied. `test/test_signal_filtfilt_nd_gust.f90` and
`tools/signal_reference/{driver_signal_filtfilt_nd_gust.f90,compare_signal_filtfilt_nd_gust_scipy.py}`
are newly authored regression and numerical comparison tools. Real-coefficient
cases use SciPy 1.17 as an independent output reference; complex-coefficient
cases use a separately written NumPy/SciPy complex least-squares reference,
accounting for SciPy 1.17's real-only observability matrix. No additional
third-party Fortran code or package dependencies were added.

## Gustafsson initial-condition filtering (2026-10-26)

`src/signal_filtfilt_gust.f90` is an original MIT-licensed implementation of
the standard Gustafsson initial-state least-squares equations. The new complex
Jacobi-SVD and QR fallback are written in Fortran without external solvers or
transplanted SciPy source. `src/signal_filtfilt.f90` forwards the optional
Gustafsson method through existing rank-1/2/3 interfaces. SciPy 1.17 serves
as an output comparison for real filters, while a separate NumPy/SciPy linear
algebra reference checks the mathematically complex coefficient case.
References: F. Gustafsson, "Determining the initial states in forward-backward
filtering," *IEEE Transactions on Signal Processing* 46(4), 988-992 (1996).

## N-dimensional direct convolution and correlation (2026-10-25)

`src/signal_convolve_nd.f90` is an original MIT-licensed Fortran 2018
implementation of the standard multidimensional convolution sum and
complex-conjugated cross-correlation. It includes independent Fortran-order
indexing, singleton broadcasting, output-layout validation, direct and FFT
method dispatch, and caller-output error semantics. The FFT option calls the
existing SciFort packed-array convolution engine in
`src/signal_fftconvolve.f90`; SciPy code is not copied.
`test/compare_signal_convolve_nd_scipy.py` uses SciPy 1.17 only to check
outputs independently against `signal.convolve`, `signal.correlate`, and
`signal.fftconvolve` for selected-axis extension cases. No new package
runtime dependencies.

## Multichannel streaming STFT and inverse STFT (2026-10-25)

`src/signal_stream_multi.f90` is an original MIT-licensed Fortran 2018
multichannel adapter built on SciFort's existing single-channel streaming
STFT, spectrogram and inverse STFT kernels. It manages synchronized channels,
common time grids, public multidimensional allocations, noncontiguous input
sections and validation without copying SciPy source. SciPy 1.17 is used
only as a numerical reference in `test/compare_signal_stream_multi_scipy.py`,
with `tools/reference_signal_stream_multi.f90` as the Fortran driver. There
are no additional runtime dependencies.

## Streaming forward STFT and spectrogram (2026-10-25)

`src/signal_stft_stream.f90` is a newly written MIT-licensed streaming
analysis implementation. It reuses SciFort's independent `forward_fft` and
`get_window` routines but introduces its own bounded circular input buffer,
chunk-processing protocol, detrending, FFT normalization, final-frame
handling, and status-preserving validation. No SciPy source code has been
copied. `test/compare_signal_stft_stream_scipy.py` uses installed SciPy 1.17
solely as an independent numerical comparison, with the Fortran reference
program in `tools/reference_signal_stft_stream.f90`. The implementation
introduces no new external dependencies.

## Rank-independent inverse STFT (2026-10-25)

`src/signal_istft_nd.f90` is a newly written MIT-licensed Fortran 2018
adapter for the existing SciFort `src/signal_istft.f90` overlap-add kernel.
It independently implements Fortran-order channel addressing, caller-output
validation, time-axis removal and transactional reconstruction. It contains
no translated or imported SciPy source. The legacy SciFort `istft` kernel
provides the FFT and synthesis mathematics; SciPy 1.17 is used only for the
independent numerical regression script
`test/compare_signal_istft_nd_scipy.py`.

## Assumed-rank multirate signal processing (2026-10-25)

The Fortran 2018 caller-allocated rank-1 through rank-15 multirate APIs in
`src/signal_multirate_nd_into.f90` are an independent MIT-licensed extension.
They reuse the existing SciFort one-dimensional kernels and implement array
lane addressing, portable `select rank` pack/unpack, and staged-result error
handling without copying SciPy source. Source generation is limited to the
rank-selector cases in `tools/generate_multirate_assumed_rank.py`.
SciPy 1.17 is used only as a numerical comparator in
`tools/signal_reference/verify_multirate_assumed_rank.py`.

## Multidimensional multirate filtering (2026-10-25)

Rank-2 and rank-3 `scipy.signal`-style resampling overloads are independently
written MIT-licensed Fortran wrappers generated from
`tools/generate_signal_multirate_nd.py`. They apply pre-existing SciFort
one-dimensional Fourier, FIR/polyphase, and decimation kernels to each
independent signal channel. The original kernels and their mathematical
formulas have not been replaced. SciPy 1.17 serves only as a numerical
reference in `tools/signal_reference/check_multirate_nd.py`; no source from
SciPy is copied and no external dependencies are introduced.

## Continuous-time Gramians and reduction (2026-10-25)

The continuous-time solver in `src/signal_clti.f90` uses the exact Cayley
identity for Lyapunov equations: `Ad=(I+h*A/2)*(I-h*A/2)^(-1)`,
`Bd=sqrt(h)*(I-h*A/2)^(-1)*B`, and
`Cd=sqrt(h)*C*(I-h*A/2)^(-1)`. It reuses the preexisting SciFort
`cont2discrete` solver and discrete Smith doubling algorithm. The continuous
balanced projection is obtained from SciFort's original Hermitian Jacobi
implementation and applied to the original continuous model, following
standard balanced-realization algebra. These routines and tests are original
MIT-licensed Fortran/Python; SciPy's `solve_continuous_lyapunov` and SVD are
used as independent numerical references, not as transplanted code. The
mathematical state-space Gramian and balancing results are classical.

## Balanced realization and continuous LTI objects (2026-10-25)

`src/signal_dlti.f90` now provides independently implemented complex
Hermitian plane-rotation diagonalization and square-root balanced state-space
truncation using the discrete Lyapunov Gramians already supplied by SciFort.
The standard Gramian-factor/SVD construction and twice-the-tail error bound
are described in B. C. Moore, "Principal component analysis in linear
systems: Controllability, observability, and model reduction," IEEE
Transactions on Automatic Control, 26(1), 17-32 (1981).
The SciFort eigensolver was implemented from the Hermitian Jacobi algebra,
not imported from a third-party source. SciPy/NumPy are test references only.

`src/signal_clti.f90` independently adds a continuous-time LTI object API.
Analog simulation, response functions, conversion, and state-space algebra
reuse existing SciFort kernels, with no third-party code copied. Its design
is informed by the public SciPy signal `lti`/`StateSpace` interfaces; a
separate continuous-time type prevents using discrete poles/frequencies or
sampling intervals where analog dynamics are intended. Its 100-case MIMO
reference comparison uses SciPy's matrix exponential and complex linear
solvers, not external Fortran routines. Both additions retain the MIT
license and introduce no new dependencies.

## Discrete LTI feedback and reduction (2026-10-25)

`src/signal_dlti.f90` now contains independently written Fortran MIMO
state-space feedback elimination, block-diagonal interconnection, Kalman
reachability/observability matrices, numerical minimal-realization projections,
and infinite-horizon discrete Lyapunov Gramians via Smith doubling. The
computations use standard linear systems identities. SciPy 1.17 and NumPy
serve only as numerical oracles for tests and simulations; the implementation
copies no SciPy source code and adds no external library dependency.

## Reusable chirp-Z plans (2026-10-24)

The `czt_plan` and `zoom_fft_plan` interfaces in `src/signal_czt.f90` were
independently implemented for SciFort. They adapt the mathematical Bluestein
FFT-convolution factorization already in SciFort's `czt` routine and retain
its Horner strategy for nonunit spirals. Public SciPy 1.17 `CZT`/`ZoomFFT`
results are used solely as numerical test references. No SciPy source was
copied; no dependency or new license was introduced.

### 2026-10-24: Whittaker-Henderson and remaining signal utilities

`src/signal_whittaker.f90` independently assembles the weighted forward-
finite-difference penalty and solves `(W + lambda * D^T D) x = W y` by
banded Cholesky factorization, for fixed or automatically selected `lambda`.
The log-penalty REML optimizer evaluates the profiled likelihood using the
log determinant returned by that same banded Cholesky factorization, with
the mathematical criterion described in SciPy 1.18 and Biessy (2023)
(https://arxiv.org/abs/2306.06932). The
penalized least-squares formulation appears in Eilers (2003), "A Perfect
Smoother," *Analytical Chemistry* 75:3631-3636, and SciPy's
`whittaker_henderson` documentation.
`src/signal_tf_normalize.f90` and `src/signal_abcd_normalize.f90` implement
standard transfer-function coefficient normalization and state-space shape
consistency. `src/signal_convolve.f90` gains a cost-based direct/FFT selector
and optional CPU-time comparison, and `src/signal_fftconvolve.f90` now owns
its local status constants rather than importing them cyclically.

The implementations are original MIT-licensed Fortran code. SciPy 1.17 is a
numerical oracle for `normalize`, FFT convolution, and correlation;
`test/compare_signal_missing_extra_scipy.py` constructs independent sparse
weighted least-squares references for the smoother. The API of
`abcd_normalize` is tested with omitted/mixed-type and zero-state matrices.
No source code from SciPy or any paper has been copied, and there are no
new FPM dependencies.

### 2026-10-24: additional SciPy signal algorithms

`src/signal_stft_dual.f90` independently uses the orthogonal residue classes
of a discrete window translated by a fixed hop to solve the dual-window
least-distance problem. `src/signal_fir2d.f90` constructs separable windowed
FIR designs and radially interpolated oversampled FIR profiles. Existing
SciFort spline recurrences are reused by `src/signal_symiir.f90`; SciFort's
original frequency-grid logic is exposed through `findfreqs` and the chirp-Z
sampling points through `czt_points`. `src/signal_unique_roots.f90` adds a
straightforward independent greedy clustering algorithm. These Fortran
implementations use original MIT-licensed code and copy no SciPy source.

SciPy 1.17.0 is used strictly for numerical validation by
`tools/signal_reference/compare_signal_missing_algorithms_scipy.py` (499
comparisons, 11,610 values), plus native public-interface and edge-case
regressions. The reference script notes two existing SciPy behavior issues:
`firwin_2d` ignores its `pass_zero` argument, and `symiirorder1` can behave
differently on strided real views. No external source is included and there
are no additional FPM dependencies.

# Code provenance

### 2026-10-24: normalized analog lowpass filter prototypes

`src/signal_analog_prototypes.f90` independently implements `buttap`,
`cheb1ap`, `cheb2ap`, and `besselap`, with algebraic circular and elliptic
pole constructions, imaginary-axis Chebyshev-II zeros, and numerical roots
of scaled reverse Bessel polynomials. The simultaneous Bessel-root solver
uses Aberth-Ehrlich updates initialized with the published Campos-Calderon
polynomial approximations. DC, group-delay, and -3 dB magnitude normalization
are derived from the zero/pole/gain formulas, without copying SciPy code.

`test/test_signal_analog_prototypes.f90` checks native normalization
identities, conjugate symmetry, invalid inputs, and edge cases. The driver
and `tools/signal_reference/compare_signal_analog_prototypes_scipy.py` compare
all four families against SciPy 1.17.0, including Bessel normalization modes.
Bessel roots were initially limited to order 12 due to double-precision
conditioning; this release has independently verified order 13 as well.
The trigonometric prototypes support orders through 64. MIT-licensed code, no additional external dependencies.

### 2026-10-24: analog lowpass filter-frequency transformations

`src/signal_frequency_transforms.f90` independently implements eight
continuous-time lowpass transformations, using Horner-based rational
polynomial substitution for transfer-function coefficients and algebraic
root mapping for the ZPK representation. It does not copy SciPy source.
`tools/signal_reference/compare_signal_frequency_transforms_scipy.py`
checks 328 real and complex coefficient/root cases against SciPy 1.17.0,
with independent native transfer-function identities and invalid-input
regressions in `test/test_signal_frequency_transforms.f90`. The new code
is MIT-licensed and requires no new runtime dependency.


### 2026-10-24: affine-parameter spline Hessians and displacement curvature

`src/signal_spline_nd.f90` extends its independent tensor-product B-spline
Hessian calculation to affine matrix/offset derivatives and per-pixel
second-order displacement sensitivities. The implementation uses the
analytic chain rule in bounded interpolation blocks and does not copy or
translate SciPy's implementation. The public modules expose both real and
complex overloads. Native rank-15, validation and finite-difference tests,
plus `compare_spline_parameter_hessians_scipy.py` (SciPy 1.17.0 numerical
reference only), verify the new behavior. No third-party source or runtime
dependency has been introduced.


### 2026-10-24: analytic spline Hessians and shaped displacement Jacobians

`src/signal_spline_nd.f90` extends the original cardinal B-spline derivative
calculation to analytic second derivatives. Tensor-product diagonal and
mixed Hessian terms share the support traversal; grid-constant singleton
axes retain their offset-dependent first and second correction factors.
Grid-shaped displacement gradients wrap the existing bounded-transform
implementation. All code is an independent Fortran implementation; no
SciPy source was copied. `compare_spline_hessians_scipy.py` uses
SciPy 1.17.0 solely to generate numerical references. Native and numerical
tests, public module exports, and API documentation were added without
external dependencies.


### 2026-10-24: grid-shaped spline derivatives and affine sensitivities

`src/signal_spline_nd.f90` extends the existing cached B-spline gradient
kernel to multidimensional coordinate grids and propagates its coordinate
Jacobian to affine matrix and offset parameters using the chain rule. The
affine implementation reuses bounded coordinate batches and does not use
finite differences at runtime. The public exports through `scifort_signal`
and `scifort_ndimage` are backward compatible. Native tests include real and
complex inputs, noncontiguous arrays, rank-changing affine maps, rank-15
inputs, invalid parameters, and finite-difference spot checks. The separate
SciPy reference computes numerical parameter derivatives from SciPy 1.17.0
`ndimage.map_coordinates`. No SciPy code was copied or translated, and no
new dependencies were introduced.


### 2026-10-24: analytic cached spline gradients and displacement sensitivities

`src/signal_spline_nd.f90` adds an original analytic derivative extension of
SciFort's centered cardinal B-spline basis, with boundary-coordinate fold
signs, singleton-axis factors for grid-constant mode, and in-pass gradient
accumulation using prefix/suffix products. The displacement-field extension
shares the existing batched geometric-transform engine. No SciPy code is
copied or translated. SciPy 1.17.0 is used only as an independent numerical
reference through coordinate central differences in
`tools/signal_reference/compare_spline_plan_gradients_scipy.py`. Native
real/complex, rank-15, invalid-input, and strided-output tests live in
`test/test_signal_spline_plan_gradients.f90` and
`test/test_signal_spline_displacement_gradients.f90`. No new external
library dependency is introduced.

### 2026-10-23: multichannel vectorized ndimage callbacks

`src/ndimage.f90` independently adds `vectorized_filter_multi_nd` for real
and complex batched reductions returning multiple channels per pixel. The
new kernel reuses existing footprint extraction and boundary handling,
adding a common interior-window fast path for single- and multi-result
callbacks. Standard output has a trailing channel axis; flat output
accommodates rank-15 images. The implementation is original MIT-licensed
Fortran, with no translated or imported SciPy source.
`compare_vectorized_multi_scipy.py` checks SciPy 1.17.0 (using separate
scalar references for its reordered-axis multi-output layout issue), and
`test_ndimage_vectorized_multi.f90` checks output layouts, noncontiguous
arrays, complex values, resource guards, rank-15 input, and callback batches.
No new dependencies are required.

### 2026-10-23: bounded-memory vectorized ndimage callbacks

`src/ndimage.f90` independently implements `vectorized_filter_nd` with real
and complex batched Fortran procedure interfaces, active-neighbor extraction
in C-order, mixed and selected-axis boundaries, and `valid` shape reduction.
It also re-exports the existing coordinate-grid interpolation entry point.
The design reuses SciFort's original stencil and boundary-index helpers;
no SciPy source is copied. Numerical parity is checked against SciPy 1.17.0
by `compare_vectorized_filter_scipy.py`, while
`test_ndimage_vectorized.f90` covers rank-15, noncontiguous outputs, mixed
modes, callback batching, NaNs, resource guards, and invalid arguments.
No new runtime dependencies are needed.

### 2026-10-23: geometric callback interpolation, indexed comprehension and selected-axis labeling

`src/signal_spline_nd.f90` adds an independent callback-driven wrapper around
SciFort's existing coordinate interpolation kernel; `src/ndimage_measurements.f90`
uses the existing heapsort/binary-search label index to form C-ordered callback
buckets in a single pass; and `src/ndimage_morphology.f90` reuses the existing
connectivity-offset embedding for selected-axis connected-component labeling.
All changes are original MIT-licensed SciFort code: no SciPy implementation
has been imported. The `compare_geometric_comprehension_scipy.py` and
`compare_label_axes_scipy.py` scripts check SciPy 1.17.0 output and independent
selected-axis slice partitions; two native regressions cover high-rank,
noncontiguous, callback-order and invalid-input behavior. No new runtime
dependencies or third-party source are included.

### 2026-10-23: complex generic derivatives and indexed label statistics

`src/ndimage.f90` independently adds complex overloads for generic gradient
magnitude and Laplacian callbacks. The gradient uses SciPy's documented
complex sum-of-squares convention. `src/ndimage_measurements.f90` implements
an original heapsort/binary-search label request index for region reductions;
`src/ndimage_morphology.f90` extends binary propagation and hole filling to
selected axes. These are original MIT-licensed implementations and do not
copy SciPy source. `compare_complex_derivatives_labels_scipy.py` and
`compare_propagation_axes_scipy.py` validate against SciPy 1.17.0 and
independent selected-axis slicing; two native tests cover rank-15 and error
paths. No runtime dependency has been added.


### 2026-10-23: ndimage brute-force distances and generic axis derivatives

`src/ndimage_distance.f90` adds an independently written bounded direct-search
implementation of `distance_transform_bf_nd` with Euclidean, city-block, and
chessboard metrics, sampled Euclidean distances, and nearest-background
witnesses. `src/ndimage.f90` provides `generic_gradient_magnitude_nd` and
`generic_laplace_nd`, with user-supplied Fortran axis-derivative procedures.
No SciPy implementation source is copied: SciPy 1.17.0 is used as an
independent numerical reference for 1,350 distances and derivative values,
including arbitrary distance-tie feature witnesses, mixed boundary modes,
and rank-15 inputs. The regression and reference scripts are
`test_ndimage_bruteforce_generic.f90` and
`compare_bruteforce_generic_scipy.py`. No runtime dependencies are added.

### 2026-10-23: ndimage ellipsoid Fourier filter, mixed order modes, generic 1-D callbacks

`src/ndimage.f90` independently implements normalized uniform-ball Fourier
multipliers through integer/half-integer Bessel formulas and a convergent
small-argument series; per-axis boundary extensions for footprint order
statistics; and line-buffer callbacks equivalent to the public SciPy
`generic_filter1d` contract. These are original MIT-licensed Fortran
algorithms; SciPy BSD-3-Clause source is not copied. SciPy 1.17.0 reference
results and a separate Bessel-function/coordinate-space model verify 9,134
output components in `compare_ellipsoid_order_line_scipy.py`.
`test_ndimage_ellipsoid_order_line.f90` checks callback alignment,
validation and rank-15 support. No new FPM dependencies are required.


This file records the origin and licensing of numerical implementation code.
Every imported, translated, or substantially adapted routine must be entered
before it is merged.

### 2026-10-22: ndimage mixed-mode filtering and Fourier-domain primitives

`src/ndimage.f90` adds independent, MIT-licensed real/complex Fourier-domain
Gaussian, uniform and shift multipliers, plus per-axis filtering boundary
modes and complex-valued generic-neighborhood callbacks. Formulas use the
ordinary Gaussian Fourier transform, sinc response, and complex exponential
translation; no SciPy implementation source is copied. SciPy 1.17.0 is a
numerical reference for 19,584 mixed-mode scalar entries and 402 Fourier
coefficients. Full-kernel mixed modes and complex callbacks are independently
verified where SciPy's corresponding public API rejects them. Tests:
`tools/ndimage_reference/compare_mixed_complex_scipy.py`,
`tools/ndimage_reference/compare_fourier_scipy.py`, and
`test/test_ndimage_mixed_fourier.f90`. No additional FPM dependencies.

### 2026-10-22: ndimage callback filtering and weighted shortest paths

`src/ndimage.f90` independently implements complex-safe selected-axis
correlation/convolution and real generic neighborhood callback evaluation,
including the Fortran/C storage-order translation for active footprints.
`src/ndimage_measurements.f90` adds original labeled callbacks with pixel
selection in deterministic C-order and explicit missing-region defaults.
`src/ndimage_chamfer.f90` extends taxicab line scans to positive per-axis
integer weights and implements general directed chamfer distance and feature
transforms using an original multi-source binary-heap Dijkstra algorithm.
This is an MIT-licensed implementation, not a translation or redistribution
of SciPy source. SciPy 1.17.0 is used solely for behavioral comparison of
its corresponding public APIs; arbitrary weighted shortest paths are tested
with an independent Python heapq reference because SciPy's array-valued
`distance_transform_cdt(metric=...)` is a connectivity mask, not a weighted
shortest-path specification. Reproducible comparisons reside in
`compare_callbacks_weighted_scipy.py` and `compare_weighted_stencil_random.py`.
The added tests cover 34,296 elements; no new FPM dependencies are introduced.

### 2026-10-22: ndimage discrete distances and region histograms

`src/ndimage_chamfer.f90` is an independent, MIT-licensed implementation of
separable L1 (city-block/taxicab) and L-infinity (chessboard) distance
transforms. The L1 transform uses two linear scans per axis, and the L-infinity
transform searches minimum-radius coverage using a range-min segment tree.
`src/ndimage.f90` adds independent full-footprint extrema enumeration while
preserving the existing separable rectangular extrema fast path.
`src/ndimage_measurements.f90` adds labeled counts, histogram binning, and
population standard deviation wrappers. These implementations are original
Fortran code; SciPy 1.17.0 and NumPy serve as behavioral/numerical test
references only, not translated or copied source. See
`tools/ndimage_reference/{compare_extrema_histogram_scipy.py,compare_chamfer_random_scipy.py}`.
The test set adds 21,900 scalar comparisons, and the new source has no
additional external dependencies. Site tie-breaking can differ from SciPy.

### 2026-10-22: ndimage Gaussian, neighborhood and binary processing

`src/ndimage.f90` and `src/ndimage_morphology.f90` implement original
MIT-licensed finite Gaussian kernel correlation and derivative polynomials,
separable rectangular min/max/mean filtering, Sobel and Laplacian stencils,
centered Boolean neighborhood enumeration, binary erosion/dilation and
union-find connected-component labeling. Public algorithm definitions and
boundary conventions are those in the SciPy ndimage documentation:
https://docs.scipy.org/doc/scipy/reference/ndimage.html. The SciPy 1.17
functions are used as **numerical tests**, not copied implementation sources.
See `tools/ndimage_reference` for reproducible comparisons. No new runtime
libraries or vendored third-party code have been added.

`src/signal_spline_nd.f90` extends the independent cardinal B-spline
collocation code with a 12-element edge halo per nonsingleton dimension for
`nearest` and `grid-constant` high-order interpolation, and separable
singleton-axis treatment. The halo width and extension semantics follow the
SciPy ndimage interpolation reference, specifically
https://github.com/scipy/scipy/blob/main/scipy/ndimage/_interpolation.py
(SciPy Developers; BSD-3-Clause). The Fortran padding, coefficient
storage, and singleton optimization code is independently implemented. The
repository retains SciPy BSD-3 attribution in `THIRD_PARTY_LICENSES.md`.

## Current implementation inventory

- `src/signal_spline_nd.f90` (`zoom_nd`, `rotate_nd`, and shape helpers),
  `test/test_signal_zoom_rotate_nd.f90`, and
  `tools/signal_reference/{zoom_rotate_reference_driver.f90,compare_signal_zoom_rotate_scipy.py}`:
  independently written MIT-licensed transform wrappers reusing SciFort's
  spline interpolation and affine-sampling numerical kernels. SciPy 1.17
  `ndimage.zoom` and `rotate` were consulted for documented coordinate,
  output-shape, and angle conventions and used solely as numerical references.
  No SciPy source code was copied or redistributed, and no new dependencies
  were introduced. A 1,076-case comparison reports the known short-array
  `reflect` spline-prefilter difference separately from strict parity.


- `src/signal_spline2d.f90`, `test/test_signal_spline2d.f90`, and
  `tools/signal_reference/{spline2d_reference_driver.f90,compare_signal_spline2d_scipy.py}`:
  SciPy-compatible two-dimensional cubic/quadratic spline coefficients and
  spline smoothing, adapted from `scipy.signal._spline_filters` and
  symmetric IIR initial-condition formulas in SciPy 1.17.0 (BSD-3-Clause).
  Reuses the project's existing `sepfir2d` for spline reconstruction.
  SciFort extends the complex smoothing path by applying the real-coefficient
  operator componentwise; SciPy 1.17 rejects this case. Original SciPy
  attribution and its BSD license are retained in THIRD_PARTY_LICENSES.md.
  Forty-four independently generated SciPy/NumPy numerical cases passed;
  no new external FPM dependencies were introduced.


- `src/signal_spline.f90`, `test/test_signal_spline.f90`, and
  `tools/signal_reference/{spline_reference_driver.f90,compare_signal_spline_scipy.py}`:
  modern Fortran adaptation of the cubic and quadratic B-spline coefficient
  recurrences, cubic smoothing recursion, cardinal B-spline basis functions,
  and mirrored evaluation logic from `scipy.signal._spline_filters`
  (SciPy 1.17.0; BSD-3-Clause, copyright notice in `THIRD_PARTY_LICENSES.md`).
  The Fortran API additionally supports complex quadratic evaluation and
  rejects numerically unstable cubic smoothing poles. SciPy provides the
  independent seeded 480-case reference, with maximum coefficient error
  `1.351e-15` and maximum evaluation error `2.798e-14`. No new external
  Fortran dependencies were introduced.

- `src/signal_spectral.f90` (`gauss_spline`), `src/signal_design.f90`
  (`kaiser_beta`, `kaiser_atten`, `kaiserord`, parameterized `firwin`),
  `test/test_signal_kaiser_spline.f90`, and
  `tools/signal_reference/{kaiser_spline_driver.f90,compare_signal_kaiser_spline.py}`:
  independently written MIT-licensed implementations of published Gaussian
  B-spline approximation and Kaiser-window formulas. The numerical formulas
  follow Bouma et al. (2007) and Oppenheim/Schafer. Public SciPy 1.17 API
  behavior guided validation; SciPy was used only as a numerical reference
  (2,388 seeded cases), with no SciPy source files copied, translated, linked,
  or bundled and no new runtime dependencies.


- `src/signal_windows_extra.f90`, `src/signal_spectral.f90` (delegating
  `get_window` wrapper), `test/test_signal_windows_{extra,public}.f90`, and
  `tools/signal_reference/{window_reference_driver.f90,compare_windows_extra.py}`:
  independently written, MIT-licensed formulas for twenty standard spectral
  window families. SciPy 1.17 (BSD-3-Clause) is used only to produce
  numerical references and corroborate periodic-versus-symmetric semantics;
  no source code is copied or imported. No external dependencies added.

- `src/signal_iir_extra.f90` (`gammatone`, `iircomb`),
  `test/test_signal_auditory_comb.f90`, and
  `tools/signal_reference/{signal_auditory_comb_driver.f90,compare_signal_auditory_comb.py}`:
  modern-Fortran adaptation of the published gammatone and IIR comb filter
  coefficients used by SciPy 1.17.0 in `scipy.signal._filter_design`,
  which are distributed under BSD-3-Clause. The gammatone filter implements
  Slaney's 1993 auditory-filter-bank parameterization, and the notch/peak
  comb follows Orfanidis, *Introduction to Signal Processing*, ch. 11 (1996).
  The existing SciPy attribution in THIRD_PARTY_LICENSES.md applies.
  The new code is integrated using SciFort's MIT-licensed result type,
  with no new external Fortran dependencies. A seeded 1,350-case SciPy 1.17
  comparison checks both numerator and denominator coefficients.

- `src/signal_multitaper.f90` (`multitaper_stream`),
  `test/test_signal_multitaper_stream*.f90`,
  `test/compare_signal_multitaper_stream.py`, and
  `reference/signal_multitaper_stream_reference_driver.f90`:
  independently written MIT-licensed Fortran streaming state management,
  cached DPSS windows and cross-spectral accumulation, using existing
  SciFort multitaper kernels. SciPy's BSD-licensed DPSS implementation and
  NumPy FFT are numerical references for 80 randomized streaming cases;
  neither library's implementation code is copied or linked into SciFort.
  No new external dependencies were introduced.

- `src/signal_multitaper.f90` (time-frequency cross-spectral extensions),
  `test/test_signal_multitaper_tf_cross*.f90`,
  `test/compare_signal_multitaper_tf_cross.py`, and
  `reference/signal_multitaper_tf_cross_reference_driver.f90`:
  independent SciFort MIT-licensed implementation of DPSS-based
  time-frequency cross-spectra and magnitude-squared coherence. The
  mathematical conventions are consistent with SciPy `csd` (conjugate the
  first FFT) and Thomson multitaper weighting, but no upstream Fortran or
  SciPy source code was copied. Independent numerical reference calculations
  use SciPy's DPSS windows and NumPy FFT under their own licenses.
- `src/signal_cwt_peaks.f90`, `test/test_signal_cwt_peaks.f90`, and
  `tools/signal_reference/{cwt_peaks_driver.f90,verify_cwt_place_poles_scipy.py}`:
  modern-Fortran implementation adapting the Ricker continuous-wavelet
  transform and ridge-line linking/filtering conventions of SciPy 1.17's
  `scipy.signal._peak_finding` and `_wavelets` (BSD-3-Clause). The source was
  studied for edge/offset and noise-percentile compatibility, including
  its behavior for non-integer widths. See the SciPy attribution and BSD-3
  terms in THIRD_PARTY_LICENSES.md. Scientific antecedent: Du, Kibbe, and
  Lin, *Bioinformatics* 22(17), 2059-2065 (2006), DOI
  10.1093/bioinformatics/btl355. No original `peaks`, `spectrum`, or
  `dynamics` code was copied for this update.
- `src/signal_cwt_peaks.f90`, `test/test_signal_cwt_custom.f90`,
  `tools/signal_reference/cwt_custom_driver.f90`, and
  `tools/signal_reference/compare_signal_cwt_custom_scipy.py` (2026-10-26):
  generalize the previously attributed SciPy-compatible CWT to accept
  user-defined real `(points,width)` callbacks, reusing the existing ridge
  detector and signal convolution conventions. Original SciFort Fortran
  integration and tests; the custom-callback API and reversing of kernels
  are based on SciPy 1.17 `find_peaks_cwt` / `_wavelets._cwt` semantics
  (BSD-3-Clause; license already retained in THIRD_PARTY_LICENSES.md).
  The Python comparison driver calls SciPy's `_cwt` solely as an independent
  numerical oracle. No external source was added to the Fortran library.
- `src/signal_place_poles.f90`, `test/test_signal_place_poles.f90`, and
  `tools/signal_reference/{place_poles_driver.f90,verify_cwt_place_poles_scipy.py}`:
  original analytic implementation using Ackermann's controllability-matrix
  formula for SISO assignment and a minimum-norm right inverse for
  fully actuated MIMO systems. SciPy's `place_poles` supplies independent
  target-pole/eigenvalue references. The original release predated the
  iterative YT/KNV0 extension described below. No external Fortran dependency
  was introduced by the analytic algorithms.
- `src/signal_place_poles_mimo.f90`, the updated `src/signal_place_poles.f90`,
  `test/test_signal_place_poles_mimo.f90`, and
  `tools/signal_reference/verify_mimo_place_poles_scipy.py`: original
  pivoted-constraint-nullspace construction for real MIMO feedback.
  Each requested pole uses `(A-lambda*I)*x=B*u`; repeated real poles and
  conjugate pairs share feasible eigenspaces while the transfer columns
  are chosen to be independent. Real-pole rank-one orthogonal projection
  is informed by the Kautsky-Nichols-Van Dooren (1985) formulation, but
  the original implementation preceded the iterative YT/KNV0 extension.
  SciPy 1.17 `place_poles` is the independent target-spectrum reference
  (470 reproducible MIMO cases). No new Fortran runtime dependencies.

- `src/signal_place_poles_mimo.f90`, `src/signal_place_poles.f90`,
  `test/test_signal_place_poles_iterative.f90`, and the
  `tools/signal_reference/{place_poles_methods_driver.f90,verify_iterative_pole_placement_scipy.py}`
  additions: independent Fortran implementation of the real rank-one KNV0
  determinant update and YT real rank-two singular-vector and complex-pair
  Hermitian projection methods, following the Tits--Yang scheduling method.
  The scientific algorithms and ordering were studied in SciPy 1.17's
  BSD-3-Clause `scipy.signal._ltisys` and its original references:
  J. Kautsky, N. K. Nichols and P. Van Dooren (1985); A. L. Tits and
  Y. Yang (1996), *Globally Convergent Algorithms for Robust Pole Assignment
  by State Feedback*. This is a new Fortran implementation rather than a
  line-by-line source translation. The full SciPy BSD notice is already
  reproduced in THIRD_PARTY_LICENSES.md. No new external dependencies.
  Independent SciPy feasibility comparisons exercise 666 method-specific
  cases; the returned gains are not required to match because the real
  transfer basis and initialization differ.

- `src/signal_relative_extrema.f90`, `src/signal_peak_excursion.f90`,
  `test/test_signal_relative_extrema_excursion.f90`, and
  `tools/signal_reference/{extrema_excursion_driver.f90,verify_extrema_excursion.py}`:
  original SciFort Fortran implementations. SciPy's public `argrelmax`,
  `argrelmin`, and `argrelextrema` specifications and numerical results are
  independent references (1,350 randomized comparisons).
  `peak_detect` is an independent implementation of the excursion-based
  peak/valley state machine described by Eli Billauer, also used by
  Jason Christopherson's `peaks` (MIT license, 2023). We did not import,
  copy, or translate upstream `peaks` source files. Attribution and the
  upstream license are recorded in THIRD_PARTY_LICENSES.md.
  https://github.com/jchristopherson/peaks
  https://github.com/jchristopherson/peaks/blob/main/LICENSE
  No new dependencies were introduced. `peak_detect` is not a SciPy
  `scipy.signal` API; it is deliberately labeled a SciFort extension.


- `src/signal_denoise.f90`, `test/test_signal_denoise.f90`,
  `test/test_signal_denoise_public.f90`,
  `tools/signal_reference/denoise_reference_driver.f90`, and
  `tools/signal_reference/compare_signal_denoise_scipy.py`: original SciFort
  MIT-licensed Fortran median and Wiener filtering implementations. Three-way
  median quickselect, explicit neighborhood indexing, zero-padded box moments,
  and boundary alignment were written independently; the SciPy 1.17 public
  `medfilt`, `medfilt2d`, and `wiener` functions were used strictly as numeric
  references (1,050 seeded cases). No SciPy code copied, translated or
  linked, and no new runtime dependencies were introduced.
  Reference: https://docs.scipy.org/doc/scipy/reference/signal.html


- `src/signal_istft.f90`, `test/test_signal_istft.f90`,
  `test/test_signal_istft_public.f90`, the two Python scripts in `test/`,
  and their Fortran stream drivers in `tools/signal_reference/`:
  original SciFort MIT-licensed inverse-DFT and weighted overlap-add
  reconstruction plus periodic window-bin sum checks for COLA/NOLA.
  Public `scipy.signal.istft`, `check_COLA` and `check_NOLA` in SciPy 1.17
  are independent numerical references. No source code was copied from
  SciPy, and no new third-party runtime dependencies were added.
  API reference: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.istft.html


- `src/signal_fftconvolve.f90` (kernel-spectrum reuse and multi-axis tiling),
  `test/test_signal_oaconvolve_tiled_multiaxis.f90`,
  `validation/reference_drivers/signal_oaconvolve_tiled_reference_driver.f90`,
  `test/compare_signal_oaconvolve_tiled_scipy.py`, and
  `benchmarks/bench_signal_oaconvolve_tiling.f90`: independently written
  optimization of the existing SciFort overlap-add algorithm. The FFT of an
  unchanged operand is cached within a call; the longer operand is tiled along
  eligible axes with a reusable workspace and cropped accumulation. SciPy
  1.17 was used only for 232 independent seeded comparisons; no SciPy code
  was copied, translated or linked. MIT licensed, no new dependencies.


- `src/signal_fftconvolve.f90` (`oaconvolve_nd`),
  `test/test_signal_oaconvolve_assumed_rank.f90`,
  `test/test_signal_oaconvolve_assumed_rank_public.f90`, and
  `tools/signal_reference/{oaconvolve_assumed_rank_driver.f90,verify_oaconvolve_assumed_rank.py}`:
  independently authored Fortran 2018 multidimensional blocked overlap-add
  using the existing SciFort radix-2 FFT and shape/broadcast helpers. The
  public SciPy 1.17 `scipy.signal.oaconvolve` API is used as an independent
  numerical reference for 400 seeded input cases; no SciPy implementation
  code was copied, translated, bundled, or linked. MIT-licensed SciFort
  source, no new third-party dependencies.

- `src/signal_filtfilt_nd.f90`, `test/test_signal_filtfilt_nd.f90`,
  `test/compare_signal_filtfilt_nd_scipy.py`, and
  `tools/signal_reference/filtfilt_nd_driver.f90`: original SciFort
  MIT-licensed Fortran 2018 assumed-rank adapters around the project's
  existing independent 1D forward/backward filtering routines. All lanes
  use the same documented SciPy-style padding/steady-state scheme (or,
  for `filtfilt_nd(method='gust')`, the separately implemented Gustafsson
  state solver) as their rank-1 counterparts; real/complex promotions
  are performed per lane. No source code was copied
  from SciPy. SciPy 1.17 served only as a numerical reference for 684 cases,
  with maximum absolute discrepancy 3.33e-16. No new dependencies.

- `src/signal_sos.f90`, `test/test_signal_sos_nd.f90`, and
  `test/compare_signal_sos_nd_scipy.py`: original SciFort MIT-licensed
  rank-independent SOS filtering implementation. It evaluates the standard
  cascaded transposed direct-form-II difference equations using independent
  per-channel delay states. Fortran 2018 `iso_c_binding` storage remapping
  supplies rank-one views of caller-owned contiguous or compiler-temporary
  buffers. For rank 15, state arrays are flattened into a documented packed
  rank-one representation rather than exceeding the Fortran rank limit.
  SciPy 1.17.0 `signal.sosfilt` serves only as a numerical reference for
  320 independently generated comparisons; neither SciPy nor any external
  implementation code was imported, copied, or vendored. Reference API:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.sosfilt.html


- `src/signal_lfilter.f90`, `test/test_signal_lfilter_nd.f90`, and
  `test/compare_signal_lfilter_nd_scipy.py`: independently authored
  MIT-licensed rank-independent extension of transposed direct-form-II
  filtering. An assumed-rank Fortran 2018 entry point uses standard
  `iso_c_binding` rank-one storage views for contiguous real/complex buffers,
  with compiler-managed copy-in/out for array sections. It reuses normalized
  coefficients and per-channel delay states without external dependencies.
  The SciPy 1.17 public `lfilter` API supplies numerical references only;
  no SciPy source code was copied or vendored.


- `src/signal_savgol.f90`, `test/test_signal_savgol_nd.f90`, and
  `test/compare_signal_savgol_nd_scipy.py`: original MIT-licensed
  rank-independent Fortran refactor. The assumed-rank API uses standard
  `iso_c_binding` storage remapping to traverse arbitrary-rank contiguous
  column-major arrays through a shared one-dimensional convolution kernel;
  Fortran manages temporaries for noncontiguous inputs and outputs. All
  polynomial filtering code and new tests are independently authored;
  SciPy 1.17 serves solely as a numerical reference. No source code was
  copied from SciPy and no external library is required.


- `src/signal_savgol.f90`, `test/test_signal_savgol_rank4.f90`,
  `test/compare_signal_savgol_rank4_scipy.py`, and
  `benchmarks/benchmark_savgol_multichannel.f90`: original MIT-licensed
  engineering changes to reuse the existing independent QR-based polynomial
  evaluation weights across signal channels, plus rank-four slicing and
  independent complex filtering. The SciPy public interface is used only
  to validate 400 seeded numerical cases; no SciPy implementation source
  is copied or linked. No additional Fortran dependencies are introduced.

- `src/signal_savgol.f90`, `src/signal_lfilter.f90`,
  `test/test_signal_savgol_axes.f90`, `test/test_signal_lfilter_mixed.f90`,
  `test/compare_signal_savgol_axes_scipy.py`, and
  `test/compare_signal_lfilter_mixed_scipy.py`: independently written
  MIT-licensed extension of existing SciFort algorithms. Rank-2/3
  Savitzky-Golay filtering applies existing polynomial least-squares
  weights along independent array lanes; complex inputs use separate
  real and imaginary filtering. Mixed-coefficient IIR wrappers promote
  their inputs to the already implemented complex direct-form-II routine.
  The public SciPy API is used solely for 600 Savitzky-Golay and 540
  mixed-coefficient numerical reference comparisons. No SciPy source code
  was copied and no third-party dependencies were introduced.

- `src/signal_sos.f90` and `src/signal_filtfilt.f90`: independently authored
  SciFort MIT-licensed rank-two/rank-three wrappers around the existing
  1D SOS and zero-phase filtering algorithms. The array slicing and
  one-based axes follow the Fortran language; SOS state arrays follow
  SciPy's publicly documented leading-section-dimension convention.
  `test/test_signal_sos_axes.f90`, `test/test_signal_filtfilt_axes.f90`,
  and `tools/signal_reference/{filter_axes_driver.f90,verify_filter_axes.py}`
  are newly authored regression/reference tools. 480 seeded numerical
  comparisons against SciPy 1.17 passed; no SciPy implementation code was
  copied and no third-party dependencies were added.


- `src/signal_fftconvolve.f90`, `src/signal_lfilter.f90`, and
  `test/test_signal_fftconvolve_nd.f90`, `test/test_signal_lfilter_axes.f90`,
  `test/test_signal_lfilter_rank3.f90`: independent MIT-licensed extensions
  using SciFort's existing radix-2 FFT and direct-form-II transposed filter.
  The multidimensional FFT uses separable axis-wise transforms and
  broadcasting. The filter routines apply the existing 1D difference
  equations to independent data lanes, with each lane's own initial/final
  states. SciPy 1.17 is used only for independent numerical reference values
  in `tools/signal_reference/verify_fftconvolve_nd.py`,
  `verify_lfilter_axes.py` and `verify_lfilter_rank3.py` (1,700 seeded cases).
  No SciPy code was copied or vendored; no new FPM dependencies.


- `src/signal_convolve2d.f90`, `test/test_signal_convolve2d.f90`,
  `tools/signal_reference/convolve2d_driver.f90`, and
  `tools/signal_reference/compare_signal_convolve2d_scipy.py` are original,
  MIT-licensed Fortran/Python implementations of the elementary discrete
  convolution sum and reversed-conjugated cross-correlation sum on two
  axes. Periodic and half-sample symmetric extensions are computed with
  modular indexing. Output trimming follows the public SciPy 1.17 API
  specifications; no third-party source code was copied, adapted, or
  vendored. The 1,600-case SciPy comparison is numerical verification,
  not an implementation dependency. Public API references:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.convolve2d.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.correlate2d.html


- `src/signal_lti_frequency.f90`, `test/test_signal_lti_schur.f90`,
  `tools/signal_reference/lti_schur_driver.f90`,
  `tools/signal_reference/compare_signal_lti_schur_scipy.py`, and
  `benchmarks/benchmark_signal_mimo_schur.f90`: independent MIT-licensed
  implementation of shifted implicit complex QR for Hessenberg matrices,
  unitary similarity updates to input/output matrices, deflation, a bounded
  iteration count, and triangular multi-right-hand-side solves. Standard
  mathematical background: Golub and Van Loan, *Matrix Computations*,
  fourth edition (2013), QR algorithm and Schur decomposition. No external
  LAPACK or SciPy implementation was adapted or copied. SciPy 1.17 is used
  only for independently calculated frequency-response reference values.

- `src/signal_lti_frequency.f90`, `test/test_signal_lti_mimo_plan.f90`,
  `tools/signal_reference/lti_mimo_channels_driver.f90`,
  `tools/signal_reference/compare_signal_lti_mimo_channels_scipy.py`, and
  `benchmarks/benchmark_signal_mimo_plan.f90`: independently written
  MIT-licensed reuse of the existing unitary Hessenberg reduction and
  one-based selection of MIMO input/output columns and rows. The seeded
  comparison independently evaluates the selected state-space resolvent
  with `scipy.linalg.solve`; no SciPy implementation code was copied.
  Default plans still store Hessenberg matrices. Users may now explicitly
  request the new independent complex-Schur reduction.


- `src/signal_lti_frequency.f90`, `test/test_signal_lti_hessenberg.f90`,
  and `benchmarks/benchmark_signal_mimo_hessenberg.f90`: independently
  written MIT-licensed complex Householder reduction and pivoted upper
  Hessenberg linear solver. The transform is unitary and reuses reduced
  matrices across all frequencies in one MIMO call. No LAPACK or SciPy
  implementation code was copied. A deterministic direct-vs-reduced test,
  480 seeded SciPy 1.17 complex/real MIMO references, and an optimized
  CPU-time benchmark document numerical agreement and runtime behavior.
  This path remains the default; a separate opt-in Schur implementation
  was added subsequently.


- `src/signal_lti_frequency.f90`, `test/test_signal_lti_mimo.f90`,
  `tools/signal_reference/lti_mimo_driver.f90`, and
  `tools/signal_reference/compare_signal_lti_mimo_scipy.py`:
  original MIT-licensed extension for MIMO state-space frequency responses
  and Bode plots, with a shared pivoted solve for multiple right-hand sides.
  No SciPy implementation source was copied. The 480-case seeded reference
  compares against independently computed `scipy.linalg.solve` state-space
  resolvents (22,648 output-input-frequency channel values), including
  complex systems, high orders, and static gain. The tests also compare
  magnitude and independently unwrapped channel phases. These interfaces
  are separate from SciPy's SISO convenience APIs and introduce no extra
  third-party dependencies.


- `src/signal_lti_frequency.f90`, `test/test_signal_lti_frequency.f90`,
  `tools/signal_reference/lti_frequency_driver.f90`, and
  `tools/signal_reference/compare_signal_lti_frequency_scipy.py`:
  original MIT-licensed wrappers and Bode magnitude/phase calculation.
  Reuses the existing SciFort analog/digital polynomial and zero/pole/gain
  frequency-response evaluators. The state-space path now uses a separately
  written shifted linear solver with partial pivoting, avoiding polynomial
  expansion. The new regression suite and seeded comparison against SciPy
  linear solves cover 560 systems and 5,600 frequency points. The
  discrete TF path zero-pads unequal polynomial degrees before evaluating
  the equivalent delay form, and continuous phase is unwrapped explicitly.
  No SciPy code is copied or imported. SciPy 1.17 is used only as an
  independent numerical reference (720 seeded cases, plus Fortran tests).
  Public API references:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.freqresp.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.bode.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.dfreqresp.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.dbode.html

- `src/signal_continuous_simulation.f90`,
  `test/test_signal_continuous_auto_time.f90`,
  `test/test_signal_continuous_mixed_zpk.f90`, and
  `tools/signal_reference/compare_signal_continuous_auto_time_scipy.py`:
  original MIT-licensed Fortran extension providing automatic response
  sampling from state-matrix characteristic roots (with triangular-matrix
  shortcuts), interpolation-safe optional argument forwarding, and
  mixed real/complex ZPK input and gain wrappers. Default time scales
  follow the public SciPy 1.17 API, with 184 independent numerical checks
  and no imported third-party code or new dependencies.

- `src/signal_continuous_simulation.f90`,
  `test/test_signal_continuous_simulation.f90`, and
  `tools/signal_reference/compare_signal_continuous_simulation_scipy.py`:
  original MIT-licensed Fortran implementation of classical LTI
  zero-order-held and first-order-held state propagation, step and impulse
  response formulas. Reuses the existing matrix exponential in SciFort
  without importing SciPy source or adding new dependencies. Validation
  uses SciPy 1.17 `lsim`, `step`, `impulse` and SciPy's `linalg.expm`
  where complex inputs or MIMO responses exceed the convenience routines.
  Reference: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lsim.html



- `src/signal_discrete_simulation.f90`,
  `test/test_signal_discrete_simulation.f90`, and
  `tools/signal_reference/compare_signal_discrete_simulation_scipy.py`:
  original MIT-licensed Fortran implementation of the standard discrete
  state-space recurrence and first-order input interpolation, supporting
  discrete system simulation, step responses and impulse responses.
  Transfer-function and ZPK calls reuse existing SciFort controller-canonical
  conversions. SciPy 1.17 serves as an independent numerical reference for
  498 seeded comparisons (with a complex-state recurrence reference where
  its `dlsim` loses imaginary parts). No SciPy or other external code has
  been imported and no new FPM dependencies are introduced. API references:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.dlsim.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.dstep.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.dimpulse.html

- `src/signal_cont2discrete.f90`, `test/test_signal_cont2discrete.f90`,
  and `tools/signal_reference/compare_signal_cont2discrete_scipy.py`:
  independent MIT-licensed Fortran implementation of standard LTI
  continuous-to-discrete conversions, including the classical scaled
  Padé(13) scaling-and-squaring matrix exponential, Gaussian elimination
  with partial pivoting, first-order-hold augmented matrix, generalized
  bilinear transformation and exact ZPK rational-substitution formulas.
  SciPy 1.17 is used for 390 seeded numerical comparisons only, not as
  copied source. Uses existing SciFort state-space/ZPK conversions and
  introduces no dependencies. API reference:
  https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.signal.cont2discrete.html .

- `src/signal_state_space.f90`, `test/test_signal_state_space.f90`, and
  `tools/signal_reference/compare_signal_state_space_scipy.py`: original
  MIT-licensed implementation of controller-canonical companion matrices
  and the classical Faddeev-LeVerrier recursion for characteristic
  polynomials and adjugate coefficients, with wrapper calls to SciFort's
  existing polynomial/ZPK conversion routines. Mathematical and numerical
  behavior was independently checked against the SciPy 1.17 `tf2ss`,
  `ss2tf`, `zpk2ss` and `ss2zpk` APIs (BSD-3-Clause), for 832 randomized
  cases. No external code was copied or vendored and no new FPM dependency
  is required. Public API reference:
  https://docs.scipy.org/doc/scipy/reference/signal.html .

- `src/signal_zpk_conversion.f90` independently implements classical
  transfer-function polynomial factorization, leading-coefficient
  normalization, monic polynomial expansion from roots, digital
  second-order section construction with conjugate-aware root pairing,
  minimal pairing by removal of temporary origin factors, inverse
  section factorization/gain accumulation, and `tf2sos` composition for
  analog/digital transfer functions. `keep_odd` pairing now includes a
  nearest-unassigned fallback for very unequal pole/zero orders.
  It calls SciFort's existing polynomial-root solver and uses no new
  external source code or dependencies. SciPy 1.17 is used solely for
  independent numerical comparisons (840 earlier randomized cases, plus
  792 minimal designs, 759 inverse SOS cases, and 992 valid TF-to-SOS
  comparisons including analog/digital filters).
  API references: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.tf2sos.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.tf2zpk.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.zpk2tf.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.zpk2sos.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.sos2zpk.html .

- `src/signal_czt.f90` is an independent MIT-licensed Fortran implementation
  of the classical chirp Z-transform identity (Bluestein's FFT convolution),
  direct Horner polynomial evaluation on nonunit contours, and the uniform
  frequency mapping of a zoom FFT. It reuses SciFort's own FFT kernel;
  no SciPy implementation code was copied or vendored. SciPy 1.17 provides
  an optional numerical reference only. API references:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.czt.html
  and https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.zoom_fft.html .


- `src/signal_remez.f90` independently implements symmetric and antisymmetric
  equiripple FIR design using the barycentric Remez exchange theorem, cosine
  polynomial interpolation, Chebyshev coefficient recovery, and standard
  linear-phase FIR reconstruction formulas. The numerical/API reference is
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.remez.html .
  No SciPy C/C++ or external source code was copied, translated, or vendored;
  this MIT-licensed Fortran module and its tests are original. SciPy 1.17 is
  used only as an optional numerical regression reference. There are no new
  FPM dependencies. The October 2026 high-order improvements use standard
  Fortran selected-real-kind arithmetic and an original parity-prefix
  extrema-selection implementation, without imported numerical code.


- `src/signal_fir_advanced.f90` independently implements `firwin2`, `firls`,
  and `minimum_phase` using frequency interpolation and inverse FFT, analytic
  weighted cosine integrals and Cholesky solution of normal equations, and
  standard homomorphic/modified-Hilbert phase conversion. The public SciPy
  API and the numerical formulas described in signal-processing literature
  provided the implementation basis. SciPy source was inspected to verify
  edge-case behavior, but no source files were copied or imported. The
  independently written MIT-licensed module and test drivers are compared
  against SciPy 1.17 numerically, with
  no new third-party code or external FPM dependencies. API references:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.firwin2.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.firls.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.minimum_phase.html


- `src/signal_analog_sos.f90` independently implements analog
  zero-pole-gain and second-order-section output for the four existing
  IIR prototype families. Real-root and conjugate-pair factorization,
  transfer-function product preservation and Horner-based s-plane response
  evaluation follow standard polynomial and LTI filter theory. No SciPy
  or external implementation code was copied; SciPy 1.17 is used as an
  independent reference for randomized numerical tests. The small change
  to `src/signal_analog_design.f90` exposes the internally computed roots
  and gain, with an optional root-only path that skips BA polynomial
  expansion while preserving the existing BA coefficients and behavior.
  Mathematical/API references:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.iirfilter.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.zpk2sos.html
  The implementation, deterministic tests and comparison driver are
  original MIT-licensed SciFort work without new dependencies.

- `src/signal_analog_design.f90` independently constructs normalized
  Butterworth, Chebyshev I, and Chebyshev II analog poles and transmission
  zeros, applies lowpass/highpass/bandpass/bandstop transformations to pole-
  zero-gain models, and multiplies polynomial factors to obtain real BA
  coefficients in descending powers of s. The existing `butter`, `cheby1`,
  `cheby2`, and `iirfilter` APIs now accept `analog=.true.`; no algorithms
  from SciPy were translated or copied. Public mathematical/API references:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.cheby1.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.cheby2.html
  and standard analog filter prototype theory. The new MIT-licensed source,
  Fortran tests, comparison driver and 1,200-case SciPy 1.17 harness add no
  third-party code, licenses, or FPM dependencies.

- `src/signal_analog.f90` and the additions to `src/signal_design.f90`
  independently implement analog transfer-function and zero-pole-gain
  responses and digital zero-pole-gain responses. Standard Horner polynomial
  evaluation, complex product evaluation, and iterative polynomial-root
  estimation supply the numerical methods. Mathematical/API references:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.freqs.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.freqs_zpk.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.freqz_zpk.html .
  New source, deterministic Fortran tests, and the SciPy-reference driver are
  original MIT-licensed SciFort work. No SciPy source was copied or imported,
  and no external dependencies or licensing obligations were introduced.

- `src/signal_periodic.f90` independently implements `sawtooth`, `square`,
  and `unit_impulse` from the piecewise definitions of periodic waves and the
  Kronecker delta, described in the public SciPy documentation:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.sawtooth.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.square.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.unit_impulse.html
  This original MIT-licensed Fortran code uses no upstream numerical source.

- `src/signal_sequences.f90` independently implements a binary Fibonacci
  linear-feedback shift register. The algorithm and the primitive-polynomial
  default tap convention follow the mathematical description and publicly
  referenced primitive-polynomial table listed at
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.max_len_seq.html .
  The mathematical source for the recurrence is the public overview of
  maximum-length sequences (https://en.wikipedia.org/wiki/Maximum_length_sequence).
  Feedback taps were checked against SciPy 1.17.0 outputs. No SciPy/Cython
  implementation code was imported, translated, or copied. The two modules,
  Fortran tests and optional SciPy regression driver are original MIT-licensed
  SciFort work without new FPM dependencies or third-party code obligations.


- `src/signal_hilbert.f90` independently implements FFT-domain analytic
  signals and the 2D single-orthant transform using standard positive-
  frequency weighting. The mathematical/API reference is the public SciPy
  documentation at
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.hilbert.html
  and
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.hilbert2.html .
  The module, deterministic Fortran tests, and optional 780-case SciPy 1.17
  reference comparison are original MIT-licensed SciFort work. No external
  numerical implementation was copied, and no new FPM dependencies or
  third-party license obligations were introduced.

- `src/signal_lfilter.f90` (additional `lfiltic` routines) and
  `src/signal_design.f90` (additional `group_delay` routines) contain original,
  independently written SciFort implementations of the transposed
  direct-form-II difference equation for past-sample initialization and the
  logarithmic phase-derivative formula for digital-filter group delay.
  Mathematical and API references are the public documentation at
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lfiltic.html
  and https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.group_delay.html .
  No SciPy implementation code was copied or adapted. The deterministic
  regression tests and optional 600-case SciPy 1.17 verification driver are
  original MIT-licensed SciFort work, with no additional external dependency.
  SciFort intentionally reports undefined group delay as NaN, not SciPy's zero.

- `src/signal_design.f90` also provides an independently written real SOS-to-
  transfer-function conversion (`sos2tf`). It multiplies ordinary quadratic
  polynomials and follows the coefficient-ordering convention documented at
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.sos2tf.html .
  No SciPy source code was copied. The tests compare selected outputs with
  SciPy 1.17.0 and with existing filter-design routines. No new dependency.
- `src/optimize/nonlin.f90` received a pure callback-dispatch restructuring
  to preserve the explicit `nonlin_linear_operator` interface across a nested
  internal Krylov callback. No external solver was imported or changed.

- `src/signal_iir_extra.f90` independently implements inverse-Chebyshev
  (Chebyshev-II) digital filter design, bilinear-transformed frequency
  selections, real biquad section pairing, scalar Butterworth/Chebyshev
  minimum-order design, and standard second-order notch/peak responses.
  The formulas are standard DSP theory; public numerical/API references are
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.cheby2.html
  and https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.iirfilter.html
  and Orfanidis, *Introduction to Signal Processing* (1996), Chapter 11.
  The module, regression test, and optional 1,350-case SciPy 1.17 reference
  harness are original SciFort work, without copied upstream source code.
  No additional third-party license notice or runtime dependency is needed.


- `src/signal_fourier_resample.f90` and updates to `src/signal_design.f90`
  and `src/signal_resample.f90` independently implement periodic Fourier
  resampling, analog Chebyshev-I lowpass poles transformed into digital
  low/high/bandpass/bandstop filters, and IIR decimation using existing SciFort
  SOS filtering. These are implementations of standard frequency truncation,
  Nyquist-pair handling, bilinear transforms, and filter cascades as described
  in the public SciPy documentation:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.resample.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.cheby1.html
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.decimate.html
  No third-party implementation code was copied or imported. The source,
  deterministic Fortran tests, and the optional randomized SciPy verification
  script are original MIT-licensed SciFort work. The Chebyshev design reuses
  SciFort's existing Butterworth frequency transformation and SOS construction.
  No new third-party license obligations or FPM dependencies were introduced.


- `src/signal_design.f90`: independently authored SciFort code for
  windowed-sinc FIR filter design, polynomial frequency responses,
  analog Butterworth pole transformations, bilinear mapping, and real
  second-order sections. Formulas follow standard DSP theory and the public
  SciPy `signal.firwin`, `butter`, `freqz`, and `sosfreqz` documentation
  (https://docs.scipy.org/doc/scipy/reference/signal.html); independent
  reference comparisons use SciPy 1.17.0. No `spectrum`, SciPy, or other
  third-party implementation source was copied or adapted, so no new
  third-party license notices are required. The regression and optional
  `tools/signal_reference/{design_driver.f90,verify_design.py}` are original
  SciFort code (MIT). The algorithms are not claimed to reproduce SciPy's
  floating-point rounding or SOS section order exactly.

- `src/signal_sos.f90` and `src/signal_filtfilt.f90`: independently authored
  SciFort implementations of cascaded transposed-direct-form-II filters,
  steady-state step initial conditions, and forward/backward edge-reflected
  filtering. Mathematical/API reference: the public SciPy documentation for
  `scipy.signal.sosfilt`, `sosfilt_zi`, `filtfilt`, and `sosfiltfilt`,
  https://docs.scipy.org/doc/scipy/reference/signal.html . No third-party
  implementation source was imported or copied. Verified against SciPy 1.17.0
  for 2,100 randomized real/complex cases. New tests and the optional
  `tools/signal_reference/{sos_driver.f90,verify_sos.py}` are original SciFort
  work under MIT. No changes to the third-party license inventory are needed.


- `src/signal_spectral_complex.f90`: independently authored complex-input,
  two-sided extensions to the SciFort FFT/Welch/STFT layer. The definitions
  follow standard spectral estimation formulas. Compared with SciPy 1.17.0
  across 960 randomized cases. No third-party source code was copied.
- `src/signal_fftconvolve.f90` (`fftconvolve_nd`),
  `test/test_signal_fftconvolve_assumed_rank.f90`, and
  `tools/signal_reference/{fftconvolve_assumed_rank_driver.f90,verify_fftconvolve_assumed_rank.py}`:
  independent MIT-licensed rank-independent extension based on the discrete
  convolution theorem and separability of the multidimensional DFT. Its
  column-major index arithmetic and Fortran 2018 assumed-rank interface are
  original to SciFort. The numerical semantics follow the public SciPy signal
  FFT-convolution specification; 576 seeded rank-15 comparisons to SciPy 1.17
  were performed. No SciPy source was copied or imported. The FFT kernel is
  reused from SciFort; no third-party dependency or license text was added.

- `src/signal_fftconvolve.f90`: independently authored radix-2 convolution
  and blocked overlap-add algorithms, checked against SciPy 1.17.0 in
  1,000 randomized real/complex cases. No third-party code copied. The FFT
  kernel is shared with the original SciFort spectral module.
- `tools/signal_reference/{spectral_complex_driver.f90,fftconvolve_driver.f90}`
  and corresponding `verify_*` scripts, plus new deterministic regressions,
  are original SciFort test implementations.

- `src/signal_spectral.f90`: original SciFort implementation of one-sided
  spectral analysis and window generation. Its Welch segment-averaging structure
  and modular spectral-estimator design were informed by Jason Christopherson's
  MIT-licensed `spectrum`, particularly `src/spectrum_periodogram.f90` and
  `src/spectrum_windows.f90` (accessed 2026-10-15). Attribution and the
  upstream copyright/license appear in `THIRD_PARTY_LICENSES.md`. The radix-2
  FFT, arbitrary-length direct DFT, SciPy-style normalization, detrending,
  boundary handling, result types, and numerical tests were authored for
  SciFort; no source code was copied verbatim. Optional test driver
  `tools/signal_reference/spectral_driver.f90` and randomized comparisons in
  `verify_spectral.py` are original SciFort code. The 1,640 comparisons against
  SciPy 1.17.0 were run on supported real-valued parameter combinations.


- `src/signal_peaks.f90`: original SciFort implementation of one-dimensional
  local peak detection (including plateaus), geometric contour prominence,
  and interpolated peak widths. Mathematical semantics and tests were checked
  against the public `scipy.signal.find_peaks`, `peak_prominences`, and
  `peak_widths` documentation and SciPy 1.17.0 output, but no third-party
  source code was imported, copied, or translated. Per-sample peak interval
  bounds extend the independent implementation and are checked in 2,200
  additional seeded SciPy comparison cases. The open-interval property flags
  follow SciPy's documented `(None, None)` request semantics, validated by
  2,400 independent comparisons. Their implementation, native regressions and
  SciPy reference drivers are original SciFort work, with no imported code.
- `src/signal_detrend.f90`: original SciFort real/complex rank-1/2/3 least-
  squares constant/linear detrending, with independently implemented segmented
  linear fits and one-based axes. The rank-1 real interface remains compatible.
  Verified against SciPy 1.17.0 `signal.detrend` in 720 seeded cases; no SciPy
  code was copied.
- `src/signal_convolve.f90`: original SciFort direct finite convolution and
  cross-correlation of real and complex sequences, and lag indexing, using
  standard discrete signal definitions. SciPy 1.17.0 reference results
  establish compatible modes and complex conjugation. No FFT algorithm or
  third-party code has been imported.
- `src/signal.f90`, `test/test_signal_*.f90`: original SciFort aggregator
  and deterministic self-contained tests; comparison values were computed
  separately with SciPy 1.17.0. All new files use MIT licensing.

- `src/optimize/assignment.f90`: original SciFort implementation of the
  rectangular primal-dual Hungarian shortest-augmenting-path formulation,
  with O(n+m) auxiliary memory. Public mathematical background: H. W. Kuhn
  (1955), *The Hungarian Method for the Assignment Problem*, Naval Research
  Logistics Quarterly 2, 83-97. The public behavior and regression targets
  were checked with `scipy.optimize.linear_sum_assignment` in SciPy 1.17.0,
  but no SciPy or other library implementation was copied or adapted.
  This is not a port of SciPy's modified Jonker-Volgenant routine.
- `src/optimize/isotonic.f90`: original weighted pool-adjacent-violators
  implementation. Public algorithm reference: F. M. T. A. Busing (2022),
  *Monotone Regression: A Simple and Fast O(n) PAVA Implementation*,
  Journal of Statistical Software, Code Snippets, 102(1),
  doi:10.18637/jss.v102.c01. Compatibility inputs, results and tests were
  checked against `scipy.optimize.isotonic_regression` in SciPy 1.17.0.
  No source, coefficient table, or implementation was copied or adapted.
  Tests in `test/test_optimize_assignment.f90` and
  `test/test_optimize_isotonic.f90` are original SciFort regressions with
  independently generated SciPy reference values.

- `src/optimize/nonlin.f90` and `src/optimize/root.f90` (wall-clock
  adaptive Krylov refresh): original SciFort timing-mode extension using
  the standard Fortran `system_clock` intrinsic and existing SciFort
  refresh scheduling. No imported implementation or numerical data.
  The existing SciPy-derived BSD-3-Clause notices are preserved.
  `test/test_optimize_root_krylov_wall_clock_refresh.f90` and the
  additional OpenMP benchmark mode are original SciFort additions.


- `benchmarks/benchmark_krylov_quadrature_refresh.f90`: original SciFort
  numerical benchmark. A nonlinear rank-two response adds a per-coordinate
  exponential integral approximated with a fixed midpoint quadrature rule.
  Analytic forward and transpose Jacobian products differentiate that same
  discretization; optional OpenMP directives parallelize coordinate loops.
  No third-party implementation or numerical data were imported. The CSV
  records measurements from GNU Fortran 14.2 `-O2` runs at one and four
  threads; values are not portable performance guarantees.

- `benchmarks/benchmark_krylov_timed_refresh.f90`: original SciFort
  benchmarking harness comparing existing fixed, measured, and CPU-timed
  scheduling on an analytic nonlinear rank-two system. Its synthetic CPU
  workload and validation cases were written for SciFort; no third-party
  source, algorithms, coefficients, or measurements were imported.
  The CSV records locally measured GNU Fortran 14.2 `-O2` results.

- `src/optimize/nonlin.f90` and `src/optimize/root.f90` (timed
  adaptive Krylov refresh): original SciFort optional scheduling using the
  standard Fortran `cpu_time` intrinsic and a minimum measurement window.
  No third-party code or constants were imported. The SciPy-derived
  nonlinear-solver license notices remain intact.
  `test/test_optimize_root_krylov_timed_refresh.f90` is an original test.

- `src/optimize/nonlin.f90` (measured work for every automatic Krylov
  preconditioner): original SciFort optional per-setup counter propagation,
  with no imported third-party code or constants. Counts were added to the
  existing SciPy-derived root solver without altering its BSD-3-Clause
  attribution. The expanded
  `test/test_optimize_root_krylov_measured_rebuild.f90` and new
  `test/test_optimize_root_krylov_measured_setup_modes.f90` are original
  SciFort regression tests; no external implementation was adapted.

- `src/optimize/nonlin.f90` and `src/optimize/root.f90` (measured cost
  of nonsymmetric rank-two Krylov preconditioner refresh): original
  SciFort callback accounting and adaptive-refresh cost policy; no
  third-party implementation, data, constants, or algorithms imported.
  `test/test_optimize_root_krylov_measured_rebuild.f90` and the expanded
  `benchmarks/benchmark_krylov_work_refresh.f90` are original SciFort
  cases. The underlying SciPy-derived solver notices are preserved.


- `src/optimize/nonlin.f90` and `src/optimize/root.f90` (work-aware
  adaptive Krylov refresh): original SciFort control-flow refinement
  using the number of inner preconditioner applications since the latest
  successful preconditioner build, a user-selected growth factor and an
  estimated rebuild cost. No external code, constants or coefficient
  tables were imported. The existing SciPy-derived BSD notices remain
  unchanged. `test/test_optimize_root_krylov_work_adaptive.f90` is an
  original deterministic regression. These keywords are SciFort
  extensions, absent from SciPy.

- `benchmarks/benchmark_krylov_nonlinear_refresh.f90` and the extended
  `test/test_optimize_root_krylov_rank_two_nonlinear_base.f90`: original
  SciFort test and benchmark constructions using a cubic diagonal term
  and two nonlinear nonsymmetric coupling factors with elementary analytic
  derivatives. No third-party numerical code, data or coefficients were
  imported; the existing solver code and licensing are unchanged.

- `src/optimize/nonlin.f90` and `src/optimize/root.f90` (nonlinear
  nonsymmetric rank-two base-diagonal callback): original SciFort-only
  extension to the existing low-product rank-two fit. It evaluates a
  user-supplied `D(x)` at the accepted nonlinear iterate, reusing the
  established fit, inverse guards, refresh policy and exact-diagonal
  fallback. No numerical implementation code or coefficients were imported.
  The existing SciPy BSD-3-Clause attribution remains unchanged. Original
  regression: `test/test_optimize_root_krylov_rank_two_nonlinear_base.f90`.
  The new callback is a SciFort extension, not part of SciPy.

- `src/optimize/nonlin.f90` and `src/optimize/root.f90` (cost-adaptive
  transpose verification in nonsymmetric rank-two Krylov discovery):
  original SciFort extension of the existing deterministic held-out
  adjoint checks. The check count is selected from a user-supplied
  relative product cost, a normalized 2-by-2 pivot determinant and
  forward-fit discrepancy. No external numerical code or constants were
  imported; the preceding SciPy adaptation notices are preserved.
  Added `test/test_optimize_root_krylov_rank_two_transpose_adaptive.f90`
  and updated `benchmarks/benchmark_krylov_rank_two_discovery.f90`.
  This keyword is a SciFort extension, not part of SciPy.


- `src/optimize/nonlin.f90` and `src/optimize/root.f90` (optional
  transpose verification for nonsymmetric rank-two Krylov discovery):
  original SciFort-only refinement of its existing elementary two-sided
  cross approximation. Two deterministic vectors are orthogonalized
  against the left fitting subspace before evaluating held-out transpose
  Jacobian products. The resulting fit check, guarded fallback and API
  wiring are original changes, not copies of third-party code. The
  preexisting SciPy adaptation attribution is retained. Regression:
  `test/test_optimize_root_krylov_rank_two_transpose_verify.f90`.
  Benchmark: `benchmarks/benchmark_krylov_rank_two_discovery.f90`.


- `src/optimize/nonlin.f90` and `src/optimize/root.f90` (six-product
  nonsymmetric rank-two discovery): original SciFort extension using the
  elementary two-sided rank-two cross approximation
  `L = (L*Q)*(P^T*L*Q)^(-1)*(P^T*L)` for `L=J-D`, with known diagonal
  `D` and analytic transpose-Jacobian-vector actions. Two held-out
  forward products check the resulting factors before the existing
  guarded Woodbury inverse is formed. Deterministic directions, safety
  checks and fallback are original; no third-party source or numerical
  constants were imported. Existing BSD-3-Clause attribution for the
  SciPy-adapted solver remains. Tests:
  `test/test_optimize_root_krylov_rank_two_nonsymmetric.f90` and
  `test/test_optimize_root_krylov_rank_two_nonsymmetric_large.f90`.
  The keywords are SciFort extensions, not SciPy APIs.

- `src/optimize/nonlin.f90` and `src/optimize/root.f90` (four-product
  symmetric rank-two discovery): independently implemented SciFort-only
  extension based on the elementary symmetric low-rank interpolation
  identity `L = Y*(S^T Y)^(-1)*Y^T`, with `Y = (J-D)*S` and known
  uncoupled diagonal `D`. The deterministic fit and held-out probes,
  numerical guards, and exact-diagonal fallback are original code;
  no third-party implementation or coefficients were imported.
  The existing BSD-3-Clause attribution for the SciPy-adapted solver
  is retained. Regression sources:
  `test/test_optimize_root_krylov_rank_two_symmetric.f90` and
  `test/test_optimize_root_krylov_rank_two_symmetric_large.f90`.
  The public keyword is a SciFort extension, absent from SciPy.

- `src/optimize/nonlin.f90` and `src/optimize/root.f90` (streamed rank-two
  discovery): original SciFort implementation of an O(n)-storage variant
  of the existing unknown-factor `D+U*V^T` decomposition. Eight or fewer
  deterministic coordinate pilot columns identify a disjoint 2-by-2 minor;
  subsequent coordinate Jv products populate the measured diagonal,
  reconstruct the remaining low-rank factors and verify each off-diagonal
  element in streaming order. This is based on elementary 2-by-2 linear
  algebra and the Woodbury matrix identity, not imported source code.
  The unchanged SciPy-adapted solver license notices remain in force.
  New tests: `test/test_optimize_root_krylov_rank_two_stream.f90` and
  `test/test_optimize_root_krylov_rank_two_stream_large.f90`. This keyword
  is a SciFort extension and not part of SciPy 1.17.0.

- `src/optimize/nonlin.f90` and `src/optimize/root.f90` (opt-in
  rank-two diagonal sketch): original SciFort implementation of rowwise
  3-by-3 linear identification for the elementary matrix model
  `J = D + U*V^T` with known `V`. Twelve deterministic bounded sign
  vectors are searched for the best worst-row geometry; the three fitted
  products are checked by the preexisting rank-two Woodbury estimator.
  Failed identification retries the established exact-diagonal path.
  `test/test_optimize_root_krylov_rank_two_sketch.f90` and
  `benchmarks/benchmark_krylov_woodbury.f90` are original validation and
  benchmarking code. No third-party source, tables, or constants were
  copied; the existing SciPy BSD attribution remains in the module header.
  The option is a SciFort extension and not part of SciPy 1.17.0.

- `src/optimize/nonlin.f90` and `src/optimize/root.f90`: original SciFort
  automatic Jacobian-diagonal approximation and update scheduling, built
  on SciFort's existing analytic or finite-difference matrix-free Jv routine.
  Zero-entry fallback, relative clipping, and MINRES positive scaling are
  original safeguards. The associated test and benchmark
  (`test/test_optimize_root_krylov_auto_diagonal.f90`,
  `benchmarks/benchmark_krylov_auto_diagonal.f90`) are original, as is
  `benchmarks/krylov_auto_diagonal_gfortran14_o2.csv`. No external source
  implementation was imported or copied. The API is a SciFort extension,
  not a SciPy 1.17.0 keyword.

- `src/optimize/nonlin.f90` and `src/optimize/root.f90`: original SciFort
  convenience interface `inner_m_diag=`, dispatching a user-provided static
  diagonal inverse-Jacobian approximation into the existing Krylov
  preconditioner interface. Input and overflow validation are original code;
  no third-party numerical implementation was copied. The benchmark and
  regression under `benchmarks/benchmark_krylov_diagonal_preconditioner.f90`
  and `test/test_optimize_root_krylov_diagonal_preconditioner.f90` are also
  original SciFort validation code.


- `benchmarks/benchmark_krylov_restart_scaling.f90` and
  `test/test_optimize_root_krylov_restart_scaling.f90`: original SciFort
  validation code using a constructed `D*Q` nonlinear residual and exact
  analytic Jacobian-vector products. The benchmark exercises already
  implemented LGMRES restart and `tol_norm` options without importing or
  adapting any third-party implementation or changing numerical kernels.

- `src/optimize/nonlin.f90` and `src/optimize/root.f90`: original SciFort
  adaptive restart-width policy for the existing LGMRES iteration. The
  deterministic 20% weak-progress threshold and bounded doubling heuristic
  use the true linear residual already calculated between restarts. This is
  a SciFort extension, not translated SciPy implementation code; no new
  third-party code was copied. Regression:
  `test/test_optimize_root_krylov_adaptive_restart.f90`.

- `src/optimize/nonlin.f90` extends the original SciFort restarted Arnoldi
  control to LGMRES with retained augmentation vectors. It reuses the existing
  flexible Arnoldi implementation and computes true residuals between cycles;
  no third-party code was copied. Public forwarding is unchanged in
  `src/optimize/root.f90`. Tests:
  `test/test_optimize_root_krylov_lgmres_restart.f90` and
  `test/test_optimize_root_krylov_gmres_restart.f90`.

- `src/optimize/nonlin.f90` and `src/optimize/root.f90`: original SciFort
  optional restarted GMRES control. The restart-cycle logic reuses the
  already implemented one-cycle flexible Arnoldi solver, recomputes the true
  residual between cycles, and limits the total number of Arnoldi steps.
  Restarted GMRES is a standard numerical method; no code was copied from a
  third-party implementation. Tests: `test/test_optimize_root_krylov_gmres_restart.f90`.


- `src/optimize/nonlin.f90` and `src/optimize/root.f90`: original SciFort
  memory-budget dispatch for analytic candidate-image caching. This
  extension chooses between preserving the exact preorthogonalization
  Arnoldi products and evaluating an exact `J(x)*candidate` callback when
  the allocated raw-image capacity would exceed a user-provided bound.
  No third-party source code was imported. Regression:
  `test/test_optimize_root_krylov_candidate_cache_budget.f90`.

- `src/optimize/krylov_augmentation.f90` and `src/optimize/nonlin.f90`:
  original SciFort memory optimization for LGMRES Jacobian-image pruning.
  Reuses the existing image array to store the incremental orthonormal
  basis when original images are not required by Arnoldi, with a guarded
  reconstruction fallback. The direct regression is
  `test/test_optimize_krylov_image_basis_storage.f90`. No third-party
  implementation was copied or introduced by this change.


- `src/optimize/krylov_augmentation.f90`, `src/optimize/nonlin.f90`, and
  `src/optimize/root.f90`: original SciFort Fortran implementation of
  optional within-iteration reuse of the current-Jacobian LGMRES image
  pruning basis. Two-pass modified Gram-Schmidt was already used for image
  filtering; this change preserves its orthonormal basis for a single
  candidate insertion at the same nonlinear iterate. The basis is rebuilt
  whenever the Jacobian point changes. No third-party source was copied.
  Coverage: `test/test_optimize_krylov_augmentation_candidate_basis.f90`,
  `test/test_optimize_root_krylov_image_basis_reuse.f90`, and the extended
  changing-Jacobian pruning comparison.


- `src/optimize/krylov_augmentation.f90`, `src/optimize/nonlin.f90`, and
  `src/optimize/root.f90`: original SciFort Fortran implementation of an
  optional LGMRES full-history QR update for solution-direction pruning.
  It uses conventional Givens rotations to triangularize the Hessenberg
  factor after removing the oldest column, then appends the next direction
  by two-pass modified Gram-Schmidt. The mathematical method is standard QR
  updating/downdating as described in Golub and Van Loan, *Matrix
  Computations*, 4th edition (Johns Hopkins University Press, 2013).
  No implementation code from that book or any third-party library was
  copied. The default solver and third-party notices are unchanged. Coverage:
  `test/test_optimize_krylov_augmentation_fifo_qr.f90` and
  `test/test_optimize_root_krylov_fifo_qr.f90`.


- `src/optimize/elementwise_batch.f90`, `src/optimize/root_scalar.f90`,
  `src/optimize/minimize_scalar.f90`, and `src/optimize.f90`: original SciFort
  Fortran 2018 heterogeneous objective argument support. Polymorphic scalar
  boxes preserve exact intrinsic/user-derived dynamic types and shapes;
  independent trailing-axis maps and active-lane copying feed the four
  existing numeric controllers without algorithm changes. GNU Fortran 14.2
  required scalar boxing instead of directly reshaping `class(*)` arrays.
  Coverage: `test/test_optimize_elementwise_typed_args.f90`. No third-party
  numerical source was copied; original SciFort implementation.


- `src/optimize/elementwise_batch.f90`, `src/optimize/root_scalar.f90`,
  `src/optimize/minimize_scalar.f90`, and `src/optimize.f90`: original SciFort
  Fortran 2018 independently-shaped real parameter lists. The constructor
  snapshots any-rank real arguments, and a centralized trailing-axis
  broadcast path assembles `(nargs,nlanes)` input for the existing packed
  callback/active-lane solvers. Generic overloads preserve the existing
  stacked argument form. No external numerical implementation is imported.
  Coverage: `test/test_optimize_elementwise_vectorized_arg_list.f90`.
  MIT-licensed original SciFort plumbing.

- `src/optimize/elementwise_batch.f90`, `src/optimize/root_scalar.f90`,
  `src/optimize/minimize_scalar.f90`, and `src/optimize.f90`: original SciFort
  assumed-rank objective-parameter broadcasting and packed callback plumbing.
  The new parameter stack has a leading `nargs` axis, trailing NumPy-style
  singleton broadcasting, and per-evaluation-column gathering by active lane
  ID. All four existing packed numerical controllers are reused and no
  numerical algorithms or external source were copied. New tests in
  `test/test_optimize_elementwise_vectorized_args.f90` cover the four solvers,
  rank 0/5/14/15, full and trailing-axis broadcasting, invalid shapes,
  noncontiguous sections, statuses, and work counts. This extension is
  MIT-licensed original SciFort code.

- `src/optimize/elementwise_batch.f90`, `root_scalar.f90`, and
  `minimize_scalar.f90`: original SciFort Fortran 2018 assumed-rank
  vectorized dispatchers (`*_vectorized_into`) that flatten independent
  grids, broadcast optional numerical arguments, and restore typed result
  shapes without copying external source. They reuse the existing packed
  active-lane controllers and preserve the older function overloads.
  `test/test_optimize_elementwise_vectorized_assumed_rank.f90` covers rank 0,
  rank 5, maximum rank 14/15, invalid shape, broadcasting and noncontiguous
  arrays. GNU Fortran 14.2 requires `intent(inout)` rather than `intent(out)`
  for default-initialized derived-type assumed-rank result dummies; each
  dispatcher explicitly resets results on entry.

- `src/optimize/root_scalar.f90` and `src/optimize/minimize_scalar.f90`:
  original SciFort rank-three and rank-four wrappers for the already attributed
  packed vector-callback root and minimum solvers. These adapters use the
  existing SciFort `elementwise_broadcast_optional` helper and rank-one
  controllers, without importing external implementation code. Test coverage
  in `test/test_optimize_elementwise_vectorized_rank34.f90` compares direct
  packed rank-one results with three- and four-dimensional grids, including
  per-lane result values, statuses, work counts, optional broadcast inputs,
  callbacks, invalid shapes, and nonfinite-input isolation.

- `src/optimize/elementwise_batch.f90`, `src/optimize/root_scalar.f90`, and
  `src/optimize/minimize_scalar.f90`: original SciFort rank-two adapters for
  the existing MIT/BSD-3-Clause-attributed packed vector-callback controllers.
  The wrappers flatten rank-two grids in Fortran order, use independent
  SciFort trailing-axis broadcasting, and reshape status/results without
  changing the original solver algorithms. No external source was imported.
  `test/test_optimize_elementwise_vectorized_rank2.f90` compares all four
  solvers to their rank-one counterparts, with bad arguments, callback
  termination, empty grids, and checked work counts.

- `src/optimize/root_scalar.f90` `bracket_root_vectorized` and
  `src/optimize/minimize_scalar.f90` `bracket_minimum_vectorized`: original
  SciFort packed active-lane orchestration, reusing the already documented
  scalar bracketing algorithms adapted from SciPy 1.17 `_bracket` (BSD-3-Clause).
  Vector callback batching, per-lane state, and bounds broadcasting are
  independent Fortran implementations. Numeric statuses/iterations/evaluations
  are checked against SciPy 1.17 and the existing scalar SciFort procedures.
  Test sources and reference generation are
  `test/test_optimize_elementwise_vectorized_bracket_root.f90`,
  `test/test_optimize_elementwise_vectorized_bracket_minimum.f90`, and
  `tools/generate_optimize_elementwise_vectorized_brackets_reference.py`.


- `src/optimize/minimize_scalar.f90` `find_minimum_vectorized`: original SciFort
  packed active-lane controller based on the already attributed SciPy 1.17
  Chandrupatla scalar minimization formulas. It evaluates initial endpoint
  triples in one vector callback and then only active trial points, keeping
  independent states, scalar stop callbacks, and per-point `nfev`. Validation
  compares against SciPy 1.17 and scalar SciFort in
  `test/test_optimize_elementwise_vectorized_minimum.f90` with reference
  constants generated by
  `tools/generate_optimize_elementwise_vectorized_minimum_reference.py`.
  No SciPy array-controller implementation was imported.

- `src/optimize/root_scalar.f90` `find_root_vectorized`: original SciFort
  vector-callback active-set controller, reusing the already attributed
  SciPy 1.17 Chandrupatla scalar formulas and the shared SciFort broadcaster.
  Packs valid endpoint samples and pending iteration points for one callback
  per active iteration; statuses, iteration/evaluation counts, and scalar
  refinement callbacks remain independent. The numerical behavior is checked
  against `scipy.optimize.elementwise.find_root` by
  `tools/generate_optimize_elementwise_vectorized_root_reference.py` and
  `test/test_optimize_elementwise_vectorized_root.f90` with generated
  `test/optimize_elementwise_vectorized_root_reference.f90`. No SciPy array
  controller source was imported. The separate vectorized minimum refinement
  routine uses the same packed-callback design.

- `src/optimize/elementwise_batch.f90`,
  `src/optimize/root_scalar.f90`, and `src/optimize/minimize_scalar.f90`:
  original SciFort Fortran 2018 assumed-rank broadcasting and flattened
  independent-problem controllers. Public generic interfaces return correctly
  shaped results through rank 15 for `bracket_root`/`bracket_minimum` and
  through rank 14 for `find_root`/`find_minimum` (one leading endpoint
  dimension). All batch ranks reuse the existing SciPy-attributed scalar
  algorithms. The new broadcasting code is independently written, without
  translating SciPy's array controller. The tests include the previously
  generated SciPy 1.17 references for ranks 1-3 and
  `test/test_optimize_elementwise_highrank.f90` for ranks 1-15. The latter
  verifies scalar parity, independent evaluation counts, invalid shapes,
  singleton expansion, callbacks, and empty grids. Real function and callback
  evaluations remain scalar in these generic overloads; a separate
  rank-one `find_root_vectorized` and `find_minimum_vectorized` interfaces
  now support vector evaluation.


- `src/optimize/root_scalar.f90` and `src/optimize/minimize_scalar.f90`:
  original SciFort rank-three independent-problem wrappers. They reuse the
  existing attributed scalar SciPy 1.17 `elementwise` algorithms, rather than
  importing SciPy's array controller. NumPy-style trailing-axis argument
  broadcasting, scalar callback ordering, per-element status/work counts,
  shape validation, and empty grids are checked with the independently generated
  `tools/generate_optimize_elementwise_rank3_reference.py`,
  `test/optimize_elementwise_rank3_reference.f90`, and
  `test/test_optimize_elementwise_rank3.f90`. A separate validation-only change
  makes scalar `find_root` reject NaN tolerances before any sign comparisons,
  avoiding a GNU Fortran 14.2 `-O3` IEEE-invalid trap.

- `src/optimize/root_scalar.f90` and `src/optimize/minimize_scalar.f90`:
  original SciFort rank-two independent-problem wrapper extensions using
  the existing attributed scalar elementwise algorithms. The NumPy-style
  trailing-axis and singleton-dimension broadcast conventions, result shape,
  per-problem statuses and evaluation counts were checked against SciPy
  1.17.0 (BSD-3-Clause) by
  `tools/generate_optimize_elementwise_rank2_reference.py` and
  `test/test_optimize_elementwise_rank2.f90`. No SciPy vectorized controller
  implementation was imported; callbacks remain scalar. The wrappers also
  normalize temporary-array lower bounds before indexed broadcasting.

- `src/optimize/root_scalar.f90` and `src/optimize/minimize_scalar.f90`:
  original SciFort rank-one independent-problem wrappers for the four SciPy
  `elementwise` bracketing and refinement functions. These forward scalar or
  rank-one broadcast arguments to the pre-existing attributed scalar kernels,
  without copying SciPy's vectorized controller. The API and output/status/work
  counts were checked against SciPy 1.17.0 (BSD-3-Clause) using
  `tools/generate_optimize_elementwise_arrays_reference.py`,
  `test/optimize_elementwise_arrays_reference.f90`, and
  `test/test_optimize_elementwise_arrays.f90`. The rank-one through rank-three wrappers now share the generalized batch
  engine. Vectorized function callbacks are not yet translated.

Original SciFort implementation code is MIT licensed. The unreleased additions
also contain the explicitly identified SciPy BSD-3-Clause adaptations and Sobol
data below. Mixed-source files retain both MIT and BSD-3-Clause notices;
third-party notices and complete license texts are recorded in
`THIRD_PARTY_LICENSES.md` and `LICENSES/`. The Landau research coefficient-table
reference is separately attributed under CC BY 4.0.

- `src/optimize/minimize.f90` `updated_constraint_hessian`: original SciFort
  dense per-output BFGS/SR1 Hessian-update implementation. The quasi-Newton
  secant formulas, auto-scaling, Powell damping, and SR1 denominator threshold
  were adapted from standard Hessian-update strategies and compared with
  SciPy 1.17.0 `scipy.optimize._hessian_update_strategy` (BSD-3-Clause).
  Each scalar constraint Hessian is weighted with its *current* multiplier,
  avoiding Hessian changes caused solely by multiplier changes. This differs
  from SciPy's per-block weighted Hessian strategy, so iteration counts need
  not agree. `test/test_optimize_trust_constr_constraint_updates.f90` and
  `tools/generate_optimize_trust_constr_constraint_updates_reference.py`
  provide regression and reference cases. No external solver was imported.

- `src/optimize/minimize.f90` `trust_constr_iteration_callback` and its SQP
  dispatch/termination hooks: original SciFort Fortran procedure interface
  modeled after SciPy 1.17.0 `scipy.optimize._trustregion_constr` callback
  invocation and status-3 early termination conventions (SciPy BSD-3-Clause).
  SciFort forwards only physical `x` and cumulative iteration number; it does
  not copy SciPy's Python callback implementation or support full `state` yet.
  `test/test_optimize_trust_constr_callback.f90` covers both SQP and barrier
  forms, with expected behavior compared in
  `tools/generate_optimize_trust_constr_callback_reference.py`.

- `src/optimize/minimize.f90` `colored_jacobian_block`: original SciFort
  implementation of a greedy coloring of the column-interference graph
  induced by a binary constraint-Jacobian sparsity pattern. Columns without
  overlapping nonzero rows may be perturbed together, with each output row's
  finite difference attributed to its unique influencing column. The concept
  is standard sparse numerical differentiation; the API and public numerical
  reference behavior were checked against SciPy 1.17.0
  `scipy.optimize.NonlinearConstraint(finite_diff_jac_sparsity=...)`.
  `test/optimize_trust_constr_sparse_jac_reference.f90` is generated by
  `tools/generate_optimize_trust_constr_sparse_jac_reference.py`.
  No SciPy sparse-difference implementation was copied or translated, and
  this feature does not import SciPy's sparse projector/LSMR implementation.

- `src/optimize/minimize.f90` nonlinear-constraint Jacobian finite-difference
  relative steps: original SciFort Fortran implementation of the signed
  probe convention and bounds-aware adjustments described by SciPy 1.17.0
  `scipy/optimize/_numdiff.py` (`_compute_absolute_step`, `approx_derivative`),
  BSD-3-Clause. The reused numerical-step helper also serves objective
  Hessian/gradient approximations. Dense and colored constraint probes are
  verified against generated SciPy numerical references in
  `tools/generate_optimize_trust_constr_constraint_rel_step_reference.py`
  and `test/test_optimize_trust_constr_constraint_rel_step.f90`.
  Signed bounded **objective-gradient** steps use the same original SciFort
  `bounded_signed_fd_step` helper at barrier initialization and in the
  equality-SQP kernel. Probe directions and gradients are independently
  checked against SciPy 1.17 `scipy.optimize.minimize` in
  `tools/generate_optimize_trust_constr_objective_signed_bounds_reference.py`
  and `test/test_optimize_trust_constr_objective_signed_bounds.f90`.

- `src/optimize/root_scalar.f90` and `src/optimize.f90`: the bisection,
  `brentq`, `brenth`, and Ridders numerical kernels are adapted from
  roots-fortran 1.0.0 by Jacob Williams
  (<https://github.com/jacobwilliams/roots-fortran/tree/1.0.0>), BSD-3-Clause.
  SciFort restructures the object-oriented/functional upstream routines into a
  SciPy-style result API, adds explicit status values, and replaces
  multiplication-based sign tests with IEEE-trap-safe comparisons. The
  `root_scalar` dispatch/result semantics, secant/Newton/Halley iteration, and
  TOMS748 organization are adapted from SciPy 1.17.0 `scipy/optimize/_root_scalar.py`
  and `scipy/optimize/_zeros_py.py` under SciPy's BSD-3-Clause license. SciFort's
  numerical-derivative Newton path delegates to `scifort_differentiate` and
  accounts for each derivative invocation as one public function call to match
  SciPy's `root_scalar` reporting convention. The scalar `bracket_root`
  endpoint-expansion rule, finite-limit interpolation, status meanings, and
  work-count convention are adapted from SciPy 1.17.0
  `scipy/optimize/_bracket.py` (BSD-3-Clause). Array broadcasting is not
  translated; `tools/generate_optimize_bracket_root_reference.py` and
  `test/optimize_bracket_root_reference.f90` preserve SciPy numerical references.
  Complete upstream notices are in `THIRD_PARTY_LICENSES.md`.
- `src/optimize/root_scalar.f90` `find_root` and its public umbrella export:
  scalar adaptation of Chandrupatla's inverse-quadratic/bisection root solver
  from SciPy 1.17.0 `scipy/optimize/_chandrupatla.py`, specifically
  `_chandrupatla` (SciPy BSD-3-Clause). The iteration update, scaled
  function and interval tolerances, bracket ordering, callback semantics,
  and status/work-count conventions are adapted to Fortran. SciFort-specific
  validation and IEEE-trap-safe nonfinite checks are included. The original
  SciPy copyright and license notice is retained in `root_scalar.f90` and
  `THIRD_PARTY_LICENSES.md`. Generated numerical references from
  `tools/generate_optimize_find_root_reference.py` are included in
  `test/optimize_find_root_reference.f90` and tested by
  `test/test_optimize_find_root.f90`. The numerical kernel is scalar;
  rank-one convenience wrappers evaluate each element independently rather
  than translating SciPy's vectorized controller.
- `src/optimize/minimize_scalar.f90` `bracket_minimum` and its public export:
  scalar Fortran adaptation of the moving-triple downhill expansion, finite
  search-limit interpolation, result ordering, iteration accounting, and
  status conventions in SciPy 1.17.0 `scipy/optimize/_bracket.py`,
  `scipy.optimize.elementwise.bracket_minimum` (BSD-3-Clause). Numerical
  comparisons are preserved by
  `tools/generate_optimize_bracket_minimum_reference.py`,
  `test/optimize_bracket_minimum_reference.f90`, and
  `test/test_optimize_bracket_minimum.f90`. This does not translate SciPy's
  vectorized/broadcasting controller, and invalid inputs are rejected before
  objective callbacks as in SciFort's existing scalar bracketing API.
  SciPy's license and copyright notices are retained in the source and in
  `THIRD_PARTY_LICENSES.md`.
- `src/optimize/minimize_scalar.f90` `find_minimum` and the umbrella export:
  scalar adaptation of Chandrupatla's three-point quadratic-fit/sectioning
  method from SciPy 1.17.0 `scipy/optimize/_chandrupatla.py`, specifically
  `_chandrupatla_minimize` (SciPy BSD-3-Clause). The iteration update,
  wider-half-bracket selection, abscissa/function stopping tests, callback
  reporting, and results are adapted with an explicit Fortran scalar
  interface. SciFort-specific checks reject malformed settings before
  calling the objective and avoid unsafe comparisons of IEEE NaNs.
  `tools/generate_optimize_find_minimum_reference.py` regenerates the
  SciPy reference constants in `test/optimize_find_minimum_reference.f90`,
  exercised by `test/test_optimize_find_minimum.f90`. Rank-one wrappers reuse
  the scalar kernel; SciPy's array execution controller has not been
  translated. The original license is retained in
  `THIRD_PARTY_LICENSES.md`.
- `src/optimize/root.f90` and its exports in `src/optimize.f90`: original SciFort
  adapter code over the public `hybrd`, `hybrj`, `lmdif`, and `lmder` interfaces
  of the community-maintained modern MINPACK project
  (<https://github.com/fortran-lang/minpack>). The external dependency is not
  vendored. SciFort maps SciPy 1.17.0 `root(method='hybr')`, `root(method='lm')`,
  and `fsolve` defaults, result/status fields, callback conventions, and public
  function/Jacobian evaluation accounting onto MINPACK. Generated references in
  `test/optimize_root_reference.f90` are public numerical outputs from SciPy
  1.17.0. No MINPACK implementation source is copied or adapted into SciFort.
- `src/optimize/nonlin.f90`, the inexact-Newton additions in
  `src/optimize/root.f90`, and their exports in `src/optimize.f90`: the iteration
  controller, default tolerance logic, Jacobian autoscaling, inverse-Jacobian
  Broyden updates, Broyden rank reduction, Anderson mixing, diagonal Broyden,
  linear/exciting mixing, and matrix-free Krylov/Newton behavior are adapted from
  SciPy 1.17.0 `scipy/optimize/_nonlin.py` and the corresponding dispatch
  conventions in `scipy/optimize/_root.py`, under SciPy's BSD-3-Clause license.
  Krylov includes SciPy-style finite-difference Jacobian-vector products,
  nonlinear forcing, LGMRES recycling, and GMRES inner solves. The dense rank-one
  Broyden Jacobian update was additionally cross-checked against
  nlesolver-fortran 1.1.0 by Jacob Williams
  (<https://github.com/jacobwilliams/nlesolver-fortran/tree/1.1.0>),
  BSD-3-Clause. Ivan Pribec's BSD-3-Clause `kiss` project
  (<https://github.com/ivan-pi/kiss>) was consulted as a modern-Fortran
  GMRES/LGMRES API and implementation reference; no `kiss` source is copied and
  SciFort does not add it as a dependency. Generated references in
  `test/optimize_root_reference.f90` are public numerical outputs from SciPy
  1.17.0. The permanent `test/test_optimize_root_nonlin.f90` regression covers
  the dense methods, Krylov/GMRES/LGMRES options, Broyden rank reductions,
  line-search choices, generic tolerances, direct wrappers, and iteration-limit
  behavior. Complete copied/adapted upstream notices are in
  `THIRD_PARTY_LICENSES.md`. The opt-in `jacobian_diff_scheme='3-point'`
  centered directional difference, default cube-root-epsilon step, and
  representability/overflow guards are SciFort-authored numerical extensions
  implementing the standard symmetric finite-difference formula; they are not
  copied from SciPy or another third-party implementation. The existing
  SciPy-derived two-point default and license notices are retained.
- `src/optimize/wolfe.f90` and the Wolfe-search path in
  `src/optimize/nonlin.f90`: adapted from SciPy 1.17.0
  `scipy/optimize/_linesearch.py` and `scipy/optimize/_dcsrch.py` under SciPy's
  BSD-3-Clause license. SciPy's `_dcsrch.py` is itself a Python translation of
  the historical MINPACK-2 DCSRCH/DCSTEP strong-Wolfe line-search algorithm by
  Jorge J. More, David J. Thuente, Brett M. Averick, and Richard G. Carter.
  Jacob Williams' unofficial MINPACK-2 mirror
  (<https://github.com/jacobwilliams/MINPACK-2>) was consulted to cross-check the
  historical `csrch/dcsrch.f` and `csrch/dcstep.f` provenance; SciFort adapts the
  SciPy BSD-licensed translation rather than vendoring the mirror.
- `src/optimize/df_sane.f90` and the `df-sane` dispatch in
  `src/optimize/root.f90`: the spectral residual iteration, relative/absolute
  convergence rule, spectral-coefficient safeguards, default forcing sequence,
  and Cruz/Cheng nonmonotone line searches are adapted from SciPy 1.17.0
  `scipy/optimize/_spectral.py` and `scipy/optimize/_linesearch.py` under SciPy's
  BSD-3-Clause license. SciFort implements the same real-vector algorithm with
  explicit Fortran work arrays and result bookkeeping; SciPy's complex-vector
  packing and Python callback/printing machinery are not translated. Generated
  `df-sane` references in `test/optimize_root_reference.f90` are public numerical
  outputs from SciPy 1.17.0. Complete SciPy notices are in
  `THIRD_PARTY_LICENSES.md`.
- `src/optimize/minimize.f90`, its exports in `src/optimize.f90`, and the
  scalar-Wolfe-2/zoom additions in `src/optimize/wolfe.f90`: Nelder-Mead
  (including bound clipping/reflection), modified Powell (including finite and
  one-sided bounded line searches), BFGS, Polak-Ribiere+ conjugate gradient, Newton-CG inner
  conjugate-gradient iteration, finite-difference gradient/Hessian-product
  bookkeeping, strong-Wolfe fallback/zoom, the dogleg trust-region controller and
  Cauchy/Newton dog-leg subproblem, the trust-ncg trust-region
  controller/Steihaug CG subproblem (including negative-curvature boundary
  handling), the trust-exact nearly-exact damping iteration (including
  indefinite-Hessian hard cases), and the trust-krylov trust-region controller,
  hot-start behavior, and GLTR forcing sequence, method defaults, termination
  rules, and result/work-count semantics are adapted from SciPy 1.17.0
  `scipy/optimize/_optimize.py`, `_minimize.py`, `_linesearch.py`,
  `_trustregion.py`, `_trustregion_dogleg.py`, `_trustregion_ncg.py`,
  `_trustregion_exact.py`, and `_trustregion_krylov.py` under SciPy's
  BSD-3-Clause license. The native Fortran Krylov subproblem implementation also
  follows the GLTR state machine, reduced tridiagonal trust-region formulation,
  and forcing rules documented and implemented by Felix Lenders, Christian
  Kirches, and Andreas Potschka in TRLIB (MIT); SciFort stores an explicit
  orthonormal Krylov basis instead of copying TRLIB's vector-free C code. Powell reuses SciFort's separately attributed
  scalar Brent implementation. Generated references in
  `test/optimize_minimize_reference.f90` are public numerical outputs from SciPy
  1.17.0 and are regenerated by `tools/generate_optimize_minimize_reference.py`.
  Complete SciPy notices are in `THIRD_PARTY_LICENSES.md`.
- The independent `nonlinear_constraint` block interface in
  `src/optimize/minimize.f90` and `src/optimize.f90` is a new SciFort Fortran
  adapter to the separate `NonlinearConstraint` conventions in SciPy 1.17.0
  `scipy/optimize/_constraints.py` (BSD-3-Clause). Its bound arrays, callback
  accounting, and multiplier-weighted Hessian aggregation are mapped onto
  SciFort's existing combined raw constraint API and canonicalizer; no new
  external numerical implementation is imported. Finite-difference probes
  use the previously attributed SciFort bounded step helper. The new reference
  generator tests direct SciPy 1.17.0 optimization results for independently
  supplied equality/inequality blocks.
- The `trust_constr` path in `src/optimize/minimize.f90` follows the
  Byrd-Omojokun trust-region SQP structure and numerical constants in SciPy
  1.17.0 `scipy/optimize/_trustregion_constr/equality_constrained_sqp.py`, the
  projected-CG subproblem and stop codes in
  `scipy/optimize/_trustregion_constr/qp_subproblem.py`, the result-state and
  stopping/accounting conventions in
  `scipy/optimize/_trustregion_constr/minimize_trustregion_constr.py`, plus the
  positive-slack logarithmic-barrier sequence in
  `scipy/optimize/_trustregion_constr/tr_interior_point.py`, and the per-component
  lower/upper constraint canonicalization in
  `scipy/optimize/_trustregion_constr/canonical_constraint.py`. SciFort uses dense
  symmetric-eigendecomposition pseudoinverses for row/null-space projections and
  reuses its equality SQP kernel for each barrier subproblem; nonfixed SciFort
  `c(x) >= 0` inequalities become `c(x)-s=0` with positive slacks, and fixed
  bounds become equalities. The implementation independently preserves SciPy's
  default `0.1` initial barrier parameter/tolerance and `0.2` decay schedule,
  and follows SciPy's equality-SQP merit function initialization/update for
  `initial_constr_penalty`. Analytic constraint Hessians follow SciPy's
  `NonlinearConstraint.hess(x, v)` weighted-Hessian convention and are added to
  the analytic objective Hessian to form the Lagrangian Hessian. For barrier
  subproblems, SciFort maps equality-form multipliers back to the original
  nonlinear equality/inequality constraint vector; linear bound/slack rows add
  no curvature. The pure equality-SQP path also follows SciPy's second-order
  correction test and minimum-norm correction based on the current Jacobian
  pseudoinverse. Nonlinear-inequality barrier subproblems now enable the same SOC
  test and clip the correction against SciPy's `0.995` fractional slack boundary
  before evaluating the corrected point. Barrier normal/tangential steps now use
  SciPy-style scaled slack coordinates and the fractional slack trust box, with
  trust-radius carry between barrier stages; bound-only barrier SOC is enabled.
  Inner SQP multipliers are retained for original-problem optimality checks, and
  missing nonlinear-constraint Hessians use a separate BFGS approximation when
  all nonlinear blocks omit exact Hessian callbacks and no numerical Hessian
  scheme is requested. The overdetermined-equality
  extension reuses the existing SciFort dense Gram-eigendecomposition
  pseudoinverse, including its minimum-norm constraint multipliers and
  null-space projections. It removes dimension-only rejection of redundant
  equality rows without importing a new numerical implementation. Tests compare
  with SciPy 1.17.0 `factorization_method='SVDFactorization'` results, including
  an inconsistent system that must terminate with status 4; the SciFort
  projector is not an implementation of SciPy's direct SVD algorithm. The
  reference data are generated by
  `tools/generate_optimize_trust_constr_redundant_reference.py`. In a mixed collection,
  SciFort retains exact per-block Hessians and reconstructs each missing one
  from finite differences of its own Jacobian, using the existing bounded
  numerical-derivative helper. This mixed-case dense approximation is a
  SciFort extension, not SciPy's Hessian-update strategy; no additional
  third-party implementation is imported. The opt-in per-block
  `%hess_fd_scheme='2-point'` / `'3-point'` follows SciPy 1.17.0
  `NonlinearConstraint.hess` numerical-differentiation conventions, using
  the existing bounded finite-difference helpers to reconstruct the full
  weighted Hessian from an analytic Jacobian. Unlike SciPy's directional
  LinearOperator implementation, SciFort symmetrizes the dense matrix for
  its SQP kernel; the scheme is not an exact copy of SciPy's operator path.
  The reference values come from
  `tools/generate_optimize_trust_constr_constraint_hess_fd_reference.py`.
  The optional SR1 objective Hessian
  update follows SciPy 1.17.0 `scipy/optimize/_hessian_update_strategy.py`
  (`SR1` and `FullHessianUpdateStrategy._auto_scale`, BSD-3-Clause), including
  its `1.0e-8` denominator skip rule, auto-scaling, and potentially indefinite
  rank-one correction. SciFort keeps this objective update separate from
  supplied exact or approximated constraint curvature. SciPy 1.17.0 reference
  values are generated by `tools/generate_optimize_trust_constr_sr1_reference.py`
  and checked in `test/test_optimize_trust_constr_sr1.f90`.
  The optional numerical objective-Hessian mode differentiates an exact
  objective-gradient callback using SciPy 1.17.0's
  `scipy/optimize/_numdiff.py` finite-difference step conventions
  (BSD-3-Clause), with SciFort's existing dense symmetric assembly and
  bounded stencil helpers. No new third-party source is copied or
  vendored. It uses `njev` rather than `nhev` to count derivative probes.
  Generated reference values are in
  `test/optimize_trust_constr_hess_fd_reference.f90` and their generator is
  `tools/generate_optimize_trust_constr_hess_fd_reference.py`.
  Barrier-stage restarts reuse cached
  objective, gradient, and analytic objective-Hessian values at the unchanged
  point, following SciPy's `ScalarFunction`/barrier initial-value handoff, and SQP
  iteration counts include the initial stopping-state update used by SciPy. The
  representative active-inequality Maratos regression therefore matches SciPy
  1.17.0 iteration and objective/gradient/Hessian callback counts exactly. The
  equality and barrier paths also report the full Lagrangian gradient, cumulative
  projected-CG iteration count, and final projected-CG stop condition, and apply
  SciPy's post-solve status `4` when a nominally successful termination remains
  infeasible. Constraint violation, optimality, barrier parameter, and barrier
  tolerance diagnostics are also reported. Final constraint values and
  Jacobians are returned in a flattened/row-stacked Fortran representation;
  evaluation-count arrays use the current SciFort API's combined user-constraint
  callback as one block and append bounds as a second zero-cost block. Optional
  `constraint_lower`/`constraint_upper` vectors canonicalize raw callback
  components into equality and nonnegative inequality residuals while preserving
  raw final values/Jacobians in that first result block; analytic constraint
  Hessian multipliers are mapped back through the same signs before the user
  callback is invoked. Bounds remain a later result block, following SciPy's
  bounds-after-user-constraints order. Cached constraint values and
  Jacobians are reused across unchanged barrier restarts and outer stopping
  checks, and analytic constraint Hessians include SciPy's initial
  zero-multiplier setup call so covered callback counts agree exactly. Following
  SciPy 1.17 `scipy/optimize/_constraints.py`,
  `canonical_constraint.py`, and `tr_interior_point.py`, the current bound-vector
  API accepts per-component `constraint_keep_feasible` and per-variable
  `bounds_keep_feasible` masks. Equality components are unaffected, infeasible
  retained inequalities are rejected at the initial point, interval keep flags are
  duplicated onto both canonical inequality sides, and retained barrier slacks are
  reset to the current constraint gap before a trial objective is evaluated. The
  strict variable-bound finite-difference behavior follows SciPy 1.17
  `scipy/optimize/_constraints.py` (`strict_bounds`) together with the one-sided
  bound adjustment in `scipy/optimize/_numdiff.py`: numerical objective-gradient
  and nonlinear-constraint Jacobian probes reverse or shorten a forward step rather
  than crossing a retained variable bound. The implementation does not copy
  SciPy's sparse projection machinery or LinearOperator plumbing, and work counts
  for cases outside the covered regressions can still differ. Regression targets
  in `test/optimize_trust_constr_reference.f90` are generated from SciPy 1.17.0 by
  `tools/generate_optimize_trust_constr_reference.py`.
- `src/optimize/cobyla_prima.f90`, the public `cobyla` wrapper in
  `src/optimize/minimize.f90`, and its export in `src/optimize.f90`: the COBYLA
  trust-region/simplex algorithm, linear-programming trust-region subproblem,
  interpolation-simplex geometry updates, merit-penalty updates, history
  selection, and termination behavior are independently adapted to modern
  Fortran primarily from SciPy 1.17.0 `scipy/_lib/pyprima/cobyla`, which is a
  faithful Python translation of the documented modern-Fortran PRIMA reference
  implementation (<https://github.com/libprima/prima>). PRIMA is by Zaikun Zhang
  and contributors and is distributed under BSD-3-Clause. SciPy's pyprima COBYLA
  translation, credited upstream to Nickolai Belakovski, is also BSD-3-Clause.
  The SciFort wrapper independently maps its
  equality-first/nonnegative-inequality callback into PRIMA's `c(x) <= 0`
  convention, converts equalities to paired inequalities, adds finite bounds,
  independently eliminates fixed coordinates before simplex construction and
  restores the full vector for objective and constraint callbacks,
  mirrors SciPy's bounds-only starting-point projection, and maps `catol`,
  `f_target`, evaluation limits, status, and `maxcv` into
  `multivariate_minimize_result`. Regression targets in
  `test/optimize_cobyla_reference.f90` are generated from SciPy 1.17.0 by
  `tools/generate_optimize_cobyla_reference.py`. Complete PRIMA and SciPy notices
  are in `THIRD_PARTY_LICENSES.md`.
- The native TNC path in `src/optimize/minimize.f90` and its public export in
  `src/optimize.f90` are an independent modern-Fortran implementation of the
  truncated-Newton/active-bound algorithm described by Stephen G. Nash. SciPy
  1.17.0 `scipy/optimize/_tnc.py` and `scipy/optimize/tnc/tnc.c` were used as
  references for tolerance/scaling conventions, projected-gradient stopping,
  return-status meanings, and bound handling. The inner Newton
  system is solved with finite-difference Hessian-vector products and BFGS-preconditioned
  truncated CG. The full BFGS inverse approximation plays the same role as the legacy
  TNC diagonal/two-step self-scaled BFGS preconditioner without copying that C
  implementation. When no analytic gradient is supplied, SciFort uses bounds-aware
  forward differences and the BFGS approximation directly to avoid nested numerical
  differentiation. Regression targets are generated from SciPy 1.17.0 by
  `tools/generate_optimize_minimize_constrained_reference.py`.

- The L-BFGS-B and SLSQP adapters in `src/optimize/minimize.f90`, their public
  exports in `src/optimize.f90`, and the dependency pins in `fpm.toml` are
  original SciFort interface code around established external numerical kernels.
  L-BFGS-B uses Jacob Williams' `lbfgsb` release 1.0.0
  (<https://github.com/jacobwilliams/lbfgsb/tree/1.0.0>), a BSD-3 modern-Fortran
  refactor of L-BFGS-B 3.0; SciFort calls its public reverse-communication
  `setulb` routine and does not vendor the implementation. SLSQP uses Jacob
  Williams' `slsqp` release 1.6.0
  (<https://github.com/jacobwilliams/slsqp/tree/1.6.0>), BSD-3-Clause plus the
  retained Dieter Kraft/ACM original notice; SciFort calls the public
  `slsqp_solver` class and does not vendor its numerical kernel. Bounds clipping,
  bounds-aware finite-difference steps, SciPy-style method selection/tolerances,
  result/status mapping, and objective/gradient evaluation accounting are adapted
  from SciPy 1.17.0 `scipy/optimize/_lbfgsb_py.py`, `_slsqp_py.py`, and
  `_minimize.py` under SciPy's BSD-3-Clause license. Generated numerical outputs
  in `test/optimize_minimize_constrained_reference.f90` come from SciPy 1.17.0
  and are regenerated by
  `tools/generate_optimize_minimize_constrained_reference.py`. Upstream license
  locations are recorded in `THIRD_PARTY_LICENSES.md`.
- `src/optimize/least_squares.f90`, its exports in `src/optimize.f90`,
  `test/test_optimize_least_squares*.f90`, and
  `tools/generate_optimize_least_squares_reference.py`: LM/leastsq use the
  existing community-maintained `fortran-lang/minpack` dependency through
  `lmder`/`lmdif`; the dense TRF and dogbox controllers, strict-feasibility and
  Coleman-Li scaling rules, reflected/dogleg step selection, trust-radius and
  termination logic, per-iteration callback timing/early-stop status, Jacobian-column
  `x_scale='jac'` updates, two-point/three-point/complex-step finite differences,
  `jac_sparsity` structural masking and grouped-column differencing, built-in and callable
  robust loss transformations, curve-fit dispatch and curve-fit numerical
  model-Jacobian selection,
  and result semantics are adapted from SciPy 1.17.0
  `scipy/optimize/_lsq/trf.py`, `_lsq/dogbox.py`, `_lsq/common.py`,
  `_lsq/least_squares.py`, `_numdiff.py`, and `_minpack_py.py` under SciPy's
  BSD-3-Clause
  license. Sparse-difference grouping follows the Curtis-Powell-Reid nonconflicting-column criterion; SciFort uses deterministic natural-column greedy order rather than SciPy's seeded permutation. SciFort replaces SciPy's dense SVD calls with the project's existing
  symmetric eigensolver applied to the normal-equation model and implements
  Fortran callback/result storage independently. `src/optimize/lsmr.f90` adapts
  the Fong-Saunders LSMR recurrence and stopping tests used by SciPy 1.17.0
  `scipy/sparse/linalg/_isolve/lsmr.py`; its matrix-free callback structure was
  cross-checked against Jacob Williams's modern Fortran LSMR project. The TRF
  LSMR two-dimensional subspace and regularization logic, including the
  LinearOperator-style `J*v`/`J**T*v` path, follows SciPy 1.17.0 `_lsq/trf.py`
  and `_lsq/common.py`. The dogbox scaled LSMR Gauss-Newton step and active-variable
  masking follow SciPy 1.17.0 `_lsq/dogbox.py`. SciFort exposes the TRF operator
  path as paired Fortran callbacks carrying the current parameter vector explicitly.
  Generated regression values are public numerical outputs from
  SciPy 1.17.0. Complete SciPy, MINPACK, and LSMR notices are in
  `THIRD_PARTY_LICENSES.md`.
- `src/optimize/minimize_scalar.f90` and the scalar-minimization exports in
  `src/optimize.f90`: the bounded Brent kernel is adapted from fmin 1.1.2 by
  Jacob Williams (<https://github.com/jacobwilliams/fmin/tree/1.1.2>),
  BSD-3-Clause. SciFort replaces fmin's standalone bounded-function interface
  with `minimize_result`, explicit status values, optional tolerance/iteration
  controls, and SciPy-style evaluation accounting. Downhill bracketing,
  unbounded Brent, golden-section minimization, method dispatch, defaults, and
  result/count semantics are adapted from SciPy 1.17.0
  `scipy/optimize/_optimize.py` and `scipy/optimize/_minimize.py` under SciPy's
  BSD-3-Clause license. Generated references in
  `test/optimize_minimize_scalar_reference.f90` come from SciPy 1.17.0.
  The public `bracket` wrapper reuses the same SciPy-derived downhill
  bracketing kernel and adds a Fortran result object and configurable
  `grow_limit=`/`maxiter=`. Its direct point/value/work-count comparisons are
  generated by `tools/generate_optimize_bracket_reference.py` and stored in
  `test/optimize_bracket_reference.f90`. The NaN precheck avoids invalid IEEE
  comparisons under trapping builds. Complete upstream notices are in
  `THIRD_PARTY_LICENSES.md`.
- `src/differentiate.f90`: the adaptive nested finite-difference stencil,
  Vandermonde weight construction, error estimate, and termination semantics are
  adapted from SciPy 1.17.0 `scipy/differentiate/_differentiate.py`
  (<https://github.com/scipy/scipy/blob/v1.17.0/scipy/differentiate/_differentiate.py>)
  under SciPy's BSD-3-Clause license. SciFort replaces NumPy broadcasting and
  array-API orchestration with scalar/vector Fortran callbacks and an independent
  dense partial-pivot solve for the small stencil systems. The scalar complex-step
  identity `aimag(f(x+i*h))/h` is adapted from ForDiff 0.1.1 by Seyed Ali Ghasemi
  (<https://github.com/gha3mi/fordiff>), BSD-3-Clause; SciFort's vector gradient
  and Jacobian forms are independently implemented by perturbing one coordinate at
  a time. NumDiff 1.7.0 by Jacob Williams
  (<https://github.com/jacobwilliams/NumDiff>) was reviewed for modern-Fortran
  Jacobian API organization, bounds/sparsity concepts, and finite-difference
  architecture; no NumDiff implementation code is copied or adapted in this file.
- `src/math.f90`: stable elementary helpers implemented from standard
  power-series and algebraic identities. `floor_real` implements the mathematical
  floor using real truncation and a negative-fraction correction, without an
  intermediate integer conversion; nonfinite values pass through unchanged.
- `src/stats/normal.f90`: standard normal `erfc` identity, a
  Mills-ratio asymptotic expansion for log-CDF tails, and safeguarded
  bisection for quantiles.
- `src/stats/uniform.f90`: direct implementation from the mathematical
  definition.
- `src/stats/exponential.f90`: direct implementation with stable
  `expm1` and `log1p` helpers.
- `src/stats/laplace.f90`: direct piecewise implementation from the
  mathematical definition.
- `src/stats/logistic.f90`: direct implementation using stable
  logistic and softplus identities.
- `src/stats/cauchy.f90`: direct implementation using reciprocal-angle
  identities and cancellation-resistant logarithmic tails.
- `src/c_api/api.f90`: interoperable wrapper layer over the native
  Fortran API.
- `src/special/incomplete_gamma.f90`: regularized incomplete gamma
  functions and inverses. Independently implemented from these
  specifications, not from any library source:
  - NIST Digital Library of Mathematical Functions (DLMF),
    https://dlmf.nist.gov/, equations 8.7.1 (series for `P`), 8.7.3 (series
    used for `Q` when `a < 1` and `x <= 1.5`), 8.9.2 (continued fraction for
    `Q`), 5.7.3 (series for `log(Gamma(1 + a))`), 5.11.1 (Stirling series),
    and 4.6.1 (series for `log(1 + t) - t`);
  - the modified Lentz algorithm for continued fractions: W. J. Lentz,
    Applied Optics 15(3) (1976) 668, doi:10.1364/AO.15.000668;
    I. J. Thompson and A. R. Barnett, Journal of Computational Physics 64
    (1986) 490-509, doi:10.1016/0021-9991(86)90046-X;
  - the Wilson-Hilferty approximation for inverse starting values:
    E. B. Wilson and M. M. Hilferty, "The Distribution of Chi-Square",
    Proceedings of the National Academy of Sciences 17 (1931) 684-688,
    doi:10.1073/pnas.17.12.684.

  The region selection, safeguarded logarithmic Newton iteration, and
  tolerances are original. The 29 constants `zeta(k) - 1`, `k = 2..30`, were
  generated with mpmath 1.3.0 at 40 digits and checked against the closed
  forms of `zeta(2)` and `zeta(4)`. The Stirling coefficients are the exact
  rationals `B(2k) / (2k (2k - 1))` for `k = 1..7`, checked with mpmath.
- `src/special/log_gamma.f90`: log-gamma and log-beta helpers shared by the
  special functions, moved from `incomplete_gamma.f90` and extended with
  `log_beta`, `log_gamma_difference`, and `log_gamma_ratio_scaled`, all
  independently derived from DLMF 5.7.3 and 5.11.1 (Stirling series with the
  large terms cancelled analytically).
- `src/special/incomplete_beta.f90`: regularized incomplete beta function and
  inverses. Independently implemented from DLMF 8.17.4 (symmetry), 8.17.7
  (series used for the small-parameter complement), 8.17.22 and 8.17.23
  (continued fraction, evaluated by the modified Lentz method cited above),
  and 5.11.1. The expansion in incomplete gamma functions for one large and
  one small parameter (`gamma_expansion`) was derived for this project from
  the substitution `t = exp(-v)` in the integral for `B_x(a, b)`, the
  factorization `1 - exp(-v) = v exp(-v/2) sinh(v/2) / (v/2)`, the
  Bernoulli generating function of DLMF 24.2.1, and the recurrence of DLMF
  8.8.2. No implementation of a similar expansion was consulted. The 30
  coefficients `B(2n) / (2n (2n)!)` were generated with mpmath 1.3.0 at 40
  digits and checked against `log(sinh(v/2) / (v/2))` to `1e-40` at `v = 1`.
- `src/stats/beta.f90`, `src/stats/student_t.f90`,
  `src/stats/f_distribution.f90`: direct implementations from the density
  definitions over the incomplete beta function.
- `src/stats/poisson.f90`: direct implementation using DLMF 8.4.8 with
  8.4.11, which express the finite Poisson sum through the incomplete gamma
  functions, and the kernel of `src/special/incomplete_gamma.f90`.
- `src/stats/binomial.f90`: direct implementation using DLMF 8.17.5, which
  expresses the binomial tail as a regularized incomplete beta function, and
  the kernel of `src/special/incomplete_beta.f90`. The quantile search
  (a normal starting approximation with a continuity correction, a doubling
  bracket, and bisection over integers) is original in both modules.
- `src/stats/bernoulli.f90`: exact specialization of the binomial distribution
  to one trial, implemented by delegating to the SciFort binomial kernel.
- `src/stats/geometric.f90`: direct implementation from the mathematical
  definition. Stable CDF and survival functions use the project-local
  `log1p_safe` and `expm1_safe` helpers. The analytic logarithmic quantile and
  integer correction are original SciFort code.
- `src/stats/negative_binomial.f90`: direct implementation from the mass
  definition and DLMF 8.17.5, which gives the CDF as a regularized incomplete
  beta function. The log-mass rearrangement through the existing beta kernel,
  and the normal-start/doubling/bisection quantile search, are original
  SciFort code.
- `src/special.f90`: re-export module for special functions.
- `src/stats/gamma.f90`, `src/stats/chi2.f90`: direct implementations from
  the density definitions over the incomplete gamma functions.
- `src/stats/lognormal.f90`: direct implementation using standard normal
  distribution functions over standardized log-transformed coordinates.
- `src/stats/weibull.f90`: direct implementation from the mathematical definition
  of the Weibull minimum distribution with stable logarithmic tails and `expm1_safe`.
- `src/stats/pareto.f90`: direct implementation from the Pareto Type I definition
  with stable logarithmic tails and `expm1_safe`.
- `src/stats/rayleigh.f90`: direct implementation from the Rayleigh distribution
  definition with stable logarithmic tails and `expm1_safe`.
- `src/stats/gumbel_r.f90` and `src/stats/gumbel_l.f90`: original SciFort
  direct implementations from the standard Gumbel density/CDF definitions.
  Stable logarithmic and inverse-tail rearrangements were derived for this
  project. SciPy documentation was consulted only for public parameterization
  and numerical parity; no SciPy source was copied or adapted.
- `src/stats/powerlaw.f90`: original SciFort direct implementation of the
  power-function density `a*x**(a-1)` on standardized `[0,1]`. The SciPy
  `scipy.stats.powerlaw` documentation was used to confirm parameterization and
  endpoint conventions; no implementation source was copied or translated.
- `src/stats/triang.f90`: original SciFort piecewise implementation of the
  triangular density/CDF on standardized `[0,1]`, including the `c=0` and
  `c=1` endpoint modes and cancellation-resistant inverse formulas. SciPy's
  `scipy.stats.triang` documentation was used only as an API/parity reference.
- `src/stats/genpareto.f90`: original SciFort implementation of the generalized
  Pareto density `(1+c*x)**(-1-1/c)` and its exponential `c=0` limit. Direct
  logarithmic tails and quantiles use project-local `log1p_safe`/`expm1_safe`.
  SciPy's `scipy.stats.genpareto` documentation was used only to confirm public
  parameterization/support and for numerical validation; no SciPy source code
  was copied, translated, or adapted.
- `src/stats/arcsine.f90`: original SciFort implementation from the standard
  arcsine density/CDF identities on `[0,1]`, with complementary-tail formulas
  chosen directly from the same identities. SciPy 1.17.0 was used only for
  numerical/API parity checks; no SciPy source was copied or translated.
- `src/stats/halfnorm.f90`: original SciFort implementation from the standard
  half-normal density and error-function CDF, reusing SciFort's normal tail and
  quantile kernels. SciPy 1.17.0 was used only for numerical/API parity checks;
  no SciPy source was copied or translated.
- `src/stats/halfcauchy.f90`: original SciFort implementation from the
  standard half-Cauchy density and arctangent CDF identities, with direct
  complementary-tail formulas for numerical stability. SciPy 1.17.0 was used
  only for numerical/API parity checks; no SciPy source was copied or translated.
- `src/stats/lomax.f90`: original SciFort implementation of the Pareto-II/Lomax
  density and survival identities, using project-local `log1p_safe` and
  `expm1_safe` helpers for stable tails and inverse transforms. SciPy 1.17.0 was
  used only for numerical/API parity checks; no SciPy source was copied or translated.
- `src/stats/chi.f90`: original SciFort implementation using the standard chi
  density and the regularized incomplete-gamma identities for the CDF,
  survival function, and quantiles. SciPy 1.17.0 was used as an API and
  numerical-validation reference; no SciPy implementation source was copied.
- `src/stats/maxwell.f90`: original SciFort specialization of the chi
  distribution at three degrees of freedom, matching the standard Maxwell
  parameterization used by `scipy.stats.maxwell`.
- `src/stats/cosine.f90`: original SciFort implementation of the bounded
  cosine density and distribution formulas used by `scipy.stats.cosine`;
  quantiles use deterministic bisection on the monotone closed-form CDF.
- `src/stats/semicircular.f90`: original SciFort implementation of the
  Wigner semicircle density/CDF used by `scipy.stats.semicircular`; quantiles
  use deterministic bisection on the monotone closed-form CDF.
- `src/stats/anglit.f90`: original SciFort implementation of the standard
  Anglit density, trigonometric CDF, and inverse-CDF identities. SciPy 1.17.0
  was used for numerical and API parity checks; no SciPy implementation code
  was copied or translated.
- `src/stats/moyal.f90`: original SciFort implementation of the standard Moyal
  density and error-function tail identities, with explicit log-tail handling
  and normal-quantile inversion for stable probabilities. SciPy 1.17.0 was used
  for numerical and API parity checks; no SciPy implementation code was copied
  or translated.
- `src/stats/hypsecant.f90`: original SciFort implementation from the standard
  hyperbolic-secant density `sech(x)/pi`, its arctangent CDF, and analytic
  inverse-CDF identity. Stable log tails and `log(cosh(x))` are independently
  rearranged to avoid overflow. SciPy 1.17.0 was used only for public
  parameterization and numerical parity checks; no SciPy source was copied.
- `src/stats/invgamma.f90`: original SciFort implementation from the standard
  inverse-gamma density and regularized incomplete-gamma identities. SciPy
  1.17.0 was used only for public parameterization and numerical parity checks;
  no SciPy implementation source was copied or translated.
- `src/stats/invgauss.f90`: original SciFort implementation of the
  inverse-Gaussian density and normal-CDF representation, with deterministic
  tail-aware bisection for quantiles. SciPy 1.17.0 was used only for numerical
  and API parity checks; no SciPy source was copied or translated.
- `src/stats/levy.f90`: original SciFort implementation from the standard Levy
  density, error-function CDF, and inverse-normal/error-function quantile
  identities. No third-party implementation code was copied or translated.
- `src/stats/loglaplace.f90`: original SciFort implementation from the
  piecewise power-law form of the log-Laplace density/CDF and its analytic
  inverse. SciPy 1.17.0 was used only for parameterization and parity checks.
- `src/stats/halflogistic.f90`: original SciFort implementation from the
  standard half-logistic density and hyperbolic-tangent CDF identities. Direct
  survival/log-survival and inverse-survival formulas avoid cancellation in the
  upper tail. SciPy 1.17.0 was used only for public parameterization and
  numerical parity checks; no SciPy source was copied.

- `src/random.f90`: clean implementation of the MRG32k3a combined multiple
  recursive generator from P. L'Ecuyer, "Good Parameters and Implementations
  for Combined Multiple Recursive Random Number Generators", Operations
  Research 47(1) (1999) 159-164, doi:10.1287/opre.47.1.159. The published
  moduli and recurrence coefficients were transcribed from the paper; no
  software implementation was copied or translated. The scalar-seed mapping,
  state API, rejection into equal 26-bit buckets, and two-chunk 52-bit
  binary64 conversion are original SciFort code.
- `src/random_variates.f90`: original SciFort inverse-transform sampling
  wrappers. They combine `src/random.f90` with the existing distribution PPF
  routines and do not incorporate code from SciPy or another RNG library.
- `src/stats/bradford.f90`, `src/stats/truncexpon.f90`, `src/stats/fisk.f90`,
  and `src/stats/dweibull.f90`: original SciFort implementations from the standard
  mathematical definitions of the Bradford, finite truncated-exponential,
  log-logistic/Fisk, and symmetric double-Weibull distributions. Stable tail and
  inverse formulas were independently rearranged using the project-local
  `log1p_safe`, `expm1_safe`, and `log1pexp` helpers. SciPy 1.17.0 was used only
  as an independent numerical validation reference; no SciPy source code was
  copied, translated, or adapted.
- `src/stats/alpha.f90`, `src/stats/fatiguelife.f90`,
  `src/stats/genlogistic.f90`, and `src/stats/gennorm.f90`: original SciFort
  implementations from the published mathematical definitions of the Alpha,
  Birnbaum-Saunders, generalized-logistic, and generalized-normal families.
  The generalized-normal CDF and inverse reuse SciFort's independently
  implemented regularized incomplete-gamma kernel. Stable tail and inverse
  rearrangements were derived independently. SciPy 1.17.0 was used only for
  public parameterization and numerical validation; no SciPy source code was
  copied, translated, or adapted.
- `src/likelihood.f90`: original SciFort reductions over the existing
  distribution `logpdf`/`logpmf` routines. The family-specific
  `*_loglikelihood` and `*_nnlf` interfaces contain no imported numerical
  algorithm or translated third-party implementation.
- `src/special/digamma.f90`: original SciFort positive-argument digamma
  helper, independently implemented from DLMF 5.5.2 (recurrence) and DLMF
  5.11.2 (asymptotic expansion). The Bernoulli-number coefficients are exact
  rational constants. No third-party implementation was copied or translated.
- `src/score.f90`: original SciFort analytic derivatives of the existing
  distribution log-density and log-mass formulas. Shape derivatives use the
  project-local positive-argument digamma helper. No SciPy or third-party
  implementation was copied, translated, or adapted.
- `src/fit.f90`: original SciFort deterministic projected coordinate/pattern
  search and family-specific fitting wrappers over `scifort_likelihood`.
  Equal-bound fixed parameters and integer-constrained fit coordinates follow
  the documented public fitting conventions used by SciPy, but the optimizer
  itself is independently implemented and does not copy or translate SciPy's
  default differential-evolution implementation or another optimizer.
- `src/math.f90` (`log1p_safe`, `expm1_safe`): series from DLMF 4.6.4 and
  4.2.19.

- `src/descriptive_lowrank.f90`: original SciFort implementations of standard
  descriptive-statistics definitions. The quantile interpolation is the
  Hyndman-Fan type 7 convention from R. J. Hyndman and Y. Fan, "Sample
  Quantiles in Statistical Packages", The American Statistician 50(4)
  (1996) 361-365, doi:10.1080/00031305.1996.10473566. Scaled compensated
  summation, scaled centered powers, merge sorting, and rank grouping were
  implemented independently for this project. The later `gmean`, `hmean`,
  `pmean`, `sem`, `variation`, `iqr`, `median_abs_deviation`, `moment`, `skew`,
  and `kurtosis` additions likewise use their standard mathematical definitions;
  SciPy 1.17.0 was used only to define compatible option semantics and generate
  independent regression values in `test/test_descriptive_scipy_compat.f90`
  through `tools/generate_statistics_reference.py`. The later rank-2 axis and
  fixed-rank keepdims wrappers reuse the same vector kernels; SciPy 1.17.0 was
  used only to generate independent axis, keepdims, and broadcasting references
  in `test/test_descriptive_axis_keepdims.f90`. The later native
  `nan_policy='propagate'` / `nan_policy='omit'` overloads are original SciFort
  wrappers around those kernels; SciPy 1.17.0 was used to define omission
  semantics and generate independent values in
  `test/test_descriptive_nan_policy.f90`. The subsequent unweighted rank-3 axis
  and keepdims wrappers also reuse only those vector kernels; SciPy 1.17.0 was
  used solely to generate independent references in
  `test/test_descriptive_rank3.f90`. The subsequent rank-3 weighted mean wrappers
  use original SciFort broadcast helpers and the same vector kernels; SciPy
  1.17.0 was used solely to generate independent rank-1/rank-2/rank-3 weight,
  singleton-expansion, keepdims, and NaN-policy references in
  `test/test_descriptive_rank3_weights.f90`. The subsequent unweighted rank-4
  axis and keepdims wrappers are original SciFort slice wrappers around the same
  vector kernels; SciPy 1.17.0 was used solely to generate independent rank-4
  axis, whole-array IQR, keepdims, invalid-axis, and NaN-policy references in
  `test/test_descriptive_rank4.f90`. `src/descriptive_highrank_means.f90`
  contains the subsequent rank-5 through rank-15 unweighted mean-family
  support as original SciFort code: one assumed-rank
  dispatcher flattens and indexes slices, then reuses the same vector kernels;
  only thin fixed-result-rank wrappers are rank-specific. SciPy 1.17.0 was
  used solely to generate independent high-rank axis, keepdims, NaN-policy,
  and zero-extent references in `test/test_descriptive_highrank_means.f90`.
  `src/descriptive_weighted_engine.f90` and the rank-split
  `src/descriptive_highrank_weights_4_9.f90` /
  `src/descriptive_highrank_weights_10_15.f90` are original SciFort code for
  weighted rank-4 through rank-15 mean-family reductions. The shared
  assumed-rank engine computes NumPy-style trailing joint-broadcast indices and
  reuses the existing vector kernels without materializing broadcast arrays.
  SciPy 1.17.0 was used solely to generate independent rank-4/rank-15 weighted,
  singleton-expansion, keepdims, and NaN-weight references in
  `test/test_descriptive_highrank_weights.f90`.
  `src/descriptive_highrank_stats_engine.f90` and the rank-split
  `src/descriptive_highrank_stats_5_9.f90` /
  `src/descriptive_highrank_stats_10_15.f90` are original SciFort code for
  rank-5 through rank-15 `sem`, `variation`, `iqr`,
  `median_abs_deviation`, `moment`, `skew`, and `kurtosis`. One shared
  assumed-rank slice engine performs column-major indexing and delegates each
  one-dimensional slice to the existing low-rank vector kernels; only thin
  fixed-result-rank wrappers are rank-specific. SciPy 1.17.0 was used solely
  to generate independent rank-15 axis, keepdims, whole-array IQR, NaN-policy,
  and nonlinear numerical references in `test/test_descriptive_highrank_stats.f90`.
  `src/descriptive_named_scale_engine.f90`,
  `src/descriptive_named_scale_lowrank.f90`, and
  `src/descriptive_named_scale_highrank.f90` are original SciFort adapters for
  SciPy-compatible `scale='normal'` on IQR and MAD. They use the mathematical
  conversion based on the standard-normal 75th percentile,
  `Phi^-1(0.75) = 0.6744897501960817`, and delegate to the existing numeric-scale
  kernels; SciPy 1.17.0 was used only to verify option semantics and generate
  independent regression values in `test/test_descriptive_named_scale.f90`.
  `src/descriptive_standardize_engine.f90` and the rank-split
  `src/descriptive_standardize_1_9.f90` /
  `src/descriptive_standardize_10_15.f90` are original SciFort implementations
  of `zscore`, `gzscore`, `zmap`, and `gstd`. The definitions are the standard
  arithmetic and geometric standard-score identities; the shared assumed-rank
  engine uses original scaled/compensated mean and variance calculations,
  column-major slice indexing, explicit broadcast-index mapping for `zmap`, and
  original flattened index decomposition for whole-array and multi-axis `gstd`
  reductions. SciPy 1.17.0 public documentation and source behavior were used to
  verify option, domain, broadcasting, axis/tuple-axis, and constant-comparison
  semantics; its public
  numerical outputs were used solely as independent regression values in
  `test/test_descriptive_standardize.f90`. SciPy's public `nan_policy='raise'`
  behavior is represented natively with Fortran 2018 `error stop` and verified
  in `test/test_descriptive_nan_raise.f90`. No NumPy or SciPy implementation
  source was copied, translated, or adapted.
  `src/descriptive_trimmed_engine.f90` and the rank-split
  `src/descriptive_trimmed_1_8.f90` / `src/descriptive_trimmed_9_15.f90` are
  original SciFort implementations of `tmean`, `tvar`, `tstd`, `tsem`, `tmin`,
  and `tmax`. They use the standard definitions after filtering observations by
  lower/upper value limits; the assumed-rank reducer and column-major slice
  indexing are original SciFort code. SciPy 1.17.0 public documentation and
  public numerical behavior were used to define defaults, endpoint and NaN
  semantics, and independent regression values in
  `test/test_descriptive_trimmed.f90` and `test/test_descriptive_nan_raise.f90`.
  No SciPy implementation source was copied, translated, or adapted.
  `src/descriptive_proportion_trim_engine.f90` and the rank-split
  `src/descriptive_proportion_trim_1_8.f90` /
  `src/descriptive_proportion_trim_9_15.f90` are original SciFort
  implementations of `trim_mean`, `trimboth`, and `trim1`. The implementation
  uses original stable merge-sorting, assumed-rank slice indexing, and
  fixed-rank reshape wrappers. SciPy 1.17.0 public documentation and public
  numerical behavior were used to define proportion rounding, axis, NaN,
  keepdims, and tail semantics and as independent reference values in
  `test/test_descriptive_proportion_trim.f90` and its child-process error
  regression `test/test_descriptive_proportion_trim_errors.f90`. No SciPy implementation source
  was copied, translated, or adapted.
  `src/descriptive_kstat_engine.f90` and the rank-split
  `src/descriptive_kstat_1_8.f90` / `src/descriptive_kstat_9_15.f90` are
  original SciFort implementations of `kstat` and `kstatvar`. The vector kernel
  uses the standard centered-moment identities for unbiased cumulant estimators;
  the assumed-rank reducer and fixed-rank reshape wrappers are original SciFort
  code. SciPy 1.17.0 public documentation at
  <https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.stats.kstat.html>
  and
  <https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.stats.kstatvar.html>
  was used to define supported orders, axis, keepdims, and NaN semantics, and
  SciPy public numerical outputs were used as independent regression references
  in `test/test_descriptive_kstat.f90`. No SciPy implementation source was
  copied, translated, or adapted.
  `src/descriptive_circular_engine.f90` and the rank-split
  `src/descriptive_circular_1_8.f90` / `src/descriptive_circular_9_15.f90` are
  original SciFort implementations of `circmean`, `circvar`, and `circstd`, using
  the standard mean-resultant-vector definitions of circular location and spread.
  SciPy 1.17.0 public API behavior was used to define default angular ranges,
  axis/keepdims, NaN, normalization, and empty-sample semantics, and its public
  numerical outputs provide independent references in
  `test/test_descriptive_circular.f90`. No SciPy source code is included.
  `src/descriptive_frequency_engine.f90` and the rank-split
  `src/descriptive_frequency_1_8.f90` / `src/descriptive_frequency_9_15.f90` are
  original SciFort implementations of `percentileofscore`, `scoreatpercentile`,
  `cumfreq`, and `relfreq`. The code uses independently written sorting, percentile
  interpolation, histogram binning, assumed-rank flattening, and fixed-rank reshape
  wrappers. SciPy 1.17.0 public documentation and public numerical behavior were
  consulted for argument/default conventions, tie handling, NaN ordering/policies,
  histogram-limit behavior, extrapoint counts, and regression references in
  `test/test_descriptive_frequency.f90`. No NumPy or SciPy implementation source was
  copied, translated, or adapted.
  `src/descriptive_rank_engine.f90` and the rank-split
  `src/descriptive_rank_1_8.f90` / `src/descriptive_rank_9_15.f90` are original
  SciFort implementations of multidimensional `rankdata` and `tiecorrect`. The ranking
  kernel uses independently written stable index sorting, explicit tie-group scanning,
  NaN filtering, assumed-rank slice extraction, and fixed-rank reshape wrappers. SciPy
  1.17.0 public documentation and public numerical behavior were consulted for the five
  tie methods, axis/default flattening, NaN-policy semantics, tie-correction behavior,
  and regression references. No NumPy or SciPy implementation source was copied,
  translated, or adapted.
  `src/descriptive_summary_engine.f90` is an original SciFort implementation of SciPy-style
  `mode` and `describe`. The implementation uses independently written merge sorting,
  tie-run scanning, assumed-rank slice extraction, and multi-output allocation helpers.
  SciPy 1.17.0 public documentation, public numerical behavior, and source were consulted
  to verify axis/defaults, tie handling, NaN policies, empty-mode conventions, `ddof`,
  bias, and the masked-invalid omission path. No SciPy implementation text was copied or
  translated; the sorting, reduction, allocation, and moment logic is independently written.
  Independent reference values are generated by
  `tools/generate_statistics_reference.py` and checked by
  `test/test_descriptive_summary.f90`. No SciPy source code was copied or translated.
  `src/descriptive_entropy_engine.f90` and the rank-split
  `src/descriptive_entropy_1_8.f90` / `src/descriptive_entropy_9_15.f90` are
  original SciFort implementations of `entropy` and `differential_entropy`.
  Shannon and Kullback-Leibler calculations use the existing SciFort `entr` and
  `rel_entr` elementwise kernels after independently written normalization,
  broadcasting, paired-NaN filtering, and assumed-rank slice extraction. The
  differential-entropy estimators are independently written from the published
  Vasicek (1976), Van Es (1992), Ebrahimi-Pflughoeft-Soofi (1994), and Correa
  (1995) spacing-estimator formulas cited by SciPy. SciPy 1.17.0 documentation,
  public numerical behavior, and source were consulted to verify API defaults,
  the `auto` method-selection thresholds, and edge-case semantics; no SciPy code
  was copied verbatim. Independent regression values are generated by
  `tools/generate_statistics_reference.py` and checked in
  `test/test_descriptive_entropy.f90`.
  `src/descriptive_expectile_sigmaclip.f90` is an original SciFort implementation
  of `expectile` and `sigmaclip`. The expectile calculation is derived independently
  from the empirical expectile estimating equation in Newey and Powell (1987),
  Econometrica 55, 819-847, and solves each finite nonnegative-weight problem by
  scanning the piecewise-linear estimating equation after sorting the observations.
  A separately written secant fallback supports negative case weights. Sigma clipping
  follows the public iterative definition based on the population mean and standard
  deviation. SciPy 1.17.0 documentation, public numerical behavior, and `expectile` /
  `sigmaclip` source were consulted only to verify defaults, broadcasting, iteration,
  and edge behavior; no SciPy implementation text was copied or translated.
  Independent references are generated by `tools/generate_statistics_reference.py`
  and checked in `test/test_descriptive_expectile_sigmaclip.f90`.
  `src/descriptive_lmoment_engine.f90` and the rank-split
  `src/descriptive_lmoment_1_8.f90` / `src/descriptive_lmoment_9_15.f90` are
  original SciFort implementations of sample L-moments. The vector kernel uses
  the standard probability-weighted-moment representation and L-moment linear
  combinations described by J. R. M. Hosking, "L-Moments: Analysis and Estimation
  of Distributions Using Linear Combinations of Order Statistics", JRSS B 52
  (1990), 105-124. Sorting, probability-weight recursion, assumed-rank slice
  extraction, IEEE-safe arithmetic, and fixed-rank result wrappers were written
  independently for SciFort. SciPy 1.17.0 documentation, public numerical behavior,
  and source were consulted to verify defaults, order/result-shape semantics,
  standardization, keepdims layout, and NaN handling; no SciPy implementation text
  was copied or translated. Independent references are generated by
  `tools/generate_statistics_reference.py` and checked by
  `test/test_descriptive_lmoment.f90` and its child-process error regression.
- `src/stats/power_transform_engine.f90`, the rank-split
  `src/stats/power_transform_1_8.f90` / `src/stats/power_transform_9_15.f90`, and
  `src/stats/power_transform.f90` are original SciFort implementations of the
  statistical Box-Cox and Yeo-Johnson power-transform layer. Transform formulas
  follow G. E. P. Box and D. R. Cox (1964) and I. Yeo and R. A. Johnson,
  Biometrika 87 (2000), 954-959. The profile log-likelihoods follow the published
  Jacobian/variance definitions. The probability-plot criterion uses Filliben's
  order-statistic-median approximation (Technometrics 17, 1975) with SciFort's
  normal quantile. The bracket expansion and golden-section maximizer are original
  SciFort code. SciPy 1.17.0 documentation, public numerical behavior, and source
  were consulted to verify defaults, likelihood conventions, Filliben probabilities,
  parameter-search criteria, NaN handling, and edge behavior; no SciPy source text
  was copied or translated. Independent reference values are generated by
  `tools/generate_statistics_reference.py` and checked by
  `test/test_power_transform.f90` and `test/test_power_transform_errors.f90`.
- `src/stats/probability_plot.f90` is an original SciFort implementation of
  probability-plot and PPCC diagnostics. Filliben order-statistic medians follow
  J. J. Filliben, "The Probability Plot Correlation Coefficient Test for Normality",
  Technometrics 17 (1975), 111-117. Least-squares and Pearson-correlation calculations,
  distribution dispatch, callback plumbing, bracket expansion, and golden-section
  refinement were written independently for SciFort. SciPy 1.17.0 documentation,
  public numerical behavior, and source were consulted to verify API defaults,
  Filliben probabilities, two-/three-point bracket semantics, empty-sample behavior,
  and Box-Cox/Yeo-Johnson normality-plot conventions; no SciPy implementation text
  was copied or translated. Independent references are generated by
  `tools/generate_statistics_reference.py` and checked by
  `test/test_probability_plot.f90` and `test/test_probability_plot_errors.f90`.
- `src/stats/sample_distance.f90` is an original SciFort implementation of the
  one-dimensional Wasserstein-1 and energy distances and multidimensional
  Wasserstein-1 distance. The 1-D path uses the standard identity expressing the
  distances as integrals of absolute/squared differences between empirical CDFs, with
  independently written stable sorting, weight normalization, merged-support
  integration, validation, and IEEE nonfinite handling. The multidimensional path
  formulates balanced Euclidean transport as a complete bipartite min-cost-flow
  problem and solves it with an independently written successive-shortest-augmenting-
  path residual-network algorithm using Bellman-Ford shortest paths. Pairwise Euclidean
  costs use a scaled norm to avoid avoidable intermediate overflow. SciPy 1.17.0
  documentation, public numerical behavior, and source were consulted to verify public
  argument conventions, the linear-programming formulation, examples, and edge
  behavior; no SciPy implementation text was copied or translated. Independent
  references are generated by `tools/generate_statistics_reference.py` and checked by
  `test/test_sample_distance.f90`, `test/test_sample_distance_errors.f90`,
  `test/test_sample_distance_nd.f90`, and `test/test_sample_distance_nd_errors.f90`.
- `src/empirical_distribution.f90` is an original SciFort implementation of
  empirical CDF/survival estimation, two-sample logrank survival comparison, exact
  quantile inference, and continuous order-statistic transforms. The uncensored ECDF
  follows the elementary empirical frequency definition; the right-censored path
  follows the Kaplan-Meier product-limit estimator (Kaplan and Meier, JASA 53, 1958,
  457-481). Confidence intervals use Greenwood's variance formula and its
  log-negative-log/exponential-Greenwood transform. The logrank implementation follows
  the standard observed-minus-expected death count statistic with the Mantel variance,
  retaining its sign and using the asymptotic standard-normal law for one- and two-sided
  p-values. Quantile-test p-values and confidence bounds are independently written from
  the exact binomial construction described by Conover and Thompson; order-statistic
  methods use the standard Beta(r, n-r+1) probability transform. SciPy 1.17.0
  documentation, public numerical behavior, and source were consulted to verify
  object/result layout, right-censoring conventions, confidence-interval methods,
  logrank alternative conventions, T1/T2 tie semantics, and edge behavior; no SciPy
  implementation text was copied or translated. Regression references are generated by
  `tools/generate_statistics_reference.py` and checked by
  `test/test_empirical_distribution.f90`, `test/test_logrank.f90`, and their
  child-process error regressions.
  `src/descriptive.f90` is the thin public facade that merges the low-rank
  and high-rank generic interfaces without changing public procedure names.
  NaN domain/finite checks in the weighted means, power mean, skewness, and
  kurtosis paths are sequenced
  explicitly so IEEE exception-trapping builds do not depend on compiler
  short-circuit
  evaluation. No NumPy, SciPy, or Fortran stdlib implementation was copied or
  translated.
- `src/bayes_mvs.f90`: original SciFort Fortran implementation of Bayesian
  posterior inference for a normal sample with unknown mean and variance under
  Jeffreys' prior. The Student-t posterior for the mean and inverse-gamma posterior
  for the variance follow the standard conjugate/Jeffreys derivation described by
  T. E. Oliphant, "A Bayesian perspective on estimating mean, variance, and
  standard-deviation from data" (2006). The generalized-gamma representation of
  the standard-deviation posterior follows algebraically by the square-root change
  of variables. SciPy 1.17.0 `_morestats.py` and public documentation were inspected
  for public API semantics, frozen-distribution parameterization, the n > 1000
  Gaussian-approximation threshold/formulas, and edge behavior; no SciPy source text
  was copied or translated. Regression values are generated independently by
  `tools/generate_bayes_mvs_reference.py` and checked by `test/test_bayes_mvs.f90`.
- `src/stats/invgamma.f90` and `src/stats/gengamma.f90`: the positive-affine
  quantile overflow guard was rewritten during `mvsdist` validation so IEEE-trapping
  builds do not evaluate an overflowing `huge()/scale` threshold when `scale < 1`;
  the numerical distribution formulas are unchanged.
- `src/hypothesis.f90`: original SciFort implementations from the standard
  definitions of Student/Welch t statistics, the exact beta null law of the
  Pearson correlation coefficient, the rank-transform approximation for
  Spearman correlation, and the tie-corrected normal approximation to the
  Mann-Whitney U statistic. The Mann-Whitney test traces to H. B. Mann and
  D. R. Whitney, "On a Test of Whether one of Two Random Variables is
  Stochastically Larger than the Other", Annals of Mathematical Statistics
  18(1) (1947) 50-60. SciPy 1.17.0 was used only as an independent numerical
  validation reference. Its documented `ttest_ind_from_stats` argument/default
  behavior and public numerical edge cases were also checked when adding the
  elemental summary-statistics interface; the implementation remains an original
  formulation of the standard pooled-variance and Welch formulas. No SciPy source
  code was copied or adapted.
- `src/hypothesis_extended.f90`: original Fortran implementations of standard
  rank, contingency, analysis-of-variance, and scale-test formulas. SciPy
  1.17.0 source and documentation were consulted for public argument/default
  semantics, Wilcoxon method-selection behavior, named Cressie-Read lambda
  conventions, O'Brien-transform semantics, and edge-case validation. The
  O'Brien formula follows the standard transformation cited by SciPy to
  S. E. Maxwell and H. D. Delaney, *Designing Experiments and Analyzing Data:
  A Model Comparison Perspective* (1990), with SciPy 1.17.0 outputs used as
  regression references. No SciPy implementation code or tables were copied or
  translated. Exact signed-rank probabilities use an
  independently implemented subset-sum dynamic program; small tied/zero
  signed-rank cases use exhaustive sign enumeration.
- `src/contingency_meta.f90`: original SciFort Fortran implementations of exact
  binomial testing, Barnard/Boschloo unconditional 2-by-2 tests, odds ratio,
  relative risk, contingency association, p-value combination, and FDR control.
  SciPy 1.17.0 `_binomtest.py`, `_hypotests.py`, `_odds_ratio.py`,
  `_relative_risk.py`, `contingency.py`, `_stats_py.py`, and `_morestats.py`
  were consulted for public defaults, formulas, numerical tolerances, and edge-
  case semantics. No SciPy implementation source was copied or translated.
  Barnard/Boschloo nuisance maximization is independently implemented as a
  deterministic grid search plus golden-section refinement rather than SHGO.
- `src/multiple_comparisons.f90`: original SciFort Fortran implementations of
  the Alexander-Govern approximation, Tukey HSD/Tukey-Kramer, Games-Howell,
  Dunnett's single-step treatment-control procedure, and the
  Krishnamoorthy-Thomson Poisson means E-test. SciPy 1.17.0 `_stats_py.py`,
  `_hypotests.py`, and `_multicomp.py` were consulted for public defaults,
  formulas, numerical conventions, and edge semantics. The implementation
  reuses SciFort's independently implemented studentized-range, multivariate-t,
  Poisson, and chi-square kernels; no SciPy implementation source was copied
  verbatim. Dunnett confidence limits use an independently implemented
  monotone bisection against the SciFort multivariate-t box probability rather
  than SciPy's noisy Brent minimization.
- `test/test_multiple_comparisons.f90` and
  `test/multiple_comparisons_reference.f90`: SciFort-owned regression logic and
  public SciPy 1.17.0 numerical outputs regenerated by
  `tools/generate_multiple_comparisons_reference.py`.
- `test/test_contingency_meta.f90` and
  `test/contingency_meta_reference.f90`: SciFort-owned regression logic and
  public SciPy 1.17.0 numerical outputs regenerated by
  `tools/generate_contingency_meta_reference.py`.
- `src/nonparametric_extended.f90`: original SciFort implementations of the
  published Ansari-Bradley, Mood, Baumgartner-Weiss-Schindler (BWS),
  Epps-Singleton, k-sample Anderson-Darling, and median-test formulas. The BWS
  statistic follows the rank-based definition described by M. Neuhaeuser,
  *Statistical Papers* 46 (2005), 1-29. SciPy 1.17.0 `_morestats.py`,
  `_hypotests.py`, and `_bws_test.py` were consulted for public defaults, tie
  conventions, alternative direction, finite/asymptotic or permutation-method
  selection, correction factors, critical values, and edge-case semantics; no
  SciPy implementation code was copied or translated. SciFort's BWS interface
  uses an explicit `rng_state` for randomized permutations and reuses the
  existing independent-sample permutation engine. The exact small-sample
  Ansari-Bradley distribution uses an independently implemented subset-count
  dynamic program. Epps-Singleton uses SciFort's existing symmetric eigensolver
  for the covariance pseudoinverse.
- `test/test_nonparametric_extended.f90` and
  `test/nonparametric_extended_reference.f90`: SciFort-owned regression logic
  and public SciPy 1.17.0 numerical outputs; the reference module is regenerated
  by `tools/generate_nonparametric_extended_reference.py`.
- `src/binned_statistic.f90`: original SciFort implementation of one-, two-, and
  arbitrary-dimensional binned statistics. The binning rules follow the standard
  half-open histogram convention with an inclusive final bin; equal-width edge
  construction, C-order flattened bin numbering, built-in statistic behavior, and
  edge-case semantics were checked against SciPy 1.17.0 public documentation and
  numerical outputs. SciPy 1.17.0 `_binned_statistic.py` was inspected to verify
  rightmost-edge and outlier-numbering semantics, but no SciPy implementation code
  was copied, translated, or adapted. The Fortran arbitrary-dimensional result is
  independently designed as a flattened statistic plus explicit shape/edge-count
  metadata, and custom statistics use a pure callback interface.
- `test/test_binned_statistic.f90` and `test/test_binned_statistic_reference.f90`:
  SciFort-owned regression logic and public SciPy 1.17.0 numerical outputs regenerated
  by `tools/generate_binned_statistic_reference.py`.
- `src/association_extended.f90`: SciFort Fortran implementations of Kendall
  tau, weighted Kendall tau, Spearman rho, Somers D ordinal association, Chatterjee xi
  correlation, point-biserial
  correlation, ordinary and robust
  line fitting, the Brunner-Munzel test, and Page's ordered-alternative test. Public
  defaults, tie conventions, method selection, and finite-sample formulas were
  checked against SciPy 1.17.0. The Somers D contingency-table concordance,
  discordance, and asymptotic standard-error formulas are adapted from SciPy
  1.17.0 `_stats_pythran.py` and `_hypotests.py` under the BSD-3-Clause license;
  the source-file notice and `THIRD_PARTY_LICENSES.md` retain the applicable terms.
  `spearmanrho` is independently implemented as Pearson correlation of average ranks,
  following the mathematical definition; SciPy 1.17.0 `_stats_py.py` was inspected for
  public method semantics and numerical parity. Its exact/randomized pairing enumeration
  and plus-one p-values reuse SciFort-owned permutation conventions, and its continuous
  Monte Carlo null is generated as independent random rank permutations. No SciPy
  `spearmanrho` implementation code was copied. Chatterjee xi follows the published
  rank formulas
  and asymptotic null variance from Chatterjee (2021); permutation inference uses
  SciFort-owned exact enumeration and explicit-state Fisher-Yates re-pairing. The
  exact Kendall inversion-count recurrence is adapted from SciPy's `_kendall_p_exact`
  implementation. Weighted Kendall's exchange-weight normalization, joint/tie
  accounting, default hyperbolic weigher, and decreasing-lexicographic rank semantics
  are adapted from SciPy 1.17.0 `_stats.pyx` / `_stats_py.py`; SciFort uses an
  independently written direct pair-count implementation and Fortran rank-mode
  spelling. The applicable SciPy BSD-3-Clause license is reproduced in
  `THIRD_PARTY_LICENSES.md`. Page exact probabilities are independently implemented
  by convolution of the single-row permutation distribution.
- `test/test_association_extended.f90` and
  `test/association_extended_reference.f90`: SciFort-owned regression logic and
  public SciPy 1.17.0 numerical outputs regenerated by
  `tools/generate_association_extended_reference.py`.
- `src/resampling.f90`: original SciFort explicit-state resampling code.
  Percentile/basic bootstrap intervals follow the ordinary nonparametric
  bootstrap framework, and BCa uses the standard bias-correction plus
  leave-one-out jackknife acceleration formulas described by Efron and
  Tibshirani. Exact permutation enumeration, randomized Fisher-Yates shuffles,
  Monte Carlo p-values, and simulated power were independently implemented.
  SciPy 1.17.0 documentation/source behavior was used to verify exact-switching,
  the `independent`/`samples`/`pairings` permutation definitions, plus-one
  randomized p-values, the `100*epsilon` comparison tolerance, one-sided
  bootstrap conventions, and the strict `p < significance` power definition.
  Selection from the finite 52-bit RNG grid uses an independently implemented
  rejection step to remove modulo bias. No third-party resampling source was
  copied or translated.

- `test/resampling_inference_reference.f90` and
  `tools/generate_resampling_inference_reference.py`: deterministic public
  SciPy 1.17.0 outputs for exact permutation, Monte Carlo, and power regression
  cases; these files contain numerical references, not SciPy implementation
  source.

- `src/special/elementary.f90`: original SciFort implementation. `erf` and
  `erfc` are thin wrappers around the Fortran standard intrinsics. The inverse
  error functions use the identities with the standard-normal quantile;
  `gammaln` and negative-argument `digamma` use the classical gamma/digamma
  reflection formulas (DLMF 5.5), while positive digamma delegates to the
  existing independently implemented recurrence/asymptotic kernel. Logistic,
  entropy, `exprel`, and `cosm1` routines follow their mathematical
  definitions with independently chosen stable rearrangements. Box-Cox and
  inverse Box-Cox functions follow G. E. P. Box and D. R. Cox, "An Analysis
  of Transformations", Journal of the Royal Statistical Society Series B
  26(2) (1964) 211-252, using the project-local `log1p_safe` and
  `expm1_safe` kernels. SciPy 1.17.0 was used only for numerical validation;
  no SciPy implementation code was copied, translated, or adapted.
- `src/special/error_fresnel.f90`: original SciFort implementation of real
  `dawsn`, `erfi`, `erfcx`, and `fresnel`. Dawson and Fresnel definitions and
  large-argument expansions follow NIST DLMF Chapter 7, especially 7.2, 7.5,
  and 7.12. The intermediate real integrals use independently written
  composite 16-point Gauss-Legendre quadrature; the nodes and weights are
  mathematical roots/weights of the degree-16 Legendre rule. SciPy 1.17.0
  public numerical outputs are used only for ordinary-range regression
  comparison, while mpmath 1.3.0 high-precision outputs validate large
  quadratic-phase cases. No SciPy, mpmath, or other third-party implementation
  source was copied, translated, or adapted.
- `src/special/faddeeva.f90`: original SciFort implementation of complex `wofz`,
  complex `erfcx`, the complex `erf`, `erfc`, `erfi`, and Dawson kernels, and
  complex Fresnel S/C support. It follows NIST DLMF Chapter 7: 7.2.3 for the
  Faddeeva definition and its error-function relation, 7.4.3 for reflection,
  7.5.8 for the Fresnel/error-function interrelations, 7.6.3 for the entire
  Faddeeva power series, 7.7.2 for the upper-half-plane real-line integral
  representation, and 7.9.3 for the large-argument continued fraction. The
  transition-region integral uses independently written composite 16-point
  Gauss-Legendre quadrature; its nodes and weights are mathematical constants of
  the degree-16 Legendre rule. The complex error/Dawson specifics use the
  standard identities `erfc(z)=exp(-z**2)*wofz(i*z)`,
  `erfi(z)=-i*erf(i*z)`, and `D(z)=i*sqrt(pi)/2*(exp(-z**2)-wofz(z))`. Complex
  Fresnel evaluation uses the two conjugate-sign forms of DLMF 7.5.8 and the
  defining Fresnel entire series in the well-conditioned diagonal sector; exact
  coordinate axes reuse the independently implemented real Fresnel kernel.
  SciPy 1.17.0 public numerical outputs and mpmath 1.3.0 high-precision values
  are used only for regression validation. No SciPy, mpmath, Steven G. Johnson
  Faddeeva-package, or other third-party implementation source was copied,
  translated, or adapted.
- `src/special/reductions_low.f90`, `src/special/reductions_high_rank.f90`,
  `src/special/reductions_high_rank_weighted_support.f90`,
  `src/special/reductions_high_rank_weighted_6_10.f90`,
  `src/special/reductions_high_rank_weighted_11_15.f90`, and the API-preserving
  `src/special/reductions.f90` facade: original SciFort max-shift implementations
  of real/complex log-sum-exp, softmax, and log-softmax, including rank-2 through
  rank-15 single-/multiple-dimension reductions, fixed-rank keepdims companions,
  and scalar/trailing rank-1 through rank-15 weights with NumPy-style joint
  broadcasting that can expand singleton input dimensions. Ranks 7 through 15
  use an original SciFort flattened column-major group/index engine shared by the
  fixed-rank public wrappers, while weighted ranks 6 through 15 share an original
  broadcast-map/reduction engine. These avoid copied nested-loop numerical kernels
  while preserving Fortran's compile-time result ranks.
  Complex weights are accepted with real or complex inputs. Real weighted sums use
  separate compensated positive/negative accumulators and log-domain subtraction.
  Complex sums shift by logarithmic term magnitude and use compensated real/imaginary
  accumulation; zero-weight terms are ignored before exponentiation in both paths.
  The complex `log_softmax` max-shift branch convention is adapted from SciPy 1.17.0
  `scipy/special/_logsumexp.py` (BSD-3-Clause), while the surrounding fixed-rank
  Fortran implementation is original. SciPy public outputs and mpmath high-precision
  values are used as regression references. SciPy 1.17.0 public outputs are also
  used to validate the rank-15 real/complex axis and normalization paths; no
  additional SciPy implementation source was copied, translated, or adapted.


- `src/stats/nakagami.f90`, `src/stats/powernorm.f90`,
  `src/stats/loggamma.f90`, and `src/stats/wald.f90`: original SciFort
  implementations from the documented probability definitions and standard
  relationships to the gamma, normal, and inverse-Gaussian distributions.
  SciPy 1.17.0 is used only as an independent numerical validation reference;
  no SciPy implementation source was copied, translated, or adapted.
- `src/stats/gompertz.f90`, `src/stats/invweibull.f90`,
  `src/stats/betaprime.f90`, and `src/stats/burr12.f90`: original SciFort
  implementations from the published mathematical definitions of these
  distributions. Beta-prime reuses SciFort's independently implemented
  incomplete-beta kernel; the other three use project-local stable elementary
  helpers. SciPy 1.17.0 is used only as an independent numerical validation
  reference; no SciPy implementation source was copied, translated, or adapted.

- `src/stats/genhalflogistic.f90`, `src/stats/exponpow.f90`,
  `src/stats/exponweib.f90`, and `src/stats/powerlognorm.f90`: original
  SciFort implementations from the documented probability definitions and
  standard transform relationships. The implementations use project-local
  stable elementary and normal-distribution kernels. SciPy 1.17.0 is used
  only as an independent numerical validation reference; no SciPy
  implementation source was copied, translated, or adapted.


- `src/stats/levy_l.f90` and `src/stats/weibull_max.f90`: original SciFort
  implementations using reflection identities over SciFort's independently
  implemented Levy and Weibull-min kernels. SciPy 1.17.0 is used only as an
  independent numerical validation reference; no SciPy implementation source
  was copied, translated, or adapted.
- `src/stats/rdist.f90`: original SciFort implementation using the standard
  affine relationship between the R-distribution and a symmetric beta law;
  it reuses SciFort's independently implemented incomplete-beta and inverse
  incomplete-beta kernels. No third-party implementation source was copied or
  translated.
- `src/stats/skewcauchy.f90`: original SciFort implementation from the public
  skew-Cauchy probability definition, using independently selected stable
  arctangent and logarithmic rearrangements. SciPy 1.17.0 is used only for
  numerical validation; no SciPy implementation source was copied, translated,
  or adapted.
- `test/test_reflected_shape_four.f90`: contains only independently generated
  SciPy 1.17.0 numerical outputs for left-Levy, Weibull maximum,
  R-distribution, and skew-Cauchy regression cases; it contains no SciPy
  implementation source.

- `src/stats/dgamma.f90`: original SciFort implementation of the symmetric
  reflected-gamma definition, reusing SciFort's independently implemented
  incomplete-gamma kernels. No SciPy implementation source was copied or adapted.
- `src/stats/laplace_asymmetric.f90`: original SciFort implementation from the
  public piecewise asymmetric-Laplace density and its integrated probabilities.
  No SciPy implementation source was copied or adapted.
- `src/stats/truncnorm.f90`: original SciFort finite-bound implementation of a
  normal law truncated between standardized limits `a < b`, using SciFort's
  normal log-tail kernels and an independently written bisection inverse.
- `src/stats/loguniform.f90`: original SciFort implementation of the reciprocal
  density `1/(x*log(b/a))` on standardized support `[a,b]`.
- `test/test_scipy_transform_four.f90` and
  `test/scipy_transform_reference.f90`: contain only independently generated
  SciPy 1.17.0 numerical outputs for regression validation; no SciPy
  implementation source is included. The reference module is regenerated by
  `tools/generate_scipy_transform_reference.py`.

- `src/stats/foldnorm.f90` and `src/stats/foldcauchy.f90`: original SciFort
  implementations of folded normal and folded Cauchy laws from their public
  reflection identities. Stable mixture/arctangent rearrangements and
  bisection inverses were independently selected; no SciPy source was copied.
- `src/stats/recipinvgauss.f90`: original SciFort implementation using the
  reciprocal relationship to SciFort's independently implemented
  inverse-Gaussian law; no third-party implementation source was copied.
- `src/stats/truncpareto.f90`: original SciFort implementation from the
  normalized truncated-Pareto density on standardized support `[1,c]`, with
  independently selected log-domain normalizer and inverse formulas.
- `test/test_folded_reciprocal_four.f90` and
  `test/folded_reciprocal_reference.f90`: contain only independently generated
  SciPy 1.17.0 numerical outputs for regression validation. The reference module
  is regenerated by `tools/generate_folded_reciprocal_reference.py`; no SciPy
  implementation source is included.

- `src/stats/exponnorm.f90`: original SciFort implementation of the
  exponentially modified normal law from the convolution/normal-tail
  identities, with independently derived direct log-tail rearrangements and
  safeguarded bisection inverses. SciPy 1.17.0 was used only for public
  parameterization and numerical validation; no implementation source was
  copied or adapted.
- `src/stats/johnsonsb.f90` and `src/stats/johnsonsu.f90`: original SciFort
  implementations from the standard Johnson normalizing transformations
  `a+b*logit(z)` and `a+b*asinh(z)`, reusing SciFort's normal and logistic
  kernels. No third-party implementation source was copied.
- `src/stats/trapezoid.f90`: original SciFort piecewise implementation of the
  trapezoidal density, cumulative law, and exact inverse on standardized
  `[0,1]`. SciPy documentation was used only to confirm public shape
  conventions and numerical parity.
- `test/test_normal_transform_four.f90` and
  `test/normal_transform_reference.f90`: contain only independently generated
  SciPy 1.17.0 numerical outputs for regression validation. The reference
  module is regenerated by `tools/generate_normal_transform_reference.py`; no
  SciPy implementation source is included.

- `src/stats/burr.f90`: original SciFort implementation of the Burr Type III
  law from its public density and cumulative definitions, with independently
  derived log-domain tails and inverse formulas. SciPy 1.17.0 is used only to
  confirm public parameterization and numerical parity; no SciPy implementation
  source was copied, translated, or adapted.
- `src/stats/mielke.f90`: original SciFort implementation of the Mielke
  beta-kappa/Dagum law from its public density, with independently derived stable
  probability and quantile formulas. Its documented equivalence to Burr Type III
  is used as an independent regression identity. No third-party source was copied.
- `src/stats/gibrat.f90`: original SciFort implementation using the mathematical
  identity between the Gibrat law and SciFort's independently implemented
  lognormal distribution with shape one. No third-party implementation source was
  copied or adapted.
- `src/stats/wrapcauchy.f90`: original SciFort implementation of the wrapped
  Cauchy law on `[0,2*pi]` from the public density and circular CDF identities,
  with independently selected `atan2` and stable denominator rearrangements.
  SciPy 1.17.0 is used only for numerical validation.
- `test/test_burr_mielke_gibrat_wrapcauchy.f90` and
  `test/burr_mielke_gibrat_wrapcauchy_reference.f90`: contain only independently
  generated SciPy 1.17.0 numerical outputs and SciFort-owned regression logic. The
  reference module is regenerated by
  `tools/generate_burr_mielke_gibrat_wrapcauchy_reference.py`; no SciPy
  implementation source is included.

- `src/stats/genextreme.f90`: original SciFort implementation of the
  generalized extreme-value law from its public CDF and density equations,
  following SciPy's documented sign convention for the shape parameter. Stable
  logarithmic tails and the quantile transform were derived independently.
- `src/stats/kappa3.f90`: original SciFort implementation of the kappa-3 law
  from its public mathematical definition, with independently derived log-domain
  probability and inverse formulas. SciPy 1.17.0 is used only for parameter and
  numerical parity checks.
- `src/stats/kappa4.f90`: original SciFort implementation of the kappa-4 nested
  transform, including its shape-dependent support and exponential limiting
  cases. The implementation was derived from the public mathematical definition;
  no third-party implementation source was copied or adapted.
- `src/stats/truncweibull_min.f90`: original SciFort implementation of a Weibull
  minimum distribution conditioned to finite standardized bounds `a` and `b`,
  with direct log-domain normalization, tails, and inverse formulas.
- `test/test_genextreme_kappa_truncweibull.f90` and
  `test/genextreme_kappa_truncweibull_reference.f90`: contain only independently
  generated SciPy 1.17.0 numerical outputs and SciFort-owned regression logic.
  The reference module is regenerated by
  `tools/generate_genextreme_kappa_truncweibull_reference.py`; no SciPy
  implementation source is included.

- `src/stats/gengamma.f90`: original SciFort implementation of the generalized
  gamma law from its public density and regularized incomplete-gamma
  representation. SciPy 1.17.0 source and documentation were inspected to
  confirm parameterization, sign conventions, and tail selection; no SciPy
  implementation source was copied, translated, or adapted.
- `src/stats/halfgennorm.f90`: original SciFort implementation of the positive
  half generalized-normal law from its public density and incomplete-gamma
  CDF. Quantile and upper-tail formulas were derived from the same mathematical
  representation. No third-party implementation source was copied.
- `src/stats/argus.f90`: original SciFort implementation of the ARGUS law from
  its public density and normalization identity. The incomplete-gamma
  normalization and finite-support inverse were implemented independently;
  SciPy 1.17.0 was used to confirm numerical parity and public conventions.
- `src/stats/erlang.f90`: original SciFort specialization of the existing
  SciFort gamma kernels. It mirrors SciPy's numerical acceptance of positive
  noninteger shape values but does not reproduce the Python warning layer.
- `test/test_gengamma_halfgennorm_argus_erlang.f90` and
  `test/gengamma_halfgennorm_argus_erlang_reference.f90`: contain SciFort-owned
  regression logic and independently generated SciPy 1.17.0 numerical outputs.
  The reference module is regenerated by
  `tools/generate_gengamma_halfgennorm_argus_erlang_reference.py`; no SciPy
  implementation source is included.

- `src/stats/crystalball.f90`: original SciFort implementation of the Crystal
  Ball law from its public piecewise density. Normalization, direct logarithmic
  tails, and closed-form inverse branches were independently derived. SciPy
  1.17.0 was used to confirm parameterization and numerical behavior.
- `src/stats/jf_skew_t.f90`: original SciFort implementation of the Jones-Faddy
  skew-t law using its public beta-transform representation. Incomplete-beta
  tails and inverse-beta quantiles reuse SciFort kernels; no third-party
  implementation source was copied or adapted.
- `src/stats/pearson3.f90`: original SciFort implementation of Pearson III by
  mapping the public skewness parameterization to SciFort gamma and normal
  kernels. SciPy 1.17.0 was inspected to confirm its small-skew normal
  transition and public conventions; the implementation is independently
  written.
- `src/stats/rel_breitwigner.f90`: the density, stable survival form, inverse
  solver, likelihood, score, fit, and RNG integration are SciFort-owned. The
  complex CDF antiderivative is adapted from `rel_breitwigner_gen._cdf` in
  SciPy 1.17.0, `scipy/stats/_continuous_distns.py`,
  https://github.com/scipy/scipy/blob/v1.17.0/scipy/stats/_continuous_distns.py.
  Upstream copyright is Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy
  Developers, under BSD-3-Clause. The notice is retained in the source file and
  the complete applicable license text is in `THIRD_PARTY_LICENSES.md`. The
  adaptation replaces NumPy complex operations with Fortran intrinsic complex
  arithmetic, adds explicit loc/scale and endpoint handling, and derives a
  reciprocal-argument survival form for upper-tail stability.
- `test/test_crystalball_jf_pearson_breitwigner.f90` and
  `test/crystalball_jf_pearson_breitwigner_reference.f90`: contain SciFort-owned
  regression logic and generated SciPy 1.17.0 numerical outputs. The reference
  module is regenerated by
  `tools/generate_crystalball_jf_pearson_breitwigner_reference.py`; no SciPy
  implementation source is included in the generated reference data.


- `src/stats/genexpon.f90`: original SciFort implementation of the generalized
  exponential law from its public density and survival representation. The
  principal-real Lambert-W inverse is independently implemented from the
  defining equation with safeguarded Halley/Newton-style refinement; no
  third-party Lambert-W source was copied or adapted. SciPy 1.17.0 was used to
  confirm public parameterization and numerical outputs.
- `src/stats/skewnorm.f90`: original SciFort skew-normal implementation from
  `2*phi(z)*Phi(a*z)` and Owen's T integral. The fixed Gauss-Legendre Owen-T
  quadrature and transformed extreme-lower-tail quadrature are independently
  written from the integral definitions; no SciPy implementation source was
  copied or adapted.
- `src/stats/tukeylambda.f90`: original SciFort implementation from the public
  Tukey-lambda quantile definition. CDF inversion, the small-shape expansion,
  support handling, density derivative, and direct symmetric ISF were derived
  independently.
- `src/stats/rice.f90`: original SciFort implementation of the Rice law. The
  scaled modified-Bessel series/asymptotics and the Poisson/incomplete-gamma
  representation of the noncentral-chi-square tails are independently written
  from standard mathematical identities. No third-party Bessel or noncentral
  chi-square implementation source was copied or adapted.
- `test/test_genexpon_skewnorm_tukeylambda_rice.f90` and
  `test/genexpon_skewnorm_tukeylambda_rice_reference.f90`: contain SciFort-owned
  regression logic and generated SciPy 1.17.0 numerical outputs. The reference
  module is regenerated by
  `tools/generate_genexpon_skewnorm_tukeylambda_rice_reference.py`; no SciPy
  implementation source is included in the generated data.

- `src/stats/dpareto_lognorm.f90`: original SciFort implementation of the
  double-Pareto lognormal law from its published density/CDF identities. Tail
  expressions are evaluated in log space and quantiles are inverted on
  `log(x)` to preserve the heavy tails. SciPy 1.17.0 was used only to verify
  public parameterization and numerical outputs; no SciPy source was copied.
- `src/stats/vonmises.f90` and `src/stats/vonmises_line.f90`: original SciFort
  implementations from the public von Mises density. The line CDF is computed
  by independently written adaptive quadrature and the circular CDF extends it
  periodically. Quantiles use safeguarded bisection. No third-party CDF or
  quantile source was copied or adapted.
- `src/stats/kstwobign.f90`: original SciFort implementation of the asymptotic
  two-sided Kolmogorov law using complementary Jacobi-theta and alternating
  Kolmogorov series derived from the standard mathematical identities. The
  density and score derivative are differentiated forms of those series.
- `src/special/elementary.f90` (`i0e`, `i1e`): the independently written scaled
  modified-Bessel series/asymptotic kernels previously private to Rice were
  generalized to real arguments and promoted to the shared special-function
  layer. No external Bessel implementation was copied or adapted.
- `src/special/bessel_classic.f90`: original SciFort compatibility wrappers for
  the common real order-zero/one Bessel functions. `j0`, `j1`, `y0`, and `y1`
  use the Fortran standard `bessel_j0`, `bessel_j1`, `bessel_y0`, and
  `bessel_y1` intrinsics. `i0` and `i1` reconstruct the existing independently
  implemented `i0e`/`i1e` kernels in the log domain, while `k0`, `k0e`, `k1`,
  and `k1e` exponentiate the existing SciFort log-K kernels with explicit
  overflow handling. SciPy 1.17.0 was used for public naming, domain behavior,
  and ordinary numerical regression values; mpmath 1.3.0 supplied independent
  high-precision `i0(713)` and `i1(713)` references. No SciPy, Cephes, or other
  third-party implementation source was copied, translated, or adapted.
- `src/special/bessel_order.f90`: original SciFort compatibility wrappers for
  SciPy-style `yn`, `kn`, `iv`, `ive`, `kv`, and `kve`. Integer-order Y uses the
  Fortran standard `bessel_yn` intrinsic. Integer/real-order K and scaled K reuse
  SciFort's existing `log_besselk`/`log_besselk_scaled` kernels. Positive-order I
  and scaled I reuse `log_besseli_scaled`; negative integer orders use the standard
  symmetry and negative noninteger orders use the connection identity DLMF 10.27.3,
  combined in signed log space. SciPy 1.17.0 was used only for public behavior and
  numerical regression values. No SciPy, AMOS, Cephes, or other third-party
  implementation source was copied, translated, or adapted.
- `src/special/spherical_bessel.f90`: original SciFort implementation of the real
  spherical Bessel functions `spherical_jn`, `spherical_yn`, `spherical_in`, and
  `spherical_kn` and their first derivatives. It was independently derived from
  NIST DLMF 10.49 (explicit half-integer forms and modified-function connection),
  10.51 (recurrences and derivative identities), and 10.53 (power series).
  Modified I evaluation reuses SciFort's independently implemented real-order
  `log_besseli_scaled` kernel. SciPy 1.17.0 was used only for public API behavior
  and ordinary numerical regression values; mpmath 1.3.0 supplied independent
  high-precision checks for tiny derivatives and high-order finite K values. No
  SciPy, AMOS, Cephes, or other third-party implementation source was copied,
  translated, or adapted.
- `src/special/bessel_jy.f90`: adapted from the California Institute of
  Technology/JPL MATH77 routines `DBESJN`, `DBESYN`, and `DBESPQ`, retrieved
  from `https://www.netlib.org/math/MATH77/dbesjn.f`,
  `https://www.netlib.org/math/MATH77/dbesyn.f`, and
  `https://www.netlib.org/math/MATH77/dbespq.f` on 2026-09-23. The upstream
  copyright is Copyright (c) 1996 California Institute of Technology, Pasadena,
  California, with Government Sponsored Research NAS7-03001. The routines
  record contributions by E. W. Ng, W. V. Snyder, S. Singletary, C. L. Lawson,
  S. Chan, F. T. Krogh, and other MATH77 contributors. Upstream is under the
  BSD-style three-clause MATH77 license at `https://netlib.org/math/license.html`;
  the complete notice and disclaimer are retained in the adapted source and in
  `THIRD_PARTY_LICENSES.md`. SciFort converted fixed-form sequence routines to
  a free-form Fortran 2018 module with pure elemental scalar interfaces, replaced
  MATH77 machine/error services with Fortran intrinsics and IEEE return values,
  replaced the very-small-x path with the equivalent standard J power series,
  and also uses that same defining series for exact integer orders when `x` is
  small relative to the order so platform `BESSEL_JN` recurrences cannot overflow
  internally in a mathematically tiny regime. It retained/adapted the middle
  Miller/Neumann normalization and large-x P/Q
  asymptotics, and added negative-order reflection plus SciPy-style domain
  wrappers. SciPy 1.17.0 supplies ordinary regression/API references; mpmath
  1.3.0 independently checks extreme reflection cases where SciPy/AMOS can
  overflow or underflow prematurely.
- `src/special/gamma_family.f90`: original SciFort implementation of real and
  complex `gamma`, `rgamma`, and `loggamma`, together with real `gammasgn` and
  `multigammaln`. Real Gamma magnitudes use SciFort's independently implemented
  `gammaln` reflection kernel and real sign intervals. Complex `loggamma` is an
  independent implementation from NIST DLMF 5.5.3 recurrence/reflection and
  5.11.1 Stirling expansion, with explicit sheet correction so analytic
  continuation from the positive axis has the single negative-real branch cut.
  Complex `gamma`/`rgamma` are formed from that branch with guarded exponential
  scaling. `multigammaln` uses the defining multivariate-gamma product. SciPy
  1.17.0/1.18.0 public documentation and numerical outputs were consulted for
  naming, branch/pole/signed-zero behavior, and regression references; mpmath
  1.3.0 supplies independent high-precision checks in extreme real regions. No
  SciPy, Hare, or other third-party implementation source was copied, translated,
  or adapted.
- `src/special/complex_digamma.f90`: original SciFort complex `digamma`/`psi`
  implementation from NIST DLMF 5.5.2 (recurrence), 5.5.4 (reflection), and
  5.11.2 (asymptotic expansion), following the computation guidance in DLMF
  5.21. The complex cotangent in the reflection term is evaluated with scaled
  trigonometric/hyperbolic identities to preserve finite near-pole results and
  avoid intermediate overflow. SciPy 1.17.0/1.18.0 public documentation and
  numerical outputs were consulted only for API, pole/infinity conventions, and
  regression references; no SciPy or other third-party implementation source was
  copied, translated, or adapted.
- `src/special/gamma_combinatorics.f90`: original SciFort implementation of
  real `beta`, generalized `binom`, `poch`, and integer-order `polygamma`.
  Gamma/beta/Pochhammer/binomial relations follow the standard NIST DLMF
  definitions and recurrence identities; polygamma uses the Hurwitz-zeta
  identity and recurrence/reflection relations from DLMF section 5.15 over
  SciFort's existing independent Hurwitz-zeta kernel. Exact integer and
  half-integer phase values are used to avoid avoidable reflection cancellation.
  SciPy 1.17.0 was consulted for public naming, real-axis pole conventions, and
  ordinary numerical references; mpmath 1.3.0 supplied independent
  high-precision checks for difficult negative-axis polygamma values. No SciPy
  or other third-party implementation source was copied, translated, or adapted.
- `src/special/factorials.f90`: original SciFort implementation of floating
  `factorial`, `factorial2`, `factorialk`, `comb`, and `perm`. Integer
  multifactorials and integer permutations use their defining finite products
  with binary64 overflow guards. Real factorial uses the identity
  `Gamma(n+1)=Gamma(1+f)*product(f+j)` with `f` the fractional part of `n`,
  anchored by SciFort's existing `gamma` kernel and guarded against binary64
  overflow; combinations and permutations reuse the generalized `binom` and
  `poch` kernels.
  SciPy 1.17.0 was consulted for public naming, real zero-extension/domain
  conventions, repetition behavior, and ordinary numerical references. The
  implementation corresponds to SciPy's floating `exact=False` interfaces;
  Python arbitrary-precision `exact=True` behavior is not reproduced. No SciPy
  or other third-party implementation source was copied, translated, or adapted.
- `src/special/exponential_integrals.f90`: mixed original/adapted SciFort module.
  The `sici` Chebyshev kernels adapt the California Institute of Technology/JPL
  MATH77 routines `DSI`, `DCI`, and `DCPVAL`, retrieved on 2026-09-24 from
  `https://www.netlib.org/math/MATH77/dsi.f`,
  `https://www.netlib.org/math/MATH77/dci.f`, and
  `https://www.netlib.org/math/MATH77/dcpval.f`. The upstream copyright is
  Copyright (c) 1996 California Institute of Technology, Pasadena, California,
  based on Government Sponsored Research NAS7-03001; the existing complete
  BSD-style MATH77 notice is retained in this source and in
  `THIRD_PARTY_LICENSES.md`. SciFort converted the fixed-form functions to one
  free-form Fortran 2018 module, removed MATH77 diagnostics/machine services,
  exposes a pure elemental two-output `sici` wrapper, and extends finite negative
  arguments through the real even Ci branch. `exp1`, `expi`, `expn`, and `shichi`
  are original SciFort implementations from the defining series, generalized-E_n
  continued fraction, positive-Ei asymptotic expansion, and standard Ei/E1
  hyperbolic-integral identities. SciPy 1.17.0 supplies public API behavior and
  ordinary regression values; mpmath 1.3.0 supplies high-precision large-positive
  Ei checks. No SciPy, Cephes, or Fullerton implementation source was copied for
  those original procedures.
- `src/special/complex_exponential_integrals.f90`: original SciFort implementation
  of complex `exp1` and `expi`. The implementation is independently derived from
  NIST DLMF 6.2.1--6.2.7 (definitions/interrelations through the entire Ein
  function) and DLMF 6.9.1 (the E1 continued fraction). Small and near-cut
  arguments use the convergent Ein power series; the complementary region uses a
  modified-Lentz evaluation of the DLMF continued fraction with guarded complex
  exponential scaling. Exact real-axis values delegate to SciFort's existing real
  kernels, with explicit signed-zero handling for SciPy-compatible branch-cut
  semantics. SciPy 1.17.0 numerical outputs and SciPy 1.18.0 public documentation
  were used only for API/branch behavior and regression references. No SciPy or
  other third-party implementation source was copied, translated, or adapted.
- `src/special/airy.f90`: original SciFort implementation of real `airy`,
  `airye`, `itairy`, `ai_zeros`, and `bi_zeros`. The small-argument basis series is
  derived directly from the Airy differential equation `y'' = x*y`; the nonzero
  real branches use NIST DLMF 9.6.2--9.6.9 relations between Airy and Bessel
  functions, reusing SciFort's existing `jv`, `ive`, and `kve` kernels. The
  very-large-positive scaled branch uses the leading DLMF 9.7 Airy asymptotic
  forms only after the Airy exponential argument would overflow binary64. The
  `itairy` implementation independently integrates the same ODE series near zero,
  reuses SciFort's existing 16-point Gauss-Legendre constants only as an internal
  moderate-range quadrature rule, and derives its tail antiderivative recurrence
  from the Airy ODE together with the definite limits and leading asymptotics in
  DLMF 9.10.4--9.10.12. Real zero tables use the leading DLMF 9.9(iv) zero phases
  as Newton starting values and refine them with `airy`; DLMF 9.9(v) and SciPy
  1.17.0 provide independent tabular regression checks. SciPy 1.17.0 was consulted
  for public naming, tuple ordering, real-axis scaling/domain behavior, and ordinary
  numerical references. mpmath 1.3.0 generated high-precision integral references
  in the transition region where SciPy 1.17.0 `itairy` loses several digits. No
  SciPy, Zhang/Jin, AMOS, Cephes, Fullerton, or other third-party implementation
  source was copied, translated, or adapted.
- `src/special/carlson_elliptic.f90`: adapts the real `DRC`, `DRD`, `DRF`,
  and `DRJ` duplication/Taylor kernels from Jacob Williams'
  `carlson-elliptic-integrals` release 1.0.0, retrieved on 2026-09-24 from
  `https://raw.githubusercontent.com/jacobwilliams/carlson-elliptic-integrals/1.0.0/src/carlson_elliptic_module.f90`.
  That project is a modern refactoring of the public-domain SLATEC Carlson
  routines and is distributed under the BSD-3-Clause license, retrieved from
  `https://raw.githubusercontent.com/jacobwilliams/carlson-elliptic-integrals/1.0.0/LICENSE`;
  the full license is retained in the adapted source and in
  `THIRD_PARTY_LICENSES.md`. SciFort fixes the implementation to binary64,
  removes diagnostic/error-status side effects, exposes pure elemental
  SciPy-style interfaces, and adds IEEE/domain behavior plus homogeneity
  normalization. The real negative-second-argument RC and negative-pole RJ
  Cauchy principal values are SciFort wrappers based on Carlson transformation
  identities; `elliprg` is an original SciFort wrapper using the standard
  Carlson RF/RD identity. The complex RC/RF/RD duplication extensions and
  branch-preserving complex RG wrapper are original SciFort implementations
  based on the Carlson definitions and duplication theorems summarized in NIST
  DLMF Chapter 19. Complex RJ is independently implemented from its defining
  Carlson integral after the substitution `t=(u/(1-u))**2`, using SciFort's
  existing 16-point Gauss-Legendre constants in an adaptive interval scheme;
  it does not translate the legacy real RJ recurrence. SciPy 1.17.0 was
  consulted for public naming, real and complex domain/limit conventions,
  principal-value behavior, sufficient complex-RJ domain checks, and ordinary
  numerical references; mpmath 1.3.0 independently checks real
  homogeneity/extreme-scale values. No SciPy implementation source was copied
  or adapted.
- `src/special/legendre_elliptic.f90`: original SciFort implementation of
  real `ellipk`, `ellipkm1`, `ellipe`, `ellipkinc`, and `ellipeinc`. The
  complete and incomplete formulas are independently implemented from NIST
  DLMF 19.25 connection formulas using SciFort's existing Carlson RF/RD/RG
  kernels. The complete-E complementary-parameter expansion follows DLMF
  19.12.2. SciPy 1.17.0 was consulted for public naming, parameterization,
  real-domain behavior, periodic conventions, and ordinary numerical
  references; mpmath 1.3.0 was used for high-precision tiny-complementary
  and extreme-parameter checks. No SciPy, Cephes, Boost, or other
  third-party implementation source was copied, translated, or adapted.
- `src/special/jacobi_elliptic.f90`: original SciFort implementation of real
  `ellipj`. The arithmetic-geometric-mean inversion and amplitude recurrence
  are independently implemented from NIST DLMF 22.20, with SciFort period
  reduction and exact `m=0`/`m=1` endpoint formulas. SciPy 1.17.0 was
  consulted for public naming, output ordering, the real `0 <= m <= 1` domain,
  and ordinary numerical references. The SciPy/xsf v0.1.0 Cephes `ellpj`
  source was inspected to confirm the pinned SciPy behavior and its documented
  near-`m=1` approximation limits, but no Cephes/SciPy implementation source
  was copied, translated, or adapted. mpmath 1.3.0 supplies independent
  high-precision checks for large arguments and parameters near one.
- `src/special/bessel_derivatives.f90`: original SciFort implementation of the
  real-axis Bessel derivatives `jvp`, `yvp`, `ivp`, and `kvp`. The finite
  shifted-order sums are independently implemented from NIST DLMF 10.6.7
  (ordinary cylindrical J/Y derivatives) and DLMF 10.29.5 (modified I/K
  derivatives). SciFort uses normalized binomial weights and a scaled finite
  accumulator so neither raw binomial coefficients nor an avoidable partial
  sum need overflow. SciPy 1.17.0 was consulted for public naming, default
  derivative-order behavior, real-axis domain parity, and numerical reference
  values. No SciPy or third-party implementation source was copied, translated,
  or adapted.
- `test/test_dpareto_vonmises_kstwobign.f90` and
  `test/dpareto_vonmises_kstwobign_reference.f90`: contain SciFort-owned
  regression logic and generated SciPy 1.17.0 numerical outputs. The reference
  module is regenerated by
  `tools/generate_dpareto_vonmises_kstwobign_reference.py`; no SciPy
  implementation source is included in the generated data.

- `src/stats/irwinhall.f90`: original SciFort implementation of the Irwin-Hall
  law. The density uses the standard cardinal-B-spline recurrence and the CDF
  the corresponding accumulated B-spline identity, chosen independently to
  avoid cancellation in the alternating polynomial representation. SciPy
  1.17.0 was used only to confirm public parameterization and numerical output.
- `src/stats/ksone.f90`: original SciFort implementation of the one-sided
  finite-sample Kolmogorov-Smirnov law from the classical
  Birnbaum-Tingey/Smirnov finite-sample sum. The logarithmically scaled survival
  and differentiated density/score sums were independently written; no SciPy
  implementation source was copied or adapted.
- `src/stats/ncx2.f90`: original SciFort noncentral chi-square implementation
  from its Poisson mixture of central chi-square laws. Mixture centering,
  direct incomplete-gamma tails, analytic score accumulation, and safeguarded
  quantile inversion were independently implemented.
- `src/stats/ncf.f90`: original SciFort noncentral F implementation from the
  Poisson mixture of beta-transformed central F components. Direct
  incomplete-beta tails, analytic mixture scores, and quantile inversion were
  independently implemented.
- `test/test_irwinhall_ksone_ncx2_ncf.f90` and
  `test/irwinhall_ksone_ncx2_ncf_reference.f90`: contain SciFort-owned
  regression logic and generated SciPy 1.17.0 numerical outputs. The reference
  module is regenerated by `tools/generate_irwinhall_ksone_ncx2_ncf_reference.py`;
  no SciPy implementation source is included in the generated data.

- `src/stats/randint.f90`: original SciFort implementation of the discrete
  uniform law on `low, ..., high-1`, including direct discrete tails and exact
  integer quantile correction. No SciPy implementation source was copied or
  translated.
- `src/stats/planck.f90`: original SciFort implementation of the Planck
  discrete exponential law from its geometric mass formula, using stable
  `expm1`/`log1p` identities for tails and normalization.
- `src/stats/dlaplace.f90`: original SciFort implementation of the symmetric
  discrete Laplace law from its closed-form mass and half-line tail sums.
- `src/stats/logser.f90`: original SciFort implementation of the logarithmic
  series law from its defining mass function. Lower and upper tail summations,
  recurrence termination, and integer quantile search were independently
  implemented.
- `test/test_randint_planck_dlaplace_logser.f90` and
  `test/randint_planck_dlaplace_logser_reference.f90`: contain SciFort-owned
  regression logic and generated SciPy 1.17.0 numerical outputs. The reference
  module is regenerated by
  `tools/generate_randint_planck_dlaplace_logser_reference.py`; no SciPy
  implementation source is included in the generated data.

- `src/stats/betabinom.f90`: original SciFort implementation of the
  beta-binomial law from the beta-binomial mass formula. Finite lower and upper
  tails are summed independently in log space and integer quantiles use binary
  search; no SciPy implementation source was copied or translated.
- `src/stats/hypergeom.f90`: original SciFort implementation of the
  hypergeometric law from products of binomial coefficients, evaluated through
  log-gamma identities with independently accumulated finite tails.
- `src/stats/nhypergeom.f90`: original SciFort implementation of the negative
  hypergeometric law from its combinatorial mass formula, with independent
  finite-tail log summation and integer quantile search.
- `src/stats/boltzmann.f90`: original SciFort implementation of the finite
  truncated-geometric/Boltzmann law. Stable `expm1`/`log1p` identities are used
  for normalization and both tails, including the small-rate regime.
- `test/test_betabinom_hypergeom_nhypergeom_boltzmann.f90` and
  `test/betabinom_hypergeom_nhypergeom_boltzmann_reference.f90`: contain
  SciFort-owned regression logic and generated public SciPy 1.17.0 numerical
  outputs. The generator additionally computes two stress tails directly from
  exact/high-precision defining formulas where SciPy 1.17.0 is less accurate;
  no SciPy implementation source is included.


- `src/special/zeta.f90`: original SciFort implementation of real Riemann and
  Hurwitz zeta functions. The positive-q Hurwitz kernel and its exponent
  derivative use the Euler-Maclaurin representation in DLMF 25.11.43; the
  Riemann continuation uses the reflection formula DLMF 25.4.2 with reduced
  `sin(pi*s/2)` phase evaluation. The SciPy-style two-argument real wrapper
  additionally handles negative noninteger `q` for integer `s` by an exact
  recurrence split into positive-q Hurwitz tails. `zetac` uses the `q=2`
  Hurwitz identity and binary exponent scaling for very large positive `s`.
  SciPy 1.17.0 supplied API/domain and ordinary numerical references; mpmath
  supplied independent high-precision checks. No SciPy source was copied or
  adapted.
- `src/special/spence.f90`: original SciFort real/complex Spence/dilogarithm
  implementation. The unchanged real path evaluates `Li_2(1-x)` using DLMF
  25.12.1, 25.12.3, 25.12.4, and 25.12.6. The complex extension uses the same
  principal dilogarithm definition and DLMF 25.12.4 inversion, plus the
  `s -> 2` logarithmic expansion associated with DLMF 25.12.12 in a narrow
  annulus around the unit circle. The 24 Bernoulli coefficients are the exact
  coefficients `-B_(2k)/[(2k)(2k+1)(2k)!]` rounded to binary64. SciPy 1.17.0
  supplied API/branch-cut compatibility references and mpmath 1.3.0 supplied
  independent high-precision checks. No SciPy or other implementation source
  was copied or adapted.

- `src/special/bessel_zeros.f90`: original SciFort implementation of positive
  real integer-order Bessel J/J'/Y/Y' zero tables and the magnitude-ordered
  `jnjnp_zeros` J/J' aggregate, plus real `y0_zeros`, `y1_zeros`, and
  `y1p_zeros` companion-value helpers. The design follows the simple-positive-zero,
  interlacing, and order behavior documented in NIST DLMF section 10.21 and uses
  bracketed safeguarded Newton refinement over SciFort's existing `jv`, `yv`,
  `jvp`, and `yvp` kernels. The aggregate builds a finite candidate grid,
  certifies that the requested cutoff precedes the first omitted order/root,
  and coalesces the exact identity `J_0'=-J_1` before deterministic tie ordering.
  The Y helper companion values follow the derivative identities in DLMF 10.6,
  including `Y_0'=-Y_1`.
  SciPy 1.17.0 documentation and numerical outputs were consulted for the public
  names, metadata conventions, 1200-entry compatibility limit, positive-interval
  convention, negative-integer-order behavior, and reference tables. No SciPy,
  Zhang/Jin, or other Bessel-zero implementation source was copied or adapted.
- `src/stats/betanbinom.f90`: original SciFort implementation of the
  beta-negative-binomial mass from binomial/beta identities, with a
  probability-ratio lower-tail recurrence and discrete quantile search.
- `src/stats/yulesimon.f90`: original SciFort implementation from
  `P(X=k)=alpha B(k,alpha+1)` and the exact beta-function survival identity.
- `src/stats/zipf.f90`: original SciFort implementation from the zeta-law mass
  and Hurwitz-zeta tail identity.
- `src/stats/zipfian.f90`: original SciFort implementation of the finite
  generalized-harmonic power law with independently accumulated finite tails.
- `test/test_betanbinom_yulesimon_zipf_zipfian.f90` and
  `test/betanbinom_yulesimon_zipf_zipfian_reference.f90`: SciFort-owned test
  logic and generated public SciPy 1.17.0 numerical outputs. The reference
  module is regenerated by
  `tools/generate_betanbinom_yulesimon_zipf_zipfian_reference.py`; no SciPy
  implementation source is included.

- `src/special/bessel_k.f90`: original SciFort implementation of real
  modified-Bessel `K_nu(x)` in logarithmic form from its standard integral
  representation, with composite Gauss-Legendre quadrature after logarithmic
  scaling and a standard large-argument asymptotic expansion. The mathematical
  specifications were taken from NIST DLMF sections 10.32 (integral
  representations) and 10.40 (asymptotic expansions). Adaptive panel sizing,
  truncation rules, and the order/argument log-derivative implementations are
  original SciFort code; no SciPy or other numerical-library source was copied
  or adapted.
- `src/stats/geninvgauss.f90`: original SciFort implementation of the
  generalized inverse-Gaussian density from its defining Bessel-K formula.
  Direct lower/upper tails use a logarithmic-coordinate integral with localized
  monotone-tail truncation and log-sum-exp quadrature. SciPy 1.17.0 public
  behavior was used for API and ordinary-value parity only.
- `src/stats/norminvgauss.f90`: original SciFort implementation of the normal
  inverse-Gaussian density from its standard Bessel-K formula. Direct tails use
  the substitution `z=sinh(y)`, localized log-domain quadrature, and direct
  log-tail inversion. SciPy 1.17.0 public behavior was used for API and
  ordinary-value parity only.
- `src/stats/skellam.f90`: original SciFort implementation of the difference of
  independent Poisson counts. The PMF uses the standard modified-Bessel-I
  identity and the tails use the standard noncentral-chi-square representation
  over SciFort's `ncx2` kernel. No SciPy source was copied or adapted.
- `test/test_gig_nig_skellam.f90` and
  `test/gig_nig_skellam_reference.f90`: SciFort-owned regression logic and
  generated public SciPy 1.17.0 numerical outputs. Additional deep-tail and
  wide-order Bessel-K constants were evaluated independently at 50-80 digit
  precision from the defining integral formulas. The SciPy reference module is
  regenerated by `tools/generate_gig_nig_skellam_reference.py`; no SciPy
  implementation source is included.

- `src/stats/genhyperbolic.f90`: original SciFort implementation of the
  generalized-hyperbolic density from its standard modified-Bessel-K formula.
  Direct CDF/SF evaluation uses the `z=sinh(y)` transformation and log-domain
  Gauss-Legendre tail integration. The exact `p<0, abs(b)=a` limiting
  normalizer is handled analytically. SciPy 1.17.0 public outputs were used only
  for API/numerical parity tests; no SciPy implementation source was copied or
  adapted.
- `src/stats/nchypergeom_fisher.f90`: original SciFort implementation of
  Fisher's noncentral hypergeometric law from its defining finite-support PMF,
  with log-sum-exp normalization and direct lower/upper sums. No code from
  SciPy's or BiasedUrn's implementations was copied, translated, or adapted.
- `test/test_genhyperbolic_fisher.f90` and
  `test/genhyperbolic_fisher_reference.f90`: SciFort-owned regression logic and
  generated public SciPy 1.17.0 numerical outputs. The reference module is
  regenerated by `tools/generate_genhyperbolic_fisher_reference.py`; no SciPy
  implementation source is included.

- `src/stats/nct.f90`: original SciFort implementation of the noncentral
  Student-t distribution from the representation `(Z+nc)/sqrt(2U/df)` with
  `U~Gamma(df/2,1)`. CDF, SF, and density are expectations in gamma-probability
  coordinates evaluated with a double-exponential endpoint map and adaptive
  spacing for small `df`. No SciPy/Boost implementation source was copied or
  adapted; SciPy 1.17.0 public outputs were used for ordinary parity only.
- `src/stats/gausshyper.f90`: original SciFort implementation of the Gauss
  hypergeometric density from its defining beta-weighted tilted density.
  Normalization and direct tails are integrated after a beta-quantile transform,
  avoiding direct endpoint-singular beta-density quadrature. No SciPy source was
  copied or adapted.
- `test/test_nct_gausshyper.f90` and `test/nct_gausshyper_reference.f90`:
  SciFort-owned regression logic and generated public SciPy 1.17.0 numerical
  outputs. Additional difficult-tail constants are independent defining-integral
  evaluations. The SciPy reference module is regenerated by
  `tools/generate_nct_gausshyper_reference.py`; no SciPy implementation source
  is included.

- `src/stats/landau.f90`: original SciFort implementation of the Landau
  distribution. The central density/CDF use the classical exponentially damped
  characteristic-function representation; the right tail uses a scaled form of
  Landau's nonoscillatory defining integral. The left-tail PDF/CDF asymptotic
  coefficients are the classical Koelbig-Schorr coefficients as tabulated in
  Takuma Yoshimura, "Numerical Evaluation and High Precision Approximation
  Formula for Landau Distribution" (2024), CC BY 4.0. No SciPy, Boost, or
  Yoshimura implementation source code was copied or adapted. SciPy 1.17.0
  public numerical outputs are used only for generated regression references.
- `test/test_landau.f90` and `test/landau_reference.f90`: SciFort-owned
  regression logic and generated public SciPy 1.17.0 numerical outputs. The
  reference module is regenerated by `tools/generate_landau_reference.py`.

- `src/stats/nchypergeom_wallenius.f90`: original SciFort implementation of
  Wallenius' noncentral hypergeometric distribution from the defining sequential
  biased-urn process. PMF/tails are obtained by exact finite-state Markov
  propagation, and the odds derivative is propagated through the same recursion.
  No SciPy or BiasedUrn source code was copied, translated, or adapted.
- `src/stats/poisson_binom.f90`: original SciFort implementation of the
  Poisson-binomial distribution from independent Bernoulli convolution. Direct
  probability/log-tail dynamic programs and leave-one-out convolutions provide
  the distribution and analytic vector score. No SciPy implementation source
  was copied or adapted.
- `test/test_wallenius_poisson_binom.f90` and
  `test/wallenius_poisson_binom_reference.f90`: SciFort-owned regression logic
  and generated public SciPy 1.17.0 numerical outputs. The reference module is
  regenerated by `tools/generate_wallenius_poisson_binom_reference.py`; only
  numerical reference values enter SciFort.

- `src/stats/kstwo.f90`: the finite-sample two-sided Kolmogorov-Smirnov
  regime selection, Durbin/Marsaglia-Tsang-Wang matrix CDF algorithm, and
  Pelz-Good large-sample expansion are translated/adapted from
  `scipy/stats/_ksstats.py` in SciPy 1.17.0, copyright Enthought, Inc. and the
  SciPy Developers, under BSD-3-Clause. SciFort omits SciPy's Pomeranz branch
  and instead extends the exact Durbin/MTW branch through the small-`n` region
  where both are accurate, while reusing SciFort's independent `ksone` kernel
  for the exact one-sided-tail regimes. The retained SciPy license is in
  `THIRD_PARTY_LICENSES.md`.
- `test/test_kstwo.f90` and `test/kstwo_reference.f90`: SciFort-owned
  regression logic and generated public SciPy 1.17.0 numerical outputs. The
  reference module is regenerated by `tools/generate_kstwo_reference.py`.

- `src/stats/levy_stable.f90`: the default S1 parameterization, S1-to-S0
  conversion, difficult-input rounding, analytic special cases, and
  Nolan/Zolotarev piecewise density/CDF formulas are translated/adapted from
  SciPy 1.17.0 `scipy/stats/_levy_stable/__init__.py` and
  `scipy/stats/_levy_stable/c_src/levyst.c`, copyright Enthought, Inc. and the
  SciPy Developers, under BSD-3-Clause. SciFort supplies its own fixed
  Gauss-Legendre and adaptive quadrature, tail inversion, likelihood, fitting,
  score, RNG, and ABI integration. The retained SciPy license is in
  `THIRD_PARTY_LICENSES.md`.
- `src/stats/studentized_range.f90`: original SciFort implementation from the
  defining normal-range/chi-square scale-mixture representation. The finite-df
  integrals correspond to the formulas described by Batista et al. (2017), and
  the infinite-df range integrals to Lund & Lund (1983). SciPy 1.17.0 was
  inspected for API semantics, the `df >= 100000` asymptotic switch, and
  numerical regression values; no SciPy implementation source is copied or
  adapted in this module.
- `test/test_levy_stable_studentized_range.f90` and
  `test/levy_stable_studentized_range_reference.f90`: SciFort-owned regression
  logic and generated public SciPy 1.17.0 numerical outputs. The reference
  module is regenerated by
  `tools/generate_levy_stable_studentized_range_reference.py`.

- `src/linalg.f90`: original SciFort compact dense-linear-algebra implementation.
  The user-supplied `Beliavsky/fortran-lapack` archive at commit
  `9b43beb7644479af46a56e743f61dbc8ae63ac7a` (BSD-3-Clause, copyright
  Federico Perini) was inspected for Cholesky interface/failure conventions and
  as a numerical-design reference. No fortran-lapack implementation source was
  copied, translated, or adapted; SciFort uses its own unblocked Cholesky,
  triangular solve, cyclic Jacobi eigensolver, and PSD wrapper.
- `src/stats/multivariate_integration.f90`: SciFort integration plumbing and
  randomized tent-transformed Halton rule around a conditional multivariate-box
  probability transform. The conditional normal/t transformation and the
  pivoted/permuted Cholesky logic are adapted from SciPy 1.17.0
  `scipy/stats/_qmvnt.py` and `scipy/stats/_qmvnt_cy.pyx` under SciPy's
  BSD-3-Clause license. SciFort does not copy SciPy's FFT/CBC lattice generator;
  the Halton sequence, randomized shifts, batching, stopping rules, singular-row
  handling, and explicit-state integration interface are SciFort-owned.
- `src/stats/multivariate_normal.f90`: original SciFort implementation of the
  multivariate-normal density, entropy, sampling, MLE, marginal-extraction, and
  public CDF/log-CDF wiring. The scalar/diagonal-covariance log-CDF path uses
  independence and logarithmic differences of the normal CDF/SF, derived from
  the error-function definition in DLMF 7.2.1 (https://dlmf.nist.gov/7.2.E1).
  This path and its identity-based regression tests are original SciFort code.
  SciPy 1.17.0 `scipy/stats/_multivariate.py` and
  `scipy/stats/_covariance.py` were inspected for dense-covariance API semantics,
  lower-triangle behavior, `_PSD` numerical rank cutoff, singular-support
  tolerance, fitting conventions, lower-limit semantics, and CDF defaults. The
  multidimensional box-probability kernel is the separately documented adapted
  code in `src/stats/multivariate_integration.f90`.
- `src/stats/multivariate_t.f90`: original SciFort implementation of the
  multivariate Student-t density/log-density, entropy, explicit-state sampling,
  marginal extraction, and public CDF wiring from the standard defining
  formulas. SciPy 1.17.0 was inspected for shape-matrix semantics, singular-matrix
  behavior, entropy asymptotics, default degrees of freedom, CDF point limits,
  and public method conventions. Its multidimensional CDF uses the separately
  documented adapted integration kernel.
- `src/stats/dirichlet.f90`, `src/stats/multinomial.f90`, and
  `src/stats/dirichlet_multinomial.f90`: original SciFort implementations from
  the standard defining probability formulas. SciPy 1.17.0 was inspected for
  parameter validation, the Dirichlet `K-1` coordinate convention, multinomial
  final-probability adjustment, entropy definition, and public method semantics.
- `test/test_multivariate.f90` and `test/multivariate_reference.f90`: SciFort-owned
  regression logic and generated public SciPy 1.17.0 numerical outputs. The
  reference module is regenerated by `tools/generate_multivariate_reference.py`.
- `src/stats/multivariate_hypergeom.f90`: original Fortran implementation of the
  standard multivariate-hypergeometric mass and moments. SciPy 1.17.0
  `scipy/stats/_multivariate.py` was inspected for parameter/support semantics and
  its sequential univariate-hypergeometric sampling decomposition; the translated
  sampling control flow is treated as a SciPy BSD-3-Clause adaptation. The
  retained SciPy license is in `THIRD_PARTY_LICENSES.md`.
- `src/stats/normal_inverse_gamma.f90`: original SciFort implementation from the
  standard normal/inverse-gamma joint density and moments. SciPy 1.17.0
  `scipy/stats/_multivariate.py` was inspected for defaults, existence conditions,
  support semantics, and the inverse-gamma/conditional-normal sampling definition.
- `src/stats/matrix_distribution_helpers.f90`, `src/stats/matrix_normal.f90`,
  `src/stats/wishart.f90`, and `src/stats/invwishart.f90`: original SciFort
  implementations of standard SPD, multivariate-gamma/digamma, matrix-normal,
  Wishart, and inverse-Wishart formulas. SciPy 1.17.0 `_multivariate.py` was
  inspected for lower-triangle semantics, parameter validation, moment/entropy
  conventions, and public API behavior. Wishart RNG uses the standard Bartlett
  construction; inverse-Wishart RNG uses a Wishart-inversion identity rather than
  SciPy's direct Axen implementation.
- `src/stats/matrix_t.f90`: SciFort density implementation from the standard
  two-sided matrix-t determinant formula. SciPy 1.17.0 `_multivariate.py` was
  inspected for the new public API, defaults, and the equivalent inverse-Wishart/
  matrix-normal mixture; the sampling decomposition is treated as a SciPy
  BSD-3-Clause adaptation. The retained SciPy license is in
  `THIRD_PARTY_LICENSES.md`.
- `test/test_matrix_distributions.f90` and
  `test/test_matrix_distributions_reference.f90`: SciFort-owned regression logic
  and generated public SciPy 1.17.0 numerical outputs. The reference module is
  regenerated by `tools/generate_matrix_distributions_reference.py`.
- `src/special/bessel_i.f90`: original SciFort implementation of the positive-term
  series for real-order modified Bessel I and its large-argument asymptotic
  expansion, based on DLMF 10.25.2 and 10.40.1. SciPy `special.ive` is used only
  for numerical regression values.
- `src/directional_stats.f90`: original SciFort implementation of directional
  mean and mean-resultant-length summaries from the standard normalized-vector
  definition. SciPy 1.17.0 `directional_stats` public behavior was used to define
  normalization, observation-axis, and zero-resultant semantics and to generate
  independent numerical regression values; no SciPy implementation code was copied
  or translated.
- `src/stats/uniform_direction.f90`: original normalized-Gaussian unit-sphere
  sampler. SciPy 1.17.0 `_multivariate.py` was inspected for public shape and
  parameter semantics; no SciPy implementation source is copied.
- `src/stats/random_matrix_helpers.f90`, `src/stats/ortho_group.f90`,
  `src/stats/special_ortho_group.f90`, and `src/stats/unitary_group.f90`: original
  SciFort Gaussian-QR/Haar implementations using re-orthogonalized modified
  Gram-Schmidt, based on the standard construction described by F. Mezzadri,
  Notices AMS 54 (2007), 592-604. SciPy 1.17.0 was used as an API/behavioral
  reference only.
- `src/stats/vonmises_fisher.f90`: original SciFort density, entropy, fit, and
  explicit-state sampling implementation. The dimension-specific sampling strategy
  follows J. Wenzel, "Numerically stable sampling of the von Mises Fisher
  distribution on S2", in dimension three and A. T. A. Wood, Communications in Statistics - Simulation and Computation 23 (1994),
  157-164, in dimensions four and above. SciPy 1.17.0 `_multivariate.py` was
  inspected for parameter validation and public behavior; no SciPy source code is
  copied.
- `src/stats/random_correlation.f90`: original Fortran implementation of the
  prescribed-eigenvalue correlation-matrix construction of P. I. Davies and
  N. J. Higham, BIT 40 (2000), 640-651. SciPy 1.17.0 `_multivariate.py` was
  inspected for tolerance defaults and validation semantics; no SciPy source code
  is copied.
- `src/stats/random_table.f90`: original fixed-margin contingency-table mass, mean,
  and exact row-conditional multivariate-hypergeometric sampler. SciPy 1.17.0 was
  inspected for public support semantics and numerical references; the sampler is
  intentionally independent of SciPy's compiled Boyett/Patefield implementations.
- `test/test_directional_random_matrix.f90` and
  `test/test_directional_random_matrix_reference.f90`: SciFort-owned regression
  logic and generated public SciPy 1.17.0 numerical outputs. The reference module
  is regenerated by `tools/generate_directional_random_matrix_reference.py`.

- `src/stats/gaussian_kde.f90`: original SciFort implementation of the
  standard Gaussian-mixture KDE formulas. SciPy 1.17.0 `scipy/stats/_kde.py` was
  inspected for public API behavior, data orientation, effective-sample-size and
  weighted-covariance conventions, bandwidth factors, integration semantics, and
  marginal behavior. No SciPy implementation source was copied or translated.
- `test/test_gaussian_kde.f90` and `test/gaussian_kde_reference.f90`: SciFort-owned
  regression logic and generated public SciPy 1.17.0 numerical outputs. The
  reference module is regenerated by `tools/generate_gaussian_kde_reference.py`.
- `src/stats/qmc.f90`: original Fortran implementation of standard QMC formulas and
  constructions. SciPy 1.17.0 `scipy/stats/_qmc.py`, `_qmc_cy.pyx`, and `_sobol.pyx`
  were inspected for public semantics, discrepancy formulas, Halton/Sobol indexing,
  LMS+digital-shift structure, strength-two OA-LHS constraints, QMC distribution
  transforms, and Poisson-disk behavior. The discrepancy equations, sequence engines,
  scrambling code, multinomial/normal transforms, Bridson sampler, Fisher-Yates
  permutations, Prim MST, and orthogonal-array construction were independently
  written in Fortran; no SciPy implementation source was copied or translated.
  Randomized paths intentionally use SciFort's explicit MRG32k3a state instead of
  NumPy RNGs. The `random-cd` proposal loop and efficient centered-discrepancy
  exchange update are adapted from SciPy 1.17.0 `scipy/stats/_qmc.py` under
  BSD-3-Clause, with one-based indexing, SciFort explicit RNG state, and an exact
  full-CD verification before accepting a predicted improvement. The underlying
  exchange formulas are attributed there to Ruichen Jin, Wei Chen, and Agus
  Sudjianto, "An efficient algorithm for constructing optimal design of computer
  experiments," Journal of Statistical Planning and Inference 134(1), 268-287
  (2005), doi:10.1016/j.jspi.2004.02.014. SciPy's BSD notice is retained in the
  source header and `THIRD_PARTY_LICENSES.md`.
- `src/stats/sensitivity_analysis.f90`: SciPy 1.17.0
  `scipy/stats/_sensitivity_analysis.py` was inspected and adapted under
  BSD-3-Clause for the public `sobol_indices` semantics, Saltelli 2010 estimator,
  mean normalization, A/B/AB sampling layout, and constant-output handling. The
  Fortran translation uses SciFort's existing Sobol engine and explicit MRG32k3a
  RNG state instead of NumPy's generator, represents independent marginal
  distributions with inverse-CDF callbacks, and keeps scalar-output indices as
  rank-two `(1,d)` arrays. SciPy's BSD notice is retained in the source header and
  `THIRD_PARTY_LICENSES.md`. The estimator is attributed to Saltelli et al.,
  "Variance based sensitivity analysis of model output. Design and estimator for
  the total sensitivity index," Computer Physics Communications 181(2), 259-270
  (2010), doi:10.1016/j.cpc.2009.09.018.
- `test/test_qmc.f90` and `test/qmc_reference.f90`: SciFort-owned regression logic
  with numerical references generated by `tools/generate_qmc_reference.py` using
  public SciPy 1.17.0 `sobol_indices` results for scalar and multi-output
  precomputed evaluations.
- `src/stats/qmc_lloyd.f90`: SciPy 1.17.0 `scipy/stats/_qmc.py` was inspected
  and adapted for the public Lloyd relaxation semantics: finite Voronoi-region
  centroids, the exponential decay sequence, unit-cube acceptance, and the
  city-block minimum-spacing stopping criterion. This portion is carried under
  SciPy's BSD-3-Clause terms. The N-dimensional Bowyer-Watson Delaunay
  triangulation, simplex circumcenter solver, facet bookkeeping, and finite
  Voronoi-vertex reconstruction are original SciFort code written from standard
  computational-geometry definitions; no Qhull source is copied or translated.
  Unlike SciPy's Qhull call with `QJ`, SciFort uses deterministic unjoggled
  geometry, so degenerate/cospherical point sets may choose a different valid
  triangulation or be rejected. The source header retains the SciPy copyright
  notice and points to `THIRD_PARTY_LICENSES.md`.
- `src/stats/qmc_sobol_data.f90`: generated Fortran representation of the Joe-Kuo
  search-criterion-6 Sobol direction-number initialization table distributed in
  SciPy 1.17.0 as `_sobol_direction_numbers.npz`. This data file is carried under
  SciPy's BSD-3-Clause terms; the original direction-number source is the UNSW
  Joe-Kuo table cited in SciPy's `_sobol.pyx`.
- `test/test_qmc.f90` and `test/qmc_reference.f90`: SciFort-owned regression logic
  and generated public SciPy 1.17.0 numerical outputs, regenerated by
  `tools/generate_qmc_reference.py`.

- `src/multiscale_graphcorr.f90`: adaptation of SciPy 1.17.0
  `scipy/stats/_mgc.py` and the MGC-specific distance-transform/local-correlation
  routines in `scipy/stats/_stats.pyx`, carried under SciPy's BSD-3-Clause terms.
  The Fortran port preserves column-wise dense distance ranks, MGC centering,
  cumulative local covariance/correlation maps, beta-thresholded largest
  4-connected significant region, smoothed optimal-scale selection, automatic
  two-sample reduction, and plus-one permutation p-values. It uses SciFort's
  explicit MRG32k3a state instead of NumPy random states and provides a native
  callback interface for custom distance matrices; process-parallel `workers` are
  not translated. SciPy's BSD notice is retained in the source header and
  `THIRD_PARTY_LICENSES.md`.
- `test/test_multiscale_graphcorr.f90` and
  `test/multiscale_graphcorr_reference.f90`: SciFort-owned regression logic and
  public SciPy 1.17.0 numerical outputs, regenerated by
  `tools/generate_multiscale_graphcorr_reference.py`.

- `src/distribution_infrastructure.f90`: native first-class continuous-distribution
  callback adapter plus generic truncation wrapper. SciPy 1.17.0
  `scipy/stats/_distribution_infrastructure.py` was inspected and adapted under
  BSD-3-Clause for the public `make_distribution`/`truncate` semantics, clipped
  support rules, normalized PDF/CDF/survival relations, and lower-/upper-tail inverse
  transforms. The pure-callback object design, frozen parameter vector, status-code
  API, fallback methods, and explicit-MRG32k3a inverse-transform sampling are original
  SciFort adaptations for Fortran. SciPy's BSD notice is retained in the source
  header and `THIRD_PARTY_LICENSES.md`.
- `test/test_distribution_infrastructure.f90` and
  `test/test_distribution_infrastructure_reference.f90`: SciFort-owned regression
  logic with public SciPy 1.17.0 numerical outputs generated by
  `tools/generate_distribution_infrastructure_reference.py` for truncated Normal,
  Rayleigh, and legacy-wrapped Uniform cases.

- `src/goodness_of_fit.f90`: SciFort implementation of goodness-of-fit and
  normality tests. The generalized Monte Carlo `goodness_of_fit` control flow,
  parameter-refitting semantics, AD/KS/CVM/Filliben statistic definitions, and
  plus-one p-value convention were implemented from the SciPy 1.17.0 public
  documentation and inspected BSD-3-Clause source in `scipy/stats/_fit.py`.
  The native Fortran implementation currently supports `norm`, `expon`, and
  `uniform` with explicit RNG state and independently written closed-form fits.
  The Shapiro-Wilk W/p-value kernel is adapted from SciPy 1.17.0
  `scipy/stats/_ansari_swilk_statistics.pyx` (Royston AS R94 / AS 181 lineage);
  the finite-sample Cramer-von Mises correction and Anderson-Darling critical
  value tables/public fitting conventions were adapted from SciPy 1.17.0
  `scipy/stats/_hypotests.py` and `_morestats.py`. The KS lattice probability,
  moment tests, optimization code, and remaining Fortran implementation were
  independently written. SciPy's BSD-3-Clause notice is retained in the source
  header and `THIRD_PARTY_LICENSES.md`.
- `test/test_goodness_of_fit.f90` and `test/goodness_of_fit_reference.f90`:
  SciFort-owned regression logic and public SciPy 1.17.0 numerical outputs; the
  reference module is regenerated by `tools/generate_goodness_of_fit_reference.py`.

Test reference data in `test/gamma_reference.f90`,
`test/beta_reference.f90`, `test/discrete_reference.f90`,
`test/discrete_extended_reference.f90`, and `test/continuous_reference.f90`
is generated by the corresponding scripts under `tools/` with mpmath, which
is BSD licensed; only computed numerical values, not mpmath source, enter the
repository. `test/extended_continuous_reference.f90`, `test/closed_form_reference.f90`,
`test/shape_family_reference.f90`, `test/transform_family_reference.f90`,
`test/positive_shape_reference.f90`, `test/scipy_transform_reference.f90`, and
`test/folded_reciprocal_reference.f90`, `test/normal_transform_reference.f90`, and
`test/burr_mielke_gibrat_wrapcauchy_reference.f90`, and
`test/genextreme_kappa_truncweibull_reference.f90`, and
`test/gengamma_halfgennorm_argus_erlang_reference.f90`, and
`test/crystalball_jf_pearson_breitwigner_reference.f90`, and
`test/genexpon_skewnorm_tukeylambda_rice_reference.f90`, and
`test/dpareto_vonmises_kstwobign_reference.f90`, and
`test/irwinhall_ksone_ncx2_ncf_reference.f90`, and
`test/randint_planck_dlaplace_logser_reference.f90`, and
`test/betabinom_hypergeom_nhypergeom_boltzmann_reference.f90`, and
`test/betanbinom_yulesimon_zipf_zipfian_reference.f90`, and
`test/gig_nig_skellam_reference.f90`, and
`test/genhyperbolic_fisher_reference.f90`,
`test/nct_gausshyper_reference.f90`, `test/landau_reference.f90`,
`test/wallenius_poisson_binom_reference.f90`, `test/kstwo_reference.f90`,
`test/levy_stable_studentized_range_reference.f90`,
`test/multivariate_reference.f90`, `test/test_matrix_distributions_reference.f90`,
`test/test_directional_random_matrix_reference.f90`, `test/gaussian_kde_reference.f90`,
and `test/qmc_reference.f90` are generated from public SciPy 1.17.0
numerical outputs by their corresponding
generator scripts under `tools/`;
only numerical reference values enter SciFort, and no SciPy implementation
source is copied or adapted. `test/test_transform_shape_four.f90` likewise
contains only independently generated SciPy 1.17.0 numerical outputs for the
generalized half-logistic, exponential-power, exponentiated-Weibull, and
power-lognormal regression cases; it contains no SciPy implementation source.

The documented SciPy source adaptations are relativistic Breit-Wigner,
finite-sample `kstwo`, Levy-stable, multivariate integration, multivariate-
hypergeometric sampling, matrix-t sampling, the Kendall exact recurrence, and
the identified goodness-of-fit kernels. The Sobol initialization data are also
redistributed from SciPy. These entries must not be described as wholly original
MIT-only implementations. The remaining entries state their individual origins.

- `src/integrate/quadrature.f90`: SciFort's one-dimensional sampled-data and
  fixed-rule integration layer. The irregular-spacing Simpson panel/correction
  formulas and the built-in Newton-Cotes coefficient table are adapted from
  SciPy 1.17.0 `scipy/integrate/_quadrature.py` under SciPy's BSD-3-Clause
  license. The Fortran interfaces, cumulative assembly, Romberg table,
  Gauss-Legendre delegation, normalized high-order Newton-Cotes solve, status
  handling, and residual refinement are SciFort implementation code.
- `src/integrate/quadpack_wrapper.f90`: original SciFort facade over the external
  `jacobwilliams/quadpack` repository, following its `master` branch in `fpm.toml`
  as of 2026-09-30. The released 2.1.x manifests use spaced `[[ test ]]` array-table
  headers that are rejected by some current FPM parsers; current `master` uses
  parser-compatible `[[test]]` headers. The wrapper calls the dependency's
  double-precision `dqags`, `dqagi`, `dqagp`, `dqawo`, `dqawf`, `dqaws`, and
  `dqawc` public routines; no QUADPACK source is copied into SciFort. Quadpack is
  BSD-3-Clause, copyright 2021-2022 Jacob Williams; its LICENSE also records the
  original SLATEC QUADPACK as public-domain work of the United States government.
  Source: https://github.com/jacobwilliams/quadpack/tree/master.
- `src/integrate/nested.f90`: original SciFort nesting/orchestration code built on
  the public `quad` facade. SciPy integration documentation and its public
  `_quadpack_py.py` dispatch layer were consulted for the `dblquad`/`tplquad`/
  `nquad` names, callback coordinate ordering, and QUADPACK mode mapping. No SciPy
  implementation body or QUADPACK numerical source was copied, translated, or
  adapted. Permanent regressions use independently recorded SciPy 1.17.0
  numerical outputs for representative nested integrals and weighted modes.
- `src/integrate/quad_vec.f90`: original SciFort implementation. SciPy 1.17.0
  documentation was consulted for the public `quad_vec` name, default tolerances,
  and `2`/`max` error-norm vocabulary. SciFort deliberately implements the vector
  result as independent calls through its existing scalar `quad` facade rather
  than adapting SciPy's `_quad_vec.py` interval-priority algorithm. No SciPy
  implementation source was copied, translated, or adapted.
- `src/integrate/qmc_quad.f90`: original SciFort implementation of replicated
  quasi-Monte Carlo integration using SciFort's existing scrambled Halton engine
  and explicit `rng_state`. SciPy 1.17.0 documentation was consulted for the
  `qmc_quad` name, defaults `n_estimates=8` and `n_points=1024`, reversed-bound
  sign convention, and sample-standard-error interpretation. No SciPy
  `_quadrature.py` implementation source was copied, translated, or adapted.
- `src/integrate/tanhsinh.f90`: original SciFort implementation of the
  double-exponential quadrature transformation from H. Takahasi and M. Mori,
  "Double Exponential Formulas for Numerical Integration", Publ. RIMS Kyoto
  Univ. 9 (1974), 721-741, DOI 10.2977/PRIMS/1195192451. SciPy 1.17.0 public
  documentation was consulted for the `tanhsinh` name and user-facing refinement
  defaults. The finite/infinite transformations, stable endpoint evaluation,
  tail guard, status handling, and Fortran API were written independently; no
  SciPy `_tanhsinh.py` implementation source was copied, translated, or adapted.
- `src/integrate/cubature.f90`: SciFort adaptive multidimensional integration
  with BSD-licensed numerical rule material adapted from SciPy 1.17.0. The
  Genz-Malik degree-7/degree-5 nodes and weights are adapted from
  `scipy/integrate/_rules/_genz_malik.py`, and the Gauss-Kronrod 15/21 nodes and
  high-order weights are adapted from
  `scipy/integrate/_rules/_gauss_kronrod.py`; the lower Gauss rule is generated
  by SciFort's existing `roots_legendre`. The global largest-error-region
  refinement and all-coordinate bisection strategy follow SciPy 1.17.0
  `scipy/integrate/_cubature.py`. The Fortran region storage, scalar/vector
  callback overloads, status API, user-point partitioning, and finite/infinite
  coordinate transformations are SciFort implementation code. The mathematical
  rule originates with A. C. Genz and A. A. Malik, "Remarks on algorithm 006:
  An adaptive algorithm for numerical integration over an N-dimensional
  rectangular region", J. Comput. Appl. Math. 6 (1980), 295-302,
  DOI 10.1016/0771-050X(80)90039-X.

- `src/integrate/odeint.f90`: original SciFort compatibility facade over the
  native SciFort LSODA `solve_ivp` path. The public callback order, time-grid and
  result orientation, Jacobian/band conventions, scalar/vector tolerance semantics,
  tolerance defaults, and legacy
  option names were based on the SciPy 1.17.0 `scipy.integrate.odeint`
  documentation, inspected 2026-09-30 at
  https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.integrate.odeint.html.
  No SciPy wrapper implementation or ODEPACK wrapper source was copied. Numerical
  outputs from SciPy 1.17.0 were used as independent regression references for
  the documented pendulum problem and backward/repeated-time integration. For
  the later `tcrit`, `mxordn`, `mxords`, and exact `mxstep` compatibility work,
  SciPy 1.17.0 `_odepackmodule.c` was inspected for wrapper semantics and the
  public-domain DLSODA documentation at
  https://jacobwilliams.github.io/odepack/proc/dlsoda.html was inspected for
  ITASK=4/TCRIT, MXORDN/MXORDS, and MXSTEP/NST call-accounting behavior. For
  `full_output` compatibility, SciPy's public documentation was used for the
  diagnostic field names/meanings, and the public-domain DLSODA documentation
  was used for the optional-output semantics and `LENRW`/`LENIW` workspace
  formulas. Local SciPy 1.17.0 runs supplied independent regression values for
  a simple decay problem. The diagnostic capture, result storage, and workspace
  calculations are original SciFort code; no SciPy wrapper or ODEPACK source was
  copied, translated, or adapted. The final `ixpr`, `mxhnil`, and `printmessg`
  controls were implemented independently from the option semantics documented
  by SciPy 1.17.0 and public DLSODA descriptions: `ixpr` reports method switches,
  `mxhnil` caps small-step messages with zero selecting the solver default, and
  `printmessg` reports the final status. SciFort's existing terminal handling for
  an LSODA too-small-step condition was intentionally retained.

- `src/integrate/ivp.f90` and `src/integrate/ivp_complex.f90`
  (`rk23_step`, `rk45_step`): the Bogacki-Shampine
  3(2) and Dormand-Prince 5(4) Runge-Kutta tableaux are adapted from Jacob
  Williams' `rklib` `src/rklib_variable_steps.f90`, branch `master`, retrieved
  2026-09-28 from
  https://raw.githubusercontent.com/jacobwilliams/rklib/master/src/rklib_variable_steps.f90.
  Upstream copyright is 2022-2023 Jacob Williams and the license is
  BSD-3-Clause; the complete notice is retained in `THIRD_PARTY_LICENSES.md`.
  The real and complex `solve_ivp` orchestration, adaptive storage, explicit-state
  result types, initial-step selection, `t_eval`, backward integration, RK23
  cubic-Hermite interpolation, scalar/vector event root storage, and scalar-or-vector
  `rtol`/`atol` validation and broadcasting are SciFort implementation code. The public
  tolerance semantics follow SciPy 1.17.0 documentation: scalars are broadcast and
  rank-one length-`n` inputs supply componentwise values used in
  `atol + rtol*abs(y)`. SciFort extends the same componentwise `rtol` weighting
  consistently through its implicit real/complex methods; their scalar Newton stopping
  heuristic uses `maxval(rtol)`. Componentwise tolerances are also passed through the
  SciPy-derived Radau/BDF finite-difference Jacobian perturbation logic and both LSODA
  branches. The complex
  variant applies the same real Runge-Kutta
  coefficients to complex state/stage vectors and uses magnitude-based error
  norms, consistent with SciPy's complex explicit-method behavior. Event crossing
  tests, positive-integer terminal occurrence semantics,
  same-step root ordering/truncation, and per-event result grouping are adapted
  from SciPy 1.17.0 `scipy/integrate/_ivp/ivp.py`, tag `v1.17.0`, inspected
  2026-09-29 at
  https://github.com/scipy/scipy/blob/v1.17.0/scipy/integrate/_ivp/ivp.py.
  SciPy is BSD-3-Clause; the applicable notice is retained in
  `THIRD_PARTY_LICENSES.md`. SciPy 1.17.0 real and complex numerical outputs were
  also used for regression comparison.

- `src/integrate/ivp.f90` and `src/integrate/ivp_complex.f90` (`rk45_step`
  dense-output coefficient block and `rk45_dense_state`): the RK45 quartic
  continuous-extension coefficient matrix
  `P`, the stage combination `Q = K.T.dot(P)`, and the normalized polynomial
  evaluation convention are adapted from SciPy 1.17.0
  `scipy/integrate/_ivp/rk.py`, tag `v1.17.0`, inspected 2026-09-29 at
  https://github.com/scipy/scipy/blob/v1.17.0/scipy/integrate/_ivp/rk.py. SciPy
  is BSD-3-Clause; the applicable notice is retained in `THIRD_PARTY_LICENSES.md`.
  The Fortran Horner evaluation, coefficient/result storage, event integration,
  and backward-step handling are SciFort implementation code. Regression values
  in `test/test_integrate_ivp.f90` were generated with SciPy 1.17.0.

- `src/integrate/sparse_jacobian.f90`: SciFort implementation of minimal
  one-based compressed-sparse-column containers for real and complex
  `solve_ivp` Jacobians. The public `jac_sparse_matrix=` behavior follows the
  documented SciPy 1.17.0 allowance for constant sparse Jacobians in Radau/BDF,
  but the CSC validation, Fortran derived types, constructor, and solver-state
  integration are original SciFort code. `jac_sparse=` reuses the same containers
  for state-dependent callback values, and `ivp_csc_from_pattern` builds an empty
  CSC value array directly from a logical finite-difference structure. No SciPy
  sparse-storage, sparse-LU, or coloring implementation is copied. CSC Newton
  systems are delegated to the original SciFort sparse LU implementation described
  below.

- `src/integrate/sparse_lu.f90`: original SciFort implementation of sparse
  Gaussian elimination with partial row pivoting for the CSC `solve_ivp`
  paths. It converts validated CSC input into sorted sparse rows,
  inserts elimination fill-in dynamically, stores the combined L/U factors
  sparsely, and applies recorded row interchanges during triangular solves. The
  real, complex, and real-Jacobian/complex-shift variants share the same
  algorithm. The mathematical method is standard Gaussian elimination with
  partial pivoting as described in Golub and Van Loan, *Matrix Computations*,
  4th edition (2013); no implementation code from that text, SciPy, SuperLU,
  UMFPACK, or another sparse solver was copied. SciFort's implementation is MIT
  licensed with the rest of the original project code.

- `src/integrate/radau.f90`: the fifth-order three-stage Radau IIA collocation
  method, embedded third-order error estimator, transformed real/complex Newton
  systems, step-size predictor, cubic collocation dense-output polynomial, and
  dense forward-difference Jacobian strategy are translated and adapted from
  SciPy 1.17.0 `scipy/integrate/_ivp/radau.py` and the dense numerical Jacobian
  routines in `scipy/integrate/_ivp/common.py`, tag `v1.17.0`, inspected
  2026-09-29 at
  https://github.com/scipy/scipy/blob/v1.17.0/scipy/integrate/_ivp/radau.py and
  https://github.com/scipy/scipy/blob/v1.17.0/scipy/integrate/_ivp/common.py.
  SciPy is BSD-3-Clause; the applicable notice is retained in
  `THIRD_PARTY_LICENSES.md`. SciFort supplies the Fortran callback API, compact
  partial-pivoting real/complex dense LU kernels, integration with the existing
  `solve_ivp` result/event machinery, and storage of Radau cubic dense
  polynomials. The `jac_matrix=` constant-matrix spelling and logical
  `jac_sparsity=` public semantics follow SciPy's documented `jac` and
  `jac_sparsity` options. Constant `jac_sparse_matrix=` and state-dependent
  `jac_sparse=` inputs use the native CSC container and sparse LU implementation
  described above, including both real and complex-shifted Radau Newton systems.
  When `jac_sparsity=` drives numerical differencing, SciFort independently fills
  that declared structure directly in CSC storage before sparse factorization.
  SciFort's finite-difference
  grouping helper in
  `src/integrate/jacobian_fd.f90` is independently implemented as a deterministic
  greedy grouping of columns with disjoint declared nonzero rows; it does not
  copy SciPy's sparse-matrix/grouping implementation. The optional rank-two
  `fun_vectorized` callback adapts SciPy 1.17.0's `OdeSolver` vectorized-RHS
  contract from `scipy/integrate/_ivp/base.py` and the batched finite-difference
  call pattern in `scipy/integrate/_ivp/common.py`: SciFort uses a separate
  Fortran procedure argument because callback rank is part of the interface, and
  keeps these hidden Jacobian probes out of `nfev`. Regression values and work
  counts were generated with installed SciPy 1.17.0.

- `src/integrate/bdf.f90` and `src/integrate/bdf_complex.f90`: the
  variable-order 1--5 quasi-constant-step BDF/NDF recurrence, step/order
  adaptation, Newton convergence tests, dense-output polynomial, and dense
  forward-difference Jacobian strategy are translated and
  adapted from SciPy 1.17.0 `scipy/integrate/_ivp/bdf.py` and the dense numerical
  Jacobian routines in `scipy/integrate/_ivp/common.py`, tag `v1.17.0`, inspected
  2026-09-29 at
  https://github.com/scipy/scipy/blob/v1.17.0/scipy/integrate/_ivp/bdf.py and
  https://github.com/scipy/scipy/blob/v1.17.0/scipy/integrate/_ivp/common.py.
  SciPy is BSD-3-Clause; the applicable notice is retained in
  `THIRD_PARTY_LICENSES.md`. SciFort supplies the Fortran callback API, compact
  dense LU factorization, integration with the existing `solve_ivp` result/event
  machinery, and storage of BDF dense polynomials. The complex translation uses
  complex state/Jacobian/LU arrays and follows SciPy `common.py` in perturbing
  complex states along the real axis, with the sign of `real(f)` selecting the
  finite-difference direction. Constant dense matrices are retained once in
  solver state and therefore report zero Jacobian evaluations, matching SciPy's
  array-like Jacobian accounting. Constant CSC matrices and state-dependent
  `jac_sparse=` callback values use the native sparse container and sparse LU
  implementation described above and are not assembled into a dense Newton
  matrix. With logical structural sparsity, numerical finite differences are
  written directly to CSC values having that pattern and are factored sparsely.
  Column grouping and the optional batched `fun_vectorized` finite-difference
  path use
  the independently implemented helper and SciPy vectorized-RHS conventions
  described in the Radau entry above. Regression values and work counts for real
  and complex BDF were generated with installed SciPy 1.17.0.

- `src/integrate/lsoda.f90`: SciFort implementation of the real-system LSODA
  compatibility path. The public method meaning and automatic nonstiff-Adams/stiff-BDF
  behavior were checked against SciPy 1.17.0 documentation and LLNL's ODEPACK
  description, inspected 2026-09-29 at
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.LSODA.html and
  https://computing.llnl.gov/projects/odepack. The active Adams path now carries a
  persistent scaled-derivative Nordsieck array through method orders 1--12. Its step
  rescaling, Pascal-triangle prediction, functional-corrector iteration, accepted
  `YH(:,j) += EL(j)*ACOR` update, successful-step `IALTH` delay, `TESCO`
  lower/current/upper-order comparison, `SM1` stability limits, 20-step switch
  reconsideration interval, order-5 Adams-to-BDF guard, roundoff-pollution switch
  branch, and ideal Adams/BDF step-factor comparison are adapted from the public-domain
  ODEPACK `DSTODA` source in Jacob Williams' modernization, inspected 2026-09-29 at
  https://raw.githubusercontent.com/jacobwilliams/odepack/master/src/M_da1/dstoda.inc.
  SciFort evaluates dense output directly from the retained Nordsieck polynomial; that
  Taylor evaluation follows from the stored scaled-derivative representation and was
  implemented independently rather than copied from an ODEPACK interpolation routine.

  Adams correction coefficients and the three `TESCO` error/order-selection constants
  for orders 1--12 are generated with the public-domain ODEPACK `DCFODE` recurrence,
  inspected 2026-09-29 at
  https://raw.githubusercontent.com/jacobwilliams/odepack/master/src/M_da1/dcfode.inc.
  The `CM1`/`CM2` values used by the finite-error method-switch comparison at orders
  1--5 are the corresponding `TESCO(2,q)*ELCO(q+1,q)` products. ODEPACK's public-domain
  notice was inspected at
  https://raw.githubusercontent.com/jacobwilliams/odepack/master/LICENSE and is
  reproduced in `THIRD_PARTY_LICENSES.md`.

  SciFort retains an independently implemented nonuniform endpoint-derivative/time
  history. It is used to reconstruct a Nordsieck state when persistent state is not
  available, to support conservative recovery, and for compatibility diagnostics; it
  is no longer rebuilt after every successful high-order Adams step. Reconstruction
  interpolates the retained derivative samples at their actual signed offsets and
  converts the polynomial coefficients to scaled solution derivatives. The accepted
  functional-iteration derivative is retained directly, matching `DSTODA`'s persistent
  corrector semantics, so an extra RHS evaluation at the corrected endpoint is no
  longer performed. For automatic LSODA startup, SciFort now follows DLSODA's
  `H0**(-2) = 1/(TOL*w0**2) + TOL*norm(F)**2` formula using the already evaluated
  initial RHS, the scalar tolerance clamp, and the DLSODA weighted maximum norm.
  During Adams trials, the error weights are frozen from the accepted state at the
  beginning of the step, `atol + rtol*abs(y_n)`, and the functional-corrector,
  error-test, and lower/current/upper-order comparisons use the corresponding
  tolerance-scaled maximum norm, matching DLSODA's `DMNORM` semantics. The Adams
  step-growth limiter also follows `DSTODA`'s `RMAX` states: 10000 during startup,
  normally 10 after the first selected rescaling, and 2 after a convergence or
  error-test failure until the next successful selected rescaling. These formulas
  were checked against the public-domain DLSODA/DSTODA sources cited above on
  2026-09-30. Installed SciPy 1.17.0 provides the compatibility work-count
  references in `test/test_integrate_lsoda_work_counts.f90`.

  The correction-rate estimate supplies the local Lipschitz information used for the
  Adams stability bound, with the most recent nonzero value retained as a `PDLAST`
  analogue and accepted-step secants used on guarded degenerate paths.

  Rejected Adams error-test trials do not run method-switch tests or Jacobian probes.
  Corrector-convergence failure uses `DSTODA`'s factor-of-four retry reduction. For the
  first through third error-test failures, SciFort now follows `DSTODA`'s rejected-step
  order-selection path: it compares current-order and one-order-lower retry factors,
  forms the lower-order estimate from the retained `YH(:,L)` column and `TESCO(1,q)`,
  applies the Adams `SM1` stability bound, forbids an order increase, and caps the
  retry factor at 0.2 after the second failure. A nonzero local correction-rate
  estimate from a converged rejected trial is retained as the `PDLAST` analogue used
  by that stability bound. After four consecutive error-test failures SciFort applies
  a conservative adaptation of `DSTODA`'s repeated-failure recovery: it resets to
  method order 1 and reduces the retry step by 0.1. ODEPACK begins the analogous
  special recovery at `KFLAG <= -3`; the fourth-rejection threshold is retained
  because the earlier hybrid controller's third-rejection threshold regressed an
  established SciFort dense-output test.

  ODEPACK keeps one Nordsieck `YH` array across method switches and changes only
  the active method coefficients/order before rescaling it. SciFort instead has
  separate Adams Nordsieck and SciPy-derived BDF backward-difference storage. On
  a BDF-to-Adams switch, SciFort independently converts the BDF polynomial into
  power-series/Nordsieck coefficients at the current endpoint and preserves the
  selected switch order through order 5; the transferred first-derivative column
  is reset to the exact accepted `h*f` value while higher coefficients are
  retained. On an Adams-to-BDF switch, SciFort now performs the inverse polynomial
  conversion at the `DSTODA`-selected order and ideal BDF step ratio rather than
  reinitializing order 1 at the previous step. These are adaptations of the
  method-state continuity semantics in public-domain `DSTODA`, not translations
  of ODEPACK's common-block storage layout.

  The LSODA stiff branch also adapts `DSTODA`'s chord-iteration freshness rules:
  while LSODA mode is active, the reused SciFort BDF Jacobian is reevaluated after
  20 accepted stiff steps or when the current Newton coefficient differs by more
  than 30 percent from that associated with the last refreshed Jacobian. Jacobian
  formation is deferred across the method switch until the first BDF trial, as in
  `DSTODA`'s `IPUP`/`PJAC` path. LSODA-only BDF trials use the pure BDF coefficients
  (the SciPy NDF `kappa` enhancement remains enabled for standalone `BDF`), freeze
  error weights from the accepted state, and use the ODEPACK weighted maximum norm.
  Internally generated Jacobians follow the public-domain `DPRJA` dense/banded
  forward-difference perturbation formula based on `sqrt(UROUND)`, the weighted
  RHS norm, and the current stiff step. Their perturbation calls are included in
  `nfev`, and the already evaluated predictor RHS is reused by the first Newton
  iteration. `DPRJA` was inspected 2026-09-30 at
  https://raw.githubusercontent.com/jacobwilliams/odepack/master/src/M_da1/dprja.inc.
  Standalone SciFort `BDF` leaves these ODEPACK-specific paths disabled, preserving
  its separately attributed SciPy BDF/NDF behavior and reference work counts. The
  `DSTODA` refresh conditions were checked against the public-domain source cited
  above; the matching DLSODA defaults `CCMAX=0.3` and `MSBP=20` were checked in
  Netlib `opkdmain.f`, inspected 2026-09-30 at https://netlib.org/odepack/opkdmain.f.
  The LSODA-only stiff controller now also uses an independently adapted ODEPACK
  Nordsieck recurrence rather than SciPy's BDF/NDF step controller. BDF correction
  polynomials and `TESCO` constants follow public-domain `DCFODE`; prediction is the
  in-place Pascal update of `YH`; chord iteration follows `DSTODA`'s `CRATE`/`CONIT`
  convergence test with `MAXCOR=3`; and successful/error/convergence-failure step
  changes follow the `IALTH`, lower/current/upper-order, `RMAX`, and factor-of-four
  rules in `DSTODA`. `DCFODE` and `DSTODA` were inspected 2026-09-30 at
  https://raw.githubusercontent.com/jacobwilliams/odepack/master/src/M_da1/dcfode.inc
  and https://raw.githubusercontent.com/jacobwilliams/odepack/master/src/M_da1/dstoda.inc.
  This is a clean adaptation to SciFort state rather than a copy of ODEPACK common
  blocks. Installed SciPy 1.17.0 supplies compatibility references: the scalar stiff
  cosine regression matches 116 RHS evaluations, 10 Jacobian evaluations, 10 LU
  factorizations, 58 accepted steps (59 stored points including the initial point),
  and the key order-transition times. The corresponding numerical-Jacobian run adds
  exactly one RHS call per scalar Jacobian evaluation. Representative nonstiff work
  counts also match installed SciPy 1.17.0 exactly.

  The stiff branch reuses SciFort's separately attributed `src/integrate/bdf.f90`. Its
  accepted endpoint derivative is recovered algebraically from the converged SciFort
  BDF equation, `d_corr = c*f_new - psi`, rather than by an additional RHS callback.
  `lband`/`uband` and `min_step` support remain original SciFort code implementing the
  public SciPy option semantics: packed analytic bands use `uband+i-j+1` in Fortran's
  1-based row indexing, and LSODA numerical bands perturb nonoverlapping columns in
  ODEPACK-style groups before unpacking into SciFort's dense BDF Newton matrix.
  Installed SciPy 1.17.0 numerical outputs and work counts are used as compatibility
  references. The later native Nordsieck/chord controller closes the representative
  stiff trajectory/work-count gap while retaining this finite-difference accounting.

- `src/integrate/dop853.f90` and `src/integrate/dop853_complex.f90`
  (`dop853_step`, `dop853_dense_coefficients`, `dop853_dense_state`): the DOP853
  stage tableau, eighth-order solution
  weights, embedded fifth-/third-order error estimator, three extra dense-output
  stages, interpolation coefficients, and seventh-order dense evaluation are
  adapted from Jacob Williams' modern Fortran DOP853 source,
  branch `master`, retrieved 2026-09-28 from
  https://raw.githubusercontent.com/jacobwilliams/dop853/master/src/dop853_module.F90.
  That modernization is Copyright (c) 2015-2022 Jacob Williams and retains
  the original DOP853 Copyright (c) 2004 Ernst Hairer; both redistribution
  notices are permissive and are reproduced in `THIRD_PARTY_LICENSES.md`.
  SciFort adapts the DOP853 step and dense-output numerical formulas in the real
  and complex modules; the complex version applies the same real tableau to
  complex stages. Solver orchestration, initial-step selection, result storage,
  `t_eval` selection, event-set semantics, root sorting/bookkeeping, and
  conditional storage of dense coefficients remain in `src/integrate/ivp.f90`
  and `src/integrate/ivp_complex.f90`. The upstream
  master snapshot was inspected
  directly; no commit identifier is asserted because one was not available
  from the retrieval interface used for this batch.

## Merge-regression references

`test/test_merge_regressions.f90` contains original SciFort regression checks
based on support endpoints, NaN propagation, translation invariance, equal-weight
softmax normalization, uniform-distribution quantiles, and independence of normal
coordinates. It reuses the existing scalar normal log tails as an identity check;
no third-party implementation or numerical tables were copied for these tests.
The verbatim `LICENSES/CC-BY-4.0.txt` was retrieved on 2026-09-23 from the
research repository and the canonical Creative Commons legal-code URL recorded
in `THIRD_PARTY_LICENSES.md`. Its SHA-256 is
`9ba9550ad48438d0836ddab3da480b3b69ffa0aac7b7878b5a0039e7ab429411`.

- `src/special/complex_sine_integrals.f90`: original SciFort implementation of
  complex Si/Ci/Shi/Chi. The small-argument kernels are independently derived from
  the defining series in NIST DLMF 6.2; the complementary region uses the DLMF
  interrelations with SciFort's independently implemented principal E1/Ei kernels.
  SciPy 1.18.0 documentation was consulted for public names, output ordering, and
  negative-real-axis complex branch conventions; SciPy 1.17.0 public numerical
  outputs were used only as regression references. No SciPy/mpmath implementation
  source was copied, translated, or adapted.

## Required entry format for future imports

Add a row or section containing:

- SciFort file and procedure names;
- upstream project, file path, and URL;
- exact version, tag, or commit hash;
- original authors and copyright holders;
- original license and location of retained license text;
- whether the code was copied, translated, adapted, or independently
  reimplemented from a publication;
- meaningful modernization changes;
- tests used to establish parity.

- `src/special/fresnel_zeros.f90`: original SciFort implementation of first-quadrant complex Fresnel S/C zero tables. Initial guesses follow NIST DLMF 7.13.5-7.13.7 (Kreyszig/Fettis-Caslin asymptotics) and are refined independently with Newton iterations using SciFort's complex Fresnel implementation. SciPy documentation was consulted for public API names, return ordering, and first-quadrant convention; no SciPy or Zhang-Jin implementation code was copied or translated.
- `src/special/faddeeva.f90` (`modfresnelp`, `modfresnelm`): original SciFort
  implementation of the modified Fresnel tail pairs from their documented integral
  definitions and the DLMF/SciFort Faddeeva identities. SciPy 1.18.0 documentation was
  consulted for public names and formulas, and SciPy 1.17.0 public numerical outputs
  were used only for nonnegative regression comparisons. No SciPy/specfun/Zhang-Jin
  implementation source was copied or translated. The negative-lower-limit `F` branch
  follows the mathematical integral definition even where SciPy 1.17.0 differs.

- `src/special/airy.f90` (complex `airy` overload): original SciFort implementation
  derived independently from NIST DLMF 9.2 connection formulas, 9.4 Maclaurin
  series, 9.7 asymptotic expansions, and the Airy ODE. SciPy documentation was
  consulted for the public complex-argument API and SciPy 1.17.0 public numerical
  outputs were used only for regression values. No AMOS, Cephes, SciPy, or other
  implementation source was copied, translated, or adapted.

- `src/special/airy.f90` (complex `airye` overload): original SciFort implementation
  of SciPy-compatible complex Airy scaling. The scaling definitions follow the
  public SciPy `airye` API. Moderate recessive Ai/Ai' values are continued inward
  from a DLMF 9.7 scaled-asymptotic anchor using Taylor coefficients derived
  directly from the Airy ODE `y''=z*y`; outer-plane Ai and Bi values combine
  DLMF 9.2 rotated-Ai identities with exponent differences formed before
  exponentiation. SciPy 1.17.0 public numerical outputs were used only as
  regression references. No AMOS, Cephes, SciPy, or other third-party
  implementation source was copied, translated, or adapted.
- `src/special/hypergeometric.f90` (`hyp0f1`): original SciFort implementation derived independently from the generalized-hypergeometric defining series in NIST DLMF 16.2, the differential equation in DLMF 16.8.3, and the Bessel-J connection in DLMF 10.16.9. Real large-argument evaluation reuses SciFort's existing real Bessel kernels; complex continuation uses a Taylor recurrence derived directly from the differential equation. SciPy 1.17.0 public numerical outputs were used only as regression references and for IEEE/pole API behavior. No SciPy, Boost, Cephes, AMOS, or other third-party implementation source was copied, translated, or adapted.
- `src/special/hypergeometric.f90` (`hyp1f1`): original SciFort implementation derived independently from the Kummer defining series and differential equation in NIST DLMF 13.2, Kummer's transformation, and the large-argument expansions in DLMF 13.7. Moderate arguments use a Taylor recurrence derived directly from the differential equation; large real and nonreal arguments use independently coded asymptotic recurrences. SciPy 1.17.0 public numerical outputs were used only as regression references and to characterize real/complex pole and IEEE API behavior. Selected complex references were cross-checked with mpmath to avoid treating known legacy complex-evaluation discrepancies as mathematical targets. No SciPy, Boost, specfun, Zhang-Jin, Cephes, or other third-party implementation source was copied, translated, or adapted.

- `src/special/hypergeometric.f90` (`hyp2f1`): original SciFort implementation derived independently from the Gauss hypergeometric defining series in NIST DLMF 15.2 and the hypergeometric differential equation and analytic-continuation framework in DLMF Chapter 15. Local Taylor coefficients are derived directly from the differential equation, and continuation paths avoid the regular singular point at `z=1`; no SciPy, Cephes, Boost, Zhang-Jin, or other third-party implementation source was copied, translated, or adapted. SciPy 1.18.0 documentation was consulted for the public API and branch/pole conventions. SciPy 1.17.0 public numerical outputs were used only as regression references, with difficult complex cases cross-checked against high-precision mpmath values.

- `src/special/hypergeometric.f90` (`hyperu`): original SciFort implementation derived independently from the confluent-hypergeometric differential equation, connection formula, Kummer transformation, and large-positive-argument asymptotic expansion in NIST DLMF Chapter 13. Moderate positive arguments are continued inward from an asymptotic anchor using the existing Taylor recurrence derived from the Kummer ODE; small positive arguments use the exact 1F1 connection formula, with symmetric Richardson extrapolation in `b` near integer parameters to remove the formula's artificial gamma singularities. Negative-integer `a` uses the generalized-Laguerre terminating recurrence. SciPy 1.18.0 documentation was consulted for the public real API and negative-axis branch-cut convention. SciPy 1.17.0 public numerical outputs and high-precision mpmath values were used only as regression references; no SciPy/specfun/Zhang-Jin/Cephes/Boost implementation source was copied, translated, or adapted.
- `src/special/parabolic_cylinder.f90` (`pbdv`): original SciFort implementation derived independently from NIST DLMF Chapter 12. Positive arguments use DLMF 12.7.14 together with SciFort's independently implemented Kummer `hyperu`; values and derivatives at zero follow DLMF 12.2.6-12.2.7 after converting from DLMF `U(a,z)` notation to `D_v(z)`, and negative arguments use a Taylor recurrence derived directly from the parabolic-cylinder ODE. Nonnegative integer orders use the DLMF Hermite recurrence to avoid numerical contamination of the recessive solution. SciPy 1.17/1.18 documentation was consulted for the public `pbdv` name, real parameter/argument domain, and paired value/derivative outputs. High-precision mpmath values and SciPy 1.17.0 public outputs were used only as regression references; no SciPy/specfun/Zhang-Jin/Cephes/Boost implementation source was copied, translated, or adapted.
- `src/special/parabolic_cylinder.f90` (`pbvv`, `pbdv_seq`, `pbvv_seq`, `pbdn_seq`): original SciFort implementations derived independently from NIST DLMF Chapter 12. `pbvv` uses the DLMF reciprocal-gamma values for the independent `V(a,z)` solution at zero, Taylor coefficients derived from the parabolic-cylinder differential equation, the large-positive-argument `V(a,z)` asymptotic series, and the exact `U`/`V` connection formula after translating DLMF's parameter `a` to SciPy's order `v`. The sequence wrappers reproduce SciPy's public order/shape convention but call SciFort scalar kernels; complex nonnegative-integer `pbdn_seq` uses the Hermite recurrence for the analytic `D_n(z)`. SciPy 1.17/1.18 documentation and public outputs were consulted for names, output ordering, and sequence shape. High-precision mpmath values were used as mathematical regression references, including complex `D_n(z)` cases where SciPy 1.17.0's legacy backend is nonanalytic. No SciPy/specfun/Zhang-Jin/Cephes/Boost implementation source was copied, translated, or adapted.
- `src/special/parabolic_cylinder.f90` (`pbwa`): original SciFort implementation derived independently from NIST DLMF Chapter 12, especially the `W(a,x)` origin normalization and confluent-hypergeometric representation in DLMF 12.14. The origin value and derivative are evaluated through complex `loggamma` ratios, and all nonzero arguments are obtained by a Taylor recurrence derived directly from `y'' + (x*x/4-a)*y = 0`. SciPy 1.17/1.18 documentation was consulted for the public `pbwa` name, paired value/derivative outputs, and the legacy `abs(a),abs(x) <= 5` domain. High-precision mpmath evaluation of the DLMF representation was used for mathematical regression references. No SciPy/specfun/Zhang-Jin/Cephes/Boost implementation source was copied, translated, or adapted.
- `src/special/owens_t.f90` (`owens_t`): original SciFort implementation from
  Owen's defining real integral. The substitution `t=tan(theta)` gives the
  bounded angular integral used for `abs(a) <= 1`; the reciprocal-slope
  relation is obtained from the symmetry of the corresponding bivariate-normal
  probability and is used only as a mathematical reduction. SciPy 1.17.0
  public outputs were used as regression references. No SciPy, Boost, Cephes,
  or other third-party implementation source was copied, translated, or adapted.
- `src/special/elementary.f90` (`agm`, `huber`, `pseudo_huber`, `ndtri_exp`):
  original SciFort implementations from the public mathematical definitions.
  AGM uses the classical arithmetic/geometric iteration; the two loss functions
  use algebraically stable forms of their defining expressions; `ndtri_exp`
  combines SciFort's normal inverse with a Newton solve of the standard normal
  logarithmic-tail asymptotic. SciPy 1.17.0 public behavior was consulted for
  names, domains, and IEEE endpoint conventions, and high-precision mpmath was
  used for deep-tail inverse-normal references. No SciPy implementation source
  was copied, translated, or adapted.
- `src/special/faddeeva.f90` (`voigt_profile`): original SciFort wrapper of the
  standard Voigt/Faddeeva identity using SciFort's existing independent `wofz`
  kernel, with explicit Cauchy and degenerate limits. SciPy 1.17.0 public outputs
  were used only as numerical/API regression references; no implementation
  source was copied or adapted.
- `src/stats/normal.f90` extreme-log-CDF guard: original SciFort numerical fix
  that detects when `-z**2/2` is below binary64 range before squaring `z`, and
  suppresses an otherwise avoidable overflow trap.


- `src/special/distribution_helpers.f90`: original SciFort wrapper implementations derived independently from the standard binomial, negative-binomial, Poisson, chi-square, F, gamma, beta, Student-t, and normal CDF identities in terms of regularized incomplete beta/gamma functions and the standard-normal quantile. The direct inverses reuse SciFort's independently implemented beta/gamma inverse kernels; the additional parameter inverses (`bdtrik`, `bdtrin`, `btdtria`, `btdtrib`, `chdtriv`, `fdtridfd`, `gdtrib`, `nbdtrik`, `nbdtrin`, `pdtrik`, `stdtridf`, `nrdtrimn`, and `nrdtrisd`) use independently written safeguarded positive-parameter root solving or exact normal transforms. SciPy 1.17.0 public documentation and public numerical outputs were consulted for function names, argument ordering, floor/count conventions, and regression references. No SciPy/CDFLIB/Boost/Cephes implementation source was copied, translated, or adapted. SciFort intentionally follows mathematical endpoint behavior rather than two observed SciPy 1.17 legacy quirks: `pdtri(k,1)=0` and `stdtrit(df,0)=-infinity`. For `fdtridfd`, whose equation can have two positive roots, SciFort explicitly returns the smaller root instead of reproducing a legacy large-sentinel failure result.

- `src/special/elementary.f90` (`cbrt`, `exp10`, `exp2`, `expm1`, `log1p`,
  `powm1`, `radian`, `round`, `sinc`, `sindg`, `cosdg`, `tandg`, `cotdg`,
  `softplus`, `diric`): original SciFort implementations from the public
  mathematical definitions and Fortran intrinsic operations. Stable near-one
  exponential/logarithmic cases reuse SciFort's existing independently written
  `expm1_safe`, `log1p_safe`, and `log1pexp` kernels. Degree trigonometric
  functions use independent modulo-360 reduction with exact quadrant handling,
  and `diric` uses the defining periodic sine ratio with removable-singularity
  limits. SciPy documentation and SciPy 1.17.0 public numerical outputs were
  consulted only for names, argument ordering, ties-to-even rounding semantics,
  and regression references. No SciPy/Cephes implementation source was copied,
  translated, or adapted.

- `src/special/noncentral_distribution_helpers.f90`: original SciFort wrapper
  implementations built around the project's existing independently implemented
  `scifort_ncx2`, `scifort_ncf`, and `scifort_nct` probability and quantile
  kernels. The parameter-inverse searches are original safeguarded bisection/
  bracketing code. SciPy public documentation and SciPy 1.17.0 public numerical
  outputs were consulted for names, argument ordering, domains, and regression
  references; the documented non-monotonicity of noncentral-F CDFs in their
  degrees of freedom motivated SciFort's explicit smallest-positive-root branch.
  No SciPy, CDFLIB, Boost, or Cephes implementation source was copied, translated,
  or adapted. SciFort intentionally maps `nctdtrit(df,nc,0)` to mathematical
  negative infinity instead of SciPy 1.17.0's observed positive-infinity legacy
  result.

- `src/special/orthogonal_polynomials.f90`: original SciFort implementation of
  the nonnegative-integer-degree Legendre, Jacobi, Gegenbauer, Laguerre,
  generalized-Laguerre, Hermite, Hermite-normalized, Chebyshev, and shifted
  evaluator families. The three-term recurrences and shifted interrelations were
  derived independently from NIST DLMF Chapter 18, especially section 18.9.
  SciPy 1.17.0 public documentation and outputs were consulted for API names,
  parameter ordering, normalizations, and regression references. Difficult
  complex cases were cross-checked against high-precision mpmath values. No
  SciPy, Cephes, Boost, or other third-party implementation source was copied,
  translated, or adapted. The `assoc_laguerre(x,n,k)` compatibility interface is
  an original thin wrapper that only reverses arguments and supplies SciPy's
  default `k=0` before calling the existing `eval_genlaguerre` implementation.

- `src/special/orthogonal_roots.f90`: original SciFort implementation of the
  classical Gaussian-quadrature node/weight families. The Golub-Welsch
  Jacobi-matrix construction, zeroth-moment weight identity, and classical
  recurrence coefficients were derived independently from NIST DLMF section
  3.5(vi) and Chapter 18. Closed-form Chebyshev rules use the standard trigonometric
  zeros and weights. SciPy 1.17.0 public documentation and public numerical
  outputs were consulted for API names, parameter ordering, weight normalizations,
  shifted-family conventions, and regression references. No SciPy, Cephes,
  Boost, LAPACK, or other third-party implementation source was copied,
  translated, or adapted; the general rules reuse SciFort's existing original
  dense symmetric Jacobi eigensolver. The historical `p_roots`, `j_roots`,
  `cg_roots`, `l_roots`, `la_roots`, `h_roots`, `he_roots`, `t_roots`,
  `u_roots`, `c_roots`, `s_roots`, `ps_roots`, `js_roots`, `ts_roots`, and
  `us_roots` names are original generic aliases of those canonical procedures;
  no additional numerical implementation is present.
- `src/special/statistical_legacy.f90`: original SciFort compatibility wrappers
  around the project's independently implemented `scifort_ksone`,
  `scifort_kstwobign`, and `scifort_tukeylambda` kernels. SciPy 1.17.0
  public documentation and public numerical outputs were consulted for the
  historical names, argument ordering, endpoint behavior, and regression
  references. No SciPy implementation source was copied, translated, or
  adapted.
- `src/special/lambert_w.f90`: original SciFort implementation of complex Lambert W and Wright omega derived independently from the defining equations, integer-branch structure, branch-point expansion, and large-argument asymptotics in NIST DLMF section 4.13, with Halley refinement written specifically for SciFort. SciPy 1.17.0 public documentation and public numerical outputs were used only for API conventions and regression references; no SciPy, Boost, mpmath, Cephes, or other implementation source was copied, translated, or adapted. Difficult branch-point behavior was cross-checked against the defining equations and high-precision mpmath values.

- `src/special/associated_legendre.f90`: original SciFort implementation of the
  real-degree, integer-order Ferrers function `lpmv`, integer-degree real/complex
  `legendre_p`/`assoc_legendre_p`, Legendre second-kind sequences `lqn`/`lqmn`,
  the normalized spherical functions `sph_legendre_p`/`sph_harm_y`, and the
  corresponding current SciPy-style `*_all` first-kind table APIs. The
  hypergeometric representations, integer-order relations, branch-factor
  products, spherical normalization, and recurrences were derived independently
  from NIST DLMF sections 14.3, 14.6, 14.10, and 14.30. The Legendre-Q exterior
  anchor uses DLMF 14.3.7 and its derivatives use the DLMF 14.10 identities;
  no Zhang-Jin/specfun implementation was imported. The first-kind derivative
  helpers and the ordinary-Legendre table derivatives are original
  truncated-Taylor-jet differentiations of their recurrences. The Fortran
  table layout uses native `-m:m` signed-order bounds and was designed for
  SciFort rather than copied from NumPy storage. SciPy public documentation
  and public numerical outputs were
  consulted only for API names, argument ordering, branch conventions,
  normalization, endpoint conventions, derivative ordering, and regression
  references. Difficult noninteger-degree cases were cross-checked against
  70-digit mpmath values. No SciPy, Zhang-Jin, specfun, mpmath, Boost, Cephes,
  or other third-party implementation source was copied, translated, or adapted.

- `src/special/bessel_integrals.f90`: original SciFort implementation of the
  order-zero cylindrical/modified-Bessel integral pairs `itj0y0`, `it2j0y0`,
  `iti0k0`, and `it2i0k0`. The convergent series and tail identities were
  derived independently from NIST DLMF sections 10.22 and 10.43. Moderate
  arguments use fixed Gauss-Legendre continuation over SciFort's existing
  J0/Y0/K0 kernels; oscillatory large-argument tails use recurrences derived
  from the Bessel differential identities. SciPy 1.17.0 public documentation
  and outputs were consulted only for API/domain compatibility and regression
  comparisons. High-precision mpmath quadrature/series values were used as
  independent numerical references. No SciPy, Zhang-Jin/specfun, mpmath,
  Boost, Cephes, or other implementation source was copied, translated, or
  adapted.
- `src/special/struve.f90`: original SciFort implementation of real Struve H,
  modified Struve L, and the historical order-zero integral helpers. The
  defining power series, Bessel-series representation, negative half-integer
  identities, large-argument H-Y and L-I difference asymptotics, and the
  `integral_0^infinity H_0(t)/t dt = pi/2` normalization were derived
  independently from NIST DLMF Chapter 11, especially sections 11.2, 11.4,
  11.6, and 11.7. The eight-point Gauss-Legendre constants used for continuation
  were generated from the defining Legendre polynomial and checked against
  exact moments through degree 15. SciPy 1.17.0 public documentation/outputs
  were consulted only for API/domain compatibility and regression references;
  difficult cases were independently checked with high-precision mpmath. No
  SciPy, Zhang-Jin/specfun, mpmath, Boost, Cephes, or other implementation
  source was copied, translated, or adapted.
- `src/special/spherical_bessel.f90` (`riccati_jn`, `riccati_yn`): original
  SciFort thin sequence wrappers around the existing independent spherical-J/Y
  kernels, using the defining Riccati-Bessel products and exact product-rule
  derivatives. SciPy public documentation was consulted for names, ordering, and
  sequence shape only. The x=0 singular Y behavior intentionally uses IEEE
  infinities rather than legacy finite sentinels.

- `src/special/integer_sequences.f90`: original SciFort implementation of Bernoulli numbers,
  Euler/secant numbers, and Stirling numbers of the second kind. Bernoulli evaluation uses the
  standard relation `B_(2m)=(-1)^(m+1) 2 (2m)! zeta(2m)/(2*pi)^(2m)` together with the
  `B_1=-1/2` convention (DLMF Chapter 24). Euler numbers use the secant generating-function
  recurrence `sum binomial(2m,2j) E_(2j)=0` and, only beyond the representable recurrence
  range, the standard Dirichlet-beta magnitude identity. Stirling numbers use the defining
  recurrence `{n,k}=k*{n-1,k}+{n-1,k-1}` from standard combinatorics references such as
  Graham, Knuth, and Patashnik, *Concrete Mathematics*, Chapter 6. SciPy 1.17.0/1.18 public
  documentation and public numerical outputs were consulted only for API names, conventions,
  and regression comparisons. High-order reference values were independently checked with
  mpmath and exact integer recurrence. No SciPy, Zhang-Jin/specfun, Cephes, Boost, mpmath,
  or other third-party implementation source was copied, translated, or adapted.

- `src/special/wright_bessel.f90`: original SciFort implementation of Wright's generalized
  Bessel function and its logarithm on nonnegative real arguments. The implementation uses the
  defining DLMF 10.46.1 series `Phi(a,b;x)=sum x^k/(k! Gamma(a k+b))`. Because all supported
  terms are nonnegative, it locates the continuous saddle with the digamma log derivative,
  sums outward from the largest term in scaled logarithmic form, and uses a local Laplace
  expansion with polygamma derivatives only for exceptionally large modal indices. SciPy
  1.17.0/1.18 public documentation and numerical outputs were consulted for API/domain
  conventions and regression values. High-precision checks used mpmath evaluation of the
  defining series and the identity `Phi(1,1;x)=I_0(2 sqrt(x))`. No SciPy, Boost, Cephes,
  mpmath, or other third-party implementation source was copied, translated, or adapted.

- `src/special/bessel_misc.f90`: original SciFort implementation of the Jahnke-Emden
  Lambda sequence and the weighted Bessel integral exposed by SciPy as `lmbda` and
  `besselpoly`. Lambda values use the mathematical identity
  `Lambda_v(x)=0F1(;v+1;-x^2/4)`, equivalent to the normalized Bessel-J definition. The
  `besselpoly` routine starts from the documented defining integral and uses original
  endpoint smoothing, analytic subtraction of the leading integrable singular term, and
  composite 16-point Gauss-Legendre quadrature. SciPy 1.17.0 public documentation and
  numerical outputs were consulted for API names, sequence ordering, and regression
  comparisons only. Difficult integral cases were checked independently with mpmath and
  `scipy.integrate.quad` applied to the defining integral. No SciPy, Cephes, Zhang-Jin,
  mpmath, or other implementation source was copied, translated, or adapted.

- `src/special/erf_zeros.f90`: original SciFort implementation of first-quadrant nonzero
  complex zeros of `erf(z)`. Initial guesses follow NIST DLMF 7.13.1-7.13.2 and are refined
  with Newton iterations using SciFort's existing independently implemented complex error
  function and the exact derivative `2*exp(-z^2)/sqrt(pi)`. SciPy public documentation and
  public numerical outputs were consulted only for the API name, ordering convention, and
  low-index regression references. No SciPy or Zhang-Jin/specfun implementation code was
  copied, translated, or adapted.

- `src/special/kelvin.f90`: original SciFort implementation of the order-zero
  Kelvin functions `ber`, `bei`, `ker`, and `kei`, their first derivatives, the
  combined `kelvin` evaluator, and positive root tables. The defining identities
  `ber(x)+i bei(x)=J_0(x exp(3 pi i/4))` and
  `ker(x)+i kei(x)=K_0(x exp(pi i/4))`, together with the cylindrical-Bessel
  asymptotic expansions, were derived independently from NIST DLMF sections
  10.61 and 10.67. Small/moderate arguments use independently written complex
  power series; root tables use original sign-bracketing/bisection code. SciPy
  1.17.0 public documentation and public numerical outputs were consulted only
  for API names, real-domain conventions, and regression values. High-precision
  mpmath Bessel evaluations were used as independent numerical references. No
  SciPy, Zhang-Jin/specfun, Cephes, Boost, mpmath, or other implementation source
  was copied, translated, or adapted.
- `src/special/mathieu.f90`: original SciFort implementation of the integer-order
  angular Mathieu characteristic values, Fourier coefficients, normalized
  `ce_m`/`se_m` functions, and the four modified/radial Mathieu functions.
  The four parity-class tridiagonal eigenproblems, Fourier normalizations,
  negative-parameter coefficient symmetries, and radial Bessel cross-product
  expansions were derived independently from NIST DLMF sections 28.4 and 28.20-28.23
  (<https://dlmf.nist.gov/28.4>, <https://dlmf.nist.gov/28.20>,
  <https://dlmf.nist.gov/28.23>). SciFort's existing pure symmetric Jacobi
  eigensolver supplies the finite Fourier eigenpairs. SciPy 1.17.0 public
  documentation and public numerical outputs were consulted only for API names,
  the degree-valued angular-input convention, and moderate-range regression
  references. No SciPy, Zhang-Jin/specfun, Cephes, Boost, or other implementation
  source was copied, translated, or adapted.
- `src/special/spheroidal.f90`: original SciFort implementation of prolate and
  oblate spheroidal characteristic values `pro_cv`/`obl_cv`, sequence helpers,
  first-kind angular functions `pro_ang1`/`obl_ang1`, and first-kind radial
  functions `pro_rad1`/`obl_rad1` plus their `_cv` variants, together with
  second-kind radial functions `pro_rad2`/`obl_rad2` and their `_cv` variants.
  The finite even/odd tridiagonal eigenproblems, Ferrers-function expansion, and
  normalization were derived independently from NIST DLMF sections 30.1, 30.3,
  30.8, 30.11, and 30.16 (<https://dlmf.nist.gov/30.1>, <https://dlmf.nist.gov/30.3>,
  <https://dlmf.nist.gov/30.8>, <https://dlmf.nist.gov/30.11>,
  <https://dlmf.nist.gov/30.16>). SciFort converts
  the DLMF recurrence matrix to a diagonally similar symmetric tridiagonal matrix
  for characteristic values, while the angular functions recover relative
  eigenvector coefficients with split continued-fraction recurrences and sum the
  Ferrers series with compensation. The radial first-kind implementation uses DLMF
  30.11.3 and 30.8.6 to form a scaled `a_n,k^-m` spherical-Bessel expansion, with
  an analytic oblate-origin limit. The second-kind implementation uses the same
  independently reconstructed coefficients with the DLMF spherical-Bessel-Y
  expansion at a stable anchor, then independently written adaptive
  Dormand-Prince continuation of the spheroidal differential equation where
  direct Y-series summation is ill-conditioned near the prolate boundary or
  oblate origin. The embedded RK5(4) tableau follows J. R. Dormand and
  P. J. Prince, "A family of embedded Runge-Kutta formulae", *Journal of
  Computational and Applied Mathematics* 6 (1980), 19-26,
  <https://doi.org/10.1016/0771-050X(80)90013-3>. The radial Wronskian supplies
  an independent normalization check.
  SciPy 1.17.0 public documentation and public numerical outputs were consulted
  only for API names, argument conventions, the characteristic-value shift,
  normalization convention, and regression references. SciPy/specfun source was
  inspected only to identify its public implementation's numerical-regime switch;
  no code, constants, recurrence formulas, or small-argument algorithm from that
  source was copied, translated, or adapted. No Zhang-Jin/specfun, Cephes, Boost,
  or other implementation source was copied, translated, or adapted.
- `src/special/complex_elementary.f90`: original SciFort implementation of
  complex `ndtr`, `log_ndtr`, `expm1`, `log1p`, `xlogy`, and `xlog1py`. The
  normal-CDF identity `Phi(z) = erfc(-z/sqrt(2))/2` and the scaled left-tail
  `erfcx` rearrangement follow standard error-function definitions (NIST DLMF
  chapter 7, <https://dlmf.nist.gov/7>); the `expm1`/`log1p` Cartesian
  rearrangements are independently derived algebraic identities. SciPy 1.17.0
  public documentation and numerical outputs were consulted only for API,
  principal-branch, zero-multiplier, edge-case semantics, and regression
  references. No SciPy/xsf or other implementation source was copied,
  translated, or adapted.
- `src/special/spherical_bessel_complex.f90`: original SciFort implementation of
  complex-argument `spherical_jn`, `spherical_yn`, `spherical_in`, and `spherical_kn`
  plus first derivatives. Values use the defining half-integer cylindrical-Bessel
  relations in NIST DLMF section 10.47 (<https://dlmf.nist.gov/10.47>) over SciFort's
  independently implemented complex J/Y/I/K layer; derivatives use the adjacent-order
  identities in DLMF section 10.51 (<https://dlmf.nist.gov/10.51>). High-order spherical
  J near its turning region uses backward Miller recurrence from the same DLMF recurrence,
  normalized by low-order half-integer values. Nonzero points on the real axis dispatch
  to the established real spherical kernels. SciPy 1.17.0 public
  documentation and numerical outputs were consulted only for API conventions, complex
  zero behavior, and regression references. No SciPy/AMOS or other implementation
  source was copied, translated, or adapted.
- `src/special/bessel_complex.f90`: original SciFort implementation of complex
  cylindrical and modified Bessel `jv`, `yv`, `iv`, `kv`, their SciPy-scaled
  variants, the complex `jn` alias, and arbitrary-order argument derivatives for
  finite real order. J/Y reconstruction uses the defining Hankel identities and
  SciFort's independently implemented scaled Hankel continuation. A direct J power
  series in the large-order/small-argument regime is derived from NIST DLMF 10.2.2
  (<https://dlmf.nist.gov/10.2>) to avoid cancellation. Complex I uses its defining
  power series from DLMF 10.25.2 (<https://dlmf.nist.gov/10.25>) when the order is
  large relative to `abs(z)` and otherwise uses the quarter-turn J connection from
  DLMF 10.27 (<https://dlmf.nist.gov/10.27>); K uses the corresponding Hankel
  quarter-turn connection. Half-plane selection preserves the principal branch and
  allows the SciPy
  `ive`/`kve` scale factors to be formed from scaled J/Hankel values. Complex
  derivatives use the shifted-order identities in DLMF 10.6.7 and 10.29.5
  (<https://dlmf.nist.gov/10.6>, <https://dlmf.nist.gov/10.29>). SciPy 1.17.0
  public documentation and numerical outputs were consulted only for API names,
  scaling conventions, branch behavior, edge cases, and regression references.
  No SciPy/AMOS or other implementation source was copied, translated, or adapted.
- `src/special/hankel.f90`: original SciFort implementation of SciPy-compatible
  `hankel1`, `hankel2`, `hankel1e`, `hankel2e`, `h1vp`, and `h2vp` for finite
  real order and real/complex arguments. Real-axis values use the defining
  identities `H1_v = J_v + i Y_v` and `H2_v = J_v - i Y_v` together with NIST
  DLMF 10.11.3-10.11.5 principal-branch continuation. For nonreal arguments,
  the implementation was derived independently from DLMF 10.27.8 (Hankel/K
  connection), the cylindrical and modified Bessel differential equations, and
  DLMF 10.6.7 (shifted-order derivatives). Scaled K and scaled Hankel solutions
  are analytically continued along constant-radius arcs using the embedded
  Dormand-Prince RK5(4) tableau already used elsewhere in SciFort; the tableau
  is from J. R. Dormand and P. J. Prince, "A family of embedded Runge-Kutta
  formulae", J. Comput. Appl. Math. 6 (1980), 19-26,
  <https://doi.org/10.1016/0771-050X(80)90013-3>. SciPy 1.17.0 public
  documentation and numerical outputs were consulted only for API names,
  scaling conventions, edge behavior, and regression references. No SciPy/AMOS
  or other implementation source was copied, translated, or adapted.
- `src/special/ellipsoidal_harmonics.f90`: the first-kind Lame-polynomial
  class matrices, diagonal symmetrization, coefficient rescaling, and
  class-prefactor evaluation are adapted from SciPy 1.17.0
  `scipy/special/_ellip_harm.pxd`
  (<https://github.com/scipy/scipy/blob/v1.17.0/scipy/special/_ellip_harm.pxd>),
  copyright Jaydeep P. Bardhan (2012), Matthew G. Knepley (2012), and
  Janani Padmanabhan (2014), under its BSD-2-Clause license. SciFort
  replaces the upstream LAPACK `dstevr` eigenproblem with the existing pure
  `symmetric_eigen_jacobi` solver. The second-kind and normalization
  routines are independently written from the public defining integral
  identities used by SciPy/NIST DLMF: endpoint singularities are removed by
  angle substitutions and the resulting smooth integrals use SciFort's own
  adaptive 16-point Gauss-Legendre quadrature. No SciPy QUADPACK/Cython
  integration code is copied or translated. SciPy 1.17.0 public numerical
  outputs are used only as regression references. The retained BSD-2-Clause
  notice is reproduced in `THIRD_PARTY_LICENSES.md`.
- `src/special/orthogonal_constructors.f90`: original SciFort implementation of
  classical orthogonal-polynomial coefficient constructors. The Legendre,
  Jacobi, Gegenbauer, Laguerre, Hermite, and Chebyshev three-term recurrences and
  shifted-family affine recurrences were derived independently from NIST DLMF
  Chapter 18 (<https://dlmf.nist.gov/18>). The monic Gegenbauer coefficients use
  the closed-form DLMF series ratio so the `alpha=0` limiting normalization is
  well defined. SciPy 1.17.0 public documentation and public `.coeffs` outputs
  were consulted only for API names, parameter domains, coefficient ordering,
  monic semantics, and regression references. No SciPy implementation source or
  other third-party polynomial-constructor code was copied, translated, or
  adapted.

- `src/optimize/nonlin.f90` (`inner_m_diag_probes`): original SciFort
  deterministic structured diagonal probing extension. It uses columns of
  the Walsh-Hadamard sign system and the algebraic identity
  `diag(J) = average(z * (J*z))` when cross terms cancel, with the previous
  coordinate estimator retained as the default. No third-party source code
  was copied or adapted for this extension. The public concept of probing
  sparse matrix diagonals is standard numerical linear algebra; the
  implementation, preflight checks, and tests are original to SciFort.

- `src/optimize/nonlin.f90` (seeded diagonal probing and verification):
  original SciFort extension implementing local-state Park-Miller
  pseudorandom Rademacher sign vectors, per-coordinate spot checks, and
  split-probe stability checks with exact-coordinate fallback. The core
  diagonal identity follows elementary stochastic diagonal estimation,
  and no third-party implementation code was copied or adapted.

- `src/optimize/nonlin.f90` (progressive diagonal probes): original SciFort
  extension to its existing Rademacher estimator. The staged running mean,
  independent-batch and coordinate checks, capped growth, and exact fallback
  are original implementations of elementary numerical estimation rules.
  The algorithm is not part of the public SciPy Krylov API. No third-party
  implementation was copied.

- `src/optimize/nonlin.f90` (rotating progressive coordinate verification):
  original SciFort extension of the progressive estimator. Successive exact
  verification samples are stratified with base-2 bit-reversed integer
  offsets, and a per-estimate logical mask retains previously sampled values.
  This is a local deterministic sampling policy, not an adaptation of SciPy
  or any third-party code. The finite-difference and analytic Jacobian-vector
  interfaces and the BSD-3-Clause adaptation notices in the source are
  unchanged.

- `src/optimize/nonlin.f90` (coverage-aware progressive exact checks): original
  SciFort implementation using a local path-compressed successor array to
  locate unchecked coordinates when the planned stratified sample is cached.
  The cyclic selection and exact-Jv reuse are independent modifications to
  the existing SciFort progressive estimator; no implementation code from
  other packages was imported. No additional license obligations are created.

- `src/optimize/nonlin.f90` (checked rank-two Woodbury preconditioner):
  original SciFort-only implementation of the elementary Woodbury matrix
  identity `(D + U*V^T)^(-1) = D^(-1) - D^(-1)*U*(I+V^T*D^(-1)*U)^(-1)*V^T*D^(-1)`.
  A pair of deterministic sign-vector products and the exact Jacobian diagonal
  identifies two coefficients per row in the known right-factor model. Two
  linearly independent structural checks and a guarded 2-by-2 inverse
  reject incompatible/ill-conditioned fits. No third-party implementation
  code was copied. The source module retains its existing SciPy BSD-3-Clause
  attribution; the new options and regression
  `test/test_optimize_root_krylov_rank_two.f90` are original SciFort additions.

- `src/optimize/nonlin.f90` (unknown-right rank-one coupling discovery):
  original SciFort-only extension to the checked Sherman-Morrison inverse.
  Exact coordinate Jacobian products are assembled once and reused to
  recover the diagonal and the off-diagonal rank-one factors. A maximum
  off-diagonal pivot, independent row/column observations, and an
  exhaustive check of off-diagonal entries identify admissible models.
  Guarded inversion and diagonal fallback reuse existing SciFort logic.
  This construction follows elementary rank-one matrix identities;
  no third-party implementation was imported. The original SciPy
  attribution in the surrounding nonlinear solver remains unchanged.
  New test: `test/test_optimize_root_krylov_rank_one_discover.f90`.

- `src/optimize/nonlin.f90` (unknown-right rank-two coupling discovery):
  original SciFort extension based on elementary rank-two matrix identities.
  A complete exact-coordinate Jacobian is used to find an off-diagonal
  2-by-2 pivot; the missing diagonal entries are masked during row/column
  factor recovery. Every off-diagonal entry is checked before the guarded
  2-by-2 Woodbury inverse is formed. This is an independently implemented
  mathematical construction, not third-party source translation or reuse;
  it leaves the existing SciPy-adapted nonlinear solver licensing and
  attribution unchanged. Test:
  `test/test_optimize_root_krylov_rank_two_discover.f90`.

  original SciFort extension based on the elementary matrix identity
  `(D + u*v^T)^(-1) = D^(-1) - (D^(-1)*u*v^T*D^(-1)) / (1 + v^T*D^(-1)*u)`.
  The vector `u` is recovered from the exact diagonal and a constant-vector
  Jacobian product for `v = ones`; two deterministic sign-vector products
  provide heuristic structural verification before enabling the update.
  A subsequent SciFort-only extension permits a known fixed right factor
  `v`: from `J*sign(v)` and the exact diagonal, solve for each row's `u`,
  normalize `v` to control scale, and apply the weighted Sherman-Morrison
  formula. This extension uses the same elementary identity, not adapted
  implementation code from SciPy or another project. The additional test is
  `test/test_optimize_root_krylov_weighted_rank_one.f90`. Tests and public
  options are original SciFort additions. The existing SciPy-adapted parts
  of this module retain their BSD-3-Clause attribution.

- `src/optimize/nonlin.f90` (adaptive diagonal refresh scheduling):
  original SciFort-only extension based on comparing consecutive accepted
  nonlinear residual norms. A maximum refresh interval bounds the age of a
  cached preconditioner, while weak progress can trigger an earlier update.
  This scheduling heuristic and its deterministic tests are original SciFort
  work; no third-party source code was imported. The surrounding SciPy-
  adapted Krylov nonlinear solver retains its existing BSD-3-Clause attribution.

- `benchmarks/benchmark_krylov_work_refresh.f90` and
  `test/test_optimize_root_krylov_work_refresh_nonlinear.f90`:
  original SciFort benchmarking and regression code using independently
  specified nonlinear diagonal and rank-two coupling functions. The
  analytic Jacobian/transpose products follow from direct differentiation
  of those test equations. No third-party implementation was imported;
  the surrounding Krylov solver and its existing attribution are unchanged.

- `src/signal_savgol.f90`: original SciFort polynomial least-squares
  construction of Savitzky-Golay convolution and interpolation weights.
  Reorthogonalized modified Gram-Schmidt QR is applied to a scaled
  Vandermonde design, followed by a triangular system for derivative
  weights. Polynomial least-squares filtering follows the standard
  Savitzky-Golay mathematical method; argument and edge-mode definitions
  were checked against the public SciPy `savgol_coeffs`/`savgol_filter`
  documentation (https://docs.scipy.org/doc/scipy/reference/generated/
  scipy.signal.savgol_filter.html). No SciPy implementation code was copied.

- `src/signal_lfilter.f90`: original SciFort implementation of the standard
  transposed direct-form-II digital filter difference equation and steady
  unit-step state recurrence. The mathematical equations and coefficient
  normalization are described publicly by SciPy's `lfilter` and `lfilter_zi`
  reference documentation (https://docs.scipy.org/doc/scipy/reference/
  generated/scipy.signal.lfilter.html). No third-party code was imported.

- `test/test_signal_savgol.f90`, `test/test_signal_lfilter.f90` and
  `tools/signal_reference/verify_filters.py`: original SciFort verification
  programs, with numerical comparisons to public SciPy 1.17.0 outputs.
  The Python tool compiles an independently authored test-only Fortran
  driver and does not introduce a SciPy dependency into the Fortran package.

- `src/signal_resample.f90`, `test/test_signal_resample.f90`, and
  `tools/signal_reference/{resample_driver.f90,verify_resample.py}`:
  independent SciFort implementation of the standard polyphase FIR
  upsample/filter/downsample operation and centered rational resampling.
  Kaiser-window sinc filtering uses the elementary Bessel I0 power series;
  causal and zero-phase FIR decimation use this module and existing
  SciFort `firwin` coefficients. The documented SciPy 1.17.0 output
  length, padding, FIR defaults, and delay-alignment conventions are
  numerical reference behavior, not third-party source code. No source
  code from `spectrum`, SciPy, or another package was copied for these
  routines. New SciFort code and tests use MIT licensing.


- `src/signal_irregular.f90`, `test/test_signal_irregular.f90`, and
  `tools/signal_reference/{irregular_driver.f90,verify_irregular.py}`:
  original SciFort implementations of independently derived weighted
  least-squares sinusoidal projections, based on the published
  generalized Lomb-Scargle model (Zechmeister and Kuerster 2009) and
  complex circular means for event-phase synchrony. Public SciPy 1.17
  `lombscargle` and `vectorstrength` function behavior was used as a
  numerical reference; no SciPy implementation code was copied.
  Implementation and tests are MIT licensed; no additional dependency.

- `src/signal_iir_extra.f90` (analog/band minimum-order extension),
  `test/test_signal_iir_order_bands.f90`, and
  `tools/signal_reference/{order_band_driver.f90,verify_iir_order_bands.py}`:
  independently authored Fortran/Python additions using the classical
  Butterworth and Chebyshev prototype-order equations, analog frequency
  transformations, bilinear prewarping, and golden-section bandstop edge
  minimization. Public SciPy API and numerical behavior were used as
  references (SciPy 1.17.0, BSD-3-Clause), including the order-selection
  equations; no SciPy implementation source was copied verbatim or vendored.
  These additions are original SciFort MIT-licensed work with no new FPM
  dependency or third-party license file requirements.

- `src/signal_elliptic.f90`, the elliptic design/order additions in
  `src/signal_iir_extra.f90`, `test/test_signal_elliptic.f90`, and
  `tools/signal_reference/{elliptic_driver.f90,verify_elliptic.py}`:
  original SciFort MIT-licensed implementation of the classical elliptic
  (Cauer) lowpass prototype and s-plane transformations. The complete
  elliptic integrals use the arithmetic-geometric mean, the incomplete
  integral uses Carlson RF duplication, and real Jacobi functions use the
  descending AGM. The degree equation follows standard elliptic filter
  theory (e.g., A. V. Oppenheim and R. W. Schafer, *Discrete-Time Signal
  Processing*, and S. J. Orfanidis, *Lecture Notes on Elliptic Filter
  Design*, https://www.ece.rutgers.edu/~orfanidi/ece521/notes.pdf).
  Public API and coefficient behavior were checked numerically against
  SciPy 1.17.0 `scipy.signal.ellip`/`ellipord` (BSD-3-Clause). No third-party
  implementation source was copied, bundled, or introduced as a runtime
  dependency. The tests and numerical reference generator are also original.

- `src/signal_fir_advanced.f90` (conditioned `firls` improvements),
  `test/test_signal_firls_rank.f90`, and
  `tools/signal_reference/verify_firls_conditioned.py`: original MIT-licensed
  SciFort improvements to the classical weighted cosine least-squares FIR
  design. Centered sine-integral formulas avoid subtraction cancellation,
  while a cyclic symmetric Jacobi eigensolver supplies a rank-revealing
  pseudoinverse of the normal-equation matrix if Cholesky fails. SciPy 1.17
  `scipy.signal.firls` was used as a numerical reference, including weighted
  in-band response checks for nearly rank-deficient problems. No SciPy source
  or additional third-party code was copied, and no dependency was added.

- `src/signal_partial_fraction.f90`, `test/test_signal_partial_fraction.f90`,
  and `tools/signal_reference/{probe_partial_fraction.f90,
  compare_signal_partial_fraction_scipy.py}`: independent MIT-licensed
  SciFort implementations of elementary polynomial division and
  multiplication, partial-fraction basis construction, Gaussian
  elimination with pivoting, and polynomial-root clustering. The
  mathematical interpretation and public numerical behavior were checked
  against the SciPy 1.17 signal-processing API reference for `residue`,
  `residuez`, `invres`, and `invresz`
  (https://docs.scipy.org/doc/scipy/reference/signal.html). SciPy source
  was not copied or translated, and no additional runtime dependency was
  introduced. Validation uses installed SciPy only as an optional test
  oracle, not as a Fortran library dependency.


- `src/signal_spatial.f90`, `test/test_signal_spatial.f90`,
  `test/test_signal_spatial_public.f90`, and
  `tools/signal_reference/{spatial_reference_driver.f90,
  compare_signal_spatial_scipy.py}`: original MIT-licensed rank-order and
  separable two-dimensional finite-impulse-response implementations. Algorithm
  behavior (zero-padded rank selection, odd footprint validation, symmetric
  edge reflection, one-based Fortran array indexing) was cross-checked with
  installed SciPy 1.17, and sparse-footprint corner cases additionally with
  `scipy.ndimage.generic_filter`. No SciPy or other third-party source was
  copied into the Fortran implementation and no library dependency was added.

- `src/signal_deconvolve.f90`, `src/signal_envelope.f90`,
  `test/test_signal_envelope_deconvolve*.f90`, and `reference/`:
  independent MIT-licensed SciFort implementations of polynomial long
  division and FFT-domain analytic-envelope filtering using the existing
  SciFort transform. Output behavior was tested against the installed
  SciPy 1.17 `scipy.signal.deconvolve` and `scipy.signal.envelope` interfaces.
  The Fortran code implements the underlying numerical algorithms without
  copying or translating SciPy source. Python reference test scripts require
  SciPy only at verification time. No new Fortran library dependencies or
  third-party source files were introduced.

- `src/signal_dpss.f90` and `src/signal_multitaper.f90`: original MIT-licensed
  implementations of the Slepian (discrete prolate spheroidal) tridiagonal
  eigenproblem, an implicit-shift QL eigensolver, Toeplitz concentration
  calculation, and multitaper periodogram averaging. Algorithm descriptions
  follow the published mathematical literature (Slepian, Bell System
  Technical Journal 57, 1978; Percival and Walden, *Spectral Analysis for
  Physical Applications*, 1993) and the published SciPy API specification:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.dpss.html
  SciPy source was not copied or translated. SciPy 1.17 and NumPy are
  optional test-time references only; no new Fortran dependency is added.

- `src/signal_dpss.f90` (selected eigenpairs, October 27, 2026),
  `test/test_signal_dpss_partial*.f90`,
  `reference/signal_dpss_partial_reference_driver.f90`, and
  `test/compare_signal_dpss_partial_scipy.py`: independently written
  MIT-licensed implementation of classical Sturm-sequence bisection,
  pivoted tridiagonal inverse iteration and Gram-Schmidt reorthogonalization
  for a *subset* of DPSS tridiagonal eigenvectors. This uses the standard
  mathematical eigensystem and general numerical linear algebra, not
  SciPy's implementation. SciPy 1.17 is used solely for optional numerical
  comparison. No SciPy code was translated, adapted or bundled; no new
  external Fortran dependency or additional third-party license applies.

- Extended `src/signal_multitaper.f90`, added
  `test/test_signal_multitaper_multichannel*.f90`,
  `test/compare_signal_multitaper_channels.py`, and
  `reference/signal_multitaper_channels_reference_driver.f90`: independent
  MIT-licensed implementations of weighted multitaper spectral matrices,
  complex/real CSD, and coherence. Mathematical definitions and frequency
  conventions are shared with NumPy/SciPy, but SciPy source code was not
  copied. Optional verification uses SciPy 1.17 DPSS windows and NumPy FFT
  operations; no additional Fortran runtime dependency or vendored code.

- `src/signal_multitaper.f90`, adaptive weighting and log-jackknife
  intervals (October 20, 2026): independent SciFort implementation of
  Thomson-style DPSS adaptive weighting, reweighted leave-one-taper-out
  log-spectral uncertainty, and Student-t pointwise intervals. Background
  methodological references include Thomson (1982, IEEE Proceedings),
  Percival and Walden (1993), and the published MNE-Python/nitime
  descriptions of adaptive weighting; their source code was not copied.
  `test/compare_signal_multitaper_adaptive.py` uses SciPy DPSS/Student-t
  and NumPy FFT only as test oracles. No new third-party Fortran sources
  or package dependencies were introduced.

### Multitaper spectrogram and harmonic F-test (October 20, 2026)

The spectral time-frequency and Thomson harmonic line-test procedures in
`src/signal_multitaper.f90` are independently written SciFort MIT code.
The DPSS and FFT kernels reuse SciFort's existing implementations; no source
was copied from external libraries and no new FPM dependency was introduced.
The line-test formula follows the taper-regression construction traditionally
attributed to Thomson (1982): regression coefficients use the zero-frequency
DPSS taper sums, and the pointwise statistic uses an approximate F(2,2K-2)
null distribution. Independent numerical validation uses SciPy 1.17 DPSS,
NumPy FFT and SciPy's F distribution, without copying those implementations.
The formulas are not a code-license dependency. The time-frequency PSD is a
SciFort extension rather than a newly ported SciPy function.

### Multitaper harmonic multiple-testing corrections (October 20, 2026)

`src/signal_multitaper.f90` now includes independent SciFort implementations
of Bonferroni, Benjamini–Hochberg, and Benjamini–Yekutieli adjustments for
the existing Thomson pointwise F-test. The algorithms use standard published
formulas and a separately implemented stable bottom-up mergesort, not
translated or copied SciPy source. SciPy 1.17 DPSS, NumPy FFT, SciPy F survival
functions, and `scipy.stats.false_discovery_control` were used strictly as
numerical reference oracles for 160 randomized cases. No new Fortran package
dependencies were introduced.

### Advanced SciPy window constructions (October 21, 2026)

`src/signal_windows_extra.f90` adds `chebwin`, `taylor`, and
`kaiser_bessel_derived`, adapting the frequency-domain Chebyshev synthesis,
Taylor-product coefficients and cosine synthesis, and KBD normalized cumulative
Kaiser weights from SciPy 1.17.0 `scipy/signal/windows/_windows.py`.
Original upstream: SciPy Developers, https://github.com/scipy/scipy,
BSD-3-Clause license (the complete notice appears in
`THIRD_PARTY_LICENSES.md`, section SciPy). The implementation uses native
Fortran arrays, SciFort's explicit validation/status convention, logarithmic
Chebyshev frequency scaling, a direct inverse DFT instead of SciPy's FFT, and
nonoverflowing log-I0 evaluation. No SciPy runtime dependency or upstream
Python source is included. `test/compare_signal_windows_advanced.py` uses
installed SciPy 1.17 as a numerical test oracle; 1,718 cases passed with
maximum absolute difference 4.82e-13.

- `src/signal_windows_extra.f90` (`general_cosine`) and
  `src/signal_spectral.f90` (DPSS-aware `get_window`) (October 21, 2026):
  independent MIT-licensed coefficient-vector cosine-series synthesis using
  the same elementary trigonometric kernel as the preexisting Blackman and
  flattop windows, plus an interface adapter that reuses SciFort's existing
  Slepian tridiagonal eigensolver. SciPy's public window conventions and the
  published DPSS definition served as references; neither SciPy nor
  Christopherson source was imported. Tests in
  `test/compare_signal_cosine_dpss.py` use installed SciPy 1.17 as an optional
  numerical oracle (1,092 comparisons). No third-party runtime dependency.

### 2-D spline coefficient evaluation (October 21, 2026)

`src/signal_spline2d_eval.f90` is an independently written MIT-licensed
implementation of tensor-product quadratic/cubic cardinal B-spline basis
weights and whole-sample mirror reflection. SciPy 1.17's public
`scipy.ndimage.map_coordinates(order=2/3, mode="mirror", prefilter=False)`
is an optional numerical oracle, not a code source; the Fortran evaluator
is built directly from the mathematical spline basis. The project adds no
third-party source or FPM runtime dependency. See
`tools/signal_reference/compare_signal_spline2d_eval_scipy.py`.

### Rank-independent spline filtering and coordinate evaluation (October 21, 2026)

`src/signal_spline_nd.f90` is an independent MIT implementation of the
cardinal B-spline interpolation equations with whole-sample mirror reflection
and a banded 1-D collocation solve along each axis. SciPy Developers,
`scipy.ndimage.spline_filter`, `spline_filter1d`, and `map_coordinates`,
https://github.com/scipy/scipy, BSD-3-Clause, provided the reference API
semantics and optional numerical comparison only. No upstream SciPy source
was copied, and the module introduces no additional runtime dependency.
See `tools/signal_reference/compare_signal_spline_nd_scipy.py`.

### Multidimensional spline boundary modes and geometric transforms (October 21, 2026)

`src/signal_spline_nd.f90` adds independent implementations of SciPy-style
boundary coordinate mappings, periodic spline collocation via a direct
Fourier-series solution, and affine/shift geometry using the existing
cardinal B-spline interpolation kernel. These algorithms were implemented
from mathematical spline formulas and publicly documented SciPy boundary
rules; no SciPy source was copied. Original reference: SciPy Developers,
`scipy.ndimage.map_coordinates`, `affine_transform`, and `shift` (BSD-3-Clause),
https://docs.scipy.org/doc/scipy/reference/ndimage.html. The independent
verification code calls installed SciPy 1.17.0 at test time only; no SciPy or
other third-party source code is distributed with the new implementation.
See `docs/signal_spline_transforms.md` for precision and prefilter limitations.

### Extended ndimage filters, morphology and exact EDT (October 22, 2026)

`src/ndimage.f90`, `src/ndimage_morphology.f90` and
`src/ndimage_distance.f90` include independently written MIT-licensed Fortran
implementations of full-dimensional correlation/convolution, Prewitt edge
filtering, discrete neighborhood order statistics, binary opening/closing,
background-hole filling, constrained geodesic propagation, and separable
Euclidean distance transforms. Rank-independent packing/scattering uses Fortran
2018 assumed-rank arrays. The exact squared-distance line transform uses the
lower envelope of parabolas described by P. F. Felzenszwalb and
D. P. Huttenlocher, "Distance Transforms of Sampled Functions," *Theory of
Computing* 8 (2012), pp. 415-428, DOI 10.4086/toc.2012.v008a019,
https://toc.cs.uchicago.edu/articles/v008a019/. The algorithm is implemented
here with native Fortran and anisotropic physical spacing, rather than by
copying code from SciPy or the authors' software. The neighborhood and
convolution formulas are standard mathematical definitions; SciPy Developers,
`scipy.ndimage` (BSD-3-Clause),
https://docs.scipy.org/doc/scipy/reference/ndimage.html, supplied API/
boundary semantics and optional numerical reference outputs only. No new
third-party code or FPM runtime dependency was imported.

Validation: `tools/ndimage_reference/compare_extended_filters_scipy.py`,
`compare_extended_morphology_scipy.py`, `compare_distance_scipy.py`, and
`compare_distance_random_scipy.py` add 349 SciPy 1.17 comparisons to the 173
original ndimage checks.

### Custom binary/grayscale morphology and EDT feature coordinates (October 22, 2026)

`src/ndimage_morphology.f90` now supports custom Boolean structuring elements
for binary morphology, connected components, propagation, and hit-or-miss
pattern matching. `src/ndimage.f90` adds independent rank-1-through-rank-15
implementations of grayscale dilation/erosion, opening/closing, morphological
gradient and Laplacian, and white/black top-hat operators with arbitrary
Boolean footprints. `src/ndimage_distance.f90` additionally propagates
nearest-background coordinates through the existing separable parabola-envelope
transform. All additions are original MIT-licensed Fortran implementations of
standard mathematical definitions. The SciPy Developers' public `ndimage`
API (BSD-3-Clause) provided reference semantics and test outputs only; no
upstream code was copied or redistributed. The expanded `ndimage_reference`
comparators add 296 SciPy 1.17 comparisons (including 128 randomized 2-D/3-D
feature-index maps) to the prior 522. There are no new runtime dependencies.

### Shifted/non-flat ndimage morphology, labeled measurements, watershed (October 22, 2026)

`src/ndimage.f90`, `src/ndimage_morphology.f90`,
`src/ndimage_measurements.f90`, and `src/ndimage_watershed.f90` add original
MIT-licensed Fortran implementations of origin-shifted binary morphology,
non-flat grayscale structuring elements, selected-label statistics, and
marker-based image foresting transform. The watershed path-cost model follows
A. X. Falcao, J. Stolfi, and R. de Alencar Lotufo, "The Image Foresting
Transform: Theory, Algorithms, and Applications," IEEE Transactions on
Pattern Analysis and Machine Intelligence 26 (2004) 19-29; standard
bucket-queue graph traversal and minimax edge differences are implemented
independently. The publicly available SciPy `ndimage` documentation and
SciPy's BSD-3-Clause reference behavior guided compatibility tests and
priority conventions; no SciPy source is embedded. The optional legacy
linear-index boundary behavior is provided solely for comparison with
SciPy 1.17. See `tools/ndimage_reference/compare_advanced_gray_scipy.py`,
`compare_binary_origin_scipy.py`, `compare_measurements_scipy.py`, and
`compare_watershed_scipy.py` for independent cross-checks. No additional
third-party dependencies were introduced.

### ndimage masked binary morphology, extrema, bounding boxes, rank footprints (October 22, 2026)

The additions to `src/ndimage_morphology.f90`, `src/ndimage_measurements.f90`,
and `src/ndimage.f90` independently implement per-iteration mask restriction,
convergence-until-stable binary morphology, extrema/position measurements,
selected-label bounding boxes, and order-statistic filters with arbitrary
Boolean footprints and shifted origins. The methods use original MIT-licensed
Fortran code built on pre-existing SciFort kernel and assumed-rank helpers;
no SciPy source code or other third-party code was copied. SciPy 1.17's
BSD-3-Clause API is used as a behavioral reference. The
`compare_regions_masks_rank_scipy.py` script supplies 8,526 additional
2D/3D scalar comparisons, and `test_ndimage_regions_masks.f90` covers invalid
inputs, C-order ties, zero/negative iterations, and rank-15 singleton axes.
No additional build or runtime dependencies were introduced.


### ndimage binary stencils, selected axes, numeric labels and 1D filters (October 23, 2026)

The implementations in `src/ndimage_morphology.f90` and `src/ndimage.f90`
include original MIT-licensed Fortran routines for generating radius-one
connectivity structures, iterating arbitrary structuring elements by Boolean
Minkowski addition, interpreting nonzero numeric inputs as connected-component
foreground, reducing generic derivatives over selected axes, and performing
axis-only uniform/minimum/maximum filters. No SciPy source code was copied.
The SciPy 1.17 `ndimage` public API and independently generated numerical
results were used for behavioral tests, including 600 structural/label/axis
comparisons and 42,000 1D filter comparisons. There are no new runtime
build dependencies.


### ndimage axis-selective filtering and linear-time extrema (October 23, 2026)

`src/ndimage.f90` independently implements a monotonic-deque sliding-window
minimum/maximum calculation and adds optional one-based axis subsets to
Gaussian, uniform, minimum and maximum filtering. The deque implementation
uses a standard sliding-window algorithm; it is not translated from SciPy
or any upstream implementation. Per-axis parameter expansion, boundary
handling, validation, and native rank-15 packing remain SciFort code.
Compatibility is tested against SciPy 1.17 in
`tools/ndimage_reference/compare_selected_axes_scipy.py` and
`compare_gaussian_axes_scipy.py`, and independently in the native
`test/test_ndimage_selected_axes_deque.f90` suite. The underlying SciPy
`ndimage` API is used as a reference (BSD-3-Clause); no source code is copied.


### ndimage selected-axis order statistics and sliding uniform means (October 23, 2026)

`src/ndimage.f90` now maps compact selected-axis rank/median/percentile
footprints into full-rank image coordinates and applies constant-time-per-
output rolling updates for uniform means along each line. The standard
sliding-sum algorithm and coordinate expansion were independently written
in Fortran; no SciPy source code was copied. SciPy 1.17 is used as a numerical
and API behavioral reference, including its ascending-axis order for rank
filter stencil sizes and its original axis order for origins. The
`compare_order_axes_running_scipy.py` suite validates 10,740 scalar/complex
values and `test_ndimage_order_axes_running.f90` exercises rank-15, invalid
inputs and noncontiguous output arrays. No new runtime dependencies.

### Selected-axis grayscale morphology and generic callbacks (October 23, 2026)

`src/ndimage.f90` extends the existing independently implemented stencil
traversals by mapping compact selected-axis offsets into full-rank Fortran
array coordinates. The implementation follows SciPy 1.17's public API
behavior for reordered axes, even-width grayscale dilation origin shifts,
and C-order callback samples without copying upstream implementation code.
The `tools/ndimage_reference/compare_gray_generic_axes_scipy.py` reference
checks 3,840 real/complex coefficients against SciPy and separate mixed-mode
coordinate-space constructions; `test/test_ndimage_gray_generic_axes.f90`
exercises rank-15 arrays, noncontiguous results and invalid inputs. SciPy is
used for validation and is not a runtime dependency. MIT licensing unchanged.

### ndimage complex line filtering, binary axis subsets, and value indices (October 23, 2026)

`src/ndimage.f90` extends the independently implemented padded line-callback
algorithm to complex inputs, and `src/ndimage_morphology.f90` maps compact
structuring elements onto selected image axes. `src/ndimage_measurements.f90`
implements an original heapsort-based compressed coordinate index for integer
labels. These are SciFort MIT-licensed routines; no SciPy source was copied.
`compare_complex_line_binary_axes_scipy.py` tests binary operations against
SciPy 1.17 and checks complex callback outputs by recombining real SciPy
line filters. `compare_value_indices_scipy.py` verifies label grouping and
C-ordered coordinates against SciPy's `value_indices`. This work introduces
no new dependencies.

### Periodic spline prefilter and shaped-coordinate maps (October 23, 2026)

`src/signal_spline_nd.f90` replaces the previous quadratic Fourier-series
solution for periodic spline interpolation with a linear-time causal and
anticausal recurrence based on the standard cardinal B-spline characteristic
poles, as described by M. Unser, A. Aldroubi and M. Eden, "Fast B-Spline
Transforms for Continuous Image Representation and Interpolation", IEEE
Transactions on Pattern Analysis and Machine Intelligence 13(3), 1991,
pp. 277-285. SciPy's `scipy/ndimage/src/ni_splines.c` (SciPy Developers,
BSD-3-Clause, https://github.com/scipy/scipy) supplied boundary behavior and
numerical constants used for cross-checking. The routines are implemented in
Fortran, with independent high-rank array packing, status validation, and
coordinate-grid reshaping; no SciPy binary or other dependency is included.
`tools/signal_reference/compare_spline_periodic_grid_scipy.py` verifies
15,428 periodic coefficient and real/complex mapped values against SciPy 1.17,
while `test/test_signal_spline_periodic_grid.f90` checks exact knot
reconstruction for short lines and noncontiguous coordinate grid outputs.
For `reflect` on very short arrays, SciPy's finite-precision recurrence may
not exactly reconstruct the original samples; SciFort retains its existing
exact finite collocation solve rather than reproducing that numerical error.

### Reusable spline plans and mixed-type batched filters (October 23, 2026)

`src/signal_spline_nd.f90` packages SciFort's existing independent cardinal
B-spline prefilter/evaluator into a reusable, privately owned coefficient
plan. No new spline coefficient formulas were imported. Its mathematical
semantics were compared with SciPy 1.17.0 `ndimage.map_coordinates` using
`tools/signal_reference/compare_spline_plan_scipy.py`.

`src/ndimage.f90` extends the existing original window-gather implementation
with real-to-complex and complex-to-real callback overloads. It shares the
existing boundary and origin validation, bounded batch size, and shape checks;
there is no imported SciPy algorithm or source file. The numerical references
in `tools/ndimage_reference/compare_vectorized_mixed_scipy.py` use SciPy
1.17.0 for validation only. The changes retain SciFort's MIT license and
introduce no dependencies.

### Separable spline tensor-grid evaluation (October 23, 2026)

`src/signal_spline_nd.f90` now provides `spline_axis_grid` and
`evaluate_spline_plan_tensor_nd`: an original sequence of one-dimensional
cardinal B-spline contractions over an already prepared interpolation plan.
The existing `support` routine supplies weights and boundary folding; no
third-party code is copied. This avoids allocating a full coordinate mesh
and repeated tensor-kernel calculations for regular grids. Validation uses
SciPy 1.17.0 `ndimage.map_coordinates` through
`tools/signal_reference/compare_spline_tensor_scipy.py` and native
`test/test_signal_spline_tensor.f90`, including complex coefficients and
singleton axes. The MIT license and FPM dependencies are unchanged.

### Batched reusable affine and geometric spline transforms (October 23, 2026)

`src/signal_spline_nd.f90` adds original bounded-memory coordinate generation
and plan-based affine/geometric evaluation using SciFort's existing spline
coefficients, interpolation kernels, and boundary handling. For a fixed input,
this avoids repeated prefiltering and full coordinate tables; no SciPy source
is incorporated. SciPy 1.17.0 `ndimage.affine_transform` is used only as an
independent numerical reference in
`tools/signal_reference/compare_spline_plan_transforms_scipy.py`. The
new implementation remains MIT-licensed with no new runtime dependencies.

### Batched nonlinear spline coordinate callbacks (October 24, 2026)

`src/signal_spline_nd.f90` adds `ndimage_coordinate_batch_map` and
`evaluate_spline_plan_geometric_batch_nd`, an original Fortran callback
interface which evaluates mappings on bounded blocks of coordinate columns.
It reuses SciFort's existing interpolation plan and kernels; no third-party
source was imported. `tools/signal_reference/compare_spline_plan_batch_mapping_scipy.py`
compares the numerical results to SciPy 1.17 `ndimage.map_coordinates`;
this external package is needed only for reference tests. The change retains
SciFort's MIT license and requires no additional FPM dependency.

### Dense displacement-field spline-plan transforms (October 24, 2026)

`src/signal_spline_nd.f90` adds flat and shaped displacement-field overloads
for its existing spline interpolation plan. The implementation independently
combines one-based output coordinates, optional zero-based affine maps, and
additive real-valued displacement components, then reuses the existing
batched spline interpolation kernels. No source code from SciPy was copied.
`tools/signal_reference/compare_spline_plan_displacements_scipy.py` uses
SciPy 1.17.0 `ndimage.map_coordinates` solely as an independent numerical
reference. The code and regression tests remain MIT-licensed and require
no additional runtime or FPM dependency.

### Bilinear transforms and elliptic prototypes (October 24, 2026)

`src/signal_bilinear.f90` independently implements polynomial bilinear
substitution and the standard zero/pole mapping, reusing SciFort's existing
frequency-transform result types. `src/signal_ellipap.f90` exposes the
existing SciFort AGM/Jacobi elliptic prototype algorithm through the
`analog_prototype_result` interface; it introduces no additional elliptic
algorithm dependency. The Bessel root-solver order limit is extended from 12
to 13 after independent SciPy checks. This original MIT-licensed work does
not copy SciPy source; SciPy 1.17.0 supplies numerical reference results in
`tools/signal_reference/compare_signal_bilinear_ellipap.py`. No new FPM
dependencies are introduced.

### Higher-order Bessel filter roots and analog/digital design (October 24, 2026)

`src/signal_analog_prototypes.f90` retains SciFort's original safeguarded
Aberth-Ehrlich root refinement and Campos-Calderon starting formulas but
replaces direct evaluation of high-degree reverse-Bessel-polynomial
coefficients with SciFort's `scifort_bessel_complex::kve_complex` for the
zero-equivalent function `exp(1/x)*K_(n+1/2)(1/x)`. Derivatives use the
standard modified Bessel K neighbor-order identity, equivalent to NIST DLMF
10.29.2. This algorithm follows the published strategy described by SciPy's
`scipy.signal.besselap` documentation (SciPy BSD-3-Clause; used as a
mathematical and numerical reference, not imported source). The starting
approximations are attributed to R. G. Campos and M. L. Calderon,
"Approximate closed-form formulas for the zeros of the Bessel Polynomials"
(2011), https://arxiv.org/abs/1105.0957. Original SciFort Fortran routines
and tests remain MIT-licensed; no code, files, or license notices were copied
from SciPy or the paper.

`src/signal_bessel_design.f90` independently composes the existing SciFort
analog prototype, the four existing ZPK frequency substitutions, the bilinear
transform, and ZPK-to-BA/SOS conversions. The API follows the public SciPy
Bessel filter normalization specification
(https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.bessel.html)
without copying its implementation. `scifort_signal_iir_extra` adds Bessel
dispatch to the existing IIR selector; `scifort_signal_analog_sos` adds
phase-normalized Bessel ZPK and analog SOS design using the same transforms.
Tests compare coefficients and stable
frequency responses with SciPy 1.17.0; SciPy is a test-time reference only.
No new external dependency or third-party source was added.

### Centered ShortTimeFFT and dual-window methods (October 24, 2026)

`src/signal_short_time_fft.f90` implements centered, windowed Fourier slices,
inverse overlap-add with canonical dual windows, tight-window construction,
and associated spectral calculations using original Fortran code. Numerical
semantics are compared with SciPy 1.17.0 `scipy.signal.ShortTimeFFT` and its
published time-frequency documentation at
https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.signal.ShortTimeFFT.html.
No third-party source code is included. The new module uses SciFort's existing
`forward_fft` and spectral-window implementation; non-power-of-two FFTs use
its existing direct DFT fallback. No additional runtime dependency is needed.

### Multidimensional ShortTimeFFT and REML selection (October 24, 2026)

`src/signal_short_time_fft.f90` now applies the original one-dimensional
transform independently across the orthogonal axes of rank-two and rank-three
real/complex arrays; an explicit frequency-axis/time-axis permutation
reconstructs rank-two/three signals. The documented dimensional conventions
are checked against SciPy 1.17.0 and channelwise SciPy 1-D reference results.

`src/signal_whittaker.f90` computes the profiled Gaussian REML criterion
using original binary64 Fortran code; its formula is documented by
SciPy 1.18.0 `scipy.signal.whittaker_henderson` (BSD-3-Clause) and Biessy,
*Whittaker-Henderson Smoothing Revisited* (2023). The implementation does
not copy code from SciPy. SciPy and NumPy are needed only for optional tests;
there are no new Fortran runtime or FPM dependencies.

### 2026-10-24: ShortTimeFFT rank-four/five extension

- `src/signal_short_time_fft.f90` now handles rank-four/five signal batches
  through independent Fortran multidimensional index traversal and SciFort's
  existing one-dimensional FFT kernel, including inverse reconstruction.
- The dimensional layout and numerical outputs were checked against SciPy
  1.17.0 (`scipy.signal.ShortTimeFFT`); no SciPy source code is incorporated.
- New regression sources use the same MIT project license.

## Discrete-time LTI system objects (2026-10-25)

`src/signal_dlti.f90` provides a SciPy-inspired discrete-time LTI object API,
including representation conversions, real/complex MIMO simulations, and
series/parallel state-space compositions. The implementation reuses SciFort's
existing Fortran numerical kernels; no external source code was copied.
Its behavior was compared with SciPy 1.17.0 state-space calculations and
response functions using new native Fortran and Python reference tests.

### Rank-independent periodogram, Welch, CSD and coherence (2026-10-25)

`src/signal_spectral_nd.f90` and the corresponding generator/tests are new
MIT-licensed Fortran adapter code over SciFort's existing one-dimensional
spectral methods. Their numerical formulas and frequency layout are based on
standard spectral estimators and checked against SciPy 1.17 (`scipy.signal`),
without incorporating SciPy source code. There is no new library dependency.

### Rank-independent STFT and spectrogram (2026-10-25)

`src/signal_timefreq_nd.f90` provides new MIT-licensed Fortran adapters over
SciFort's existing one-dimensional real/complex time-frequency routines.
Their shape semantics, spectrogram modes, and frequency/time coordinates
follow the published SciPy `stft` and `spectrogram` interfaces. Numerical
validation uses SciPy 1.17 as a test-only reference, without copying SciPy
source code. No new third-party dependencies are required.

### Colored-noise multitaper background and Holm significance (2026-10-26)

`src/signal_multitaper.f90` now includes original Fortran median-based local
PSD background estimation and Holm step-down adjusted p-values. This change
uses the standard running-order-statistic definition of a median and the
standard Holm (1979) sequentially rejective multiple-test formula; it does
not import third-party library source. Numerical conventions were checked
independently against SciPy 1.17 DPSS/FFT/median filtering and against
Statsmodels 0.14.6 Holm adjustments. These Python packages are optional
verification tools, not runtime or FPM dependencies. The numerical
background is descriptive rather than a new hypothesis-test calibration.


### Cost-aware FFT convolution padding (2026-10-26)

The 2/3/5/7-smooth padding search and multidimensional convolution cost
heuristic in `src/signal_fftconvolve.f90` are new MIT-licensed Fortran code.
They reuse the pre-existing SciFort radix-2 and mixed-radix transforms and
apply a conservative, locally benchmarked work-array cost estimate. The
mathematical convolution and Fourier conventions are unchanged. Numerical
verification uses independent SciPy 1.17 convolution and overlap-add results
without copying SciPy source. There are no new runtime dependencies.

### Specialized mixed-radix butterflies (2026-10-26)

`src/signal_spectral.f90` now applies algebraic radix-3/5/7 DFT butterfly
identities based on sums and differences of symmetrically indexed inputs.
This new MIT-licensed Fortran implementation was independently authored;
no source was copied from FFTW, SciPy, or other third-party FFT libraries.
Validation used independent SciPy 1.17 FFT results as a test-only reference,
plus direct Fourier sums in native tests. No new runtime dependencies were
introduced. `benchmarks/bench_signal_fft_butterflies.f90` records the
repeatable old-versus-new performance measurement protocol.

### FFT-based DPSS spectral-concentration ratios (2026-10-27)

`src/signal_dpss.f90` now evaluates the existing sinc Toeplitz concentration
Rayleigh quotient with a linear autocorrelation computed using the discrete
Fourier convolution/correlation identity (Wiener-Khinchin theorem) and a
zero-padded radix-2 Cooley-Tukey FFT. The direct calculation is retained as
an independent verification option. The FFT, compensated summation, and
method-selection code are original MIT-licensed SciFort work; no SciPy,
FFTW, or third-party Fortran source has been incorporated. A private FFT
routine avoids a circular `signal_spectral` module dependency. SciPy 1.17
DPSS results are test-only reference values. No new runtime dependencies.

### Multichannel irregular-sampling Lomb-Scargle (2026-10-27)

`src/signal_irregular_nd.f90` is new MIT-licensed Fortran glue written
for SciFort. It applies the project's existing weighted least-squares
`scifort_signal_irregular` kernel independently along lanes of an assumed-rank
Fortran array; the array-indexing and transactional-output implementation
is original. The underlying sinusoidal least-squares model is the
Lomb-Scargle model already documented for `signal_irregular.f90`.
`test/test_signal_irregular_nd.f90` and the new Fortran/Python numeric-reference
drivers were created for SciFort. SciPy 1.17.0 is consulted solely as an
independent numerical reference; no SciPy implementation source or other
third-party code was copied. No additional dependency is introduced.

### Shared-trigonometry and weighted/missing multichannel Lomb-Scargle (2026-10-27)

The updated `src/signal_irregular_nd.f90` retains its original MIT license.
It applies the same weighted sinusoidal least-squares model as the existing
rank-one SciFort implementation but caches sine and cosine sample values
for the shared time/frequency grid and reuses weighted Gram calculations
across channels. For per-channel weights and omitted samples, independent
calls use the existing scalar kernel with explicitly masked data.
This is original SciFort implementation work and adds no vendored code or
runtime dependencies. The underlying Lomb-Scargle model provenance is
recorded above for `src/signal_irregular.f90`. SciPy 1.17 is used solely as
an independent numerical reference, using one-dimensional calls after
explicit masking; no SciPy source code has been copied.

### Reusable weighted Lomb-Scargle plans (2026-10-27)

The `lombscargle_plan` type and planned rank-1..15 adapters in
`src/signal_irregular_nd.f90` are original MIT-licensed SciFort Fortran
work. They precompute a weighted, phase-rotated trigonometric design and
weighted Gram denominators, then reuse the fixed design for independent
signal records. The underlying weighted sinusoidal least-squares model is
already attributed in the earlier Lomb-Scargle provenance sections; no
third-party source was imported or adapted. The native regression,
reference driver, SciPy 1.17 comparison harness, and benchmark were
written for SciFort. SciPy's numerical output is used only for independent
reference values. No new runtime dependencies were introduced.

## Rank-independent analytic envelope adapter (2026-10-27)

`src/signal_envelope.f90` adds original MIT-licensed Fortran 2018 assumed-rank
`envelope_into` adapters around SciFort's existing FFT-based envelope kernel.
The adapters pack time-axis lanes in column-major order, call the existing
numerical routine, and scatter the fully staged results only after all
channels succeed. This is an independent SciFort interface extension; no
SciPy implementation code or third-party library was imported.

The public `scipy.signal.envelope` mathematical behavior (FFT-bin bandpass,
residual modes, optional resampling and squared magnitude) is described at
https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.envelope.html .
Verification uses the installed SciPy 1.17 package in
`reference/compare_signal_envelope_nd_scipy.py` and a matching original
Fortran reference driver. Native tests are in
`test/test_signal_envelope_rank15_into.f90`. No new runtime dependencies or
third-party licensing obligations are introduced.

## Spectral multiple coherence and conditional power (2026-10-29)

`src/signal_multiple_coherence.f90` implements the standard multiple-coherence
and Gaussian conditional-variance *identity* using the diagonal of the inverse
spectral matrix. It is original MIT-licensed SciFort code that reuses the
existing Hermitian Cholesky solver in `scifort_signal_partial_coherence` and
its documented diagonal shrinkage and trace ridge. The mathematical identity
is a direct Schur complement (see e.g. standard block-matrix inverse theory);
no third-party implementation was imported, translated, or adapted. The new
native test, output driver, NumPy comparison and public interface test are
also original SciFort work. NumPy linear solves were used only for separate
numerical verification; no runtime dependency or license obligation is added.

## Band-integrated multivariate spectral coherence (2026-10-29)

`src/signal_band_coherence.f90` and its static/time-frequency integration,
partial/multiple coherence wrappers, direct real/complex multitaper entry
points, native/public tests and reference oracle are original MIT-licensed
SciFort contributions. The quadrature computes the analytic integral of
piecewise-linear Hermitian matrix entries with no external interpolation
library. It reuses SciFort's preexisting Hermitian Cholesky precision solver.
`numpy.interp`, `numpy.trapezoid` and `numpy.linalg.inv` are used only for an
independently written test oracle; no third-party source code was copied,
translated or vendored. This is explicitly a SciFort extension, not a
standard `scipy.signal` API. No new runtime dependencies were introduced.

## Statistic-specific canonical moving-block pilot (2026-10-30)

`src/signal_canonical_block_selection.f90` now adds original SciFort
MIT-licensed delete-one-window sensitivity estimation for conditional group
coherence and squared canonical roots. `src/signal_canonical_bootstrap.f90`
forwards optional statistic, mode and frequency targeting while keeping the
existing spectral pilot and local deterministic RNG. The new source, native
and public tests, independent NumPy reference driver/oracle and the reproducible
three-channel VAR(1) demonstration were written for SciFort. NumPy/SciPy are
used only during testing and simulation, not linked into Fortran or vendored.
The autocovariance/Bartlett selection heuristic is the same original SciFort
pilot documented in `docs/signal_canonical_auto_block_length.md`; it is not
claimed optimal or calibrated for canonical coherence.

## Band-stop order objective and arbitrary-center exponential windows (2026-10-31)

`src/signal_band_stop_obj.f90` is an independent, MIT-licensed implementation
of the standard analog band-stop order formulas using SciFort's existing
AGM-based complete elliptic integrals. It adds no runtime dependency.
`src/signal_windows_extra.f90` and `src/signal_spectral.f90` add a small
independently written extension of the exponential/Poisson window equation
to a displaced center; existing tau-only callers are unchanged. The tests
in `test/` and `tools/signal_reference/` were written for SciFort and
use SciPy/NumPy only as independently executable numerical references.
No SciPy source code was copied, vendored or translated. See
`docs/signal_band_stop_objective.md` and
`docs/signal_exponential_center.md` for conventions and limitations.

## Rectilinear-grid nearest/multilinear interpolation (2026-10-31)

`src/interpolate_regular_grid.f90`, `test/test_interpolate_regular_grid.f90`,
`tools/interpolate_reference/regular_grid_driver.f90` and its randomized
NumPy/SciPy oracle were written as original MIT-licensed SciFort code. The
implementation is the standard rectilinear-cell bracketing algorithm followed
by tensor-product linear corner weights or nearest-cell selection. The interface
and boundary semantics are matched to the published SciPy 1.17
`RegularGridInterpolator`/`interpn` API using independent randomized
comparisons. None of SciPy's or Jacob Williams's `finterp` implementation
source was translated, imported, or copied. The upstream BSD-style `finterp`
package remains an optional future FPM dependency; it is not part of this
source archive's runtime requirements.

## Multidimensional PCHIP and banded B-spline collocation (2026-10-31)

`src/interpolate_regular_grid.f90` now adds an independently implemented
Fritsch–Butland PCHIP evaluator, applied from the final grid axis to the
first to match SciPy's documented multidimensional algorithm. It also
replaces the original dense B-spline collocation factorization with a
purpose-built banded LU solver (half-bandwidth degree <= 5), based on the
standard locality and total nonnegativity of B-spline bases. The native tests,
SciPy comparison scripts, and benchmark driver are original SciFort/MIT code.
SciPy was used only as a black-box numerical reference; the historical dense
SciFort version was used only as a performance baseline. No source from SciPy,
`bspline-fortran`, `finterp`, or any other external interpolation package was
copied or incorporated. This update introduces no FPM dependencies.

## PPoly/BPoly and general Hermite interpolation (2026-10-31)

`src/interpolate_piecewise.f90` provides original MIT-licensed power-basis
and Bernstein-basis piecewise polynomial implementations based on algebraic
Horner/de Casteljau evaluation, binomial basis changes, analytic calculus,
and Hermite endpoint-derivative formulas. `ppoly_from_cubic` is an original
bridge from existing SciFort cubic polynomial coefficients. The numerical
comparison programs, drivers and native tests were written for SciFort.
SciPy 1.17's public `PPoly`, `BPoly`, and `BPoly.from_derivatives` interfaces
serve as black-box numeric references. No SciPy, `bspline-fortran`, `finterp`,
or other external interpolation implementation was copied, and this update
adds no third-party licenses or dependencies.

## 2026-10-31: barycentric interpolation and optional root backend

- `src/interpolate_barycentric.f90`: independent SciFort MIT implementation
  from the public barycentric interpolation and spectral differentiation
  formulas. SciPy 1.17 is an external numerical reference; its source was
  not copied.
- `src/interpolate_piecewise.f90`: SciFort-owned extension to accept an
  optional real-polynomial root-solver callback; original solver retained.
- `optional/polyroots_backend/`: SciFort-owned glue calling the `rpoly`
  interface from **Jacob Williams, polyroots-fortran v1.2.0**. Source:
  https://github.com/jacobwilliams/polyroots-fortran ; permissive BSD-style
  license, with additional upstream algorithm-specific notices. Dependency
  is optional and version-pinned; no external source has been copied into
  SciFort. Upstream FPM manifest links system LAPACK and BLAS; no successful
  upstream installation or benchmark is asserted here.
- `tools/interpolate_reference/polyroots_rpoly_test_stub.f90` is a
  SciFort-owned **test double**, not an algorithm implementation or a
  replacement for the third-party library.

## Verified optional polyroots-fortran build (2026-10-31)

`tools/interpolate_reference/check_optional_polyroots_local.py` and
`optional/polyroots_backend/test/test_polyroots_real_stress.f90` are
SciFort-owned MIT-licensed build/test code. They compile Jacob Williams's
**separately supplied** BSD-style `polyroots-fortran` source in a temporary
directory and do not redistribute it. The upstream ZIP includes a license
with attribution to additional algorithm contributors. Anyone redistributing
the upstream source or binary must retain all applicable notices. The
SciFort optional FPM dependency remains a pinned separate upstream package.

## Global interpolation: Krogh and Floater–Hormann (2026-10-31)

`src/interpolate_global.f90` is newly authored, MIT-licensed SciFort code.
Its confluent divided differences, differentiated Newton/Horner recurrence,
and Floater–Hormann rational barycentric formulas are standard mathematical
constructions. SciPy 1.17's `KroghInterpolator` and
`FloaterHormannInterpolator` are independent black-box reference calculations;
no SciPy, Jacob Williams `polyroots-fortran`, or other external implementation
was copied, incorporated, or bundled in this increment. The optional
previously validated `polyroots-fortran` root callback remains a separate
upstream FPM dependency.

## RBF and legacy interp1d algorithms (2026-10-31)

`src/interpolate_rbf.f90` and `src/interpolate_legacy.f90` are original
MIT-licensed SciFort implementations of standard radial-basis interpolation,
augmented polynomial collocation, Gaussian elimination, and one-dimensional
interval selection. The spline modes reuse SciFort's previously authored
`src/interpolate_bspline.f90`. The `scipy.interpolate.RBFInterpolator` and
`scipy.interpolate.interp1d` classes serve solely as black-box numerical
references; no SciPy or external Fortran source code was copied. This update
introduces no new FPM or BLAS/LAPACK dependencies.

## Scattered-data interpolation (2026-10-31)

`src/interpolate_scattered.f90` implements native SciFort nearest-neighbor
search, incremental Bowyer–Watson 2D Delaunay triangulation and barycentric
interpolation. It is original MIT-licensed code created for SciFort, not a
translation or copy of SciPy's Qhull implementation. SciPy 1.17 serves as
an independent numerical test reference. `src/interpolate_rbf.f90`'s
heap-based local-neighbor selection is also original SciFort code. Neither
change incorporates any third-party library or alters the FPM dependencies.

## Cached 3D circumsphere screening and general simplices (2026-10-31)

The circumsphere prefilter, Gaussian-elimination barycentric solver and
optional external-simplex constructor in `src/interpolate_scattered.f90` are
original SciFort code based on standard computational-geometry and linear
algebra formulas, not copied or adapted from SciPy/Qhull or an upstream
Fortran root-finding/interpolation library. SciPy and Qhull are independent
black-box numerical and geometry references used by the new test scripts.
No additional third-party licenses, runtime libraries or mandatory FPM
dependencies are introduced.

## 2026-10-31 native four- and five-dimensional triangulation

The generic `delaunay_nd`, facet adjacency, cached barycentric transforms,
neighbor walk and BVH fallback in `src/interpolate_scattered.f90` are
independently written SciFort code. SciPy 1.17/Qhull are used only as black-box
reference implementations in `compare_auto_nd_scipy.py`. No Qhull or other
third-party code was copied; the project's existing MIT license applies.


## FITPACK TCK and optional perazz backend (2026-11-01)

`src/interpolate_fitpack_compat.f90` is new MIT-licensed SciFort code using
existing SciFort B-spline routines. `perazz/fitpack` is only an optional FPM
dependency (Git tag 1.0.0), not copied or vendored.
`optional/fitpack_backend/src/scifort_optional_fitpack.f90` is an original
SciFort adapter based on the documented public `fitpack_curve` API; a
synthetic interface compilation test has passed, but the actual upstream
backend is not tested in this environment. Its upstream BSD-3-Clause
license and attribution are in `LICENSES/perazz-fitpack-BSD-3-Clause.txt`.

## FITPACK-style parametric curves (2026-11-01)

The nonperiodic `splprep` and parametric `splev` implementations in
`src/interpolate_fitpack_compat.f90` are original SciFort code. They
reuse SciFort's existing B-spline and weighted least-squares solvers.
The SciPy/FITPACK `splprep` behavior is used as an independent numeric
reference; no additional SciPy, FITPACK, or perazz source is incorporated.
The full optional `perazz/fitpack` dependency still requires genuine upstream
build verification. Delaunay remains unchanged.

## Offline upstream FITPACK integration verification (2026-11-01)

The SciFort-local builder and updated SciPy reference checker are original
SciFort test harnesses, written to compile the user-supplied upstream FITPACK
2.0.0 source ZIP and check the existing SciFort adapter. No upstream FITPACK
Fortran implementation files were added to SciFort. BSD-3-Clause license
notices and the optional separation are retained.

## Optional FITPACK 2.0.0 periodic/parametric adapters (2026-11-01)

`optional/fitpack_backend/src/scifort_optional_fitpack.f90` now includes
SciFort-written wrappers for FITPACK's public periodic scalar and
open/closed parametric curve types. The wrappers and accompanying tests,
comparison scripts and offline FPM configuration helper are original MIT
SciFort code. They were compiled against a user-supplied `fitpack-main.zip`
whose manifest reports version 2.0.0, and cross-checked using independent
SciPy `BSpline` evaluations. No FITPACK source was copied into this project.
Original copyright and BSD-3-Clause terms remain upstream and are listed at
`LICENSES/perazz-fitpack-BSD-3-Clause.txt`. See
`docs/interpolate_fitpack_periodic_parametric.md` for compatibility limits.

## FITPACK bivariate surface fitting and evaluation (2026-11-01)

`src/interpolate_fitpack_bivariate.f90` and new
`bivariate_smoothing_spline` / `lsq_bivariate_spline` procedures in
`optional/fitpack_backend/src/scifort_optional_fitpack.f90` are original
MIT-licensed SciFort wrapper/evaluation code. They use SciFort's 1D B-spline
kernel and the documented public API of Federico Perini's separate BSD-3-Clause
`perazz/fitpack` project, based on Paul Dierckx's original algorithms. No
third-party FITPACK implementation files have been copied into SciFort.
The adapter's workspace sizing follows the public `surfit` interface and was
validated against the *user-supplied* FITPACK 2.0.0 source. The numerical
regressions use `scipy.interpolate.bisplev` and `LSQBivariateSpline` as
independent test references; their code was not copied. See
`docs/interpolate_fitpack_bivariate.md` and
`LICENSES/perazz-fitpack-BSD-3-Clause.txt`.

## End-to-end hierarchical clustering and cophenetic correlation (2026-11-03)

The new `fcluster` dispatch, `fclusterdata` convenience operation, and
`cophenetic_correlation` in `src/cluster_hierarchy.f90` are independent
MIT-licensed SciFort implementations. They compose pre-existing native
SciFort algorithms and use the standard Pearson sample correlation formula,
with scaling to avoid intermediate overflow. Public SciPy 1.17 API
documentation and numerical reference outputs were consulted for behavior
and verification; **no SciPy or third-party implementation code was copied**.
Test sources, including independent randomized SciPy numerical comparison
scripts, are original SciFort code.

## Threshold incomplete LU (ILUT) and reusable sparse preconditioners (2026-11-05)

`scifort_sparse_linalg.spilu` is original MIT-licensed SciFort code implementing
standard row-wise no-pivot incomplete LU with magnitude dropping and per-row
fill budgets. It reuses the original SciFort sparse-row helpers and ILU(0)
factor storage/triangular solves. SciPy/SuperLU `spilu` is an independent
reference for randomized solution and residual comparisons; no SciPy,
SuperLU, or third-party solver source code is copied or bundled. This
implementation is *not* SuperLU ILUTP and lacks its pivoting/reordering.

## Complex sparse Krylov matrix-function actions (2026-11-06)

`src/sparse_linalg_complex_funm.f90` and its native/reference tests are
original MIT-licensed SciFort implementations, adapting the algorithmic
structure of SciFort's existing real `funm_multiply_krylov` routine. The
complex Arnoldi and Lanczos recurrences and restarted projected matrix
functions are standard numerical-analysis techniques. SciPy 1.17's public
function signatures, and independent `scipy.linalg` dense matrix-function
values, were consulted for API compatibility and validation. No SciPy or
third-party source code was copied; no additional license obligation or
numerical dependency was introduced.

## Selected tridiagonal eigenvalues with Sturm bisection (2026-11-08)

The extensions to `src/linalg_banded_eigen.f90`, together with the new
selection tests, comparison driver, and example, are original MIT-licensed
SciFort code. They implement the classical Sturm-sequence count from the
symmetric tridiagonal characteristic-polynomial recurrence and bisection
using Gershgorin eigenvalue bounds. Numerical concepts follow the standard
mathematical description of tridiagonal eigenproblems in Golub and Van Loan,
*Matrix Computations*, 4th edition, and SciPy's publicly documented
`scipy.linalg.eigh_tridiagonal` / `eigvalsh_tridiagonal` selection semantics
(https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh_tridiagonal.html).
SciPy 1.17.0 is used exclusively as an independent numerical test reference.
No SciPy/LAPACK implementation source or other third-party code was copied,
and there are no new license or runtime dependency requirements.

## Selected Hermitian band eigenpairs (2026-11-08)

The selected eigenpair extensions in `src/linalg_banded_eigen.f90` are
original MIT-licensed SciFort code, combining the classical Hermitian inertia
count, bisection, pivoted banded LU and inverse iteration with two-pass
reorthogonalization. The algorithmic outline follows the public numerical
literature, including the discussion in *Comparison of eigensolvers for
symmetric band matrices* (2015, https://pmc.ncbi.nlm.nih.gov/articles/PMC4617464/)
and publicly documented SciPy `eig_banded` selection semantics
(https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eig_banded.html).
No SciPy/LAPACK or other third-party source code was copied. SciPy 1.17.0
was used only as an independent comparison reference; no additional license
or mandatory numerical library dependency was introduced.

## Real symmetric LU path in banded selected eigenvectors (2026-11-08)

The `real_input` specialization in `src/linalg_banded_eigen.f90` reuses
SciFort's existing MIT-licensed real banded LU factor and solve interfaces.
It adds no copied external source or third-party dependency. Mathematical
behavior is checked against SciPy 1.17.0 and the complex-input path.

## 2026-11-08 — Band-to-tridiagonal Givens bulge chasing

`src/linalg_banded_tridiagonal.f90`, the `method='bulge'` routing in
`src/linalg_banded_eigen.f90`, native/reference tests, and associated
examples and documentation are original MIT-licensed SciFort code. The
independently derived complex Givens similarity and its out-of-band bulge
chasing are standard numerical linear-algebra methods. The public
mathematical descriptions in Golub and Van Loan, *Matrix Computations*,
4th edition (2013), and the LAPACK DSBTRD documentation
(https://www.netlib.org/lapack/explore-html/d0/d62/dsbtrd_8f_source.html)
were consulted for general algorithmic structure. No LAPACK/Netlib/SciPy
source code was copied or adapted. SciPy 1.17.0 is used only as an
independent numerical test reference. No additional third-party code,
license obligation, or external runtime dependency was introduced.

## 2026-11-08 — Real-arithmetic banded Givens reduction and auto selection

The real specialization of `src/linalg_banded_tridiagonal.f90` and automatic
method selection in `src/linalg_banded_eigen.f90` are original MIT-licensed
SciFort work. The real rotations are derived algebraically from the existing
SciFort Hermitian Givens similarity, without importing source code from LAPACK,
SciPy, or other projects. The crossover rule is a conservative heuristic from
benchmarks performed with the supplied SciFort implementation, not an
adaptation of any external performance model. Independent SciPy 1.17.0
numerical comparisons are used only for validation. No additional third-party
code, licenses, or dependencies are introduced.

## 2026-11-08 — Generalized complex Schur factors (`qz`)

`src/linalg_qz.f90`, its native tests, example, and Python reference driver
are original MIT-licensed SciFort work. The algorithm follows algebraic
generalized-Schur identities: construct a Schur decomposition of an invertible
pencil reduction (`B**(-1) A`, `A**(-1) B`, or a nonsingular shifted pencil),
then triangularize the transformed coefficient via QR. This is **not** a
translation of SciPy's or LAPACK's QZ/`gges` source. The mathematical output
convention was checked against the public documentation of `scipy.linalg.qz`
(https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.qz.html),
with SciPy 1.17.0 used only for numerical comparison. No third-party
implementation source was copied, and no additional license notices or
runtime dependencies are required.

## 2026-11-09 — Ordered complex generalized Schur (`ordqz`)

`src/linalg_ordqz.f90`, the new native tests, Fortran example and reference
harness are original MIT-licensed SciFort implementations. Adjacent triangular
pencil blocks are exchanged using a directly derived pair of two-by-two
unitary transformations: first construct a right generalized eigenvector of
the second diagonal block, then align the resulting first columns with a
left unitary rotation. Homogeneous eigenvalues provide safe region selection.
The public mathematical behavior was checked against the documentation of
`scipy.linalg.ordqz` (https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.ordqz.html).
SciPy is used only as an independent numerical validation reference. No
SciPy/LAPACK/Netlib source code was copied or adapted, and no additional
third-party licenses, sources, or runtime dependencies were introduced.

## 2026-11-09 -- Pivoted dense LDL factorization

`src/linalg_ldl.f90`, `test/test_linalg_ldl.f90`, the Fortran example, API
notes and Python reference harness are new MIT-licensed SciFort code.
The diagonal-pivoting choice and 1x1/2x2 Schur-complement updates were
independently derived from the Bunch-Kaufman matrix factorization described
in the public mathematical documentation of LAPACK `DSYTF2`:
https://www.netlib.org/lapack/explore-html/db/d8d/group__hetf2_ga91a9fdb3b98a38e1a91689b65c366569.html .
The output convention is verified against the public SciPy documentation:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.ldl.html .
No SciPy, LAPACK, or third-party code has been copied or adapted. SciPy 1.17.0
is used only as an independent test oracle. No license notices or additional
runtime dependencies were needed.

## 2026-11-09 — Polynomial companion and Fiedler matrix constructors

The `scifort_linalg_companion` and `scifort_linalg_fiedler` modules,
`test/test_linalg_companion.f90`, the Fortran example and numerical
comparison driver are original MIT-licensed SciFort code. Polynomial
companion matrices are formed directly from the standard Frobenius
first-row/subdiagonal identity; the Fiedler companion matrix uses the
standard interleaving of coefficient and unit entries, and the ordinary
Fiedler matrix uses `abs(a(i)-a(j))`. The SciPy reference implementations
were used to verify public API conventions and mathematical outputs,
not copied or translated. Independent comparisons used SciPy 1.17.0,
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.companion.html,
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.fiedler_companion.html,
and https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.fiedler.html.
No third-party source, library, license notice or external runtime
requirement was added.

## 2026-11-09 — Cosine–sine decomposition (`cossin`)

`src/linalg_cossin.f90`, `test/test_linalg_cossin.f90`, the Fortran
example, reference driver, comparison script, and documentation are new
MIT-licensed SciFort code. The implementation derives the complete CS
decomposition from the project’s original one-sided Jacobi SVD of the
upper-left partition, matrix-block identities for the complementary
subspaces, and two-pass Gram–Schmidt completion of orthonormal bases.
The public mathematical layout and input conventions were checked against
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cossin.html.
SciPy 1.17.0 is used only as an independent numerical reference. Neither
SciPy, LAPACK, nor third-party source code was copied or incorporated.
No extra dependency or licensing obligation was introduced.

## 2026-11-09 — QR rank updates, insertion and deletion

`src/linalg_qr_updates.f90`, `test/test_linalg_qr_updates.f90`,
`reference/linalg_qr_updates_driver.f90`,
`validation/compare_linalg_qr_updates.py`,
`example/linalg_qr_updates_example.f90`, and their documentation are
original MIT-licensed SciFort contributions. They use QR projection,
two-pass Gram–Schmidt augmentation and existing SciFort Householder QR.
The update operations were independently formulated from standard
orthogonal-factor identities; no SciPy, LAPACK, or third-party source
code was copied into these modules. SciPy's public API documentation was
used only to establish behavior, array dimensions and mathematical
conventions:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.qr_update.html,
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.qr_insert.html,
and https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.qr_delete.html.
SciPy 1.17.0 was used as an independent numerical validation oracle.
No new build dependencies or external licenses were introduced.

## 2026-11-09 — Adaptive randomized matrix-rank estimation

`src/linalg_estimate_rank.f90`, `test/test_linalg_estimate_rank.f90`,
`validation/reference_drivers/driver_linalg_estimate_rank.f90`,
`test/compare_linalg_estimate_rank_scipy.py`,
`example/linalg_estimate_rank_example.f90`, and documentation are original
MIT-licensed SciFort contributions. They implement randomized range sampling
with two-pass modified Gram–Schmidt. The public SciPy documentation at
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.interpolative.estimate_rank.html
was consulted for the behavior and input conventions. SciPy 1.17.0 and NumPy
2.3.5 were used only for independent validation, and their randomized
algorithms and source code were not copied. No new dependencies or third-party
license obligations were introduced.

## 2026-11-10 — Constrained and mixed-integer differential evolution

`src/optimize/differential_evolution.f90` was extended with original MIT-licensed
SciFort Fortran code for integer variable projection, linear/nonlinear constraint
callbacks, feasibility-first trial selection, and polishing hooks. Its original
Storn–Price / Ian Morrall MIT attribution remains in the file and
`LICENSES/ian-morrall-differential-evolution-MIT.txt`. No third-party code was
added for this extension. `src/optimize/differential_evolution_polish.f90` is
an original SciFort bridge to the existing L-BFGS-B implementation already
used by `src/optimize/minimize.f90`. Public API conventions were compared to
SciPy 1.17.0's `scipy.optimize.differential_evolution`, but SciPy source code
was not copied.

## 2026-11-10 — Original DIRECT deterministic global optimizer

`src/optimize/direct.f90`, `test/optimize_direct_reference.f90`,
`validation/optimize_direct_driver.f90`,
`validation/compare_optimize_direct_scipy.py`,
`example/optimize_direct_example.f90`, and documentation are original
MIT-licensed SciFort work based on published mathematical specifications:

- D. R. Jones, C. D. Perttunen, B. E. Stuckman (1993), "Lipschitzian
  Optimization Without the Lipschitz Constant", Journal of Optimization
  Theory and Applications 79(1), 157-181.
  https://doi.org/10.1007/BF00941892
- D. R. Jones and J. R. R. A. Martins (2021), "The DIRECT algorithm:
  25 years later", Journal of Global Optimization 79, 521-566.
  https://doi.org/10.1007/s10898-020-00952-6

The reference behavior and option descriptions were checked against the
public `scipy.optimize.direct` documentation, with SciPy 1.17.0 used as an
independent validation oracle. No SciPy source code, external Fortran DIRECT
implementation, or third-party numerical library was copied into this module.
No third-party license notice or dependency was added.

## 2026-11-10 — Original bounded simplicial global search (SHGO subset)

`src/optimize/shgo.f90`, `src/optimize/shgo_minimize.f90`,
`test/optimize_shgo_reference.f90`, `validation/optimize_shgo_driver.f90`,
`validation/compare_optimize_shgo_scipy.py`, `docs/optimize_shgo.md`,
`example/optimize_shgo_example.f90` and
`scripts/validate_optimize_shgo.sh` contain original SciFort work released
under the project's MIT license. No third-party implementation source was
copied, imported, or translated. The local-polishing adapter *calls* existing
SciFort Nelder-Mead and Powell routines; their existing provenance remains
unchanged.

The independent construction is based on the published algorithmic ideas:

- S. C. Endres, C. Sandrock, W. W. Focke (2018), "A simplicial homology
  algorithm for Lipschitz optimisation", *Journal of Global Optimization*
  72(2), 181-217, https://doi.org/10.1007/s10898-018-0645-y.
- Kuhn-Freudenthal triangulations of Cartesian hypercubes (classical
  combinatorial construction). Edges are generated directly from the
  coordinatewise monotone simplex definition, not third-party code.
- Public SciPy SHGO documentation was used as a behavioral reference:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.shgo.html
  (consulted 2026-11-10), but SciPy source code was not copied.

The implementation uses structured dyadic meshes rather than SciPy's adaptive
complex and does not claim the paper's full convergence guarantees. A graph
local-minimum is a birth in the zero-dimensional sublevel-set filtration;
no full persistent homology or homology-rank stopping has been implemented.
There are no new third-party dependencies or licenses for this addition.

## 2026-11-10 — SHGO selective adaptive extension

The optional dyadic cell-refinement implementation added to
`src/optimize/shgo.f90` and its public-API forwarding in
`src/optimize/shgo_minimize.f90` are original SciFort MIT-licensed work. The
vertex cache, per-cell subdivisions, deterministic best/average-corner
ranking, and local Kuhn cell-edge graph were independently implemented.
No SciPy or other source code was copied. General SHGO terminology and API
scope were checked against
https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.shgo.html .
The adaptive structure has hanging nodes and does not reproduce SciPy's
conforming triangulation or global convergence proof. New validation inputs
and test drivers in `test/optimize_shgo_adaptive_reference.f90` and
`validation/compare_optimize_shgo_adaptive_scipy.py` are original; SciPy 1.17.0
is used solely as an independent numerical reference.

## 2026-11-10 — SHGO sampled homology-growth stabilization

The optional stopping criteria and diagnostics in `src/optimize/shgo.f90`,
forwarding through `src/optimize/shgo_minimize.f90` and `src/optimize.f90`,
and connected equal-valued plateau handling are independent original
SciFort MIT-licensed additions. The new tests, example, Fortran driver,
SciPy comparison script and validation report are also original. No upstream
SciPy SHGO code was copied or adapted.

Conceptual references: Endres, Sandrock and Focke (2018),
https://doi.org/10.1007/s10898-018-0645-y, and the public SciPy SHGO
options documentation, https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.shgo.html
(consulted 2026-11-10). SciPy uses a different `minhgrd` stopping rule;
this extension monitors sampled graph-minimum counts, requires repeated
stability and a best-value guard, and does not claim topological convergence
or global optimality. No new dependencies or license obligations arise.

## 2026-11-10 — SHGO conforming planar adaptive triangulation

`src/optimize/shgo.f90`, `src/optimize/shgo_minimize.f90`,
`test/optimize_shgo_conforming_reference.f90`,
`validation/optimize_shgo_conforming_driver.f90`,
`validation/compare_optimize_shgo_conforming_scipy.py`,
`example/optimize_shgo_conforming_example.f90`, and corresponding tests
and documentation are original SciFort MIT-licensed work. No SciPy
implementation code was copied. The construction uses elementary planar
fan triangulation, exact dyadic coordinates and standard triangle-edge
incidence tests. Conceptual background: Freudenthal (1942), and Endres,
Sandrock, and Focke (2018), DOI 10.1007/s10898-018-0645-y.
The SciPy SHGO public API and numerical results are comparison references:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.shgo.html .

## 2026-11-10 — SHGO conforming three-dimensional adaptive mesh

`src/optimize/shgo.f90`, `src/optimize/shgo_minimize.f90`,
`test/optimize_shgo_conforming_3d_reference.f90`,
`validation/optimize_shgo_conforming_3d_driver.f90`,
`validation/compare_optimize_shgo_conforming_3d_scipy.py`,
`example/optimize_shgo_conforming_3d_example.f90`,
`validation/optimize_shgo_conforming_3d_2026-11-10.md`, and the corresponding
script/doc updates are original SciFort MIT-licensed implementation work.
No SciPy SHGO source code, third-party tetrahedralization implementation,
or externally licensed mesh library was copied. The implementation uses
standard convex polygon ear clipping, octree subdivisions, exact dyadic
lattice coordinates, coning tetrahedralization and simplex-incidence counting.
Mathematical background remains the Freudenthal (1942) cube triangulation
and Endres, Sandrock, Focke (2018), DOI 10.1007/s10898-018-0645-y.
SciPy SHGO API/quality is a validation comparison only:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.shgo.html .

## 2026-11-10 — SHGO inequality/equality constrained optimization

The constraint evaluator and feasibility-first simplicial ranking in
`src/optimize/shgo.f90`, and the SLSQP bridge/constraint conversion in
`src/optimize/shgo_minimize.f90`, are original MIT-licensed SciFort code.
They reuse the public vector-bound constraint convention already established
by SciFort's differential evolution implementation and the existing
`scifort_optimize_minimize::slsqp` wrapper; they do not vendor new solver
implementations or introduce new dependencies. SciPy SHGO and
NonlinearConstraint documentation were checked for compatible mathematical
meaning of equality and inequality bounds, not copied as source. See
https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.shgo.html .
All new tests, examples, and independent SciPy validation scripts are
original. Constrained core tests run without FPM; adapter tests use a
clearly labelled test-only minimizer stand-in. A full SLSQP link/test remains
unverified in this execution environment.

## 2026-11-10 — SHGO analytic Jacobians and feasibility projection

The optional nonlinear constraint Jacobian interface, canonical SLSQP bridge,
small dense Gaussian-elimination solver, damped projected Gauss-Newton
feasibility-restoration routine, and all new tests/examples/docs are original
SciFort MIT-licensed code. They implement elementary differentiation and
constrained least-squares mathematics without adapting or vendoring any SciPy
or third-party solver implementation. The conceptual SHGO background is
Endres, Sandrock and Focke (2018), DOI 10.1007/s10898-018-0645-y. SciPy's
public constrained-optimization conventions are used for interface
compatibility and independent validation only:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.shgo.html .
No new dependencies or third-party license texts are required. The bundled
SLSQP test stand-in is not a numerical implementation; actual FPM linkage
and performance of the public SLSQP backend remain to be verified.


## 2026-11-10 — SHGO independent nonlinear constraints and mixed Jacobians

The new `shgo_nonlinear_constraint` public type, metadata validator,
constraint-value concatenation and mixed analytic/finite-difference
Jacobian handling in `src/optimize/shgo.f90`, the SLSQP bridge in
`src/optimize/shgo_minimize.f90`, and all related tests, scripts, examples
and documentation are original MIT-licensed SciFort additions. There is no
new dependency and no third-party source was copied. The reference for
constraint semantics is the public SciPy SHGO and NonlinearConstraint API:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.shgo.html
https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.NonlinearConstraint.html
(consulted 2026-11-10). Results from SciPy 1.17.0 are used as independent
numerical benchmarks; comparisons independently verify feasibility because
solver success flags do not by themselves certify constraint satisfaction.
The provided SLSQP smoke-test stand-in checks constraint and Jacobian layout;
it does **not** replace compiling and testing against SciFort's real SLSQP
dependency with Fortran Package Manager.
