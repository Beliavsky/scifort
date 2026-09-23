# Third-party licenses

Most SciFort numerical source is original MIT-licensed implementation.
`CODE_PROVENANCE.md` identifies the small amount of adapted third-party code.
Standard mathematical formulas, language intrinsics, and numerical identities
do not by themselves create a third-party source-code dependency.

## SciPy

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
from SciPy 1.17.0's distributed `_sobol_direction_numbers.npz`. The applicable SciPy
BSD-3-Clause license is reproduced below.

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
